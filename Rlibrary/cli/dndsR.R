#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(cli)
  library(dndsR)
})

# =========================
# Dual-form flag alias helper
# =========================
# Given a long flag like "--output-dir", returns c("--output-dir", "--output_dir")
.alias <- function(longflag) {
  us <- sub("--", "--", gsub("-", "_", longflag, fixed = TRUE), fixed = TRUE)
  unique(c(longflag, us))
}

# =========================
# Global/base options
# =========================
base_opts <- list(
  make_option(c("-t","--threads"), type="integer", default=8,
              help="CPUs/threads [default: %default]"),
  make_option(c("-v","--verbose"), action="store_true", default=FALSE,
              help="Verbose logging")
)

# =========================
# Helpers
# =========================

# pretty dump of parsed options when --verbose
.dump_opts <- function(cmd, o) {
  if (isTRUE(o$verbose)) {
    cli::cli_h2(paste0("Running: ", cmd))
    lst <- as.list(o)
    safe_names <- names(lst)
    names(lst) <- ifelse(nzchar(safe_names), safe_names, paste0("V", seq_along(lst)))
    utils::str(lst, vec.len = 6)
  }
}

# Build arg list by adding only user-set values
.add_if <- function(lst, name, val, conv = identity) {
  if (!is.null(val) && !(is.logical(val) && is.na(val)) && !(length(val) == 1 && is.na(val))) {
    lst[[name]] <- conv(val)
  }
  lst
}

# Tri-state resolver for paired store_true flags:
# prefer 'off' when negative flag is TRUE, else 'on' when positive flag is TRUE, else NULL
.resolve_toggle <- function(pos_flag, neg_flag) {
  if (!is.null(neg_flag) && !is.na(neg_flag) && isTRUE(neg_flag)) return(FALSE)
  if (!is.null(pos_flag) && !is.na(pos_flag) && isTRUE(pos_flag)) return(TRUE)
  NULL
}

# parse --oa key=value pairs with light type casting
.cast_scalar <- function(x) {
  x_trim <- trimws(x)
  xl <- tolower(x_trim)
  if (xl %in% c("true","false")) return(as.logical(xl))
  if (xl == "null") return(NULL)       # scalar path returns NULL
  if (xl == "na") return(NA)
  if (xl %in% c("inf","+inf")) return(Inf)
  if (xl == "-inf") return(-Inf)
  if (grepl("^\\d+$", x_trim)) return(as.integer(x_trim))
  if (grepl("^\\d*\\.?\\d+(e[+-]?\\d+)?$", x_trim, ignore.case = TRUE)) return(as.numeric(x_trim))
  x_trim
}

.parse_oa <- function(v) {
  if (is.null(v)) return(list())
  out <- list()
  for (kv in v) {
    sp <- regexpr("=", kv, fixed = TRUE)
    if (sp < 1) stop("Invalid --oa argument (expected key=value): ", kv)
    key <- substr(kv, 1, sp - 1)
    val <- substr(kv, sp + 1, nchar(kv))
    if (grepl(",", val, fixed = TRUE)) {
      parts <- strsplit(val, ",", fixed = TRUE)[[1]]
      casted <- lapply(parts, .cast_scalar)
      # For vector context, treat NULL tokens as NA
      casted[vapply(casted, is.null, logical(1))] <- list(NA)
      vec_chr <- vapply(casted, function(y) {
        if (length(y) == 0) NA_character_ else as.character(y)
      }, FUN.VALUE = character(1), USE.NAMES = FALSE)
      out[[key]] <- type.convert(vec_chr, as.is = TRUE)
    } else {
      out[[key]] <- .cast_scalar(val)
    }
  }
  out
}

.parse_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  unique(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}
.parse_csv_required <- function(x) {
  v <- .parse_csv(x); if (is.null(v)) character(0) else v
}

.parse_kv_list <- function(v) {
  if (is.null(v)) return(NULL)
  out <- list()
  for (kv in v) {
    sp <- regexpr("=", kv, fixed = TRUE)
    if (sp < 1) stop("Invalid key=value: ", kv)
    key <- tolower(trimws(substr(kv, 1, sp - 1)))
    val <- trimws(substr(kv, sp + 1, nchar(kv)))
    out[[key]] <- val
  }
  out
}
.parse_kv_list_of_csv <- function(v) {
  raw <- .parse_kv_list(v); if (is.null(raw)) return(NULL)
  lapply(raw, .parse_csv_required)
}

# normalizer for separator lists
.parse_seps <- function(x) {
  v <- .parse_csv_required(x)
  if (!length(v)) return(c(";", "|", ","))  # hard default (only used if user set this flag)
  v <- trimws(v[nzchar(v)])
  v[v %in% c("comma", "<comma>", "\\,")] <- ","
  v
}

# =========================
# split_comparisons
# =========================
opts_split <- list(
  make_option(c("-C", .alias("--comparison-file")), type="character", dest="comparison_file",
              help="TSV with columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff"),
  make_option(c("-m", .alias("--mode")), type="character", default="subgenome",
              help="subgenome|haplotype|custom [default: %default]"),
  make_option(.alias("--custom-regex"), type="character", default=NULL,
              help="Regex with exactly one capture group when --mode=custom"),
  make_option(.alias("--case-insensitive"),  action="store_true", default=TRUE,
              help="Case-insensitive label matching [default: %default]"),
  make_option(.alias("--case-sensitive"),   action="store_true", default=FALSE,
              help="Override to make matching case-sensitive")
)
run_split <- function(o){
  .dump_opts("split_comparisons", o)
  stopifnot(!is.null(o$comparison_file))
  case_insensitive <- if (isTRUE(o$case_sensitive)) FALSE else isTRUE(o$case_insensitive)
  p <- dndsR::split_comparisons(
    comparison_file  = o$comparison_file,
    mode             = o$mode,
    custom_regex     = o$custom_regex,
    case_insensitive = case_insensitive
  )
  cli::cli_alert_success("Wrote: {p}")
}

# =========================
# extract_cds
# =========================
opts_extract_cds <- list(
  make_option(.alias("--comparison-name"), type = "character", dest = "comparison_name",
              help = "Unique identifier for a comparison (e.g., 'CheAl_v_CheFo')"),
  make_option(.alias("--query-fasta"), type = "character", dest = "query_fasta",
              help = "Path to query genome FASTA"),
  make_option(.alias("--subject-fasta"), type = "character", dest = "subject_fasta",
              help = "Path to subject genome FASTA (optional in single-genome mode)"),
  make_option(.alias("--query-gff"), type = "character", dest = "query_gff",
              help = "Path to query GFF3"),
  make_option(.alias("--subject-gff"), type = "character", dest = "subject_gff",
              help = "Path to subject GFF3 (optional in single-genome mode)"),
  make_option(c("-O", .alias("--output-dir")), type = "character", dest = "output_dir", default = NULL,
              help = "Output directory; one subdir per comparison [inherits function default if omitted]"),
  make_option(.alias("--overwrite"), action = "store_true", dest = "overwrite", default = NA,
              help = "Overwrite existing outputs"),
  make_option(c("-C", .alias("--comparison-file")), type = "character", dest = "comparison_file",
              help = "TSV with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff"),
  make_option(.alias("--group-by"), type = "character", dest = "group_by", default = NULL,
              help = "CDS grouping: gene|tx"),
  make_option(.alias("--export-proteins"), action = "store_true", dest = "export_proteins", default = NA,
              help = "Also write translated proteins (*.faa)"),
  make_option(.alias("--genetic-code"), type = "integer", dest = "genetic_code", default = NULL,
              help = "NCBI genetic code for translation"),
  make_option(.alias("--keep-internal-stops"), action = "store_true", dest = "keep_internal_stops", default = NA,
              help = "Keep sequences with internal stops (otherwise drop them)")
)

run_extract_cds <- function(o) {
  .dump_opts("extract_cds", o)

  args <- list(
    comparison_name = o$comparison_name,
    query_fasta     = o$query_fasta,
    subject_fasta   = o$subject_fasta,
    query_gff       = o$query_gff,
    subject_gff     = o$subject_gff
  )
  args <- .add_if(args, "output_dir",          o$output_dir)
  args <- .add_if(args, "overwrite",           if (!is.na(o$overwrite) && isTRUE(o$overwrite)) TRUE else NULL)
  args <- .add_if(args, "verbose",             if (isTRUE(o$verbose)) TRUE else NULL)
  args <- .add_if(args, "comparison_file",     o$comparison_file)
  args <- .add_if(args, "group_by",            o$group_by)
  args <- .add_if(args, "export_proteins",     if (!is.na(o$export_proteins) && isTRUE(o$export_proteins)) TRUE else NULL)
  args <- .add_if(args, "genetic_code",        o$genetic_code, as.integer)
  args <- .add_if(args, "keep_internal_stops", if (!is.na(o$keep_internal_stops) && isTRUE(o$keep_internal_stops)) TRUE else NULL)

  res <- do.call(dndsR::extract_cds, args)
  cli::cli_alert_success("extract_cds completed.")
  invisible(res)
}

# =========================
# calculate_dnds
# =========================
opts_calculate_dnds <- list(
  make_option(c("-C", .alias("--comparison-file")), type="character", dest="comparison_file",
              help="TSV (or data.frame via R) with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff"),
  make_option(.alias("--comparison-name"), type="character", dest="comparison_name",
              help="Unique comparison ID (e.g., 'CheAl_v_CheFo')"),
  make_option(.alias("--subject-fasta"), type="character", dest="subject_fasta",
              help="Subject genome FASTA (single mode)"),
  make_option(.alias("--query-fasta"), type="character", dest="query_fasta",
              help="Query genome FASTA (single mode)"),
  make_option(.alias("--subject-gff"), type="character", dest="subject_gff",
              help="Subject GFF3 (single mode)"),
  make_option(.alias("--query-gff"), type="character", dest="query_gff",
              help="Query GFF3 (single mode)"),
  make_option(c("-O", .alias("--output-dir")), type="character", dest="output_dir", default=NULL,
              help="Output directory [inherits function default if omitted]"),
  make_option(.alias("--comp-cores"), type="integer", dest="comp_cores", default=NULL,
              help="Cores passed to orthologr::dNdS via comp_cores [default: --threads]"),
  make_option(.alias("--aligner"), type="character", dest="aligner", default=NULL,
              help="orthologr aligner"),
  make_option(.alias("--sensitivity-mode"), type="character", dest="sensitivity_mode", default=NULL,
              help="orthologr sensitivity_mode"),
  make_option(.alias("--dnds-method"), type="character", dest="dnds_method", default=NULL,
              help="orthologr dnds_est.method"),
  make_option(.alias("--oa"), type = "character", action = "append", dest = "oa", default = NULL,
              help = paste0(
                "Pass through to orthologr::dNdS as key=value (repeatable). ",
                "Examples: --oa eval=1E-6 --oa ortho_detection=RBH --oa clean_folders=TRUE ",
                "--oa include_only_chr=1,2,3"
              ))
)

run_calculate_dnds <- function(o){
  .dump_opts("calculate_dnds", o)
  extra <- .parse_oa(o$oa)
  cores <- if (!is.null(o$comp_cores)) as.integer(o$comp_cores) else as.integer(o$threads)

  args <- list(
    comparison_file = o$comparison_file,
    comparison_name = o$comparison_name,
    subject_fasta   = o$subject_fasta,
    query_fasta     = o$query_fasta,
    subject_gff     = o$subject_gff,
    query_gff       = o$query_gff
  )
  args <- .add_if(args, "output_dir",       o$output_dir)
  args <- .add_if(args, "comp_cores",       cores)
  args <- .add_if(args, "aligner",          o$aligner)
  args <- .add_if(args, "sensitivity_mode", o$sensitivity_mode)
  args <- .add_if(args, "dnds_method",      o$dnds_method)

  paths <- do.call(dndsR::calculate_dnds, c(args, extra))

  if (is.character(paths) && length(paths)) {
    cli::cli_alert_success("calculate_dnds completed. Wrote {length(paths)} file{?s}.")
    for (p in paths) cli::cli_bullets(c(v = paste0("<", p, ">")))
  } else {
    cli::cli_alert_success("calculate_dnds completed.")
  }
  invisible(paths)
}

# =========================
# append_annotations
# =========================
opts_append_annotations <- list(
  make_option(.alias("--dnds-file"), type = "character", dest = "dnds_file",
              help = "Path to a dN/dS TSV (single mode)"),
  make_option(.alias("--query-gff"), type = "character", dest = "query_gff",
              help = "Path to query GFF3 (single mode)"),
  make_option(.alias("--subject-gff"), type = "character", dest = "subject_gff",
              help = "Path to subject GFF3 (single mode)"),
  make_option(c("-o", .alias("--output-file")), type = "character", dest = "output_file", default = NULL,
              help = "Optional output path (single mode)"),
  make_option(c("-C", .alias("--comparison-file")), type = "character", dest = "comparison_file",
              help = "TSV with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff (batch mode)"),
  make_option(c("-O", .alias("--output-dir")), type = "character", dest = "output_dir", default = NULL,
              help = "Root directory containing per-comparison folders (batch mode) [inherits function default if omitted]"),
  make_option(.alias("--overwrite"), action = "store_true", dest = "overwrite", default = NA,
              help = "Recompute and overwrite existing <comp>_dnds_annot.tsv"),
  make_option(.alias("--custom"), type = "character", dest = "custom", default = NULL,
              help = "Comma patterns like 'REX{5},TOM{8}' → q_prefix/s_prefix columns")
)

run_append_annotations <- function(o) {
  .dump_opts("append_annotations", o)
  options(dndsR.threads = as.integer(o$threads))
  threads <- as.integer(o$threads); if (is.na(threads) || threads < 1L) threads <- 1L

  args <- list(
    dnds_file       = o$dnds_file,
    query_gff       = o$query_gff,
    subject_gff     = o$subject_gff
  )
  args <- .add_if(args, "output_file",    o$output_file)
  args <- .add_if(args, "comparison_file",o$comparison_file)
  args <- .add_if(args, "output_dir",     o$output_dir)
  args <- .add_if(args, "custom",         o$custom)
  args <- .add_if(args, "threads",        threads)
  args <- .add_if(args, "overwrite",      if (!is.na(o$overwrite) && isTRUE(o$overwrite)) TRUE else NULL)

  res <- do.call(dndsR::append_annotations, args)

  if (is.character(res) && length(res)) {
    cli::cli_alert_success("append_annotations completed. Wrote {length(res)} file{?s}.")
    for (p in res) cli::cli_bullets(c(v = paste0("<", p, ">")))
  } else {
    cli::cli_alert_success("append_annotations completed.")
  }
  invisible(res)
}

# =========================
# ipr_enrichment
# =========================
opts_ipr_enrichment <- list(
  make_option(.alias("--dnds-annot-file"), type="character", dest="dnds_annot_file",
              help="Path to <comp>_dnds_annot.tsv (single mode)"),
  make_option(c("-C", .alias("--comparison-file")), type="character", dest="comparison_file",
              help="TSV with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff (batch mode)"),
  make_option(c("-O", .alias("--output-dir")), type="character", dest="output_dir", default=NULL,
              help="Root directory of per-comparison folders [inherits function default if omitted]"),
  make_option(.alias("--sides"), type="character", dest="sides", default=NULL,
              help="Comma-separated among query,subject"),
  make_option(.alias("--pos-threshold"), type="double", dest="pos_threshold", default=NULL,
              help="dNdS > pos_threshold defines positive selection"),
  make_option(.alias("--max-dnds"), type="double", dest="max_dnds", default=NULL,
              help="Drop rows with dNdS >= max_dnds or NA"),
  make_option(.alias("--filter-expr"), type="character", dest="filter_expr", default=NULL,
              help='Logical expression evaluated in data (e.g., "q_seqname == s_seqname")'),
  make_option(c(.alias("--make-plots")) , action="store_true", dest="make_plots", default=NA,
              help="Write top-N bubble plots"),
  make_option(c(.alias("--no-plots"))   , action="store_true", dest="no_plots", default=NA,
              help="Disable plots (overrides --make-plots)"),
  make_option(.alias("--top-n"), type="integer", dest="top_n", default=NULL,
              help="Top-N rows for plot"),
  make_option(c(.alias("--drop-rows-without-term")), action="store_true", dest="drop_rows_without_term", default=NA,
              help="Remove rows with no IPR from both pos and background"),
  make_option(c(.alias("--keep-rows-without-term")), action="store_true", dest="keep_rows_without_term", default=NA,
              help="Override to keep rows without IPR"),
  make_option(.alias("--min-total"), type="integer", dest="min_total", default=NULL,
              help="Minimum total occurrences (pos+nonpos) required"),
  make_option(.alias("--min-pos"), type="integer", dest="min_pos", default=NULL,
              help="Minimum positive occurrences required"),
  make_option(.alias("--fdr-method"), type="character", dest="fdr_method", default=NULL,
              help='One of BH,IHW,qvalue,none'),
  make_option(.alias("--alpha"), type="double", dest="alpha", default=NULL,
              help="FDR alpha for IHW"),
  make_option(.alias("--term-sep"), type="character", dest="term_sep", default=NULL,
              help="Separator used in q_ipr/s_ipr strings"),
  make_option(.alias("--include-types"), type="character", dest="include_types", default=NULL,
              help='Comma-separated InterPro ENTRY_TYPEs to keep in pooled mode (e.g., "Domain,Homologous_superfamily")'),
  make_option(c(.alias("--stratify-by-type")), action="store_true", dest="stratify_by_type", default=NA,
              help="Run separate analyses per ENTRY_TYPE with type-specific backgrounds"),
  make_option(.alias("--types"), type="character", dest="types", default=NULL,
              help='When stratified, comma-separated ENTRY_TYPEs to run (default: infer)'),
  make_option(.alias("--adjust-scope"), type="character", dest="adjust_scope", default=NULL,
              help='When stratified: "global" or "per_type"'),
  make_option(.alias("--entries-source"), type="character", dest="entries_source", default=NULL,
              help='Where to load InterPro entries: auto|local|remote|none'),
  make_option(.alias("--entries-path"), type="character", dest="entries_path", default=NULL,
              help="Path to InterPro entry.list TSV (ENTRY_AC,ENTRY_TYPE,ENTRY_NAME)"),
  make_option(.alias("--entries-url"), type="character", dest="entries_url", default=NULL,
              help="Remote URL for current InterPro entry.list"),
  make_option(.alias("--entries-timeout-s"), type="integer", dest="entries_timeout_s", default=NULL,
              help="Timeout (seconds) for remote fetch"),
  make_option(c(.alias("--keep-unmatched")), action="store_true", dest="keep_unmatched", default=NA,
              help="Keep IPRs not found in entry.list when filtering"),
  make_option(c(.alias("--drop-unmatched")), action="store_true", dest="drop_unmatched", default=NA,
              help="Override to drop unmatched IPRs"),
  make_option(.alias("--exclude-ids"), type="character", dest="exclude_ids", default=NULL,
              help='Comma-separated IPRs to exclude globally (e.g., "IPR000123,IPR012345")'),
  make_option(.alias("--exclude-iprs"), type="character", dest="exclude_iprs", default=NULL,
              help="Deprecated alias of --exclude-ids (comma-separated)"),
  make_option(.alias("--term-trees"), type="character", dest="term_trees", default=NULL,
              help="Path to IPR parent/child edgelist TSV (two columns)"),
  make_option(c(.alias("--exclude-descendants")), action="store_true", dest="exclude_descendants", default=NA,
              help="Expand exclude IDs to descendants using --term-trees"),
  make_option(.alias("--exclude-descendants-depth"), type="double", dest="exclude_descendants_depth", default=NULL,
              help="Depth limit for descendant exclusion"),
  make_option(.alias("--exclude-descendants-limit"), type="integer", dest="exclude_descendants_limit", default=NULL,
              help="Hard cap on number of excluded IPRs")
)

run_ipr_enrichment <- function(o){
  .dump_opts("ipr_enrichment", o)

  # vectors / scalars possibly NULL
  sides_vec         <- .parse_csv(o$sides)
  include_types_vec <- .parse_csv(o$include_types)
  types_vec         <- .parse_csv(o$types)
  exclude_ids_vec   <- .parse_csv(o$exclude_ids)
  exclude_iprs_vec  <- .parse_csv(o$exclude_iprs)

  # tri-state booleans
  make_plots_val     <- .resolve_toggle(o$make_plots, o$no_plots)
  drop_rows_val      <- .resolve_toggle(o$drop_rows_without_term, o$keep_rows_without_term)
  keep_unmatched_val <- .resolve_toggle(o$keep_unmatched, o$drop_unmatched)
  stratify_val       <- if (!is.na(o$stratify_by_type) && isTRUE(o$stratify_by_type)) TRUE else if (isFALSE(o$stratify_by_type)) FALSE e
lse NULL
  exclude_desc_val   <- if (!is.na(o$exclude_descendants) && isTRUE(o$exclude_descendants)) TRUE else NULL

  args <- list(
    dnds_annot_file = o$dnds_annot_file,
    comparison_file = o$comparison_file
  )
  args <- .add_if(args, "output_dir",           o$output_dir)
  args <- .add_if(args, "sides",                sides_vec)
  args <- .add_if(args, "pos_threshold",        o$pos_threshold, as.numeric)
  args <- .add_if(args, "max_dnds",             o$max_dnds,      as.numeric)
  args <- .add_if(args, "filter_expr",          o$filter_expr)
  args <- .add_if(args, "make_plots",           make_plots_val)
  args <- .add_if(args, "top_n",                o$top_n,         as.integer)
  args <- .add_if(args, "drop_rows_without_term", drop_rows_val)
  args <- .add_if(args, "min_total",            o$min_total,     as.integer)
  args <- .add_if(args, "min_pos",              o$min_pos,       as.integer)
  args <- .add_if(args, "fdr_method",           o$fdr_method)
  args <- .add_if(args, "alpha",                o$alpha,         as.numeric)
  args <- .add_if(args, "term_sep",             o$term_sep)
  args <- .add_if(args, "include_types",        include_types_vec)
  args <- .add_if(args, "stratify_by_type",     stratify_val)
  args <- .add_if(args, "types",                types_vec)
  args <- .add_if(args, "adjust_scope",         o$adjust_scope)
  args <- .add_if(args, "entries_source",       o$entries_source)
  args <- .add_if(args, "entries_path",         o$entries_path)
  args <- .add_if(args, "entries_url",          o$entries_url)
  args <- .add_if(args, "entries_timeout_s",    o$entries_timeout_s, as.integer)
  args <- .add_if(args, "keep_unmatched",       keep_unmatched_val)
  args <- .add_if(args, "exclude_ids",          exclude_ids_vec)
  args <- .add_if(args, "exclude_iprs",         exclude_iprs_vec)
  args <- .add_if(args, "term_trees",           o$term_trees)
  args <- .add_if(args, "exclude_descendants",  exclude_desc_val)
  args <- .add_if(args, "exclude_descendants_depth",  o$exclude_descendants_depth)
  args <- .add_if(args, "exclude_descendants_limit",  o$exclude_descendants_limit, as.integer)

  res <- do.call(dndsR::ipr_enrichment, args)

  if (is.character(res) && length(res)) {
    cli::cli_alert_success("ipr_enrichment completed. Wrote {length(res)} file{?s}.")
    for (p in res) cli::cli_bullets(c(v = paste0("<", p, ">")))
  } else {
    cli::cli_alert_success("ipr_enrichment completed.")
  }
  invisible(res)
}

# =========================
# go_enrichment
# =========================
opts_go_enrichment <- list(
  make_option(.alias("--dnds-annot-file"), type="character", dest="dnds_annot_file",
              help="Path to <comp>_dnds_annot.tsv (single mode)"),
  make_option(c("-C", .alias("--comparison-file")), type="character", dest="comparison_file",
              help="TSV with comparison_name,... (batch mode)"),
  make_option(c("-O", .alias("--output-dir")), type="character", dest="output_dir", default=NULL,
              help="Root directory of per-comparison folders [inherits function default if omitted]"),

  make_option(.alias("--sides"), type="character", dest="sides", default=NULL,
              help="Comma-separated among query,subject"),
  make_option(.alias("--ontologies"), type="character", dest="ontologies", default=NULL,
              help="Comma-separated among BP,MF,CC"),
  make_option(.alias("--algorithm"), type="character", dest="algorithm", default=NULL,
              help="topGO algorithm: weight01|elim|classic|weight"),
  make_option(.alias("--statistic"), type="character", dest="statistic", default=NULL,
              help="topGO statistic: fisher|ks"),

  make_option(.alias("--pos-threshold"), type="double", dest="pos_threshold", default=NULL,
              help="dNdS > pos_threshold defines positive selection"),
  make_option(.alias("--max-dnds"), type="double", dest="max_dnds", default=NULL,
              help="Drop rows with dNdS >= max_dnds or NA"),
  make_option(.alias("--filter-expr"), type="character", dest="filter_expr", default=NULL,
              help='Logical expression evaluated in data (e.g., "q_seqname == s_seqname")'),

  make_option(c(.alias("--drop-rows-without-go")), action="store_true", dest="drop_rows_without_go", default=NA,
              help="Universe = genes with ≥1 GO for that side"),
  make_option(c(.alias("--keep-rows-without-go")), action="store_true", dest="keep_rows_without_go", default=NA,
              help="Override to keep rows without GO (disables the drop behavior)"),

  make_option(.alias("--node-min"), type="integer", dest="node_min", default=NULL,
              help="Minimum term size (topGO nodeSize)"),
  make_option(.alias("--node-max"), type="double", dest="node_max", default=NULL,
              help="Maximum term size to keep post-test"),
  make_option(.alias("--p-adjust"), type="character", dest="p_adjust", default=NULL,
              help="BH or none"),

  make_option(c(.alias("--make-plots")), action="store_true", dest="make_plots", default=NA,
              help="Write top-N bubble plots"),
  make_option(c(.alias("--no-plots")), action="store_true", dest="no_plots", default=NA,
              help="Disable plots (overrides --make-plots)"),
  make_option(.alias("--top-n"), type="integer", dest="top_n", default=NULL,
              help="Top-N rows to plot"),

  make_option(.alias("--exclude-gos"), type="character", dest="exclude_gos", default=NULL,
              help='Comma-separated GO IDs to exclude globally (e.g., "GO:0008150,GO:0009987")'),
  make_option(c(.alias("--exclude-descendants")), action="store_true", dest="exclude_descendants", default=NA,
              help="Also exclude descendants of each GO ID"),
  make_option(.alias("--exclude-descendants-depth"), type="double", dest="exclude_descendants_depth", default=NULL,
              help="Depth limit for descendant exclusion"),
  make_option(.alias("--exclude-descendants-limit"), type="integer", dest="exclude_descendants_limit", default=NULL,
              help="Hard cap on total excluded nodes"),

  make_option(c(.alias("--include-definition")), action="store_true", dest="include_definition", default=NA,
              help="Append GO term definition from GO.db"),
  make_option(c(.alias("--no-definition")), action="store_true", dest="no_definition", default=NA,
              help="Disable adding GO term definitions (overrides --include-definition)")
)

run_go_enrichment <- function(o){
  .dump_opts("go_enrichment", o)
  sides_vec <- .parse_csv(o$sides)
  onts_vec  <- .parse_csv(o$ontologies)
  exclude_vec <- .parse_csv(o$exclude_gos)

  drop_rows_val   <- .resolve_toggle(o$drop_rows_without_go, o$keep_rows_without_go)
  make_plots_val  <- .resolve_toggle(o$make_plots, o$no_plots)
  include_def_val <- .resolve_toggle(o$include_definition, o$no_definition)
  exclude_desc_val<- if (!is.na(o$exclude_descendants) && isTRUE(o$exclude_descendants)) TRUE else NULL

  args <- list(
    dnds_annot_file = o$dnds_annot_file,
    comparison_file = o$comparison_file
  )
  args <- .add_if(args, "output_dir",            o$output_dir)
  args <- .add_if(args, "sides",                 sides_vec)
  args <- .add_if(args, "ontologies",            onts_vec)
  args <- .add_if(args, "algorithm",             o$algorithm)
  args <- .add_if(args, "statistic",             o$statistic)
  args <- .add_if(args, "pos_threshold",         o$pos_threshold, as.numeric)
  args <- .add_if(args, "max_dnds",              o$max_dnds,      as.numeric)
  args <- .add_if(args, "filter_expr",           o$filter_expr)
  args <- .add_if(args, "drop_rows_without_go",  drop_rows_val)
  args <- .add_if(args, "node_min",              o$node_min,      as.integer)
  args <- .add_if(args, "node_max",              o$node_max)
  args <- .add_if(args, "p_adjust",              o$p_adjust)
  args <- .add_if(args, "make_plots",            make_plots_val)
  args <- .add_if(args, "top_n",                 o$top_n,         as.integer)
  args <- .add_if(args, "exclude_gos",           exclude_vec)
  args <- .add_if(args, "exclude_descendants",   exclude_desc_val)
  args <- .add_if(args, "exclude_descendants_depth", o$exclude_descendants_depth)
  args <- .add_if(args, "exclude_descendants_limit", o$exclude_descendants_limit, as.integer)
  args <- .add_if(args, "include_definition",    include_def_val)

  res <- do.call(dndsR::go_enrichment, args)

  if (is.character(res) && length(res)) {
    cli::cli_alert_success("go_enrichment completed. Wrote {length(res)} file{?s}.")
    for (p in res) cli::cli_bullets(c(v = paste0("<", p, ">")))
  } else {
    cli::cli_alert_success("go_enrichment completed.")
  }
  invisible(res)
}

# =========================
# term_enrichment (generic)
# =========================
opts_term_enrichment <- list(
  make_option(.alias("--dnds-annot-file"), type="character", dest="dnds_annot_file",
              help="Path to <comp>_dnds_annot.tsv (single mode)"),
  make_option(c("-C", .alias("--comparison-file")), type="character", dest="comparison_file",
              help="TSV with comparison_name,... (batch mode)"),
  make_option(c("-O", .alias("--output-dir")), type="character", dest="output_dir", default=NULL,
              help="Root directory of per-comparison folders [inherits function default if omitted]"),

  make_option(.alias("--terms"), type="character", dest="terms", default=NULL,
              help='Comma list of term types to test (e.g., "kegg,pfam"). Default: auto-detect'),
  make_option(.alias("--exclude-terms"), type="character", dest="exclude_terms", default=NULL,
              help='Comma list of types to exclude entirely'),

  make_option(.alias("--sides"), type="character", dest="sides", default=NULL,
              help="Comma-separated among query,subject"),
  make_option(.alias("--pos-threshold"), type="double", dest="pos_threshold", default=NULL,
              help="dNdS > pos_threshold defines positive selection"),
  make_option(.alias("--max-dnds"), type="double", dest="max_dnds", default=NULL,
              help="Drop rows with dNdS >= max_dnds or NA"),
  make_option(.alias("--filter-expr"), type="character", dest="filter_expr", default=NULL,
              help='Logical expression evaluated in data (e.g., "q_seqname == s_seqname")'),

  make_option(c(.alias("--make-plots")), action="store_true", dest="make_plots", default=NA,
              help="Write top-N bubble plots"),
  make_option(c(.alias("--no-plots")), action="store_true", dest="no_plots", default=NA,
              help="Disable plots (overrides --make-plots)"),
  make_option(.alias("--top-n"), type="integer", dest="top_n", default=NULL,
              help="Top-N rows for plot"),

  make_option(c(.alias("--drop-rows-without-term")), action="store_true", dest="drop_rows_without_term", default=NA,
              help="Annotation-aware universe (drop rows without the tested term type)"),
  make_option(c(.alias("--keep-rows-without-term")), action="store_true", dest="keep_rows_without_term", default=NA,
              help="Override to keep rows without the term type"),

  make_option(.alias("--min-total"), type="integer", dest="min_total", default=NULL,
              help="Minimum total occurrences required"),
  make_option(.alias("--min-pos"), type="integer", dest="min_pos", default=NULL,
              help="Minimum positive occurrences required"),
  make_option(.alias("--fdr-method"), type="character", dest="fdr_method", default=NULL,
              help="BH|IHW|qvalue|none"),
  make_option(.alias("--alpha"), type="double", dest="alpha", default=NULL,
              help="FDR level for IHW weighting"),

  make_option(.alias("--term-sep"), type="character", dest="term_sep", default=NULL,
              help="Default separator for term strings"),
  make_option(.alias("--term-seps"), type="character", dest="term_seps", default=NULL,
              help='Candidate separators for auto-detect per column (comma-separated). Use "<comma>" for ",".'),

  make_option(.alias("--term-blocklist"), type="character", dest="term_blocklist", default=NULL,
              help="Comma list of suffixes to ignore as term families"),

  # Exclusions / trees / metadata
  make_option(.alias("--exclude-ids"), type="character", dest="exclude_ids", default=NULL,
              help='Global comma list of IDs to exclude (applies to all types)'),
  make_option(.alias("--eid"), type="character", action="append", dest="eid", default=NULL,
              help='Per-type exclude IDs as key=value (repeatable). Example: --eid pfam=PF00001,PF00002'),
  make_option(.alias("--term-trees-path"), type="character", dest="term_trees_path", default=NULL,
              help="Single edgelist path (applies to all types)"),
  make_option(.alias("--tt"), type="character", action="append", dest="tt", default=NULL,
              help='Per-type tree path as key=value (repeatable). Example: --tt pfam=/path/pfam_tree.tsv'),
  make_option(c(.alias("--exclude-descendants")), action="store_true", dest="exclude_descendants", default=NA,
              help="Expand exclude IDs using term_trees"),
  make_option(.alias("--exclude-descendants-depth"), type="double", dest="exclude_descendants_depth", default=NULL,
              help="Depth limit for descendant exclusion"),
  make_option(.alias("--exclude-descendants-limit"), type="integer", dest="exclude_descendants_limit", default=NULL,
              help="Hard cap on excluded nodes per type"),

  make_option(.alias("--term-metadata-path"), type="character", dest="term_metadata_path", default=NULL,
              help="Single metadata TSV (applies to all types; needs columns: ID, NAME [optional DEFINITION])"),
  make_option(.alias("--tm"), type="character", action="append", dest="tm", default=NULL,
              help='Per-type metadata path as key=value (repeatable). Example: --tm kegg=/path/kegg_meta.tsv'),
  make_option(c(.alias("--keep-unmatched-ids")) , action="store_true", dest="keep_unmatched_ids", default=NA,
              help="Keep terms lacking metadata"),
  make_option(c(.alias("--drop-unmatched-ids")) , action="store_true", dest="drop_unmatched_ids", default=NA,
              help="Drop terms lacking metadata (overrides keep)")
)

run_term_enrichment <- function(o){
  .dump_opts("term_enrichment", o)
  terms_vec          <- .parse_csv(o$terms)
  exclude_terms_vec  <- .parse_csv(o$exclude_terms)
  sides_vec          <- .parse_csv(o$sides)
  term_seps_vec      <- if (is.null(o$term_seps)) NULL else .parse_seps(o$term_seps)
  term_blocklist_vec <- if (is.null(o$term_blocklist)) NULL else .parse_csv_required(o$term_blocklist)
  exclude_ids_global <- .parse_csv(o$exclude_ids)
  exclude_ids_per    <- .parse_kv_list_of_csv(o$eid)
  trees_per          <- .parse_kv_list(o$tt)
  meta_per           <- .parse_kv_list(o$tm)

  term_trees_val <- if (!is.null(o$term_trees_path)) o$term_trees_path else NULL
  term_meta_val  <- if (!is.null(o$term_metadata_path)) o$term_metadata_path else NULL

  exclude_ids_final <- NULL
  if (!is.null(exclude_ids_per)) {
    exclude_ids_final <- exclude_ids_per
    if (!is.null(exclude_ids_global)) {
      for (k in names(exclude_ids_final)) {
        exclude_ids_final[[k]] <- unique(c(exclude_ids_global, exclude_ids_final[[k]]))
      }
    }
  } else if (!is.null(exclude_ids_global)) {
    exclude_ids_final <- exclude_ids_global
  }

  term_trees_final <- if (!is.null(trees_per)) trees_per else term_trees_val
  term_meta_final  <- if (!is.null(meta_per))  meta_per  else term_meta_val

  make_plots_val         <- .resolve_toggle(o$make_plots, o$no_plots)
  drop_rows_val          <- .resolve_toggle(o$drop_rows_without_term, o$keep_rows_without_term)
  keep_unmatched_ids_val <- .resolve_toggle(o$keep_unmatched_ids, o$drop_unmatched_ids)
  exclude_desc_val       <- if (!is.na(o$exclude_descendants) && isTRUE(o$exclude_descendants)) TRUE else NULL

  args <- list(
    dnds_annot_file = o$dnds_annot_file,
    comparison_file = o$comparison_file
  )
  args <- .add_if(args, "output_dir",              o$output_dir)
  args <- .add_if(args, "terms",                   terms_vec)
  args <- .add_if(args, "exclude_terms",           exclude_terms_vec)
  args <- .add_if(args, "sides",                   sides_vec)
  args <- .add_if(args, "pos_threshold",           o$pos_threshold, as.numeric)
  args <- .add_if(args, "max_dnds",                o$max_dnds,      as.numeric)
  args <- .add_if(args, "filter_expr",             o$filter_expr)
  args <- .add_if(args, "make_plots",              make_plots_val)
  args <- .add_if(args, "top_n",                   o$top_n,         as.integer)
  args <- .add_if(args, "drop_rows_without_term",  drop_rows_val)
  args <- .add_if(args, "min_total",               o$min_total,     as.integer)
  args <- .add_if(args, "min_pos",                 o$min_pos,       as.integer)
  args <- .add_if(args, "fdr_method",              o$fdr_method)
  args <- .add_if(args, "alpha",                   o$alpha,         as.numeric)
  args <- .add_if(args, "term_sep",                o$term_sep)
  args <- .add_if(args, "term_seps",               term_seps_vec)
  args <- .add_if(args, "term_blocklist",          term_blocklist_vec)
  args <- .add_if(args, "exclude_ids",             exclude_ids_final)
  args <- .add_if(args, "term_trees",              term_trees_final)
  args <- .add_if(args, "exclude_descendants",     exclude_desc_val)
  args <- .add_if(args, "exclude_descendants_depth", o$exclude_descendants_depth)
  args <- .add_if(args, "exclude_descendants_limit", o$exclude_descendants_limit, as.integer)
  args <- .add_if(args, "term_metadata",           term_meta_final)
  args <- .add_if(args, "keep_unmatched_ids",      keep_unmatched_ids_val)

  res <- do.call(dndsR::term_enrichment, args)

  if (is.character(res) && length(res)) {
    cli::cli_alert_success("term_enrichment completed. Wrote {length(res)} file{?s}.")
    for (p in res) cli::cli_bullets(c(v = paste0("<", p, ">")))
  } else {
    cli::cli_alert_success("term_enrichment completed.")
  }
  invisible(res)
}

# =========================
# dnds_ideogram
# =========================
opts_dnds_ideogram <- list(
  make_option(.alias("--dnds-annot"), type="character", dest="dnds_annot",
              help="Single prebuilt table with q_/s_ coords (single mode; function mainly intended for batch)"),
  make_option(c("-C", .alias("--comparison-file")), type="character", dest="comparison_file",
              help="TSV: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff (batch mode)"),
  make_option(c("-O", .alias("--output-dir")), type="character", dest="output_dir", default=NULL,
              help="Root directory of per-comparison folders [inherits function default if omitted]"),

  make_option(.alias("--sides"), type="character", dest="sides", default=NULL,
              help="Comma-separated among query,subject"),
  make_option(.alias("--window-size"), type="integer", dest="window_size", default=NULL,
              help="Window size (bp) for binning genes"),
  make_option(.alias("--max-dnds"), type="double", dest="max_dnds", default=NULL,
              help="Drop rows with dNdS >= max_dnds or NA"),
  make_option(.alias("--filter-expr"), type="character", dest="filter_expr", default=NULL,
              help='Logical filter evaluated on the augmented table (e.g., "q_chr == s_chr")'),

  make_option(c(.alias("--make-png")), action="store_true", dest="make_png", default=NA,
              help="Also write PNG via convertSVG"),
  make_option(c(.alias("--no-png")), action="store_true", dest="no_png", default=NA,
              help="Disable PNG export (overrides --make-png)"),
  make_option(.alias("--overwrite"), action="store_true", dest="overwrite", default=NA,
              help="Overwrite existing outputs"),
  make_option(.alias("--keep-intermediate"), action="store_true", dest="keep_intermediate", default=NA,
              help="Keep <comp>_ideogram_input.tsv"),

  # Chromosome label normalization (opt-in; raw by default)
  make_option(.alias("--chr-strip-leading-chr0"), action="store_true", dest="chr_strip_leading_chr0", default=NA,
              help='Strip leading "chr" and optional "0" (e.g., chr01 -> 1)'),
  make_option(.alias("--chr-strip-leading"), type="character", dest="chr_strip_leading", default=NULL,
              help='Regex to strip from start (e.g., "Cf0")'),
  make_option(.alias("--chr-strip-trailing"), type="character", dest="chr_strip_trailing", default=NULL,
              help='Regex to strip from end (e.g., "_random")'),
  make_option(c(.alias("--chr-case-insensitive")), action="store_true", dest="chr_case_insensitive", default=NA,
              help="Apply normalization regex case-insensitively"),
  make_option(c(.alias("--chr-case-sensitive")), action="store_true", dest="chr_case_sensitive", default=NA,
              help="Override to make normalization case-sensitive"),

  # GFF indexing
  make_option(.alias("--restrict-gff-to-gene"), action="store_true", dest="restrict_gff_to_gene", default=NA,
              help='Use only rows with type=="gene" when indexing coordinates')
)

run_dnds_ideogram <- function(o){
  .dump_opts("dnds_ideogram", o)
  sides_vec <- .parse_csv(o$sides)

  make_png_val          <- .resolve_toggle(o$make_png, o$no_png)
  chr_case_ins_val      <- .resolve_toggle(o$chr_case_insensitive, o$chr_case_sensitive)
  overwrite_val         <- if (!is.na(o$overwrite) && isTRUE(o$overwrite)) TRUE else NULL
  keep_intermediate_val <- if (!is.na(o$keep_intermediate) && isTRUE(o$keep_intermediate)) TRUE else NULL
  chr_strip_chr0_val    <- if (!is.na(o$chr_strip_leading_chr0) && isTRUE(o$chr_strip_leading_chr0)) TRUE else NULL
  restrict_gff_val      <- if (!is.na(o$restrict_gff_to_gene) && isTRUE(o$restrict_gff_to_gene)) TRUE else NULL

  args <- list(
    dnds_annot      = o$dnds_annot,
    comparison_file = o$comparison_file
  )
  args <- .add_if(args, "output_dir",             o$output_dir)
  args <- .add_if(args, "sides",                  sides_vec)
  args <- .add_if(args, "window_size",            o$window_size, as.integer)
  args <- .add_if(args, "max_dnds",               o$max_dnds,    as.numeric)
  args <- .add_if(args, "filter_expr",            o$filter_expr)
  args <- .add_if(args, "make_png",               make_png_val)
  args <- .add_if(args, "overwrite",              overwrite_val)
  args <- .add_if(args, "keep_intermediate",      keep_intermediate_val)
  args <- .add_if(args, "chr_strip_leading_chr0", chr_strip_chr0_val)
  args <- .add_if(args, "chr_strip_leading",      o$chr_strip_leading)
  args <- .add_if(args, "chr_strip_trailing",     o$chr_strip_trailing)
  args <- .add_if(args, "chr_case_insensitive",   chr_case_ins_val)
  args <- .add_if(args, "restrict_gff_to_gene",   restrict_gff_val)
  args <- .add_if(args, "verbose",                if (isTRUE(o$verbose)) TRUE else NULL)

  res <- do.call(dndsR::dnds_ideogram, args)

  if (is.character(res) && length(res)) {
    cli::cli_alert_success("dnds_ideogram completed. Wrote {length(res)} file{?s}.")
    for (p in res) cli::cli_bullets(c(v = paste0("<", p, ">")))
  } else {
    cli::cli_alert_success("dnds_ideogram completed.")
  }
  invisible(res)
}

# =========================
# Subcommands registry
# =========================
subcommands <- list(
  split_comparisons = list(
    opts = c(base_opts, opts_split),      fun = run_split,
    help = "Split genomes by label (subgenome or haplotype) and emit a new comparison file"
  ),
  extract_cds    = list(opts = c(base_opts, opts_extract_cds), fun = run_extract_cds,
                        help = "Extract CDS (and optionally proteins) from genome+GFF"),
  calculate_dnds = list(opts = c(base_opts, opts_calculate_dnds), fun = run_calculate_dnds,
                        help = "Run dN/dS with orthologr from pre-extracted CDS"),
  append_annotations = list(opts = c(base_opts, opts_append_annotations), fun = run_append_annotations,
                            help = "Append query/subject GFF annotations to dN/dS results"),
  ipr_enrichment = list(opts = c(base_opts, opts_ipr_enrichment), fun = run_ipr_enrichment,
                        help = "InterPro (IPR) enrichment (Fisher) for q/s IPR terms in dNdS results"),
  go_enrichment = list(opts = c(base_opts, opts_go_enrichment), fun = run_go_enrichment,
                       help = "GO enrichment with topGO for query/subject columns"),
  term_enrichment = list(opts = c(base_opts, opts_term_enrichment), fun = run_term_enrichment,
                         help = "Generic term enrichment (Fisher) for q_*/s_* term columns"),
  dnds_ideogram = list(
    opts = c(base_opts, opts_dnds_ideogram),
    fun  = run_dnds_ideogram,
    help = "Genome-wide dN/dS ideograms (RIdeogram) for query/subject"
  )
)

usage <- function(){
  cat("\nUsage: dndsR <subcommand> [options]\n\nSubcommands:\n")
  for (nm in names(subcommands))
    cat(sprintf("  %-14s %s\n", nm, subcommands[[nm]]$help))
  cat("\nUse: dndsR <subcommand> --help   for subcommand options\n\n")
}

# =========================
# Entry point
# =========================
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) == 0 || argv[1] %in% c("-h","--help")) { usage(); quit(status=0) }
cmd <- argv[1]; args <- argv[-1]
if (!cmd %in% names(subcommands)) { cli::cli_alert_danger(paste("Unknown:", cmd)); usage(); quit(status=1) }
spec <- subcommands[[cmd]]
parser <- optparse::OptionParser(usage=paste0("%prog ", cmd, " [options]"), option_list = spec$opts)
opts <- optparse::parse_args(parser, args=args, positional_arguments = FALSE)

# ---- Verbosity bridge for all subcommands ----
if (isTRUE(opts$verbose) || identical(tolower(Sys.getenv("DNDSR_VERBOSE", "0")), "1")) {
  options(dndsR.verbose = TRUE, dndsR.cli.verbose = TRUE)
  if (!isTRUE(opts$verbose) && identical(tolower(Sys.getenv("DNDSR_VERBOSE", "0")), "1")) {
    opts$verbose <- TRUE
  }
}

# ---- Threads bridge for all subcommands ----
threads <- if (!is.null(opts$threads)) as.integer(opts$threads) else 1L
if (is.na(threads) || threads < 1L) threads <- 1L
options(dndsR.threads = threads)

tryCatch({
  spec$fun(opts)
  if (isTRUE(opts$verbose)) cli::cli_alert_success("Done.")
}, error = function(e){
  cli::cli_alert_danger(conditionMessage(e)); quit(status=1)
})
