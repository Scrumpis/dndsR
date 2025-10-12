#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(cli)
  library(dndsR)
})

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
  if (!length(v)) return(c(";", "|", ","))  # hard default
  v <- trimws(v[nzchar(v)])
  v[v %in% c("comma", "<comma>", "\\,")] <- ","
  v
}

# =========================
# split_comparisons
# =========================
opts_split <- list(
  make_option(c("-C","--comparison-file"), type="character", dest="comparison_file",
              help="TSV with columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff"),
  make_option(c("-m","--mode"), type="character", default="subgenome",
              help="subgenome|haplotype|custom [default: %default]"),
  make_option(c("-r","--custom-regex"), type="character", default=NULL,
              help="Regex with exactly one capture group when --mode=custom"),
  make_option(c("--case-insensitive"), action="store_true", default=TRUE,
              help="Case-insensitive label matching [default: %default]"),
  make_option(c("--case-sensitive"),   action="store_true", default=FALSE,
              help="Override to make matching case-sensitive")
)
run_split <- function(o){
  .dump_opts("split_comparisons", o)
  stopifnot(!is.null(o$comparison_file))
  case_insensitive <- if (isTRUE(o$`case-sensitive`)) FALSE else isTRUE(o$`case-insensitive`)
  p <- dndsR::split_comparisons(
    comparison_file = o$comparison_file,
    mode            = o$mode,
    custom_regex    = o$`custom-regex`,
    case_insensitive= case_insensitive
  )
  cli_alert_success("Wrote: {p}")
}

# =========================
# extract_cds
# =========================
opts_extract_cds <- list(
  make_option(c("--comparison-name"), type = "character", dest = "comparison_name",
              help = "Unique identifier for a comparison (e.g., 'CheAl_v_CheFo')"),
  make_option(c("--query-fasta"), type = "character", dest = "query_fasta",
              help = "Path to query genome FASTA"),
  make_option(c("--subject-fasta"), type = "character", dest = "subject_fasta",
              help = "Path to subject genome FASTA (optional in single-genome mode)"),
  make_option(c("--query-gff"), type = "character", dest = "query_gff",
              help = "Path to query GFF3"),
  make_option(c("--subject-gff"), type = "character", dest = "subject_gff",
              help = "Path to subject GFF3 (optional in single-genome mode)"),
  make_option(c("-O","--output-dir"), type = "character", default = ".",
              help = "Output directory; one subdir per comparison [default: %default]"),
  make_option(c("--overwrite"), action = "store_true", default = FALSE,
              help = "Overwrite existing outputs"),
  make_option(c("-C","--comparison-file"), type = "character", dest = "comparison_file",
              help = "TSV with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff"),
  make_option(c("--group-by"), type = "character", default = "gene",
              help = "CDS grouping: gene|tx [default: %default]"),
  make_option(c("--export-proteins"), action = "store_true", default = FALSE,
              help = "Also write translated proteins (*.faa)"),
  make_option(c("--genetic-code"), type = "integer", default = 1L,
              help = "NCBI genetic code for translation [default: %default]"),
  make_option(c("--keep-internal-stops"), action = "store_true", default = FALSE,
              help = "Keep sequences with internal stops (otherwise drop them)")
)

run_extract_cds <- function(o) {
  .dump_opts("extract_cds", o)
  res <- dndsR::extract_cds(
    comparison_name      = o$comparison_name,
    query_fasta          = o$query_fasta,
    subject_fasta        = o$subject_fasta,
    query_gff            = o$query_gff,
    subject_gff          = o$subject_gff,
    output_dir           = o$`output-dir`,
    overwrite            = isTRUE(o$overwrite),
    verbose              = isTRUE(o$verbose),
    comparison_file      = o$comparison_file,
    group_by             = o$`group-by`,
    export_proteins      = isTRUE(o$`export-proteins`),
    genetic_code         = as.integer(o$`genetic-code`),
    keep_internal_stops  = isTRUE(o$`keep-internal-stops`)
  )
  cli::cli_alert_success("extract_cds completed.")
  invisible(res)
}

# =========================
# calculate_dnds
# =========================
opts_calculate_dnds <- list(
  make_option(c("-C","--comparison-file"), type="character", dest="comparison_file",
              help="TSV (or data.frame via R) with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff"),
  make_option(c("--comparison-name"), type="character", dest="comparison_name",
              help="Unique comparison ID (e.g., 'CheAl_v_CheFo')"),
  make_option(c("--subject-fasta"), type="character", dest="subject_fasta",
              help="Subject genome FASTA (single mode)"),
  make_option(c("--query-fasta"), type="character", dest="query_fasta",
              help="Query genome FASTA (single mode)"),
  make_option(c("--subject-gff"), type="character", dest="subject_gff",
              help="Subject GFF3 (single mode)"),
  make_option(c("--query-gff"), type="character", dest="query_gff",
              help="Query GFF3 (single mode)"),
  make_option(c("-O","--output-dir"), type="character", default=getwd(),
              help="Output directory [default: %default]"),
  make_option(c("--comp-cores"), type="integer", default=NULL,
              help="Cores passed to orthologr::dNdS via comp_cores [default: --threads]"),
  make_option(c("--aligner"), type="character", default="diamond",
              help="orthologr aligner [default: %default]"),
  make_option(c("--sensitivity-mode"), type="character", default="fast",
              help="orthologr sensitivity_mode [default: %default]"),
  make_option(c("--dnds-method"), type="character", default="Comeron",
              help="orthologr dnds_est.method [default: %default]"),
  make_option(c("--oa"), type = "character", action = "append", default = NULL,
              help = paste0(
                "Pass through to orthologr::dNdS as key=value (repeatable). ",
                "Examples: --oa eval=1E-6 --oa ortho_detection=RBH --oa clean_folders=TRUE ",
                "--oa include_only_chr=1,2,3"
              ))
)

run_calculate_dnds <- function(o){
  .dump_opts("calculate_dnds", o)
  extra <- .parse_oa(o$oa)
  cores <- if (!is.null(o$`comp-cores`)) as.integer(o$`comp-cores`) else as.integer(o$threads)
  base <- list(
    comparison_file = o$comparison_file,
    comparison_name = o$comparison_name,
    subject_fasta   = o$subject_fasta,
    query_fasta     = o$query_fasta,
    subject_gff     = o$subject_gff,
    query_gff       = o$query_gff,
    output_dir      = o$`output-dir`,
    comp_cores      = cores,
    aligner         = o$aligner,
    sensitivity_mode= o$`sensitivity-mode`,
    dnds_method     = o$`dnds-method`
  )
  paths <- do.call(dndsR::calculate_dnds, c(base, extra))

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
  make_option(c("--dnds-file"), type = "character", dest = "dnds_file",
              help = "Path to a dN/dS TSV (single mode)"),
  make_option(c("--query-gff"), type = "character", dest = "query_gff",
              help = "Path to query GFF3 (single mode)"),
  make_option(c("--subject-gff"), type = "character", dest = "subject_gff",
              help = "Path to subject GFF3 (single mode)"),
  make_option(c("-o","--output-file"), type = "character", dest = "output_file", default = NULL,
              help = "Optional output path (single mode)"),
  make_option(c("-C","--comparison-file"), type = "character", dest = "comparison_file",
              help = "TSV with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff (batch mode)"),
  make_option(c("-O","--output-dir"), type = "character", dest = "output_dir", default = getwd(),
              help = "Root directory containing per-comparison folders (batch mode) [default: %default]"),
  make_option(c("--overwrite"), action = "store_true", default = FALSE,
              help = "Recompute and overwrite existing <comp>_dnds_annot.tsv [default: %default]"),
  make_option(c("--custom"), type = "character", default = NULL,
              help = "Comma patterns like 'REX{5},TOM{8}' → q_prefix/s_prefix columns")
)

run_append_annotations <- function(o) {
  .dump_opts("append_annotations", o)

  # advertise globally for helpers
  options(dndsR.threads = as.integer(o$threads))

  threads <- as.integer(o$threads)
  if (is.na(threads) || threads < 1L) threads <- 1L

  res <- dndsR::append_annotations(
    dnds_file       = o$dnds_file,
    query_gff       = o$query_gff,
    subject_gff     = o$subject_gff,
    output_file     = o$output_file,
    comparison_file = o$comparison_file,
    output_dir      = o$output_dir,
    custom          = o$custom,
    threads         = threads,          # use base_opts' -t/--threads
    overwrite       = isTRUE(o$overwrite)
  )

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
  make_option(c("--dnds-annot-file"), type="character", dest="dnds_annot_file",
              help="Path to <comp>_dnds_annot.tsv (single mode)"),
  make_option(c("-C","--comparison-file"), type="character", dest="comparison_file",
              help="TSV with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff (batch mode)"),
  make_option(c("-O","--output-dir"), type="character", dest="output_dir", default=getwd(),
              help="Root directory of per-comparison folders [default: %default]"),
  make_option(c("--sides"), type="character", default="query,subject",
              help="Comma-separated among query,subject [default: %default]"),
  make_option(c("--pos-threshold"), type="double", default=1,
              help="dNdS > pos_threshold defines positive selection [default: %default]"),
  make_option(c("--max-dnds"), type="double", default=10,
              help="Drop rows with dNdS >= max_dnds or NA [default: %default]"),
  make_option(c("--filter-expr"), type="character", default=NULL,
              help='Logical expression evaluated in data (e.g., "q_seqname == s_seqname")'),
  make_option(c("--make-plots"), action="store_true", default=TRUE,
              help="Write top-N bubble plots [default: %default]"),
  make_option(c("--no-plots"), action="store_true", default=FALSE,
              help="Disable plots (overrides --make-plots)"),
  make_option(c("--top-n"), type="integer", default=20,
              help="Top-N rows for plot [default: %default]"),
  make_option(c("--drop-rows-without-term"), action="store_true", default=TRUE,
              help="Remove rows with no IPR from both pos and background [default: %default]"),
  make_option(c("--keep-rows-without-term"), action="store_true", default=FALSE,
              help="Override to keep rows without IPR"),
  make_option(c("--min-total"), type="integer", default=0,
              help="Minimum total occurrences (pos+nonpos) required [default: %default]"),
  make_option(c("--min-pos"), type="integer", default=0,
              help="Minimum positive occurrences required [default: %default]"),
  make_option(c("--fdr-method"), type="character", default="BH",
              help='One of BH,IHW,qvalue,none [default: %default]'),
  make_option(c("--alpha"), type="double", default=0.05,
              help="FDR alpha for IHW [default: %default]"),
  make_option(c("--term-sep"), type="character", default=";",
              help="Separator used in q_ipr/s_ipr strings [default: %default]"),
  make_option(c("--include-types"), type="character", default=NULL,
              help='Comma-separated InterPro ENTRY_TYPEs to keep in pooled mode (e.g., "Domain,Homologous_superfamily")'),
  make_option(c("--stratify-by-type"), action="store_true", default=FALSE,
              help="Run separate analyses per ENTRY_TYPE with type-specific backgrounds"),
  make_option(c("--types"), type="character", default=NULL,
              help='When stratified, comma-separated ENTRY_TYPEs to run (default: infer)'),
  make_option(c("--adjust-scope"), type="character", default="global",
              help='When stratified: "global" or "per_type" [default: %default]'),
  make_option(c("--entries-source"), type="character", default="auto",
              help='Where to load InterPro entries: auto|local|remote|none [default: %default]'),
  make_option(c("--entries-path"), type="character", default=NULL,
              help="Path to InterPro entry.list TSV (ENTRY_AC,ENTRY_TYPE,ENTRY_NAME)"),
  make_option(c("--entries-url"), type="character",
              default="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list",
              help="Remote URL for current InterPro entry.list [default: %default]"),
  make_option(c("--entries-timeout-s"), type="integer", default=20,
              help="Timeout (seconds) for remote fetch [default: %default]"),
  make_option(c("--keep-unmatched"), action="store_true", default=TRUE,
              help="Keep IPRs not found in entry.list when filtering [default: %default]"),
  make_option(c("--drop-unmatched"), action="store_true", default=FALSE,
              help="Override to drop unmatched IPRs"),
  make_option(c("--exclude-ids"), type="character", default=NULL,
              help='Comma-separated IPRs to exclude globally (e.g., "IPR000123,IPR012345")'),
  make_option(c("--exclude-iprs"), type="character", default=NULL,
              help="Deprecated alias of --exclude-ids (comma-separated)"),
  make_option(c("--term-trees"), type="character", default=NULL,
              help="Path to IPR parent/child edgelist TSV (two columns)"),
  make_option(c("--exclude-descendants"), action="store_true", default=FALSE,
              help="Expand exclude IDs to descendants using --term-trees"),
  make_option(c("--exclude-descendants-depth"), type="double", default=Inf,
              help="Depth limit for descendant exclusion [default: Inf]"),
  make_option(c("--exclude-descendants-limit"), type="integer", default=5000,
              help="Hard cap on number of excluded IPRs [default: %default]")
)

run_ipr_enrichment <- function(o){
  .dump_opts("ipr_enrichment", o)
  sides_vec <- .parse_csv_required(o$sides)
  include_types_vec <- .parse_csv(o$`include-types`)
  types_vec <- .parse_csv(o$types)
  exclude_ids_vec <- .parse_csv(o$`exclude-ids`)
  exclude_iprs_vec <- .parse_csv(o$`exclude-iprs`)

  make_plots <- if (isTRUE(o$`no-plots`)) FALSE else isTRUE(o$`make-plots`)
  drop_rows_without_term <- if (isTRUE(o$`keep-rows-without-term`)) FALSE else isTRUE(o$`drop-rows-without-term`)
  keep_unmatched <- if (isTRUE(o$`drop-unmatched`)) FALSE else isTRUE(o$`keep-unmatched`)

  res <- dndsR::ipr_enrichment(
    dnds_annot_file = o$dnds_annot_file,
    comparison_file = o$comparison_file,
    output_dir      = o$output_dir,
    sides           = sides_vec,
    pos_threshold   = as.numeric(o$`pos-threshold`),
    max_dnds        = as.numeric(o$`max-dnds`),
    filter_expr     = o$`filter-expr`,
    make_plots      = make_plots,
    top_n           = as.integer(o$`top-n`),
    drop_rows_without_term = drop_rows_without_term,
    min_total       = as.integer(o$`min-total`),
    min_pos         = as.integer(o$`min-pos`),
    fdr_method      = o$`fdr-method`,
    alpha           = as.numeric(o$alpha),
    term_sep        = o$`term-sep`,
    include_types   = include_types_vec,
    stratify_by_type= isTRUE(o$`stratify-by-type`),
    types           = types_vec,
    adjust_scope    = o$`adjust-scope`,
    entries_source  = o$`entries-source`,
    entries_path    = o$`entries-path`,
    entries_url     = o$`entries-url`,
    entries_timeout_s = as.integer(o$`entries-timeout-s`),
    keep_unmatched  = keep_unmatched,
    exclude_ids     = exclude_ids_vec,
    exclude_iprs    = exclude_iprs_vec,
    term_trees      = o$`term-trees`,
    exclude_descendants = isTRUE(o$`exclude-descendants`),
    exclude_descendants_depth  = o$`exclude-descendants-depth`,
    exclude_descendants_limit  = as.integer(o$`exclude-descendants-limit`)
  )

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
  make_option(c("--dnds-annot-file"), type="character", dest="dnds_annot_file",
              help="Path to <comp>_dnds_annot.tsv (single mode)"),
  make_option(c("-C","--comparison-file"), type="character", dest="comparison_file",
              help="TSV with comparison_name,... (batch mode)"),
  make_option(c("-O","--output-dir"), type="character", dest="output_dir", default=getwd(),
              help="Root directory of per-comparison folders [default: %default]"),

  make_option(c("--sides"), type="character", default="query,subject",
              help="Comma-separated among query,subject [default: %default]"),
  make_option(c("--ontologies"), type="character", default="BP,MF,CC",
              help="Comma-separated among BP,MF,CC [default: %default]"),
  make_option(c("--algorithm"), type="character", default="weight01",
              help="topGO algorithm: weight01|elim|classic|weight [default: %default]"),
  make_option(c("--statistic"), type="character", default="fisher",
              help="topGO statistic: fisher|ks [default: %default]"),

  make_option(c("--pos-threshold"), type="double", default=1,
              help="dNdS > pos_threshold defines positive selection [default: %default]"),
  make_option(c("--max-dnds"), type="double", default=10,
              help="Drop rows with dNdS >= max_dnds or NA [default: %default]"),
  make_option(c("--filter-expr"), type="character", default=NULL,
              help='Logical expression evaluated in data (e.g., "q_seqname == s_seqname")'),

  make_option(c("--drop-rows-without-go"), action="store_true", default=TRUE,
              help="Universe = genes with ≥1 GO for that side [default: %default]"),
  make_option(c("--keep-rows-without-go"), action="store_true", default=FALSE,
              help="Override to keep rows without GO (disables the drop behavior)"),

  make_option(c("--node-min"), type="integer", default=10,
              help="Minimum term size (topGO nodeSize) [default: %default]"),
  make_option(c("--node-max"), type="double", default=Inf,
              help="Maximum term size to keep post-test [default: Inf]"),
  make_option(c("--p-adjust"), type="character", default="BH",
              help="BH or none [default: %default]"),

  make_option(c("--make-plots"), action="store_true", default=TRUE,
              help="Write top-N bubble plots [default: %default]"),
  make_option(c("--no-plots"), action="store_true", default=FALSE,
              help="Disable plots (overrides --make-plots)"),
  make_option(c("--top-n"), type="integer", default=20,
              help="Top-N rows to plot [default: %default]"),

  make_option(c("--exclude-gos"), type="character", default=NULL,
              help='Comma-separated GO IDs to exclude globally (e.g., "GO:0008150,GO:0009987")'),
  make_option(c("--exclude-descendants"), action="store_true", default=FALSE,
              help="Also exclude descendants of each GO ID"),
  make_option(c("--exclude-descendants-depth"), type="double", default=Inf,
              help="Depth limit for descendant exclusion [default: Inf]"),
  make_option(c("--exclude-descendants-limit"), type="integer", default=5000,
              help="Hard cap on total excluded nodes [default: %default]"),

  make_option(c("--include-definition"), action="store_true", default=TRUE,
              help="Append GO term definition from GO.db [default: %default]"),
  make_option(c("--no-definition"), action="store_true", default=FALSE,
              help="Disable adding GO term definitions (overrides --include-definition)")
)

run_go_enrichment <- function(o){
  .dump_opts("go_enrichment", o)
  sides_vec <- .parse_csv_required(o$sides)
  onts_vec  <- .parse_csv_required(o$ontologies)
  exclude_vec <- .parse_csv(o$`exclude-gos`)

  drop_rows_without_go <- if (isTRUE(o$`keep-rows-without-go`)) FALSE else isTRUE(o$`drop-rows-without-go`)
  make_plots <- if (isTRUE(o$`no-plots`)) FALSE else isTRUE(o$`make-plots`)
  include_definition <- if (isTRUE(o$`no-definition`)) FALSE else isTRUE(o$`include-definition`)

  res <- dndsR::go_enrichment(
    dnds_annot_file = o$dnds_annot_file,
    comparison_file = o$comparison_file,
    output_dir      = o$output_dir,
    sides           = sides_vec,
    ontologies      = onts_vec,
    algorithm       = o$algorithm,
    statistic       = o$statistic,
    pos_threshold   = as.numeric(o$`pos-threshold`),
    max_dnds        = as.numeric(o$`max-dnds`),
    filter_expr     = o$`filter-expr`,
    drop_rows_without_go = drop_rows_without_go,
    node_min        = as.integer(o$`node-min`),
    node_max        = o$`node-max`,
    p_adjust        = o$`p-adjust`,
    make_plots      = make_plots,
    top_n           = as.integer(o$`top-n`),
    exclude_gos     = exclude_vec,
    exclude_descendants       = isTRUE(o$`exclude-descendants`),
    exclude_descendants_depth = o$`exclude-descendants-depth`,
    exclude_descendants_limit = as.integer(o$`exclude-descendants-limit`),
    include_definition = include_definition
  )

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
  make_option(c("--dnds-annot-file"), type="character", dest="dnds_annot_file",
              help="Path to <comp>_dnds_annot.tsv (single mode)"),
  make_option(c("-C","--comparison-file"), type="character", dest="comparison_file",
              help="TSV with comparison_name,... (batch mode)"),
  make_option(c("-O","--output-dir"), type="character", dest="output_dir", default=getwd(),
              help="Root directory of per-comparison folders [default: %default]"),

  make_option(c("--terms"), type="character", default=NULL,
              help='Comma list of term types to test (e.g., "kegg,pfam"). Default: auto-detect'),
  make_option(c("--exclude-terms"), type="character", default="ipr,go",
              help='Comma list of types to exclude entirely [default: %default]'),

  make_option(c("--sides"), type="character", default="query,subject",
              help="Comma-separated among query,subject [default: %default]"),
  make_option(c("--pos-threshold"), type="double", default=1,
              help="dNdS > pos_threshold defines positive selection [default: %default]"),
  make_option(c("--max-dnds"), type="double", default=10,
              help="Drop rows with dNdS >= max_dnds or NA [default: %default]"),
  make_option(c("--filter-expr"), type="character", default=NULL,
              help='Logical expression evaluated in data (e.g., "q_seqname == s_seqname")'),

  make_option(c("--make-plots"), action="store_true", default=TRUE,
              help="Write top-N bubble plots [default: %default]"),
  make_option(c("--no-plots"), action="store_true", default=FALSE,
              help="Disable plots (overrides --make-plots)"),
  make_option(c("--top-n"), type="integer", default=20,
              help="Top-N rows for plot [default: %default]"),

  make_option(c("--drop-rows-without-term"), action="store_true", default=TRUE,
              help="Annotation-aware universe (drop rows without the tested term type) [default: %default]"),
  make_option(c("--keep-rows-without-term"), action="store_true", default=FALSE,
              help="Override to keep rows without the term type"),

  make_option(c("--min-total"), type="integer", default=0,
              help="Minimum total occurrences required [default: %default]"),
  make_option(c("--min-pos"), type="integer", default=0,
              help="Minimum positive occurrences required [default: %default]"),
  make_option(c("--fdr-method"), type="character", default="BH",
              help="BH|IHW|qvalue|none [default: %default]"),
  make_option(c("--alpha"), type="double", default=0.05,
              help="FDR level for IHW weighting [default: %default]"),

  make_option(c("--term-sep"), type="character", default=";",
              help="Default separator for term strings [default: %default]"),
  make_option(c("--term-seps"), type="character", default=";,|,<comma>",
              help='Candidate separators for auto-detect per column (comma-separated). Use "<comma>" for ",". Default resolves to c(";", "|", ",").'),

  make_option(c("--term-blocklist"), type="character",
              default="attributes,attribute,attr,notes,note,description,product,name,id,gene_id,transcript_id,parent,dbxref,source,target,type,seqname,seqid,star
t,end,strand,phase,biotype,class,len",
              help="Comma list of suffixes to ignore as term families [default: a built-in list]"),

  # Exclusions / trees / metadata
  make_option(c("--exclude-ids"), type="character", default=NULL,
              help='Global comma list of IDs to exclude (applies to all types)'),
  make_option(c("--eid"), type="character", action="append", default=NULL,
              help='Per-type exclude IDs as key=value (repeatable). Example: --eid pfam=PF00001,PF00002'),
  make_option(c("--term-trees-path"), type="character", default=NULL,
              help="Single edgelist path (applies to all types)"),
  make_option(c("--tt"), type="character", action="append", default=NULL,
              help='Per-type tree path as key=value (repeatable). Example: --tt pfam=/path/pfam_tree.tsv'),
  make_option(c("--exclude-descendants"), action="store_true", default=FALSE,
              help="Expand exclude IDs using term_trees"),
  make_option(c("--exclude-descendants-depth"), type="double", default=Inf,
              help="Depth limit for descendant exclusion [default: Inf]"),
  make_option(c("--exclude-descendants-limit"), type="integer", default=5000,
              help="Hard cap on excluded nodes per type [default: %default]"),

  make_option(c("--term-metadata-path"), type="character", default=NULL,
              help="Single metadata TSV (applies to all types; needs columns: ID, NAME [optional DEFINITION])"),
  make_option(c("--tm"), type="character", action="append", default=NULL,
              help='Per-type metadata path as key=value (repeatable). Example: --tm kegg=/path/kegg_meta.tsv'),
  make_option(c("--keep-unmatched-ids"), action="store_true", default=TRUE,
              help="Keep terms lacking metadata [default: %default]"),
  make_option(c("--drop-unmatched-ids"), action="store_true", default=FALSE,
              help="Drop terms lacking metadata (overrides keep)")
)

run_term_enrichment <- function(o){
  .dump_opts("term_enrichment", o)
  terms_vec          <- .parse_csv(o$terms)
  exclude_terms_vec  <- .parse_csv(o$`exclude-terms`)
  sides_vec          <- .parse_csv_required(o$sides)
  term_seps_vec      <- .parse_seps(o$`term-seps`)
  term_blocklist_vec <- .parse_csv_required(o$`term-blocklist`)
  exclude_ids_global <- .parse_csv(o$`exclude-ids`)
  exclude_ids_per    <- .parse_kv_list_of_csv(o$eid)
  trees_per          <- .parse_kv_list(o$tt)
  meta_per           <- .parse_kv_list(o$tm)

  term_trees_val <- if (!is.null(o$`term-trees-path`)) o$`term-trees-path` else NULL
  term_meta_val  <- if (!is.null(o$`term-metadata-path`)) o$`term-metadata-path` else NULL

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

  make_plots <- if (isTRUE(o$`no-plots`)) FALSE else isTRUE(o$`make-plots`)
  drop_rows_without_term <- if (isTRUE(o$`keep-rows-without-term`)) FALSE else isTRUE(o$`drop-rows-without-term`)
  keep_unmatched_ids <- if (isTRUE(o$`drop-unmatched-ids`)) FALSE else isTRUE(o$`keep-unmatched-ids`)

  res <- dndsR::term_enrichment(
    dnds_annot_file          = o$`dnds-annot_file`,
    comparison_file          = o$comparison_file,
    output_dir               = o$output_dir,
    terms                    = terms_vec,
    exclude_terms            = exclude_terms_vec,
    sides                    = sides_vec,
    pos_threshold            = as.numeric(o$`pos-threshold`),
    max_dnds                 = as.numeric(o$`max-dnds`),
    filter_expr              = o$`filter-expr`,
    make_plots               = make_plots,
    top_n                    = as.integer(o$`top-n`),
    drop_rows_without_term   = drop_rows_without_term,
    min_total                = as.integer(o$`min-total`),
    min_pos                  = as.integer(o$`min-pos`),
    fdr_method               = o$`fdr-method`,
    alpha                    = as.numeric(o$alpha),
    term_sep                 = o$`term-sep`,
    term_seps                = term_seps_vec,
    term_blocklist           = term_blocklist_vec,
    exclude_ids              = exclude_ids_final,
    term_trees               = term_trees_final,
    exclude_descendants      = isTRUE(o$`exclude-descendants`),
    exclude_descendants_depth= o$`exclude-descendants-depth`,
    exclude_descendants_limit= as.integer(o$`exclude-descendants-limit`),
    term_metadata            = term_meta_final,
    keep_unmatched_ids       = keep_unmatched_ids
  )

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
  make_option(c("--dnds-merged-file"), type="character", dest="dnds_merged_file",
              help="Single prebuilt table with q_/s_ coords (single mode; function mainly intended for batch)"),
  make_option(c("-C","--comparison-file"), type="character", dest="comparison_file",
              help="TSV: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff (batch mode)"),
  make_option(c("-O","--output-dir"), type="character", dest="output_dir", default=getwd(),
              help="Root directory of per-comparison folders [default: %default]"),

  make_option(c("--sides"), type="character", default="query,subject",
              help="Comma-separated among query,subject [default: %default]"),
  make_option(c("--window-size"), type="integer", default=300000,
              help="Window size (bp) for binning genes [default: %default]"),
  make_option(c("--max-dnds"), type="double", default=10,
              help="Drop rows with dNdS >= max_dnds or NA [default: %default]"),
  make_option(c("--filter-expr"), type="character", default=NULL,
              help='Logical filter evaluated on the augmented table (e.g., "q_chr == s_chr")'),

  make_option(c("--make-png"), action="store_true", default=TRUE,
              help="Also write PNG via convertSVG [default: %default]"),
  make_option(c("--no-png"), action="store_true", default=FALSE,
              help="Disable PNG export (overrides --make-png)"),
  make_option(c("--overwrite"), action="store_true", default=FALSE,
              help="Overwrite existing outputs"),
  make_option(c("--keep-intermediate"), action="store_true", default=FALSE,
              help="Keep <comp>_ideogram_input.tsv"),

  # Chromosome label normalization (opt-in; raw by default)
  make_option(c("--chr-strip-leading-chr0"), action="store_true", default=FALSE,
              help='Strip leading "chr" and optional "0" (e.g., chr01 -> 1)'),
  make_option(c("--chr-strip-leading"), type="character", default=NULL,
              help='Regex to strip from start (e.g., "Cf0")'),
  make_option(c("--chr-strip-trailing"), type="character", default=NULL,
              help='Regex to strip from end (e.g., "_random")'),
  make_option(c("--chr-case-insensitive"), action="store_true", default=TRUE,
              help="Apply normalization regex case-insensitively [default: %default]"),
  make_option(c("--chr-case-sensitive"), action="store_true", default=FALSE,
              help="Override to make normalization case-sensitive"),

  # GFF indexing
  make_option(c("--restrict-gff-to-gene"), action="store_true", default=FALSE,
              help='Use only rows with type=="gene" when indexing coordinates')
)

run_dnds_ideogram <- function(o){
  .dump_opts("dnds_ideogram", o)
  sides_vec <- .parse_csv_required(o$sides)
  make_png <- if (isTRUE(o$`no-png`)) FALSE else isTRUE(o$`make-png`)
  chr_case_insensitive <- if (isTRUE(o$`chr-case-sensitive`)) FALSE else isTRUE(o$`chr-case-insensitive`)

  res <- dndsR::dnds_ideogram(
    dnds_merged_file        = o$`dnds-merged-file`,
    comparison_file         = o$`comparison-file`,
    output_dir              = o$`output-dir`,
    sides                   = sides_vec,
    window_size             = as.integer(o$`window-size`),
    max_dnds                = as.numeric(o$`max-dnds`),
    filter_expr             = o$`filter-expr`,
    make_png                = make_png,
    overwrite               = isTRUE(o$overwrite),
    keep_intermediate       = isTRUE(o$`keep-intermediate`),
    chr_strip_leading_chr0  = isTRUE(o$`chr-strip-leading-chr0`),
    chr_strip_leading       = o$`chr-strip-leading`,
    chr_strip_trailing      = o$`chr-strip-trailing`,
    chr_case_insensitive    = chr_case_insensitive,
    restrict_gff_to_gene    = isTRUE(o$`restrict-gff-to-gene`),
    verbose                 = isTRUE(o$verbose)
  )

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
  # keep flags consistent if env var turned it on
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
