#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(cli); library(dndsR)
})

base_opts <- list(
  make_option(c("-t","--threads"), type="integer", default=8,
              help="CPUs/threads [default: %default]"),
  make_option(c("-v","--verbose"), action="store_true", default=FALSE,
              help="Verbose logging")
)

# parse --oa key=value pairs with light type casting
.cast_scalar <- function(x) {
  if (x %in% c("TRUE","FALSE")) return(as.logical(x))
  if (x == "NULL") return(NULL)
  if (grepl("^\\d+$", x)) return(as.integer(x))
  if (grepl("^\\d*\\.?\\d+(e[+-]?\\d+)?$", x, ignore.case = TRUE)) return(as.numeric(x))
  x
}

.parse_oa <- function(v) {
  if (is.null(v)) return(list())
  out <- list()
  for (kv in v) {
    # split only on the first '='
    sp <- regexpr("=", kv, fixed = TRUE)
    if (sp < 1) stop("Invalid --oa argument (expected key=value): ", kv)
    key <- substr(kv, 1, sp - 1)
    val <- substr(kv, sp + 1, nchar(kv))
    # support comma-separated vectors: a,b,c
    if (grepl(",", val, fixed = TRUE)) {
      parts <- strsplit(val, ",", fixed = TRUE)[[1]]
      out[[key]] <- vapply(parts, .cast_scalar, FUN.VALUE = NA, USE.NAMES = FALSE)
    } else {
      out[[key]] <- .cast_scalar(val)
    }
  }
  out
}

# ---- split_by_label ----
opts_split <- list(
  make_option(c("-C","--comparison-file"), type="character", dest="comparison_file",
              help="TSV with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff"),
  make_option(c("-m","--mode"), type="character", default="subgenome",
              help="subgenome|haplotype|custom [default: %default]"),
  make_option(c("-r","--custom-regex"), type="character", default=NULL,
              help="Regex with exactly one capture group when --mode=custom"),
  make_option(NULL, "--case-insensitive", action="store_true", default=TRUE,
              help="Case-insensitive label matching [default: %default]")
)
run_split <- function(o){
  stopifnot(!is.null(o$comparison_file))
  p <- dndsR::split_comparisons_by_label(
    comparison_file = o$comparison_file,
    mode            = o$mode,
    custom_regex    = o$`custom-regex`,
    case_insensitive= o$`case-insensitive`
  )
  cli_alert_success("Wrote: {p}")
}

# ---- extract_cds ----
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
  # note: --verbose is global via base_opts
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
  # Let the R function enforce its own required-arg rules.
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

# ---- calculate_dnds ----
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
  make_option(c("--comp-cores"), type="integer", default=4,
              help="Cores passed to orthologr::dNdS via comp_cores [default: %default]"),
  make_option(c("--aligner"), type="character", default="diamond",
              help="orthologr aligner [default: %default]"),
  make_option(c("--sensitivity-mode"), type="character", default="fast",
              help="orthologr sensitivity_mode [default: %default]"),
  make_option(c("--dnds-method"), type="character", default="Comeron",
              help="orthologr dnds_est.method [default: %default]"),
  make_option(c("--oa"), type = "character", action = "append", default = NULL,
              help = "Pass through to orthologr::dNdS as key=value (repeatable). Examples: --oa eval=1E-6 --oa ortho_detection=RBH --oa clean_folders=TRUE")

  # NOTE: your function accepts ...; if you want CLI passthrough later,
  # we can add a --orthologr-args JSON or key=val repeatable flag.
)

run_calculate_dnds <- function(o){
  extra <- .parse_oa(o$oa)
  base <- list(
    comparison_file = o$comparison_file,
    comparison_name = o$comparison_name,
    subject_fasta   = o$subject_fasta,
    query_fasta     = o$query_fasta,
    subject_gff     = o$subject_gff,
    query_gff       = o$query_gff,
    output_dir      = o$`output-dir`,
    comp_cores      = as.integer(o$`comp-cores`),
    aligner         = o$aligner,
    sensitivity_mode= o$`sensitivity-mode`,
    dnds_method     = o$`dnds-method`
  )
  # Forward everything, letting your R function handle defaults + merges
  paths <- do.call(dndsR::calculate_dnds, c(base, extra))

  if (is.character(paths) && length(paths)) {
    cli::cli_alert_success("calculate_dnds completed. Wrote {length(paths)} file{?s}.")
    for (p in paths) cli::cli_bullets(c(v = paste0("<", p, ">")))
  } else {
    cli::cli_alert_success("calculate_dnds completed.")
  }
  invisible(paths)
}

# ---- append_annotations ----
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
              help = "Whitespace-delimited file with: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff (batch mode)"),
  make_option(c("-O","--output-dir"), type = "character", dest = "output_dir", default = getwd(),
              help = "Root directory containing per-comparison folders (batch mode) [default: %default]"),
  make_option(c("--custom"), type = "character", default = NULL,
              help = "Comma patterns like 'REX{5},TOM{8}' → q_prefix/s_prefix columns")
  # --verbose comes from base_opts
)

run_append_annotations <- function(o) {
  res <- dndsR::append_annotations(
    dnds_file       = o$dnds_file,
    query_gff       = o$query_gff,
    subject_gff     = o$subject_gff,
    output_file     = o$output_file,
    comparison_file = o$comparison_file,
    output_dir      = o$output_dir,
    custom          = o$custom
  )
  if (is.character(res) && length(res)) {
    cli::cli_alert_success("append_annotations completed. Wrote {length(res)} file{?s}.")
    for (p in res) cli::cli_bullets(c(v = paste0("<", p, ">")))
  } else {
    cli::cli_alert_success("append_annotations completed.")
  }
  invisible(res)
}

# ---- enrichment ----
opts_enrich <- list(
  make_option(c("-i","--in"), type="character", dest="infile", help="Annotated table"),
  make_option(c("-p","--pval"), type="double", default=0.05,
              help="BH-adjusted p-value cutoff [default: %default]"),
  make_option(c("-o","--out"), type="character", default="enrichment.tsv",
              help="Output TSV [default: %default]")
)
run_enrich <- function(o){
  stopifnot(!is.null(o$infile))
  dndsR::term_enrichment(
    input = o$infile, alpha = o$pval,
    out   = o$out, threads = dndsR_cli_threads(o$threads, 8L)
  )
  cli_alert_success("Wrote: {o$out}")
}

subcommands <- list(
  split_by_label = list(opts = c(base_opts, opts_split),      fun = run_split,
                        help = "Split genomes by label and emit a new comparison file"),
  extract_cds    = list(opts = c(base_opts, opts_extract_cds), fun = run_extract_cds,
                        help = "Extract CDS (and optionally proteins) from genome+GFF"),
  calculate_dnds = list(opts = c(base_opts, opts_calculate_dnds), fun = run_calculate_dnds,
                        help = "Run dN/dS with orthologr from pre-extracted CDS"),
  append_annotations = list(opts = c(base_opts, opts_append_annotations), fun = run_append_annotations,
                            help = "Append query/subject GFF annotations to dN/dS results"),
  enrichment     = list(opts = c(base_opts, opts_enrich),     fun = run_enrich,
                        help = "Run enrichment analysis on annotated table")
)

usage <- function(){
  cat("\nUsage: dndsR.R <subcommand> [options]\n\nSubcommands:\n")
  for (nm in names(subcommands))
    cat(sprintf("  %-14s %s\n", nm, subcommands[[nm]]$help))
  cat("\nUse: dndsR.R <subcommand> --help   for subcommand options\n\n")
}

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) == 0 || argv[1] %in% c("-h","--help")) { usage(); quit(status=0) }
cmd <- argv[1]; args <- argv[-1]
if (!cmd %in% names(subcommands)) { cli::cli_alert_danger(paste("Unknown:", cmd)); usage(); quit(status=1) }
spec <- subcommands[[cmd]]
parser <- optparse::OptionParser(usage=paste0("%prog ", cmd, " [options]"), option_list = spec$opts)
opts <- optparse::parse_args(parser, args=args, positional_arguments = FALSE)

tryCatch({
  spec$fun(opts)
  if (isTRUE(opts$verbose)) cli::cli_alert_success("Done.")
}, error = function(e){
  cli::cli_alert_danger(conditionMessage(e)); quit(status=1)
})
