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

# ---- compute_dnds (example mapping) ----
opts_compute <- list(
  make_option(c("-a","--orthologs"), type="character", help="Ortholog pairs TSV"),
  make_option(c("-c","--cds"), type="character", help="CDS FASTA"),
  make_option(c("-m","--method"), type="character", default="yn00",
              help="dN/dS method (yn00|codeml|...) [default: %default]"),
  make_option(c("-o","--out"), type="character", default="dnds.tsv",
              help="Output TSV [default: %default]")
)
run_compute <- function(o){
  need <- c("orthologs","cds"); miss <- need[!nzchar(unlist(o[need]))]
  if (length(miss)) stop("Missing: ", paste(miss, collapse=", "))
  dndsR::compute_dnds(
    orthologs = o$orthologs,
    cds       = o$cds,
    method    = o$method,
    out       = o$out,
    threads   = dndsR_cli_threads(o$threads, 8L)
  )
  cli_alert_success("Wrote: {o$out}")
}

# ---- append_gff ----
opts_append <- list(
  make_option(c("-i","--in"), type="character", dest="infile", help="Input dN/dS TSV/CSV"),
  make_option(c("-g","--gff"), type="character", help="Annotations GFF3"),
  make_option(c("-o","--out"), type="character", default="dnds.annot.tsv",
              help="Output file [default: %default]"),
  make_option(NULL, "--tsv", action="store_true", default=FALSE,
              help="Write TSV instead of CSV")
)
run_append <- function(o){
  stopifnot(!is.null(o$infile), !is.null(o$gff))
  fmt <- if (o$tsv) "tsv" else "csv"
  dndsR::append_gff_annotations(
    input  = o$infile, gff = o$gff,
    out    = o$out, format = fmt,
    threads= dndsR_cli_threads(o$threads, 8L)
  )
  cli_alert_success("Wrote: {o$out}")
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
  split_by_label = list(opts = c(base_opts, opts_split),   fun = run_split,
                        help = "Split genomes by label and emit a new comparison file"),
  compute_dnds   = list(opts = c(base_opts, opts_compute), fun = run_compute,
                        help = "Compute dN/dS from orthologs + CDS"),
  append_gff     = list(opts = c(base_opts, opts_append),  fun = run_append,
                        help = "Append GFF annotations to dN/dS table"),
  enrichment     = list(opts = c(base_opts, opts_enrich),  fun = run_enrich,
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
