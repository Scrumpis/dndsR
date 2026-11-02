#' Copy the bundled example to a folder
#'
#' Copies two small synthetic genomes, their annotations, a comparison file, and the
#' designed truth used by [verify_example()]. The sequences have no biological meaning.
#'
#' @param output_dir Folder to create. It must not already contain the example.
#' @return The folder path, invisibly.
#' @examples
#' folder <- dndsr_example(file.path(tempdir(), "dndsr_example"))
#' list.files(folder)
#' @export
dndsr_example <- function(output_dir = "dndsr_example") {
  source_dir <- system.file("extdata", "example", package = "dndsR", mustWork = TRUE)
  if (file.exists(file.path(output_dir, "comparison_file.txt"))) {
    stop("The example is already in ", output_dir, ". Choose another folder or remove it.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(list.files(source_dir, full.names = TRUE), output_dir, recursive = TRUE)
  if (!all(ok)) stop("Could not copy the example to ", output_dir, call. = FALSE)
  message("Example copied to ", output_dir)
  invisible(output_dir)
}

#' Check example results against independently established values
#'
#' Run this after `extract_cds` and `calculate_dnds` on the bundled example (see
#' [dndsr_example()]). It confirms that the installation works as expected. It cannot
#' show that results for your own data are correct.
#'
#' Checks:
#' \itemize{
#'   \item CDS match the sequences designed before the genomes were generated.
#'   \item Invalid CDS were excluded.
#'   \item Each designed class of gene has dN/dS in its expected range.
#'   \item dN and dS are close to an independent Nei-Gojobori calculation. The methods
#'     differ, so this uses a tolerance rather than exact equality.
#'   \item dN and dS match a stored reference run, which detects a change in behavior
#'     between versions.
#' }
#'
#' @param output_dir Folder containing the example and its results.
#' @param comparison Name of the comparison to check.
#' @return A data frame with one row per check, invisibly. It stops with an error if any
#'   check fails.
#' @export
verify_example <- function(output_dir = ".", comparison = "example_AvB") {
  path <- function(...) file.path(output_dir, ...)
  dnds_file <- path(comparison, paste0(comparison, "_dnds.tsv"))
  if (!file.exists(dnds_file)) {
    stop("No results at ", dnds_file, ". Run extract_cds and calculate_dnds on the example first.", call. = FALSE)
  }
  result <- utils::read.delim(dnds_file)
  result$pair <- sub("^spA_", "", result$query_id)
  value <- function(pair, column) result[[column]][match(pair, result$pair)]

  checks <- list()
  check <- function(name, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = isTRUE(passed), detail = detail)
  }
  read_fasta <- function(file) {
    lines <- readLines(file, warn = FALSE)
    starts <- grep("^>", lines)
    ends <- c(starts[-1] - 1L, length(lines))
    seqs <- mapply(function(s, e) paste(lines[seq.int(s + 1L, length.out = e - s)], collapse = ""), starts, ends)
    stats::setNames(toupper(seqs), sub("^>(\\S+).*", "\\1", lines[starts]))
  }

  # Designed sequences. g30's longest valid isoform is 45 nt longer by design.
  for (species in c("spA", "spB")) {
    got <- read_fasta(path(comparison, paste0(species, "_CDS.fasta")))
    want <- read_fasta(path("truth", paste0(species, "_longest_isoform_cds.fasta")))
    same <- !grepl("g30$", names(got))
    check(
      sprintf("%s CDS match the designed sequences", species),
      all(names(got) %in% names(want)) && identical(unname(got[same]), unname(want[names(got)[same]])) &&
        identical(unname(nchar(got[!same])), unname(nchar(want[names(got)[!same]]) + 45L)),
      sprintf("%d sequences", length(got))
    )
  }

  excluded <- c("g31", "g32", "g33", "g34")
  check("invalid and unmatched genes are excluded", !any(excluded %in% result$pair), paste(excluded, collapse = ", "))
  check("multi-exon minus-strand gene g28 is conserved", value("g28", "dNdS") < 0.2,
        sprintf("dN/dS = %.3f; it was 1.6 before the fix in 0.3.0", value("g28", "dNdS")))

  designed <- utils::read.delim(path("truth", "designed_substitutions.tsv"))
  ratios <- function(class) value(designed$pair[designed$class == class], "dNdS")
  in_range <- function(class, low, high) {
    x <- ratios(class)
    check(sprintf("%s genes: dN/dS in %g to %g", class, low, high), all(!is.na(x) & x >= low & x <= high),
          sprintf("observed %.3f to %.3f", min(x), max(x)))
  }
  in_range("conserved", 0, 0.2)
  in_range("neutral_like", 0.7, 1.3)
  in_range("excess_nonsynonymous", 2, Inf)
  in_range("structure", 0, 0.2)
  in_range("multi_isoform", 0, 0.2)
  check("identical sequences: dN = dS = 0, dN/dS undefined",
        value("g01", "dN") == 0 && value("g01", "dS") == 0 && is.na(value("g01", "dNdS")))
  check("synonymous changes only: dN = 0, dS > 0", value("g02", "dN") == 0 && value("g02", "dS") > 0)

  # Independent reference. Method differences are bigger when dS is small.
  ng86 <- utils::read.delim(path("truth", "ng86_reference.tsv"), na.strings = c("NA", ""))
  shared <- setdiff(intersect(result$pair, ng86$pair), c("g01", "g02", "g03", "g30"))
  dn_gap <- abs(value(shared, "dN") - ng86$dN[match(shared, ng86$pair)])
  ds_gap <- abs(value(shared, "dS") - ng86$dS[match(shared, ng86$pair)])
  check("dN close to independent Nei-Gojobori values",
        all(dn_gap < 0.01 + 0.15 * ng86$dN[match(shared, ng86$pair)]), sprintf("largest gap %.4f", max(dn_gap)))
  check("dS close to independent Nei-Gojobori values",
        all(ds_gap < 0.015 + 0.25 * ng86$dS[match(shared, ng86$pair)]), sprintf("largest gap %.4f", max(ds_gap)))

  reference <- utils::read.delim(path("truth", "reference_run.tsv"))
  gap <- max(abs(value(reference$pair, "dN") - reference$dN), abs(value(reference$pair, "dS") - reference$dS))
  check("dN and dS match the reference run", !is.na(gap) && gap < 0.001, sprintf("largest gap %.5f", gap))

  out <- do.call(rbind, checks)
  for (i in seq_len(nrow(out))) {
    message(sprintf("[%s] %s%s", if (out$passed[i]) "PASS" else "FAIL", out$check[i],
                    if (nzchar(out$detail[i])) paste0(" (", out$detail[i], ")") else ""))
  }
  if (!all(out$passed)) stop(sum(!out$passed), " of ", nrow(out), " example checks failed.", call. = FALSE)
  message(sprintf("All %d example checks passed.", nrow(out)))
  invisible(out)
}
