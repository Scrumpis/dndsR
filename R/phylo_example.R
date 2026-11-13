#' Copy the bundled phylogenetic example to a folder
#'
#' Copies seven synthetic codon alignments of up to eight taxa, a species tree, a
#' foreground file, and the design and direct HyPhy results used by
#' [verify_phylo_example()]. The sequences were simulated and have no biological meaning.
#' Four orthogroups are testable; three are designed to be skipped.
#'
#' @param output_dir Folder to create. It must not already contain the example.
#' @return The folder path, invisibly.
#' @examples
#' folder <- phylo_example(file.path(tempdir(), "dndsr_phylo_example"))
#' list.files(folder, recursive = TRUE)
#' @export
phylo_example <- function(output_dir = "dndsr_phylo_example") {
  source_dir <- system.file("extdata", "phylo_example", package = "dndsR", mustWork = TRUE)
  if (file.exists(file.path(output_dir, "foreground.tsv"))) {
    stop("The phylogenetic example is already in ", output_dir, ". Choose another folder or remove it.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(list.files(source_dir, full.names = TRUE), output_dir, recursive = TRUE)
  if (!all(ok)) stop("Could not copy the example to ", output_dir, call. = FALSE)
  message("Phylogenetic example copied to ", output_dir)
  invisible(output_dir)
}

#' Check phylogenetic example results against the design and direct HyPhy runs
#'
#' Run this after [phylo_busted()] on the bundled example (see [phylo_example()]), with its
#' default settings. It confirms that dndsR prepares inputs as designed, labels the
#' intended branches, and reproduces HyPhy run directly on independently prepared inputs.
#' It checks the software, not the biology: it is not evidence that BUSTED detects
#' selection correctly in your data.
#'
#' @param output_dir Folder containing the example and its results.
#' @param run_name The \code{run_name} given to [phylo_busted()].
#' @return A data frame with one row per check, invisibly. It stops if any check fails.
#' @export
verify_phylo_example <- function(output_dir = ".", run_name = "busted") {
  path <- function(...) file.path(output_dir, ...)
  results_file <- path(run_name, paste0(run_name, "_results.tsv"))
  if (!file.exists(results_file)) {
    stop("No results at ", results_file, ". Run phylo_busted on the example first.", call. = FALSE)
  }
  res <- utils::read.delim(results_file, colClasses = c(orthogroup = "character"), na.strings = "NA")
  branches <- utils::read.delim(path(run_name, paste0(run_name, "_branches.tsv")), na.strings = "NA")
  prov <- jsonlite::read_json(path(run_name, paste0(run_name, "_provenance.json")))
  design <- utils::read.delim(path("truth", "design.tsv"), colClasses = "character")
  direct <- utils::read.delim(path("truth", "direct_hyphy.tsv"), colClasses = c(orthogroup = "character", hyphy_version = "character"))

  checks <- list()
  check <- function(name, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = isTRUE(passed), detail = detail)
  }
  row <- function(og) res[match(og, res$orthogroup), ]

  got <- ifelse(res$status == "tested", "tested", paste0(res$status, ":", res$reason))
  want <- design$expected_status[match(res$orthogroup, design$orthogroup)]
  check("every orthogroup has its designed status", setequal(res$orthogroup, design$orthogroup) && identical(got, want),
        paste(sprintf("%s=%s", res$orthogroup, got), collapse = ", "))
  check("terminal stop codons masked in OG0001", identical(as.integer(row("OG0001")$terminal_stops_masked), 8L),
        sprintf("%s masked", row("OG0001")$terminal_stops_masked))
  fg <- function(og) sort(branches$branch[branches$orthogroup == og & branches$foreground])
  check("clade hypothesis labels spA, spB, and their stem", identical(row("OG0001")$n_fg_branches, 3L) &&
          identical(fg("OG0001")[fg("OG0001") %in% c("spA", "spB")], c("spA", "spB")))
  merged <- branches$species_tree_branches[branches$orthogroup == "OG0004" & branches$foreground]
  check("with spB missing, spA's branch includes the clade's stem, and is foreground",
        identical(fg("OG0004"), "spA") && length(merged) == 1L && grepl("mrca(spA,spB)", merged, fixed = TRUE) && grepl("spA", merged),
        paste("species-tree branches:", merged))

  # Agreement with HyPhy run directly. Values come from the JSON, at full precision.
  same_setup <- identical(prov$random_seed, 1L) && identical(prov$hyphy_cpu, 1L) &&
    identical(prov$settings$srv, "Yes") && identical(prov$settings$multiple_hits, "None") &&
    identical(prov$settings$error_sink, "No") && identical(prov$settings$rates, 3L) && identical(prov$settings$syn_rates, 3L) &&
    identical(prov$settings$starting_points, 1L) && identical(prov$settings$grid_size, 250L)
  check("run used the reference settings (seed 1, 1 HyPhy thread, HyPhy defaults)", same_setup,
        if (same_setup) "" else "rerun phylo_busted on the example with its defaults")
  gaps <- vapply(direct$orthogroup, function(og) {
    j <- tryCatch(jsonlite::read_json(path(run_name, row(og)$json)), error = function(e) NULL)
    if (is.null(j)) return(Inf)
    d <- direct[direct$orthogroup == og, ]
    con <- j$fits[["Constrained model"]][["Log Likelihood"]]
    max(abs(c(j[["test results"]]$LRT - d$lrt, j[["test results"]][["p-value"]] - d$p_value,
              j$fits[["Unconstrained model"]][["Log Likelihood"]] - d$lnl_unconstrained,
              if (is.na(d$lnl_constrained)) if (is.null(con)) 0 else Inf else con - d$lnl_constrained)))
  }, numeric(1))
  check("LRT, p-value, and log-likelihoods equal direct HyPhy runs (tolerance 1e-6)", all(gaps < 1e-6),
        sprintf("largest difference %.3g; HyPhy %s here, %s in the reference", max(gaps), prov$hyphy_version, direct$hyphy_version[1]))
  check("tested branches equal the reference count", identical(as.integer(row(direct$orthogroup)$n_fg_branches), as.integer(direct$n_tested_branches)))

  tested <- res$status == "tested"
  bh <- stats::p.adjust(res$p_value[tested], "BH")
  check("p_adj is Benjamini-Hochberg over the tested orthogroups", isTRUE(all.equal(res$p_adj[tested], bh, tolerance = 1e-12)) &&
          all(res$family_size == sum(tested)) && all(is.na(res$p_adj[!tested])), sprintf("family of %d", sum(tested)))
  md5_ok <- all(vapply(which(tested), function(i) identical(unname(tools::md5sum(path(run_name, res$json[i]))), res$json_md5[i]), logical(1)))
  check("provenance records HyPhy's version and each JSON's hash", md5_ok && !is.null(prov$hyphy_version) && !is.null(prov$random_seed))

  out <- do.call(rbind, checks)
  for (i in seq_len(nrow(out))) {
    message(sprintf("[%s] %s%s", if (out$passed[i]) "PASS" else "FAIL", out$check[i],
                    if (nzchar(out$detail[i])) paste0(" (", out$detail[i], ")") else ""))
  }
  if (!all(out$passed)) stop(sum(!out$passed), " of ", nrow(out), " phylogenetic example checks failed.", call. = FALSE)
  message(sprintf("All %d phylogenetic example checks passed.", nrow(out)))
  invisible(out)
}
