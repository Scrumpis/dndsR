#' Copy a bundled phylogenetic example to a folder
#'
#' Two synthetic examples, simulated with no biological meaning:
#' \itemize{
#'   \item \code{"genomes"} (default): six genome FASTA and GFF3 files, an orthogroup table
#'     in dndsR's long format and in OrthoFinder's format, a species tree, a foreground,
#'     three gene trees, and the design. Fourteen orthogroups; eight can be tested, and six
#'     are designed to be excluded (a paralog, invalid gene models, missing taxa). Run
#'     [phylo_prepare()] and then [phylo_busted()].
#'   \item \code{"alignments"}: seven codon alignments of up to eight taxa, for
#'     [phylo_busted()] directly. Four can be tested; three are designed to be excluded.
#' }
#' Each includes the direct HyPhy results used by [verify_phylo_example()].
#'
#' @param output_dir Folder to create. It must not already contain the example.
#' @param type "genomes" or "alignments".
#' @return The folder path, invisibly.
#' @examples
#' folder <- phylo_example(file.path(tempdir(), "dndsr_phylo_example"))
#' list.files(folder, recursive = TRUE)
#' @export
phylo_example <- function(output_dir = "dndsr_phylo_example", type = c("genomes", "alignments")) {
  type <- match.arg(type)
  source_dir <- system.file("extdata", if (type == "genomes") "phylo_genomes" else "phylo_example", package = "dndsR", mustWork = TRUE)
  if (file.exists(file.path(output_dir, "foreground.tsv"))) {
    stop("A phylogenetic example is already in ", output_dir, ". Choose another folder or remove it.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(list.files(source_dir, full.names = TRUE), output_dir, recursive = TRUE)
  if (!all(ok)) stop("Could not copy the example to ", output_dir, call. = FALSE)
  message(sprintf("Phylogenetic example (%s) copied to %s", type, output_dir))
  invisible(output_dir)
}

#' Check phylogenetic example results against the design and direct HyPhy runs
#'
#' Run this in the example folder (see [phylo_example()]) after the documented commands,
#' with default settings. It confirms that dndsR builds inputs as designed, excludes what it
#' should for the stated reasons, labels the intended branches, and reproduces HyPhy, which
#' it runs directly here on inputs prepared without dndsR (\code{truth/direct/}), with a
#' command written out in full rather than dndsR's own. HyPhy's results differ slightly
#' between processors with and without AVX-512, so stored results from both kinds are
#' reported for comparison, not checked. It checks the software, not the biology: it
#' is not evidence that BUSTED detects selection correctly in your data.
#'
#' @param output_dir Folder containing the example and its results.
#' @param run_name The \code{run_name} given to [phylo_busted()].
#' @param prepared_name The \code{run_name} given to [phylo_prepare()] (genomes example).
#' @param threads HyPhy runs at once for the direct comparison.
#' @return A data frame with one row per check, invisibly. It stops if any check fails.
#' @export
verify_phylo_example <- function(output_dir = ".", run_name = "busted", prepared_name = "prepared", threads = 2L) {
  path <- function(...) file.path(output_dir, ...)
  genomes <- file.exists(path("genomes.tsv"))
  results_file <- path(run_name, paste0(run_name, "_results.tsv"))
  if (!file.exists(results_file)) {
    stop("No results at ", results_file, ". Run the example's commands first.", call. = FALSE)
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
  fg <- function(og) sort(branches$branch[branches$orthogroup == og & branches$foreground])
  status <- ifelse(res$status == "tested", "tested", paste0(res$status, ":", res$reason))
  want_col <- if (genomes) "expected_busted" else "expected_status"
  check("every orthogroup has its designed outcome",
        setequal(res$orthogroup, design$orthogroup) && identical(status, design[[want_col]][match(res$orthogroup, design$orthogroup)]),
        paste(sprintf("%s=%s", res$orthogroup, status), collapse = ", "))

  if (genomes) {
    pre_dir <- path(prepared_name)
    pre <- utils::read.delim(file.path(pre_dir, paste0(prepared_name, "_summary.tsv")), colClasses = c(orthogroup = "character"), na.strings = "NA")
    excluded <- design[design$expected_basis != "", ]
    check("exclusions have their designed basis and stage",
          identical(res$exclusion_basis[match(excluded$orthogroup, res$orthogroup)], excluded$expected_basis) &&
            all(res$stage[match(excluded$orthogroup, res$orthogroup)] == "prepare"),
          paste(sprintf("%s=%s", excluded$orthogroup, res$exclusion_basis[match(excluded$orthogroup, res$orthogroup)]), collapse = ", "))
    read_fa <- function(file) {
      lines <- readLines(file, warn = FALSE)
      starts <- grep("^>", lines)
      ends <- c(starts[-1] - 1L, length(lines))
      seqs <- mapply(function(s, e) paste(lines[seq.int(s + 1L, length.out = e - s)], collapse = ""), starts, ends)
      stats::setNames(toupper(seqs), sub("^>([^|[:space:]]+).*", "\\1", lines[starts]))
    }
    # Every prepared alignment, without gaps, is the designed CDS without its stop codon.
    truth <- read_fa(path("truth", "designed_cds.fasta"))
    members <- utils::read.delim(path("orthogroup_table.tsv"), colClasses = "character")
    prepared <- pre$orthogroup[pre$status == "prepared"]
    same <- vapply(prepared, function(og) {
      a <- readLines(file.path(pre_dir, "alignments", paste0(og, ".fasta")))
      ids <- sub("^>[^|]*\\|", "", a[c(TRUE, FALSE)])
      seqs <- gsub("-", "", a[c(FALSE, TRUE)], fixed = TRUE)
      want <- truth[ids]
      identical(unname(seqs), unname(substr(want, 1L, nchar(want) - 3L))) &&
        setequal(ids, members$gene_id[members$orthogroup == og])
    }, logical(1))
    check("each prepared alignment holds exactly the designed coding sequences, stop codons removed", all(same),
          sprintf("%d alignments", length(same)))
    og10 <- read_fa(file.path(pre_dir, "alignments", "OG10.fasta"))
    true10 <- read_fa(path("truth", "OG10_true_codon_alignment.fasta"))
    check("codon insertion and deletion are aligned as designed (OG10)",
          identical(unname(substring(og10[names(true10)], 4L)), unname(true10)))
    merged <- branches$species_tree_branches[branches$orthogroup == "OG09" & branches$foreground]
    check("with spB missing, spA's branch includes the clade's stem, and is foreground",
          identical(fg("OG09"), "spA") && length(merged) == 1L && grepl("mrca(spA,spB)", merged, fixed = TRUE),
          paste("species-tree branches:", merged))
    check("foreground of complete orthogroups is spA, spB, and their stem",
          all(vapply(c("OG01", "OG02", "OG03"), function(og) length(fg(og)) == 3L && all(c("spA", "spB") %in% fg(og)), logical(1))))
    direct <- direct[direct$topology == "species", ]
  } else {
    check("terminal stop codons masked in OG0001", identical(as.integer(row("OG0001")$terminal_stops_masked), 8L),
          sprintf("%s masked", row("OG0001")$terminal_stops_masked))
    check("clade hypothesis labels spA, spB, and their stem", identical(row("OG0001")$n_fg_branches, 3L) &&
            identical(fg("OG0001")[fg("OG0001") %in% c("spA", "spB")], c("spA", "spB")))
    merged <- branches$species_tree_branches[branches$orthogroup == "OG0004" & branches$foreground]
    check("with spB missing, spA's branch includes the clade's stem, and is foreground",
          identical(fg("OG0004"), "spA") && length(merged) == 1L && grepl("mrca(spA,spB)", merged, fixed = TRUE) && grepl("spA", merged),
          paste("species-tree branches:", merged))
  }

  # Agreement with HyPhy run directly, here, on inputs prepared without dndsR
  # (truth/direct/). Values come from the JSON, at full precision.
  same_setup <- identical(prov$random_seed, 1L) && identical(prov$hyphy_cpu, 1L) && !length(prov$extra_seeds) &&
    identical(prov$settings$srv, "Yes") && identical(prov$settings$multiple_hits, "None") &&
    identical(prov$settings$error_sink, "No") && identical(prov$settings$rates, 3L) && identical(prov$settings$syn_rates, 3L) &&
    identical(prov$settings$starting_points, 1L) && identical(prov$settings$grid_size, 250L)
  check("run used the reference settings (seed 1, 1 HyPhy thread, HyPhy defaults)", same_setup,
        if (same_setup) "" else "rerun with the defaults")
  runs <- if (genomes) unique(direct$run) else unique(direct$orthogroup)
  og_of <- if (genomes) stats::setNames(direct$orthogroup, direct$run)[runs] else stats::setNames(runs, runs)
  here <- .verify_direct_hyphy(path("truth", "direct"), runs, prov$hyphy_path, threads)
  values <- function(j) {
    con <- j$fits[["Constrained model"]][["Log Likelihood"]]
    c(j[["test results"]]$LRT, j[["test results"]][["p-value"]], j$fits[["Unconstrained model"]][["Log Likelihood"]],
      if (is.null(con)) NA_real_ else con)
  }
  gaps <- vapply(runs, function(r) {
    j <- tryCatch(jsonlite::read_json(path(run_name, row(og_of[[r]])$json)), error = function(e) NULL)
    if (is.null(j) || is.null(here[[r]])) return(Inf)
    a <- values(j); b <- values(here[[r]])
    if (!identical(is.na(a), is.na(b))) return(Inf)
    max(abs(a - b), na.rm = TRUE)
  }, numeric(1))
  check("LRT, p-value, and log-likelihoods equal HyPhy run directly here, on inputs prepared without dndsR", all(gaps < 1e-9),
        sprintf("%d orthogroups; largest difference %.3g", length(gaps), max(gaps)))
  n_ref <- if (genomes) direct$tested_branches else direct$n_tested_branches
  check("tested branches equal the reference count", identical(as.integer(row(og_of)$n_fg_branches), as.integer(n_ref[match(runs, if (genomes) direct$run else direct$orthogroup)])))
  # Stored results from two machine classes: reported, not a check. HyPhy's arithmetic
  # differs slightly between processors with and without AVX-512.
  match_class <- vapply(unique(direct$cpu_class), function(cl) {
    d <- direct[direct$cpu_class == cl, ]
    all(vapply(runs, function(r) {
      x <- d[(if (genomes) d$run else d$orthogroup) == r, ]
      !is.null(here[[r]]) && isTRUE(all.equal(values(here[[r]])[1], x$lrt, tolerance = 0))
    }, logical(1)))
  }, logical(1))
  notes <- sprintf("HyPhy %s here; results here %s", prov$hyphy_version,
                   if (any(match_class)) paste0("equal the stored results from a ", names(match_class)[match_class][1], " machine")
                   else "differ from both stored results (avx512, avx2): see truth/direct_hyphy.tsv")

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
  message("[NOTE] ", notes)
  if (!all(out$passed)) stop(sum(!out$passed), " of ", nrow(out), " phylogenetic example checks failed.", call. = FALSE)
  message(sprintf("All %d phylogenetic example checks passed.", nrow(out)))
  invisible(out)
}

#' Run HyPhy BUSTED directly on stored example inputs, independently of dndsR's runner
#'
#' The command is written out in full here on purpose: it must not share code with the
#' runner it checks.
#' @keywords internal
.verify_direct_hyphy <- function(direct_dir, runs, hyphy_path, threads = 2L) {
  bin <- if (!is.null(hyphy_path) && file.exists(hyphy_path)) hyphy_path else Sys.which("hyphy")
  if (!nzchar(bin)) stop("HyPhy is needed to check the example.", call. = FALSE)
  work <- tempfile("phylo-verify-")
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  one <- function(r) {
    d <- file.path(work, r)
    dir.create(d, recursive = TRUE)
    file.copy(file.path(direct_dir, r, c("alignment.fas", "tree.nwk")), d)
    old <- setwd(d)
    on.exit(setwd(old))
    system2(bin, c("CPU=1", "ENV=RANDOM_SEED=1", "busted", "--alignment", "alignment.fas", "--tree", "tree.nwk",
                   "--branches", "FG", "--code", "Universal", "--srv", "Yes", "--rates", "3", "--syn-rates", "3",
                   "--multiple-hits", "None", "--error-sink", "No", "--grid-size", "250", "--starting-points", "1",
                   "--output", "busted.json"), stdout = "hyphy.log", stderr = "hyphy.log")
    tryCatch(jsonlite::read_json("busted.json"), error = function(e) NULL)
  }
  out <- if (.Platform$OS.type == "windows" || threads <= 1L) lapply(runs, one) else parallel::mclapply(runs, one, mc.cores = threads)
  stats::setNames(out, runs)
}
