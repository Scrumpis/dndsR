# phylo_busted with real HyPhy on the bundled examples, against HyPhy run directly, on this
# machine, on inputs prepared without dndsR (made by tools/dev/validation/phylo-direct-hyphy*.sh
# and stored in each example's truth/direct/), with a command written out in full
# (.verify_direct_hyphy). Stored values are not used: HyPhy's arithmetic differs slightly
# between processors with and without AVX-512. Skipped without HyPhy (and MAFFT for the
# genome example); slow, so not in fast mode.
hyphy_or_skip <- function() {
  skip_on_os("windows")
  skip_if_not_installed("ape")
  skip_if(fast_tests(), "slow: runs HyPhy")
  bin <- tryCatch(.hyphy_path(), error = function(e) NULL)
  skip_if(is.null(bin), "HyPhy not found")
  bin
}
json_values <- function(dir, res, og) {
  values(jsonlite::read_json(file.path(dir, res$json[res$orthogroup == og])))
}
values <- function(j) c(lrt = j[["test results"]]$LRT, p = j[["test results"]][["p-value"]],
                        lnl = j$fits[["Unconstrained model"]][["Log Likelihood"]])
direct_here <- function(ex, runs, bin) lapply(.verify_direct_hyphy(file.path(ex, "truth", "direct"), runs, bin, test_cores()), values)

test_that("the alignment example reproduces direct HyPhy runs, and the example check passes", {
  bin <- hyphy_or_skip()
  work <- withr::local_tempdir()
  ex <- file.path(work, "ex")
  phylo_example(ex, type = "alignments")
  run <- function(...) suppressMessages(phylo_busted(alignment_dir = file.path(ex, "alignments"), species_tree_file = file.path(ex, "species_tree.nwk"),
                                                     foreground_file = file.path(ex, "foreground.tsv"), output_dir = ex, ...))
  res <- run(threads = test_cores())
  here <- direct_here(ex, paste0("OG000", 1:4), bin)
  for (og in names(here)) expect_identical(json_values(file.path(ex, "busted"), res, og), here[[og]], label = og)
  expect_silent(suppressMessages(verify_phylo_example(ex)))

  # A rerun with the same inputs, settings, and HyPhy reuses results, even with a HyPhy
  # that would fail; a changed seed reruns.
  broken <- file.path(work, "hyphy-broken")
  writeLines(c("#!/bin/sh", sprintf('if [ "$1" = --version ]; then "%s" --version; exit 0; fi', bin), "exit 1"), broken)
  Sys.chmod(broken, "0755")
  expect_identical(run(hyphy = broken)$p_value, res$p_value)
  expect_warning(run(hyphy = broken, seed = 2), "4 orthogroup\\(s\\) failed")
})

test_that("extra starts add diagnostics without changing the primary result", {
  bin <- hyphy_or_skip()
  work <- withr::local_tempdir()
  ex <- file.path(work, "ex")
  phylo_example(ex, type = "alignments")
  dir.create(file.path(ex, "one"))
  file.copy(file.path(ex, "alignments", "OG0001.fasta"), file.path(ex, "one"))
  res <- suppressMessages(phylo_busted(alignment_dir = file.path(ex, "one"), species_tree_file = file.path(ex, "species_tree.nwk"),
                                       foreground_file = file.path(ex, "foreground.tsv"), output_dir = ex, extra_seeds = 2:3, threads = 1))
  expect_identical(res$lrt, direct_here(ex, "OG0001", bin)$OG0001[["lrt"]])
  expect_identical(res$n_starts, 3L)
  expect_identical(res$start_seeds, "1,2,3")
  expect_true(all(file.exists(file.path(ex, "busted", "hyphy", "OG0001", c("start_seed2", "start_seed3"), "busted.json"))))
  expect_true(res$lrt_min_starts <= res$lrt && res$lrt <= res$lrt_max_starts)
  expect_true(res$lnl_unconstrained_best >= res$lnl_unconstrained)
  prov <- jsonlite::read_json(file.path(ex, "busted", "busted_provenance.json"), simplifyVector = TRUE)
  expect_identical(prov$extra_seeds, 2:3)
  # One tested orthogroup: BH leaves its p-value as it is, and the column sits beside p_best_fits.
  expect_identical(res$p_best_fits_adj, res$p_best_fits)
  expect_identical(match("p_best_fits_adj", names(res)), match("p_best_fits", names(res)) + 1L)
})

test_that("genomes to BUSTED: preparation, results, gene trees, and direct HyPhy agree", {
  bin <- hyphy_or_skip()
  skip_if(is.null(tryCatch(.mafft_path(), error = function(e) NULL)), "MAFFT not found")
  work <- withr::local_tempdir()
  ex <- file.path(work, "ex")
  phylo_example(ex)
  old <- setwd(ex)
  on.exit(setwd(old), add = TRUE)
  suppressMessages(phylo_prepare("orthogroup_table.tsv", "species_tree.nwk", "foreground.tsv", genome_file = "genomes.tsv", threads = test_cores()))
  res <- suppressMessages(phylo_busted(prepared_dir = "prepared", threads = test_cores()))
  expect_silent(suppressMessages(verify_phylo_example()))
  expect_identical(sum(res$status == "tested"), 8L)
  expect_true(all(res$stage[res$status == "excluded"] == "prepare"))

  # Gene trees, given to HyPhy, against direct runs on hand-written gene trees.
  suppressMessages(phylo_prepare("orthogroup_table.tsv", "species_tree.nwk", "foreground.tsv", genome_file = "genomes.tsv",
                                 gene_tree_dir = "gene_trees", run_name = "prepared_gt"))
  gt <- suppressMessages(phylo_busted(prepared_dir = "prepared_gt", run_name = "busted_gt", threads = test_cores()))
  here <- direct_here(".", c("OG04_gene", "OG05_gene"), bin)
  for (og in c("OG04", "OG05")) expect_identical(json_values("busted_gt", gt, og), here[[paste0(og, "_gene")]], label = og)
  st <- function(og) gt[gt$orthogroup == og, ]
  expect_identical(st("OG06")$reason, "foreground_not_in_gene_tree")
  expect_identical(st("OG01")$reason, "gene_tree_missing")
  expect_identical(c(st("OG04")$gt_missing_species_splits, st("OG05")$gt_missing_species_splits), c(1L, 2L))
  expect_identical(c(st("OG04")$gt_discordant_near_foreground, st("OG05")$gt_discordant_near_foreground), c(FALSE, TRUE))

  # The prepared hypothesis cannot be changed between steps.
  writeLines(c("scope\ttaxa", "stem\tspA,spB"), "prepared/foreground.tsv")
  expect_error(suppressMessages(phylo_busted(prepared_dir = "prepared", run_name = "tampered")), "differs from the one phylo_prepare used")
})
