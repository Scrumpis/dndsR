# phylo_busted with real HyPhy on the bundled example, against HyPhy run directly on
# independently prepared inputs (tools/dev/validation/phylo-direct-hyphy.sh, stored in
# truth/direct_hyphy.tsv). Skipped without HyPhy; slow, so not in fast mode.
hyphy_or_skip <- function() {
  skip_on_os("windows")
  skip_if_not_installed("ape")
  skip_if(fast_tests(), "slow: runs HyPhy")
  bin <- tryCatch(.hyphy_path(), error = function(e) NULL)
  skip_if(is.null(bin), "HyPhy not found")
  bin
}

test_that("the example reproduces direct HyPhy runs, and the example check passes", {
  bin <- hyphy_or_skip()
  work <- withr::local_tempdir()
  phylo_example(file.path(work, "ex"))
  ex <- file.path(work, "ex")
  res <- suppressMessages(phylo_busted(file.path(ex, "alignments"), file.path(ex, "species_tree.nwk"),
                                       file.path(ex, "foreground.tsv"), output_dir = ex, threads = test_cores()))
  direct <- utils::read.delim(file.path(ex, "truth", "direct_hyphy.tsv"))
  for (og in direct$orthogroup) {
    j <- jsonlite::read_json(file.path(ex, "busted", res$json[res$orthogroup == og]))
    d <- direct[direct$orthogroup == og, ]
    expect_equal(j[["test results"]]$LRT, d$lrt, tolerance = 1e-6, label = paste(og, "LRT"))
    expect_equal(j$fits[["Unconstrained model"]][["Log Likelihood"]], d$lnl_unconstrained, tolerance = 1e-6, label = paste(og, "lnL"))
  }
  expect_silent(suppressMessages(verify_phylo_example(ex)))

  # A rerun with the same inputs, settings, and HyPhy reuses results, even with a HyPhy
  # that would fail; a changed seed reruns.
  broken <- file.path(work, "hyphy-broken")
  writeLines(c("#!/bin/sh", sprintf('if [ "$1" = --version ]; then "%s" --version; exit 0; fi', bin), "exit 1"), broken)
  Sys.chmod(broken, "0755")
  again <- suppressMessages(phylo_busted(file.path(ex, "alignments"), file.path(ex, "species_tree.nwk"),
                                         file.path(ex, "foreground.tsv"), output_dir = ex, hyphy = broken))
  expect_identical(again$p_value, res$p_value)
  expect_warning(suppressMessages(phylo_busted(file.path(ex, "alignments"), file.path(ex, "species_tree.nwk"),
                                               file.path(ex, "foreground.tsv"), output_dir = ex, hyphy = broken, seed = 2)),
                 "4 orthogroup\\(s\\) failed")
})
