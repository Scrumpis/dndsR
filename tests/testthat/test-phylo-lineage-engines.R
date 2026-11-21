# phylo_branch_omega, phylo_absrel, and phylo_relax against codeml and HyPhy run directly,
# on the alignment example's OG0001 (all taxa) and OG0004 (spB missing, so the foreground
# collapses onto spA's branch). Direct inputs were prepared without dndsR
# (tools/dev/validation/phylo-direct-lineage.sh, stored in truth/direct/); commands and
# control files are written out here. HyPhy is compared exactly, as in test-phylo-busted.R.
# codeml is compared within 1e-3 log-likelihood units and 0.1% in omega, since the order of
# sequences and branches can move its optimizer slightly. Slow; skipped without the programs.
skip_on_os("windows")
skip_if_not_installed("ape")
skip_if(fast_tests(), "slow: runs codeml and HyPhy")

ex <- withr::local_tempdir(.local_envir = teardown_env())
suppressMessages(phylo_example(file.path(ex, "ex"), type = "alignments"))
ex <- file.path(ex, "ex")
direct <- file.path(ex, "truth", "direct")
ogs <- c("OG0001", "OG0004")
run_in <- function(dir, bin, args) {
  old <- setwd(dir); on.exit(setwd(old))
  system2(bin, args, stdout = "run.log", stderr = "run.log", stdin = "/dev/null")
}

test_that("phylo_branch_omega reproduces codeml run directly", {
  bin <- tryCatch(dndsR:::.codeml_path(), error = function(e) NULL)
  skip_if(is.null(bin), "codeml not found")
  res <- suppressMessages(phylo_branch_omega(alignment_dir = file.path(ex, "alignments"), species_tree_file = file.path(ex, "species_tree.nwk"),
                                             foreground_file = file.path(ex, "foreground.tsv"), output_dir = ex, threads = test_cores()))
  ctl <- function(model, omega, fix, tree) c("seqfile = alignment.phy", paste("treefile =", tree), "outfile = out.mlc", "noisy = 0",
    "verbose = 1", "runmode = 0", "seqtype = 1", "CodonFreq = 2", "clock = 0", "aaDist = 0", paste("model =", model), "NSsites = 0",
    "icode = 0", "fix_kappa = 0", "kappa = 2", paste("fix_omega =", fix), paste("omega =", omega), "fix_alpha = 1", "alpha = 0",
    "Malpha = 0", "ncatG = 1", "getSE = 0", "RateAncestor = 0", "Small_Diff = .5e-6", "cleandata = 0", "method = 0")
  best <- function(og, model, tree, starts) {
    fits <- lapply(starts, function(w) {
      d <- withr::local_tempdir(); file.copy(file.path(direct, og, c("alignment.phy", tree)), d)
      writeLines(ctl(model, w, 0, tree), file.path(d, "c.ctl")); run_in(d, bin, "c.ctl")
      x <- readLines(file.path(d, "out.mlc"))
      l <- grep("^lnL\\(", x, value = TRUE)
      w_line <- x[grep("^w ratios as node labels", x) + 1]
      list(lnl = as.numeric(sub("^.*\\):\\s*(-?[0-9.]+).*$", "\\1", l)),
           omega = as.numeric(strsplit(trimws(sub("^.*[=:]\\s*", "", grep("^(omega \\(dN/dS\\)|w \\(dN/dS\\) for branches)", x, value = TRUE))), "\\s+")[[1]]),
           tips = if (length(w_line)) w_line else NA_character_)
    })
    fits[[which.max(vapply(fits, `[[`, numeric(1), "lnl"))]]
  }
  for (og in ogs) {
    r <- res[res$orthogroup == og, ]
    m0 <- best(og, 0, "tree_codeml.nwk", c(0.4, 0.1, 2))
    fr <- best(og, 1, "tree_codeml.nwk", c(0.4, 0.1, 2))
    tw <- best(og, 2, "tree_codeml_fg.nwk", c(0.4, 0.1, 2))
    expect_lt(abs(r$lnl_m0 - m0$lnl), 1e-3, label = paste(og, "M0"))
    expect_equal(r$omega_m0, m0$omega, tolerance = 1e-3, label = paste(og, "omega M0"))
    expect_lt(abs(r$lnl_free - fr$lnl), 1e-3, label = paste(og, "free"))
    expect_lt(abs(r$lnl_two_ratio - tw$lnl), 1e-3, label = paste(og, "two-ratio"))
    expect_equal(c(r$omega_bg, r$omega_fg), tw$omega, tolerance = 1e-3, label = paste(og, "omegas"))
    # Terminal-branch omega from codeml's own labelled tree, against dndsR's branch table.
    b <- utils::read.delim(file.path(ex, "branch_omega", "branch_omega_branches.tsv"))
    b <- b[b$orthogroup == og & b$terminal, ]
    for (tip in b$branch) {
      direct_w <- as.numeric(sub(paste0(".*\\b", tip, " #([0-9.eE+-]+).*"), "\\1", fr$tips, perl = TRUE))
      expect_equal(b$omega_free[b$branch == tip], direct_w, tolerance = 1e-3, label = paste(og, tip))
    }
  }
  # Files kept: codeml's control file and main output per model and start, one alignment and
  # one tree per model per orthogroup; none of codeml's auxiliary files.
  kept <- list.files(file.path(ex, "branch_omega", "paml", "OG0001"), recursive = TRUE)
  expect_false(any(basename(kept) %in% c("rst", "rst1", "rub", "lnf", "codeml.log")))
  expect_lte(length(kept), 25L)
  # keep_runs = FALSE leaves no codeml folder for a finished orthogroup.
  lean <- suppressMessages(phylo_branch_omega(alignment_dir = file.path(ex, "alignments"), species_tree_file = file.path(ex, "species_tree.nwk"),
                                              output_dir = ex, run_name = "lean", models = "M0", keep_runs = FALSE, threads = test_cores()))
  expect_false(dir.exists(file.path(ex, "lean", "paml", "OG0001")))
  expect_equal(lean$lnl_m0[lean$orthogroup == "OG0001"], res$lnl_m0[res$orthogroup == "OG0001"], tolerance = 1e-6)
  # OG0004 lacks spB: the foreground is spA's branch, merged with the clade's stem.
  b4 <- utils::read.delim(file.path(ex, "branch_omega", "branch_omega_branches.tsv"))
  b4 <- b4[b4$orthogroup == "OG0004" & b4$foreground, ]
  expect_identical(b4$branch, "spA")
  expect_match(b4$species_tree_branches, "+", fixed = TRUE)
})

hyphy_bin <- function() {
  bin <- tryCatch(dndsR:::.hyphy_path(), error = function(e) NULL)
  skip_if(is.null(bin), "HyPhy not found")
  bin
}

test_that("phylo_absrel reproduces aBSREL run directly", {
  bin <- hyphy_bin()
  res <- suppressMessages(phylo_absrel(alignment_dir = file.path(ex, "alignments"), species_tree_file = file.path(ex, "species_tree.nwk"),
                                       output_dir = ex, threads = test_cores()))
  b <- utils::read.delim(file.path(ex, "absrel", "absrel_branches.tsv"))
  for (og in ogs) {
    d <- withr::local_tempdir(); file.copy(file.path(direct, og, c("alignment.fas", "tree_absrel.nwk")), d)
    run_in(d, bin, c("CPU=1", "ENV=RANDOM_SEED=1", "absrel", "--alignment", "alignment.fas", "--tree", "tree_absrel.nwk",
                     "--branches", "FG", "--code", "Universal", "--srv", "No", "--multiple-hits", "None", "--output", "a.json"))
    j <- jsonlite::read_json(file.path(d, "a.json"))
    expect_identical(res$lnl_adaptive[res$orthogroup == og], j$fits[["Full adaptive model"]][["Log Likelihood"]], label = og)
    ba <- j[["branch attributes"]][["0"]]
    tips <- b$branch[b$orthogroup == og & b$terminal]
    # dndsR's own JSON is HyPhy's output, so compared exactly; the TSV, written to 15
    # significant digits, within 1e-12.
    mine <- jsonlite::read_json(file.path(ex, "absrel", res$json[res$orthogroup == og]))[["branch attributes"]][["0"]]
    for (tip in tips) {
      row <- b[b$orthogroup == og & b$branch == tip, ]
      for (k in c("Uncorrected P-value", "Corrected P-value", "LRT", "Baseline MG94xREV omega ratio")) {
        expect_identical(mine[[tip]][[k]], ba[[tip]][[k]], label = paste(og, tip, k))
      }
      expect_equal(c(row$p_raw, row$p_holm, row$baseline_omega),
                   as.numeric(c(ba[[tip]][["Uncorrected P-value"]], ba[[tip]][["Corrected P-value"]], ba[[tip]][["Baseline MG94xREV omega ratio"]])),
                   tolerance = 1e-12, label = paste(og, tip))
    }
    expect_identical(res$n_branches_tested[res$orthogroup == og], length(ba))
  }
})

test_that("phylo_relax reproduces RELAX run directly", {
  bin <- hyphy_bin()
  res <- suppressMessages(phylo_relax(alignment_dir = file.path(ex, "alignments"), species_tree_file = file.path(ex, "species_tree.nwk"),
                                      foreground_file = file.path(ex, "foreground.tsv"), output_dir = ex, threads = test_cores(), models = "Minimal"))
  for (og in ogs) {
    d <- withr::local_tempdir(); file.copy(file.path(direct, og, c("alignment.fas", "tree_relax.nwk")), d)
    run_in(d, bin, c("CPU=1", "ENV=RANDOM_SEED=1", "relax", "--alignment", "alignment.fas", "--tree", "tree_relax.nwk",
                     "--test", "TEST", "--reference", "REFERENCE", "--code", "Universal", "--models", "Minimal", "--rates", "3",
                     "--srv", "No", "--multiple-hits", "None", "--grid-size", "250", "--starting-points", "1", "--output", "r.json"))
    tr <- jsonlite::read_json(file.path(d, "r.json"))[["test results"]]
    r <- res[res$orthogroup == og, ]
    expect_identical(c(r$k, r$lrt, r$p_value), c(tr[["relaxation or intensification parameter"]], tr$LRT, tr[["p-value"]]), label = og)
    expect_identical(r$direction, if (r$k > 1) "intensification" else "relaxation")
  }
})
