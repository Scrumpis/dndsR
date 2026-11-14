# Per-pair reliability flags. Expected counts are hand-derived, or from
# tools/dev/validation/ng86.py, an independent Python NG86 that shares no code with dndsR.

test_that("NG86 sites and differences match hand counts", {
  # TTA (Leu): position 1 has 1 synonymous change of 3 (CTA); position 2 has 1 sense change
  # (TCA, Ser; TAA and TGA are stops); position 3 has 1 of 3 (TTG).
  expect_equal(.ng86_sites()[["TTA"]], 2 / 3)
  # TTA to CTC: via CTA, two synonymous steps; via TTC (Phe), two nonsynonymous. Mean 1 and 1.
  expect_equal(.ng86_codon_diffs("TTA", "CTC"), c(1, 1))
  # TAC (Tyr) to TGG (Trp): the pathway through TAG (stop) is excluded; via TGC (Cys), 2 nonsynonymous.
  expect_equal(.ng86_codon_diffs("TAC", "TGG"), c(0, 2))
  # Gaps, ambiguity, and stops are not compared.
  k <- .ng86_counts("TTA---NNNTAA", "CTCGGGAAATAG")
  expect_identical(k$codons, 1L); expect_equal(k$syn_differences, 1); expect_equal(k$nonsyn_differences, 1)
})

test_that("NG86 counts agree with the independent ng86.py on the example", {
  truth <- system.file("extdata", "example", "truth", package = "dndsR")
  read <- function(f) { l <- readLines(f); h <- grep("^>", l); stats::setNames(vapply(seq_along(h), function(i)
    paste(l[(h[i] + 1):(c(h[-1] - 1, length(l))[i])], collapse = ""), ""), sub("^>", "", l[h])) }
  a <- read(file.path(truth, "spA_longest_isoform_cds.fasta")); b <- read(file.path(truth, "spB_longest_isoform_cds.fasta"))
  # ng86.py output (Sd, Nd, S), printed to 1 and 4 decimals.
  ref <- list(g02 = c(30, 0, 226.0833), g03 = c(0, 30, 242.6667), g04 = c(40, 8, 231.0833), g09 = c(40, 8, 230.0833))
  for (g in names(ref)) {
    k <- .ng86_counts(a[[paste0("spA_", g)]], b[[paste0("spB_", g)]])
    expect_equal(c(k$syn_differences, k$nonsyn_differences, k$syn_sites), ref[[g]], tolerance = 1e-6, label = g)
  }
})

test_that("flags mark the documented conditions and change no value", {
  dir <- withr::local_tempdir()
  set.seed(1)
  long <- paste(sample(c("GCT", "GGT", "CTG", "AAA", "GAT"), 150, TRUE), collapse = "")
  mut <- function(s, at, to) { substr(s, at, at + 2) <- to; s }
  pairs <- list(
    ok = c(long, Reduce(function(s, k) mut(s, k, "GCC"), seq(1, 60, by = 3) * 3 - 2 + 0, long)),
    short = c("GCTGGTCTG", "GCCGGTCTG"),
    none = c(long, long)
  )
  files <- vapply(names(pairs), function(n) { f <- file.path(dir, paste0(n, "_pal2nal.aln"))
    writeLines(c(paste0(">q_", n), pairs[[n]][1], paste0(">s_", n), pairs[[n]][2]), f); f }, "")
  res <- data.frame(query_id = c("q_ok", "q_short", "q_none", "q_missing"), subject_id = c("s_ok", "s_short", "s_none", "s_missing"),
                    dS = c(1e-6, 0.5, 0, 0.1), dNdS = c(50, 0.1, NA, 0.2))
  before <- res
  f <- .dnds_pair_flags(res, files)
  expect_identical(res, before)
  has <- function(i, x) grepl(paste0("(^|;)", x, "($|;)"), f$flags[i])
  expect_true(has(1, "ds_below_observed")); expect_true(has(1, "extreme_ratio")); expect_true(has(1, "low_ds"))
  expect_false(has(1, "few_codons"))
  expect_true(has(2, "few_codons")); expect_true(has(2, "few_synonymous_differences"))
  expect_true(has(3, "few_synonymous_differences")); expect_false(has(3, "ds_below_observed"))
  expect_true(has(4, "no_alignment"))
  expect_identical(f$codons_compared[2], 3)
})
