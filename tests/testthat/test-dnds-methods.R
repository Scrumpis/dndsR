# dN/dS methods other than Comeron (KaKs_Calculator). Expected values are hand-derived, or
# from PAML yn00 4.10.10, which shares no code with dndsR or KaKs_Calculator.

test_that("only supported methods are accepted", {
  expect_error(calculate_dnds(comparison_name = "x", query_seq = "a", subject_seq = "b", dnds_method = "GY"),
               "dnds_method must be one of")
  expect_error(calculate_dnds(comparison_name = "x", query_seq = "a", subject_seq = "b", dnds_method = c("NG", "YN")),
               "dnds_method must be one of")
  expect_error(calculate_dnds(comparison_name = "x", query_seq = "a", subject_seq = "b", dnds_method = "NG",
                              kaks_calculator = file.path(tempdir(), "missing")), "KaKs_Calculator is needed")
})

test_that("zero changes are detected from the alignment, by hand", {
  z <- .pair_zero_changes
  # TTA (Leu) to CTA (Leu): synonymous only.
  expect_identical(z("TTAGGG", "CTAGGG"), list(no_n = TRUE, no_s = FALSE))
  # GGG (Gly) to AGG (Arg): one position, nonsynonymous only.
  expect_identical(z("TTAGGG", "TTAAGG"), list(no_n = FALSE, no_s = TRUE))
  # GGG to AAG (Lys): both pathways (via AGG Arg, via GAG Glu) change the amino acid at each step.
  expect_identical(z("TTAGGG", "TTAAAG"), list(no_n = FALSE, no_s = TRUE))
  # TTA to CTC, both Leu: one pathway passes through TTC (Phe), so dN need not be 0, even
  # though the proteins are identical. (dndsR once reported 0 here; NG86 gives dN > 0.)
  expect_identical(z("TTAGGG", "CTCGGG"), list(no_n = FALSE, no_s = FALSE))
  # AGA to CGG, both Arg: both pathways (via CGA, via AGG) stay Arg.
  expect_identical(z("AGA", "CGG"), list(no_n = TRUE, no_s = FALSE))
  # TAC (Tyr) to TGG (Trp): one pathway passes through the stop TAG, so neither is known.
  expect_identical(z("TAC", "TGG"), list(no_n = FALSE, no_s = FALSE))
  # Identical; gaps, ambiguity, and stops are skipped, as KaKs_Calculator skips them.
  expect_identical(z("TTA---GGGTAA", "TTACCCNGGTAG"), list(no_n = TRUE, no_s = TRUE))
  # Site-class methods: only identical compared codons give 0.
  for (m in c("LWL", "LPB", "MLWL")) {
    expect_identical(z("TTAGGG", "CTAGGG", m), list(no_n = FALSE, no_s = FALSE))
    expect_identical(z("TTA---GGG", "TTACCCNGG", m), list(no_n = TRUE, no_s = TRUE))
  }
  expect_identical(z("TTAGGG", "CTAGGG", "YN"), list(no_n = TRUE, no_s = FALSE))
})

test_that("KaKs_Calculator output is read with zeros, ratios, and names preserved", {
  skip_on_os("windows")   # the stand-in below is a /bin/sh script
  dir <- withr::local_tempdir()
  # A fake KaKs_Calculator: fixed output in its format, with "NA" for 0 and "-0".
  fake <- file.path(dir, "KaKs_Calculator")
  out <- c("Sequence\tMethod\tKa\tKs\tKa/Ks",
           "p1\tNG\tNA\tNA\tNA", "p2\tNG\tNA\t0.15\t-0", "p3\tNG\t0.05\tNA\tNA", "p4\tNG\t0.01\t0.2\t0.05", "p5\tNG\tNA\tNA\tNA",
           "p6\tNG\t0.1\t99\t0.00101")
  writeLines(c("#!/bin/sh", 'if [ "$1" = -h ]; then echo "Program: KaKs_Calculator Toolbox"; echo "Version: 2.0, June. 2009"; exit 0; fi',
               'while [ $# -gt 0 ]; do [ "$1" = -o ] && o=$2; shift; done',
               paste0("cat > \"$o\" <<'EOF'\n", paste(out, collapse = "\n"), "\nEOF")), fake)
  Sys.chmod(fake, "755")
  expect_identical(.kaks_find(fake)$version, "2.0")
  pairs <- list(c("TTAGGG", "TTAGGG"),   # identical: 0 and 0
                c("TTAGGG", "CTAGGG"),   # synonymous only: dN 0
                c("TTAGGG", "TTAAGG"),   # nonsynonymous at one position: dS 0
                c("TTAGGG", "CTAAGG"),
                c("TTAGGG", "CCCAAA"),   # NA kept: changes of both kinds
                c("TTAGGG", "CTAAGG"))   # KaKs_Calculator's cap, 99: NA
  names <- list(c("gene-a-1", "gene-b-1"), c("q2", "s2"), c("q3", "s3"), c("q4", "s4"), c("q5", "s5"), c("q6", "s6"))
  files <- vapply(seq_along(pairs), function(k) {
    f <- file.path(dir, sprintf("q%d_pal2nal.aln", k))
    writeLines(c(paste0(">", names[[k]][1]), pairs[[k]][1], paste0(">", names[[k]][2]), pairs[[k]][2]), f); f
  }, character(1))
  est <- .kaks_estimate(files, "NG", fake, dir)
  expect_identical(est$query_id, c("gene-a-1", "q2", "q3", "q4", "q5", "q6"))
  expect_identical(est$subject_id, c("gene-b-1", "s2", "s3", "s4", "s5", "s6"))
  expect_identical(est$dN, c(0, 0, 0.05, 0.01, NA, 0.1))
  expect_identical(est$dS, c(0, 0.15, 0, 0.2, NA, NA))
  expect_equal(est$dNdS, c(NA, 0, NA, 0.05, NA, NA))
  expect_identical(attr(est, "capped"), 1L)
  res <- data.frame(query_id = c("q4", "gene-a-1", "qX"), subject_id = c("s4", "gene-b-1", "sX"), dN = 9, dS = 9, dNdS = 9, other = 1:3)
  r <- .dnds_replace_values(res, est)
  expect_identical(r$dN, c(0.01, 0, NA)); expect_identical(r$other, 1:3); expect_identical(attr(r, "unmatched"), 1L)
})

test_that("an existing result made with another method is reported", {
  dir <- withr::local_tempdir()
  q <- file.path(dir, "q.fa"); s <- file.path(dir, "s.fa"); writeLines(">a\nATG", q); writeLines(">b\nATG", s)
  p <- file.path(dir, "prov.json")
  jsonlite::write_json(.dnds_provenance(q, s, list(dnds_est.method = "Comeron"), "NG"), p, auto_unbox = TRUE)
  expect_null(.provenance_mismatch(p, q, s, "NG"))
  expect_match(.provenance_mismatch(p, q, s, "Comeron"), "dnds_method NG, not Comeron")
  # A record from before dnds_method was recorded.
  jsonlite::write_json(list(query_md5 = unname(tools::md5sum(q)), subject_md5 = unname(tools::md5sum(s)),
                            parameters = list(dnds_est.method = "Comeron")), p, auto_unbox = TRUE)
  expect_null(.provenance_mismatch(p, q, s, "Comeron"))
  expect_match(.provenance_mismatch(p, q, s, "YN"), "Comeron, not YN")
})

test_that("NG and YN agree with PAML yn00 on example pairs", {
  bin <- Sys.which("KaKs_Calculator")
  skip_if(!nzchar(bin), "KaKs_Calculator not installed")
  dir <- withr::local_tempdir()
  truth <- system.file("extdata", "example", "truth", package = "dndsR")
  read <- function(f) { l <- readLines(f); h <- grep("^>", l); stats::setNames(vapply(seq_along(h), function(i)
    paste(l[(h[i] + 1):(c(h[-1] - 1, length(l))[i])], collapse = ""), ""), sub("^>", "", l[h])) }
  a <- read(file.path(truth, "spA_longest_isoform_cds.fasta")); b <- read(file.path(truth, "spB_longest_isoform_cds.fasta"))
  genes <- c("g02", "g03", "g10")
  files <- vapply(genes, function(g) {
    x <- a[[paste0("spA_", g)]]; y <- b[[paste0("spB_", g)]]
    x <- substr(x, 1, nchar(x) - 3); y <- substr(y, 1, nchar(y) - 3)   # terminal stops, as yn00 was given
    f <- file.path(dir, paste0(g, "_pal2nal.aln")); writeLines(c(paste0(">spA_", g), x, paste0(">spB_", g), y), f); f
  }, character(1))
  # PAML yn00 4.10.10 (icode 0, weighting 0, commonf3x4 0), rounded to 4 decimals by yn00.
  paml <- list(NG = rbind(c(0, 0.1454), c(0.0471, 0), c(0.0120, 0.1952)),
               YN = rbind(c(0, 0.1496), c(0.0484, 0), c(0.0119, 0.2100)))
  for (m in names(paml)) {
    est <- .kaks_estimate(files, m, bin, dir)
    expect_identical(est$dN[1], 0); expect_identical(est$dS[2], 0)
    # YN: within yn00's rounding. NG: KaKs_Calculator counts about 1.5% fewer synonymous
    # sites than PAML's NG86 from the same differences, so its dS reads about that much higher
    # and its dN slightly lower (docs/validation/estimators.md).
    rel <- function(x, y) abs(x - y) / pmax(y, 1e-3)
    expect_lt(max(rel(est$dN[2:3], paml[[m]][2:3, 1])), if (m == "NG") 0.01 else 0.005, label = paste(m, "dN"))
    expect_lt(max(rel(est$dS[c(1, 3)], paml[[m]][c(1, 3), 2])), if (m == "NG") 0.03 else 0.001, label = paste(m, "dS"))
    expect_true(is.na(est$dNdS[2]))
  }
})
