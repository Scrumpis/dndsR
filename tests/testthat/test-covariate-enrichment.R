# Expected values here come from pbinom, phyper, and hand arithmetic, not from the code under test.
test_that("the Poisson-binomial tail matches the binomial when probabilities are equal", {
  for (k in c(1L, 3L, 7L)) {
    expect_equal(dndsR:::.poisson_binomial_upper(rep(0.2, 15), k), stats::pbinom(k - 1L, 15, 0.2, lower.tail = FALSE), tolerance = 1e-12)
  }
})

test_that("the Poisson-binomial tail matches hand calculation with unequal probabilities", {
  p <- c(0.1, 0.5, 0.9)
  # P(all three) = 0.045, P(exactly two) = .1*.5*.1 + .1*.5*.9 + .9*.5*.9 = 0.455
  expect_equal(dndsR:::.poisson_binomial_upper(p, 3L), 0.045, tolerance = 1e-12)
  expect_equal(dndsR:::.poisson_binomial_upper(p, 2L), 0.045 + 0.455, tolerance = 1e-12)
  expect_equal(dndsR:::.poisson_binomial_upper(p, 0L), 1)
  expect_equal(dndsR:::.poisson_binomial_upper(p, 4L), 0)
})

# A term of short genes. Positive genes are short, whatever their term.
confounded <- function() {
  set.seed(11)
  n <- 3000
  len <- round(exp(stats::runif(n, log(50), log(3000))))
  ds <- exp(stats::runif(n, log(0.01), log(1)))
  positive <- stats::runif(n) < stats::plogis(2.5 - 1.2 * log(len))
  term <- ifelse(len < 120 & stats::runif(n) < 0.5, "IPR000001", ifelse(stats::runif(n) < 0.05, "IPR000002", "IPR000003"))
  data.frame(query_id = paste0("g", seq_len(n)), subject_id = paste0("s", seq_len(n)), dN = ifelse(positive, 2, 0.2) * ds, dS = ds,
             dNdS = ifelse(positive, 2, 0.2), q_len = len, s_len = len, q_ipr = term, s_ipr = NA, stringsAsFactors = FALSE)
}

test_that("a term confounded with length passes Fisher's test but not the adjusted test", {
  dir <- withr::local_tempdir(); f <- file.path(dir, "N_dnds_annot.tsv")
  d <- confounded()
  utils::write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  out <- suppressMessages(covariate_enrichment(dnds_annot_file = f, terms = "ipr", sides = "query"))
  got <- utils::read.delim(out)
  row <- got[got$term == "IPR000001", ]
  in_term <- d$q_ipr == "IPR000001"; k <- sum(d$dNdS[in_term] > 1); K <- sum(d$dNdS > 1)
  expect_equal(row$pos_count, k)
  expect_equal(row$total_count, sum(in_term))
  expect_equal(row$pos_total, K)
  expect_equal(row$fisher_p, stats::phyper(k - 1L, K, nrow(d) - K, sum(in_term), lower.tail = FALSE), tolerance = 1e-10)
  expect_lt(row$fisher_p, 1e-4)
  expect_gt(row$p_value, 0.01)
  expect_equal(row$n_tested, 3L)
})

test_that("a term with a real excess beyond length and dS stays significant", {
  dir <- withr::local_tempdir(); f <- file.path(dir, "N_dnds_annot.tsv")
  d <- confounded()
  extra <- which(d$q_ipr == "IPR000002" & d$q_len > 500)[1:25]
  d$dNdS[extra] <- 2; d$dN[extra] <- 2 * d$dS[extra]
  utils::write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  got <- utils::read.delim(suppressMessages(covariate_enrichment(dnds_annot_file = f, terms = "ipr", sides = "query")))
  expect_lt(got$p_adj[got$term == "IPR000002"], 0.001)
})

test_that("with no link to length or dS the fit is close to the overall rate, and BH matches by hand", {
  dir <- withr::local_tempdir(); f <- file.path(dir, "N_dnds_annot.tsv")
  set.seed(5); n <- 4000
  terms <- sample(sprintf("IPR%06d", 1:5), n, TRUE)
  d <- data.frame(query_id = paste0("g", 1:n), subject_id = paste0("s", 1:n), dN = 0.1, dS = exp(stats::runif(n, -4, 0)),
                  dNdS = ifelse(stats::runif(n) < 0.1, 2, 0.2), q_len = round(exp(stats::runif(n, 4, 8))), s_len = 100, q_ipr = terms, stringsAsFactors = FALSE)
  utils::write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
  got <- utils::read.delim(suppressMessages(covariate_enrichment(dnds_annot_file = f, terms = "ipr", sides = "query", min_pos = 0)))
  expect_equal(sum(got$expected), sum(d$dNdS > 1), tolerance = 1e-6)
  o <- order(got$p_value); manual <- got$p_value[o] * nrow(got) / seq_len(nrow(got)); manual <- rev(cummin(rev(manual)))
  expect_equal(got$p_adj[o], pmin(1, manual), tolerance = 1e-10)
})

test_that("GO terms reach their ancestors, each ontology has its own universe, and BH spans them", {
  skip_if_not_installed("GO.db")
  dir <- withr::local_tempdir(); f <- file.path(dir, "N_dnds_annot.tsv")
  set.seed(3); n <- 600
  # GO:0005576 is cellular component, and GO:0008150 the root of biological process.
  go <- ifelse(seq_len(n) <= 300, "GO:0005576", ifelse(seq_len(n) <= 500, "GO:0006355", NA))
  d <- data.frame(query_id = paste0("g", 1:n), subject_id = paste0("s", 1:n), dN = 0.1, dS = exp(stats::runif(n, -4, 0)),
                  dNdS = ifelse(stats::runif(n) < 0.2, 2, 0.2), q_len = round(exp(stats::runif(n, 4, 8))), s_len = 100, q_go = go, stringsAsFactors = FALSE)
  utils::write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  got <- utils::read.delim(suppressMessages(covariate_enrichment(dnds_annot_file = f, terms = "go", sides = "query", min_pos = 0)))
  cc <- got[got$ontology == "CC", ]; bp <- got[got$ontology == "BP", ]
  expect_true(all(cc$universe == 300L))
  expect_true(all(bp$universe == 200L))
  expect_true("GO:0005575" %in% cc$term) # the root of cellular component, an ancestor
  expect_equal(got$n_tested[1], nrow(got))
  expect_equal(unique(cc$total_count[cc$term == "GO:0005575"]), 300L)
})

test_that("a floor of a tenth of the positive rate limits the tail, and few positives bring a notice", {
  dir <- withr::local_tempdir(); f <- file.path(dir, "N_dnds_annot.tsv")
  set.seed(9); n <- 2000
  # All 6 positive genes are the 6 longest, so a term of the shortest genes has a fit near zero.
  len <- round(exp(seq(log(60), log(6000), length.out = n)))
  ds <- exp(stats::runif(n, -4, 0))
  positive <- rank(len) > n - 6
  term <- ifelse(rank(len) <= 40, "IPR000001", "IPR000002")
  positive[1:2] <- TRUE # 2 positive in the short term; 8 positive overall
  d <- data.frame(query_id = paste0("g", 1:n), subject_id = paste0("s", 1:n), dN = 0.1, dS = ds, dNdS = ifelse(positive, 2, 0.2),
                  q_len = len, s_len = len, q_ipr = term, stringsAsFactors = FALSE)
  utils::write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_message(out <- covariate_enrichment(dnds_annot_file = f, terms = "ipr", sides = "query", verbose = FALSE), "Power is low")
  row <- utils::read.delim(out)
  row <- row[row$term == "IPR000001", ]
  # Rate 8/2000 = 0.004, floor 0.0004. The expected count over 40 genes cannot be below 40 * 0.0004.
  expect_gte(row$expected, 40 * 0.0004)
  expect_gte(row$p_value, stats::pbinom(1, 40, 0.0004, lower.tail = FALSE) * (1 - 1e-9))
})

test_that("missing inputs fail with a clear message", {
  dir <- withr::local_tempdir(); f <- file.path(dir, "N_dnds_annot.tsv")
  utils::write.table(data.frame(query_id = "a", dS = 1, dNdS = 1), f, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(covariate_enrichment(dnds_annot_file = f, sides = "query"), "q_len")
  expect_error(covariate_enrichment(), "either comparison_file")
  expect_error(covariate_enrichment(dnds_annot_file = file.path(dir, "no.tsv")), "not found")
})

test_that("a quantile threshold equals the same number given directly, and is recorded", {
  dir <- withr::local_tempdir(); f <- file.path(dir, "N_dnds_annot.tsv")
  d <- confounded()
  set.seed(3); d$dNdS <- exp(stats::rnorm(nrow(d), -1, 1)); d$dNdS[d$q_len < 120] <- d$dNdS[d$q_len < 120] * 3; d <- d[d$dNdS < 10, ]
  utils::write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  by_q <- utils::read.delim(suppressMessages(covariate_enrichment(dnds_annot_file = f, terms = "ipr", sides = "query", pos_threshold = "q0.95")))
  cutoff <- stats::quantile(d$dNdS, 0.95, names = FALSE, type = 7)
  by_n <- utils::read.delim(suppressMessages(covariate_enrichment(dnds_annot_file = f, terms = "ipr", sides = "query", pos_threshold = cutoff)))
  expect_equal(by_q, by_n)
  expect_equal(unique(by_q$threshold), cutoff)
  expect_equal(unique(by_q$pos_total), sum(d$dNdS > cutoff))
  expect_equal(sum(d$dNdS > cutoff), ceiling(0.05 * nrow(d)))
  spelled <- utils::read.delim(suppressMessages(covariate_enrichment(dnds_annot_file = f, terms = "ipr", sides = "query", pos_threshold = "quantile:0.95")))
  expect_equal(spelled, by_q)
})

test_that("a threshold that cannot be used stops, rather than finding no positive genes", {
  dir <- withr::local_tempdir(); f <- file.path(dir, "N_dnds_annot.tsv")
  utils::write.table(confounded(), f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  for (bad in list("q95", "high", c(1, 2), NA_real_, "q1.5")) {
    expect_error(covariate_enrichment(dnds_annot_file = f, terms = "ipr", sides = "query", pos_threshold = bad), "pos_threshold must be")
  }
  expect_error(suppressMessages(term_enrichment(dnds_annot_file = f, pos_threshold = "q0.95")), "pos_threshold must be one number")
  expect_error(suppressMessages(ipr_enrichment(dnds_annot_file = f, pos_threshold = "q0.95")), "pos_threshold must be one number")
})
