# The notice must reach the user and the results, without changing any table.
small_table <- function(dir) {
  set.seed(21); n <- 400
  dnds <- ifelse(stats::runif(n) < 0.1, 2, 0.2)
  go <- sample(c("GO:0006468", "GO:0006508", "GO:0055085"), n, TRUE)
  d <- data.frame(query_id = paste0("q", 1:n), subject_id = paste0("s", 1:n), dN = dnds * 0.1, dS = 0.1, dNdS = dnds, q_len = 300, s_len = 300,
                  q_pfam = sample(c("PF00001", "PF00002"), n, TRUE), s_pfam = sample(c("PF00001", "PF00002"), n, TRUE), q_go = go, s_go = go, stringsAsFactors = FALSE)
  f <- file.path(dir, "N_dnds_annot.tsv")
  utils::write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
  f
}

test_that("term_enrichment warns and records the notice beside its tables", {
  dir <- withr::local_tempdir(); f <- small_table(dir)
  expect_warning(suppressMessages(term_enrichment(dnds_annot_file = f, sides = "query", make_plots = FALSE, threads = 1, verbose = FALSE)), "does not adjust for gene length or dS")
  meta <- jsonlite::read_json(file.path(dir, "N_term_enrichment_notice.json"))
  expect_equal(meta$command, "term_enrichment")
  expect_false(meta$adjusted_for_length_and_dS)
  expect_match(meta$calibration_notice, "about 70%")
  expect_equal(meta$settings$pos_threshold, 1)
  expect_equal(meta$dndsR_version, as.character(utils::packageVersion("dndsR")))
  tab <- utils::read.delim(file.path(dir, "N_q_PFAM_enrichment.tsv"))
  expect_false(any(grepl("notice", names(tab))))
})

test_that("go_enrichment warns and records the notice", {
  skip_if_not_installed("topGO"); skip_if_not_installed("GO.db")
  dir <- withr::local_tempdir(); f <- small_table(dir)
  expect_warning(suppressMessages(go_enrichment(dnds_annot_file = f, sides = "query", ontologies = "BP", make_plots = FALSE, threads = 1, node_min = 1)), "does not adjust for gene length or dS")
  meta <- jsonlite::read_json(file.path(dir, "N_go_enrichment_notice.json"))
  expect_equal(meta$settings$algorithm, "weight01")
  expect_equal(meta$settings$p_adjust, "BH")
  expect_match(meta$calibration_notice, "43%")
})

test_that("two identical runs write identical notices", {
  dir <- withr::local_tempdir(); f <- small_table(dir)
  run <- function() suppressWarnings(suppressMessages(term_enrichment(dnds_annot_file = f, sides = "query", make_plots = FALSE, threads = 1, verbose = FALSE)))
  run(); a <- readLines(file.path(dir, "N_term_enrichment_notice.json"))
  run(); expect_identical(readLines(file.path(dir, "N_term_enrichment_notice.json")), a)
})
