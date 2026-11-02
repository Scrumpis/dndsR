# Constructed cases whose answers follow directly from the hypergeometric distribution.
upper_tail <- function(a, term, positive, universe) {
  stats::phyper(a - 1L, term, universe - term, positive, lower.tail = FALSE)
}

annotated_table <- function(path, dnds, ipr = NA_character_, go = NA_character_) {
  ids <- sprintf("g%03d", seq_along(dnds))
  d <- data.frame(query_id = ids, subject_id = paste0("s", ids), dN = dnds * 0.1, dS = 0.1, dNdS = dnds,
                  q_ipr = ipr, s_ipr = ipr, q_go = go, s_go = go)
  utils::write.table(d, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

domains <- function(n) {
  entries <- utils::read.delim(system.file("extdata", "interpro_entry.list.tsv", package = "dndsR"))
  entries$ENTRY_AC[entries$ENTRY_TYPE == "Domain"][seq_len(n)]
}

run_ipr <- function(dir, ...) {
  suppressMessages(ipr_enrichment(
    dnds_annot_file = file.path(dir, "t_dnds_annot.tsv"), sides = "query", make_plots = FALSE,
    entries_source = "local", entries_path = system.file("extdata", "interpro_entry.list.tsv", package = "dndsR"),
    auto_detect_release = FALSE, max_prop = 1, ...
  ))
  utils::read.delim(file.path(dir, "t_q_IPR_Domain_enrichment.tsv"))
}

test_that("InterPro Fisher tests and BH adjustment match hand calculation", {
  dir <- withr::local_tempdir()
  ids <- domains(2)
  # 40 annotated genes: 10 of 12 positive genes carry term 1, against 2 of 28 others.
  term <- c(rep(ids[1], 10), rep(ids[2], 2), rep(ids[1], 2), rep(ids[2], 26))
  dnds <- c(rep(2, 12), rep(0.1, 28))
  annotated_table(file.path(dir, "t_dnds_annot.tsv"), dnds, ipr = term)
  got <- run_ipr(dir, fdr_method = "BH", tree_source = "none", min_total = 1, min_pos = 1)
  row <- got[got$IPR == ids[1], ]
  p1 <- upper_tail(10L, 12L, 12L, 40L)
  p2 <- upper_tail(2L, 28L, 12L, 40L)
  expect_equal(row$pos_count, 10L)
  expect_equal(row$total_count, 12L)
  expect_equal(row$p_value, p1, tolerance = 1e-10)
  expect_equal(got$p_value[got$IPR == ids[2]], p2, tolerance = 1e-10)
  expect_equal(row$p_adj, stats::p.adjust(c(p1, p2), "BH")[1], tolerance = 1e-10)
  expect_equal(row$enrichment, (10 / 12) / (2 / 28), tolerance = 1e-10)
})

test_that("local InterPro files need no release lookup and no internet access", {
  dir <- withr::local_tempdir()
  ids <- domains(2)
  term <- c(rep(ids[1], 10), rep(ids[2], 2), rep(ids[1], 2), rep(ids[2], 26))
  annotated_table(file.path(dir, "t_dnds_annot.tsv"), c(rep(2, 12), rep(0.1, 28)), ipr = term)
  withr::local_envvar(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9", ftp_proxy = "http://127.0.0.1:9")
  # auto_detect_release stays at its default, TRUE.
  ipr_enrichment(dnds_annot_file = file.path(dir, "t_dnds_annot.tsv"), sides = "query", make_plots = FALSE,
                 max_prop = 1, min_total = 1, min_pos = 1, tree_source = "none", entries_source = "local",
                 entries_path = system.file("extdata", "interpro_entry.list.tsv", package = "dndsR"))
  got <- utils::read.delim(file.path(dir, "t_q_IPR_Domain_enrichment.tsv"))
  expect_equal(got$p_value[got$IPR == ids[1]], upper_tail(10L, 12L, 12L, 40L), tolerance = 1e-10)
})

test_that("InterPro enrichment works with no internet access", {
  dir <- withr::local_tempdir()
  ids <- domains(2)
  term <- c(rep(ids[1], 10), rep(ids[2], 2), rep(ids[1], 2), rep(ids[2], 26))
  annotated_table(file.path(dir, "t_dnds_annot.tsv"), c(rep(2, 12), rep(0.1, 28)), ipr = term)
  # A proxy that refuses connections stands in for a compute node without internet.
  withr::local_envvar(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9", ftp_proxy = "http://127.0.0.1:9")
  expect_message(
    ipr_enrichment(dnds_annot_file = file.path(dir, "t_dnds_annot.tsv"), sides = "query", make_plots = FALSE,
                   max_prop = 1, min_total = 1, min_pos = 1, entries_timeout_s = 2, tree_timeout_s = 2),
    "bundled"
  )
  got <- utils::read.delim(file.path(dir, "t_q_IPR_Domain_enrichment.tsv"))
  expect_equal(got$p_value[got$IPR == ids[1]], upper_tail(10L, 12L, 12L, 40L), tolerance = 1e-10)
  prov <- jsonlite::read_json(file.path(dir, "t_interpro_provenance.json"))
  expect_match(paste(unlist(prov), collapse = " "), "bundled")
})

test_that("the parent-child method removes genes explained by a significant child", {
  dir <- withr::local_tempdir()
  ids <- domains(3)
  parent <- ids[1]; child <- ids[2]; lonely <- ids[3]
  # Genes 1-10 carry child and parent and are positive; 11-20 carry only the parent and
  # are not. 21-25 carry an unrelated term and are positive. 75 others have no term.
  term <- c(rep(paste(parent, child, sep = ";"), 10), rep(parent, 10), rep(lonely, 5), rep(NA_character_, 75))
  dnds <- c(rep(2, 10), rep(0.1, 10), rep(2, 5), rep(0.1, 75))
  annotated_table(file.path(dir, "t_dnds_annot.tsv"), dnds, ipr = term)
  got <- run_ipr(
    dir, method = "parent_child", fdr_method = "none", drop_rows_without_term = FALSE,
    tree_source = "local", term_trees = data.frame(parent = parent, child = child),
    min_total = 2, min_pos = 2
  )
  expect_true(child %in% got$IPR)
  expect_equal(got$p_value[got$IPR == child], upper_tail(10L, 10L, 15L, 100L), tolerance = 1e-10)
  # Its positive genes all belong to the child, so nothing is left to test.
  expect_false(parent %in% got$IPR)
  # A term outside the tree is still tested.
  expect_true(lonely %in% got$IPR)
  expect_equal(got$p_value[got$IPR == lonely], upper_tail(5L, 5L, 15L, 100L), tolerance = 1e-10)
})

test_that("terms of unknown type are not tested in any stratum", {
  dir <- withr::local_tempdir()
  ids <- domains(1)
  term <- c(rep(paste(ids, "IPR999999", sep = ";"), 10), rep(ids, 10), rep(NA_character_, 10))
  annotated_table(file.path(dir, "t_dnds_annot.tsv"), c(rep(2, 10), rep(0.1, 20)), ipr = term)
  got <- run_ipr(dir, fdr_method = "none", tree_source = "none", min_total = 1, min_pos = 1)
  expect_false("IPR999999" %in% got$IPR)
})

test_that("GO tests count descendants and match hand calculation", {
  skip_if_not_installed("topGO")
  skip_if_not_installed("GO.db")
  dir <- withr::local_tempdir()
  # GO:0006952 (defense response) is a child of GO:0006950 (response to stress).
  go <- c(rep("GO:0006952", 10), rep("GO:0006950", 10), rep("GO:0008152", 20))
  dnds <- c(rep(2, 12), rep(0.1, 28)) # positives: all 10 children and 2 parent-only genes
  annotated_table(file.path(dir, "t_dnds_annot.tsv"), dnds, go = go)
  run <- function(...) {
    suppressMessages(go_enrichment(
      dnds_annot_file = file.path(dir, "t_dnds_annot.tsv"), sides = "query", ontologies = "BP",
      make_plots = FALSE, node_min = 1, ...
    ))
    utils::read.delim(file.path(dir, sprintf("t_q_GO_BP_topGO_%s_fisher.tsv", list(...)$algorithm)), quote = "")
  }
  got <- run(algorithm = "classic", p_adjust = "BH")
  child <- got[got$GO_ID == "GO:0006952", ]
  parent <- got[got$GO_ID == "GO:0006950", ]
  expect_equal(child$annotated, 10L)
  expect_equal(child$significant, 10L)
  expect_equal(parent$annotated, 20L) # includes the child's genes
  expect_equal(parent$significant, 12L)
  expect_equal(child$raw_p, upper_tail(10L, 10L, 12L, 40L), tolerance = 1e-8)
  expect_equal(parent$raw_p, upper_tail(12L, 20L, 12L, 40L), tolerance = 1e-8)
  expect_equal(child$enrichment, Inf)
  expect_equal(parent$enrichment, (12 / 12) / (8 / 28), tolerance = 1e-10)
  expect_equal(got$p_adj, stats::p.adjust(got$raw_p, "BH"), tolerance = 1e-8)
  expect_true(all(got$p_adjust_method == "BH"))

  # Every row's counts belong to its own term.
  expect_true(all(got$significant <= got$annotated))
  root <- got[got$GO_ID == "GO:0008150", ]
  expect_equal(root$annotated, 40L)
  expect_equal(root$significant, 12L)
  expect_equal(root$raw_p, 1)

  # weight01 tests a term without its significant child's genes, so the parent weakens.
  weighted <- run(algorithm = "weight01", p_adjust = "none")
  expect_equal(weighted$raw_p[weighted$GO_ID == "GO:0006952"], child$raw_p, tolerance = 1e-8)
  expect_gt(weighted$raw_p[weighted$GO_ID == "GO:0006950"], parent$raw_p)
  expect_true(all(weighted$p_adj == weighted$raw_p))
  expect_true(all(weighted$p_adjust_method == "none"))
})

test_that("ideogram windows average the right genes", {
  got <- .ideogram_windows(
    dnds = c(0.5, 2, 1, 0.1),
    chr = c("a", "a", "a", "b"),
    start = c(10, 90, 150, 5),
    window_size = 100
  )
  first <- got[got$Chr == "a" & got$Start == 0, ]
  expect_equal(first$n_genes, 2L)
  expect_equal(first$avg_ln_dNdS, mean(log(c(1.5, 3))))
  expect_equal(first$avg_ln_neg_dNdS, log(1.5))
  expect_equal(first$avg_ln_pos_dNdS, log(3))
  second <- got[got$Chr == "a" & got$Start == 100, ]
  expect_equal(second$n_genes, 1L) # exactly 1 is neither
  expect_equal(second$neg_genes, 0L)
  expect_equal(second$pos_genes, 0L)
  expect_true(is.na(second$avg_ln_neg_dNdS))
  expect_true(is.na(second$avg_ln_pos_dNdS))
  expect_equal(got$End[got$Chr == "b"], 100)
})
