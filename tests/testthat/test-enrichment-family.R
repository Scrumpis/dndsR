# The original defect: min_pos removed terms before correction, so a term's own outcome
# decided whether it counted as a test. With dN/dS assigned at random, over half of runs on a
# real genome pair reported a significant term. See docs/validation/enrichment.md.
upper_tail <- function(a, term, positive, universe) {
  stats::phyper(a - 1L, term, universe - term, positive, lower.tail = FALSE)
}

# Written out here rather than calling p.adjust, to be independent of the code under test.
bh_by_hand <- function(p) {
  n <- length(p)
  ranks <- rank(p, ties.method = "first")
  adjusted <- p * n / ranks
  for (i in seq_len(n)) adjusted[i] <- min(1, adjusted[p >= p[i]])
  adjusted
}

entries <- function(type) {
  e <- utils::read.delim(system.file("extdata", "interpro_entry.list.tsv", package = "dndsR"))
  e$ENTRY_AC[e$ENTRY_TYPE == type]
}

# 10 positive genes among 183 annotated. Three terms contain positive genes, and thirty do not.
design <- function() {
  ids <- entries("Domain")[1:33]
  term <- c(rep(ids[1], 2), rep(ids[2], 3), rep(ids[3], 5), # positive genes
            rep(ids[1], 2), rep(ids[2], 1), rep(ids[3], 20), rep(ids[4:33], each = 5), rep(NA, 17))
  list(ids = ids, term = term, dnds = c(rep(2, 10), rep(0.1, 190)),
       p = c(upper_tail(2L, 4L, 10L, 183L), upper_tail(3L, 4L, 10L, 183L), upper_tail(5L, 25L, 10L, 183L), rep(1, 30)))
}

write_table <- function(path, dnds, term) {
  ids <- sprintf("g%03d", seq_along(dnds))
  utils::write.table(
    data.frame(query_id = ids, subject_id = paste0("s", ids), dN = dnds * 0.1, dS = 0.1, dNdS = dnds, q_ipr = term, s_ipr = term),
    path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  path
}

ipr <- function(file, out = "t_q_IPR_Domain_enrichment.tsv", ...) {
  suppressMessages(ipr_enrichment(
    dnds_annot_file = file, sides = "query", make_plots = FALSE, threads = 1, auto_detect_release = FALSE, max_prop = 1,
    entries_source = "local", entries_path = system.file("extdata", "interpro_entry.list.tsv", package = "dndsR"), ...
  ))
  utils::read.delim(file.path(dirname(file), out), quote = "")
}

term <- function(file, ...) {
  suppressMessages(term_enrichment(dnds_annot_file = file, terms = "ipr", exclude_terms = NULL, sides = "query",
                                   make_plots = FALSE, threads = 1, verbose = FALSE, ...))
  utils::read.delim(file.path(dirname(file), "t_q_IPR_enrichment.tsv"), quote = "")
}

for (fn in c("ipr", "term")) local({
  run <- function(...) {
    out <- get(fn, mode = "function")(...)
    if (fn == "term") out$IPR <- out$term
    out
  }
  extra <- if (fn == "ipr") list(tree_source = "none") else list()
  call_run <- function(file, ...) do.call(run, c(list(file), extra, list(...)))

  test_that(paste(fn, "corrects across every term of adequate size, whatever its outcome"), {
    d <- design()
    file <- write_table(file.path(withr::local_tempdir(), "t_dnds_annot.tsv"), d$dnds, d$term)
    got <- call_run(file, fdr_method = "BH", min_total = 2, min_pos = 2)
    expect_setequal(got$IPR, d$ids[1:3]) # min_pos still decides what is reported
    expect_equal(unique(got$n_tested), 33L)
    expect_equal(unique(got$p_adjust_method), "BH")
    expected <- stats::setNames(bh_by_hand(d$p), d$ids)
    expect_equal(got$p_adj, unname(expected[got$IPR]), tolerance = 1e-10)
    expect_equal(got$p_adj, stats::p.adjust(d$p, "BH")[match(got$IPR, d$ids)], tolerance = 1e-10)
    # The filtered family of the old code would have given a far smaller value.
    expect_gt(got$p_adj[got$IPR == d$ids[1]], stats::p.adjust(d$p[1:3], "BH")[1] * 2)
  })

  test_that(paste(fn, "reporting filters do not change p-values"), {
    d <- design()
    file <- write_table(file.path(withr::local_tempdir(), "t_dnds_annot.tsv"), d$dnds, d$term)
    one <- call_run(file, min_total = 2, min_pos = 1)
    three <- call_run(file, min_total = 2, min_pos = 3)
    shared <- match(three$IPR, one$IPR)
    expect_lt(nrow(three), nrow(one))
    expect_equal(three$p_adj, one$p_adj[shared], tolerance = 1e-12)
    expect_equal(three$n_tested, one$n_tested[shared])
  })

  test_that(paste(fn, "min_total sets the family, without using the outcome"), {
    d <- design()
    file <- write_table(file.path(withr::local_tempdir(), "t_dnds_annot.tsv"), d$dnds, d$term)
    got <- call_run(file, min_total = 5, min_pos = 2)
    # Terms 1 and 2 have 4 genes, and go; the 31 larger ones stay, including 30 with no positives.
    expect_equal(got$IPR, d$ids[3])
    expect_equal(got$n_tested, 31L)
    expect_equal(got$p_adj, bh_by_hand(c(d$p[3], rep(1, 30)))[1], tolerance = 1e-10)
  })

  test_that(paste(fn, "results do not depend on row order"), {
    d <- design()
    dir <- withr::local_tempdir()
    a <- call_run(write_table(file.path(dir, "t_dnds_annot.tsv"), d$dnds, d$term))
    withr::with_seed(1, o <- sample(length(d$dnds)))
    dir2 <- withr::local_tempdir()
    b <- call_run(write_table(file.path(dir2, "t_dnds_annot.tsv"), d$dnds[o], d$term[o]))
    expect_equal(b[, c("IPR", "p_value", "p_adj")], a[, c("IPR", "p_value", "p_adj")], tolerance = 1e-12)
  })
})

test_that("the other corrections receive the whole family", {
  d <- design()
  file <- write_table(file.path(withr::local_tempdir(), "t_dnds_annot.tsv"), d$dnds, d$term)
  by <- ipr(file, tree_source = "none", fdr_method = "BY")
  expect_equal(by$p_adj, stats::p.adjust(d$p, "BY")[match(by$IPR, d$ids)], tolerance = 1e-10)
  none <- ipr(file, tree_source = "none", fdr_method = "none")
  expect_equal(none$p_adj, none$p_value)
  expect_equal(unique(none$p_adjust_method), "none")
})

test_that("global correction spans InterPro types, before any reporting filter", {
  d <- design()
  fam <- entries("Family")[1:5]
  extra <- rep(fam, each = 4)
  file <- write_table(file.path(withr::local_tempdir(), "t_dnds_annot.tsv"), c(d$dnds, rep(0.1, 20)), c(d$term, extra))
  got <- ipr(file, tree_source = "none", fdr_method = "BH", adjust_scope = "global", min_total = 2, min_pos = 2)
  # Domain universe unchanged by the Family genes, which carry no Domain term. Family has no positives,
  # so it has no rows to report, but its five tests still count.
  expect_false(file.exists(file.path(dirname(file), "t_q_IPR_Family_enrichment.tsv")))
  expect_equal(unique(got$n_tested), 38L)
  expected <- bh_by_hand(c(d$p, rep(1, 5)))
  expect_equal(got$p_adj, expected[match(got$IPR, d$ids)], tolerance = 1e-10)
})

test_that("pooled analysis corrects across every eligible term", {
  d <- design()
  file <- write_table(file.path(withr::local_tempdir(), "t_dnds_annot.tsv"), d$dnds, d$term)
  got <- ipr(file, out = "t_q_IPR_enrichment.tsv", tree_source = "none", stratify_by_type = FALSE)
  expect_equal(unique(got$n_tested), 33L)
  expect_equal(got$p_adj, stats::p.adjust(d$p, "BH")[match(got$IPR, d$ids)], tolerance = 1e-10)
})

test_that("max_prop limits the family by size alone", {
  d <- design()
  file <- write_table(file.path(withr::local_tempdir(), "t_dnds_annot.tsv"), d$dnds, d$term)
  suppressMessages(ipr_enrichment(
    dnds_annot_file = file, sides = "query", make_plots = FALSE, threads = 1, auto_detect_release = FALSE, max_prop = 0.1,
    entries_source = "local", entries_path = system.file("extdata", "interpro_entry.list.tsv", package = "dndsR"),
    tree_source = "none"
  ))
  got <- utils::read.delim(file.path(dirname(file), "t_q_IPR_Domain_enrichment.tsv"), quote = "")
  # Term 3 has 25 of 183 genes (13.7%), so it is not tested at all.
  expect_false(d$ids[3] %in% got$IPR)
  expect_equal(unique(got$n_tested), 32L)
})

test_that("the parent-child method keeps eligible terms in the family", {
  ids <- entries("Domain")[1:13]
  parent <- ids[1]; child <- ids[2]; lonely <- ids[3]; others <- ids[4:13]
  # As in test-enrichment-analytic.R, with ten extra terms that contain no positive genes.
  term <- c(rep(paste(parent, child, sep = ";"), 10), rep(parent, 10), rep(lonely, 5), rep(others, each = 5), rep(NA, 25))
  dnds <- c(rep(2, 10), rep(0.1, 10), rep(2, 5), rep(0.1, 75))
  file <- write_table(file.path(withr::local_tempdir(), "t_dnds_annot.tsv"), dnds, term)
  got <- ipr(file, method = "parent_child", fdr_method = "BH", drop_rows_without_term = FALSE,
             tree_source = "local", term_trees = data.frame(parent = parent, child = child), min_total = 2, min_pos = 2)
  expect_setequal(got$IPR, c(child, lonely)) # the parent has no positive genes left
  expect_equal(unique(got$n_tested), 13L)
  p <- c(upper_tail(10L, 10L, 15L, 100L), upper_tail(5L, 5L, 15L, 100L), rep(1, 11))
  expect_equal(got$p_adj[match(c(child, lonely), got$IPR)], bh_by_hand(p)[1:2], tolerance = 1e-10)
})

test_that("go_enrichment gives the same results with or without definitions", {
  skip_if_not_installed("topGO")
  skip_if_not_installed("GO.db")
  dir <- withr::local_tempdir()
  go <- c(rep("GO:0006952", 10), rep("GO:0006950", 10), rep("GO:0008152", 20))
  dnds <- c(rep(2, 12), rep(0.1, 28))
  ids <- sprintf("g%03d", seq_along(dnds))
  file <- file.path(dir, "t_dnds_annot.tsv")
  utils::write.table(data.frame(query_id = ids, subject_id = paste0("s", ids), dN = dnds * 0.1, dS = 0.1, dNdS = dnds, q_go = go, s_go = go),
                     file, sep = "\t", quote = FALSE, row.names = FALSE)
  run <- function(definition) {
    suppressMessages(go_enrichment(dnds_annot_file = file, sides = "query", ontologies = "BP", make_plots = FALSE,
                                   node_min = 1, threads = 1, include_definition = definition))
    utils::read.delim(file.path(dir, "t_q_GO_BP_topGO_weight01_fisher.tsv"), quote = "")
  }
  with_definition <- run(TRUE)
  without <- run(FALSE)
  expect_true("term_def" %in% names(with_definition))
  expect_true(any(nzchar(with_definition$term_def)))
  expect_false("term_def" %in% names(without))
  expect_gt(nrow(without), 0L)
  expect_identical(without, with_definition[, setdiff(names(with_definition), "term_def")])
})

test_that("ipr_enrichment corrects across entry types by default, and per type on request", {
  d <- design()
  fam <- entries("Family")[1:2]
  # Two Family terms among 20 further genes, of which 3 are positive.
  term <- c(d$term, rep(fam[1], 10), rep(fam[2], 10))
  dnds <- c(d$dnds, rep(2, 3), rep(0.1, 17))
  file <- write_table(file.path(withr::local_tempdir(), "t_dnds_annot.tsv"), dnds, term)
  domain_p <- c(upper_tail(2L, 4L, 10L, 183L), upper_tail(3L, 4L, 10L, 183L), upper_tail(5L, 25L, 10L, 183L), rep(1, 30))
  family_p <- c(upper_tail(3L, 10L, 3L, 20L), 1)
  read_both <- function(...) {
    domain <- ipr(file, tree_source = "none", min_total = 2, min_pos = 1, ...)
    family <- utils::read.delim(file.path(dirname(file), "t_q_IPR_Family_enrichment.tsv"), quote = "")
    list(domain = domain, family = family)
  }
  default <- read_both()
  expect_equal(unique(c(default$domain$n_tested, default$family$n_tested)), 35L)
  expect_equal(unique(c(default$domain$adjust_scope, default$family$adjust_scope)), "global")
  all_adj <- bh_by_hand(c(domain_p, family_p))
  expect_equal(default$domain$p_adj, all_adj[match(default$domain$IPR, d$ids)], tolerance = 1e-10)
  expect_equal(default$family$p_adj[default$family$IPR == fam[1]], all_adj[34], tolerance = 1e-10)

  per_type <- read_both(adjust_scope = "per_type")
  expect_equal(unique(per_type$domain$n_tested), 33L)
  expect_equal(unique(per_type$family$n_tested), 2L)
  expect_equal(unique(per_type$domain$adjust_scope), "per_type")
  expect_equal(per_type$domain$p_adj, bh_by_hand(domain_p)[match(per_type$domain$IPR, d$ids)], tolerance = 1e-10)
  expect_equal(per_type$family$p_adj[per_type$family$IPR == fam[1]], bh_by_hand(family_p)[1], tolerance = 1e-10)
})

test_that("parent-child correction spans entry types by default, and per type on request", {
  dom <- entries("Domain")[1:3]; parent <- dom[1]; child <- dom[2]; lonely <- dom[3]
  fam <- entries("Family")[1:2]
  term <- c(rep(paste(parent, child, sep = ";"), 10), rep(parent, 10), rep(lonely, 5), rep(NA, 75),
            rep(paste(fam[1], fam[2], sep = ";"), 10), rep(fam[1], 10))
  dnds <- c(rep(2, 10), rep(0.1, 10), rep(2, 5), rep(0.1, 75), rep(2, 2), rep(0.1, 18))
  file <- write_table(file.path(withr::local_tempdir(), "t_dnds_annot.tsv"), dnds, term)
  tree <- data.frame(parent = c(parent, fam[1]), child = c(child, fam[2]))
  run <- function(...) ipr(file, method = "parent_child", fdr_method = "BH", drop_rows_without_term = FALSE,
                           tree_source = "local", term_trees = tree, min_total = 2, min_pos = 2, ...)
  default <- run()
  family_global <- utils::read.delim(file.path(dirname(file), "t_q_IPR_Family_enrichment.tsv"), quote = "")
  per_type <- run(adjust_scope = "per_type")
  family_per_type <- utils::read.delim(file.path(dirname(file), "t_q_IPR_Family_enrichment.tsv"), quote = "")
  expect_equal(unique(default$adjust_scope), "global")
  expect_equal(unique(per_type$adjust_scope), "per_type")
  expect_equal(unique(per_type$n_tested), 3L) # the Domain terms
  expect_equal(unique(default$n_tested), 5L) # with the Family terms
  expect_equal(unique(family_per_type$n_tested), 2L)
  expect_equal(unique(family_global$n_tested), 5L)
  # Same tests either way: only the correction differs, and it is BH on all five.
  expect_equal(default$p_value, per_type$p_value[match(default$IPR, per_type$IPR)])
  expect_equal(family_global$p_value, family_per_type$p_value[match(family_global$IPR, family_per_type$IPR)])
  # 17 positive genes among 120 rows. The Domain parent has none left once its child claims them.
  child_p <- upper_tail(10L, 10L, 17L, 120L); lonely_p <- upper_tail(5L, 5L, 17L, 120L)
  expect_equal(default$p_value[match(c(child, lonely), default$IPR)], c(child_p, lonely_p), tolerance = 1e-10)
  raw <- c(child_p, lonely_p, 1, upper_tail(2L, 10L, 17L, 120L), upper_tail(2L, 20L, 17L, 120L))
  expect_equal(default$p_adj[match(c(child, lonely), default$IPR)], bh_by_hand(raw)[1:2], tolerance = 1e-10)
})

go_table <- function(dir, dnds, go) {
  ids <- sprintf("g%03d", seq_along(dnds))
  file <- file.path(dir, "t_dnds_annot.tsv")
  utils::write.table(data.frame(query_id = ids, subject_id = paste0("s", ids), dN = dnds * 0.1, dS = 0.1, dNdS = dnds, q_go = go, s_go = go),
                     file, sep = "\t", quote = FALSE, row.names = FALSE)
  file
}

go_tables <- function(file, ontologies = c("BP", "MF", "CC"), algorithm = "classic", ...) {
  suppressMessages(go_enrichment(dnds_annot_file = file, sides = "query", ontologies = ontologies, make_plots = FALSE,
                                 node_min = 1, threads = 1, algorithm = algorithm, ...))
  rows <- lapply(ontologies, function(o) {
    path <- file.path(dirname(file), sprintf("t_q_GO_%s_topGO_%s_fisher.tsv", o, algorithm))
    if (file.exists(path)) utils::read.delim(path, quote = "") else NULL
  })
  do.call(rbind, rows)
}

go_design <- function() {
  # Defense response (BP), DNA binding (MF), and nucleus (CC), on separate genes.
  list(go = c(rep("GO:0006952", 20), rep("GO:0008152", 20), rep("GO:0003677", 20), rep("GO:0003824", 20),
              rep("GO:0005634", 20), rep("GO:0005737", 20)),
       dnds = c(rep(2, 8), rep(0.1, 32), rep(2, 6), rep(0.1, 34), rep(2, 4), rep(0.1, 36)))
}

test_that("go_enrichment corrects across all ontologies by default", {
  skip_if_not_installed("topGO")
  skip_if_not_installed("GO.db")
  d <- go_design()
  file <- go_table(withr::local_tempdir(), d$dnds, d$go)
  got <- go_tables(file)
  expect_equal(unique(got$p_adjust_method), "BH")
  expect_equal(unique(got$adjust_scope), "global")
  expect_equal(unique(got$n_tested), nrow(got))
  expect_equal(got$p_adj, bh_by_hand(got$raw_p), tolerance = 1e-10)
  # A raw p-value from outside topGO: 8 of 8 positive BP genes are in a 20-gene term, of 40.
  defense <- got[got$GO_ID == "GO:0006952", ]
  expect_equal(defense$raw_p, upper_tail(8L, 20L, 8L, 40L), tolerance = 1e-8)
  expect_gt(defense$p_adj, defense$raw_p)
})

test_that("go_enrichment can correct per ontology, or not at all", {
  skip_if_not_installed("topGO")
  skip_if_not_installed("GO.db")
  d <- go_design()
  file <- go_table(withr::local_tempdir(), d$dnds, d$go)
  per <- go_tables(file, adjust_scope = "per_ontology")
  expect_equal(unique(per$adjust_scope), "per_ontology")
  for (o in c("BP", "MF", "CC")) {
    one <- per[per$ontology == o, ]
    expect_equal(unique(one$n_tested), nrow(one))
    expect_equal(one$p_adj, bh_by_hand(one$raw_p), tolerance = 1e-10)
  }
  none <- go_tables(file, p_adjust = "none")
  expect_equal(unique(none$p_adjust_method), "none")
  expect_equal(unique(none$adjust_scope), "none")
  expect_equal(none$p_adj, none$raw_p)
  # Asking for one ontology makes it the family.
  single <- go_tables(file, ontologies = "BP")
  expect_equal(single$p_adj, bh_by_hand(single$raw_p), tolerance = 1e-10)
  expect_equal(unique(single$n_tested), nrow(single))
})

test_that("every GO algorithm has p_adj corrected from its own raw p-values", {
  skip_if_not_installed("topGO")
  skip_if_not_installed("GO.db")
  d <- go_design()
  file <- go_table(withr::local_tempdir(), d$dnds, d$go)
  for (algorithm in c("weight01", "elim", "classic", "weight")) {
    got <- go_tables(file, algorithm = algorithm)
    expect_equal(got$p_adj, bh_by_hand(got$raw_p), tolerance = 1e-10, info = algorithm)
    expect_true(all(got$p_adj >= got$raw_p - 1e-12), info = algorithm)
  }
})

test_that("an ontology with no positive genes writes nothing but stays in the GO family", {
  skip_if_not_installed("topGO")
  skip_if_not_installed("GO.db")
  d <- go_design()
  dnds <- d$dnds; dnds[81:120] <- 0.1 # CC now has no positive gene
  dir <- withr::local_tempdir()
  file <- go_table(dir, dnds, d$go)
  got <- go_tables(file)
  expect_false(file.exists(file.path(dir, "t_q_GO_CC_topGO_classic_fisher.tsv")))
  expect_gt(unique(got$n_tested), nrow(got))
  extra <- unique(got$n_tested) - nrow(got)
  expect_equal(got$p_adj, bh_by_hand(c(got$raw_p, rep(1, extra)))[seq_len(nrow(got))], tolerance = 1e-10)
})

test_that("GO correction does not depend on row order, and reports the same term rank", {
  skip_if_not_installed("topGO")
  skip_if_not_installed("GO.db")
  d <- go_design()
  a <- go_tables(go_table(withr::local_tempdir(), d$dnds, d$go))
  withr::with_seed(1, o <- sample(length(d$dnds)))
  b <- go_tables(go_table(withr::local_tempdir(), d$dnds[o], d$go[o]))
  cols <- c("GO_ID", "annotated", "significant", "raw_p", "p_adj", "n_tested")
  a <- a[order(a$GO_ID), cols]; b <- b[order(b$GO_ID), cols]
  rownames(a) <- rownames(b) <- NULL
  expect_equal(b, a, tolerance = 1e-12)
})
