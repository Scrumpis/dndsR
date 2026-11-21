# Lineage-level phylogenetic analyses: parsing, branch identity, statistics, and summaries,
# without external programs. Expected values are read from codeml's own output
# (tests/testthat/fixtures/codeml, from codeml 4.10.10 run by hand on the alignment
# example's OG0001) or computed by hand.
skip_if_not_installed("ape")

fixture <- function(f) test_path("fixtures", "codeml", f)
tips8 <- c("spA", "spB", "spC", "spD", "spE", "spF", "spG", "spH")

test_that("codeml's one-ratio, free-ratio, and two-ratio outputs are read", {
  m0 <- dndsR:::.codeml_parse(fixture("model0.mlc"), tips8)
  expect_equal(m0$lnl, -2595.409595); expect_identical(m0$np, 15L)
  expect_equal(m0$kappa, 2.03177); expect_equal(m0$omega, 0.69368)
  m1 <- dndsR:::.codeml_parse(fixture("model1.mlc"), tips8)
  expect_equal(m1$lnl, -2576.086117); expect_identical(m1$np, 27L)
  expect_equal(nrow(m1$branches), 13L)
  # From the branch table: 11..1 is spA (t 0.102, dN/dS 3.3842, S*dS 3.2).
  a <- m1$branches[vapply(m1$branches$tips, identical, logical(1), "spA"), ]
  expect_equal(c(a$t, a$omega, a$SdS), c(0.102, 3.3842, 3.2))
  # 10..11 is the stem of spA and spB.
  ab <- m1$branches[vapply(m1$branches$tips, function(x) setequal(x, c("spA", "spB")), logical(1)), ]
  expect_equal(ab$omega, 1.2808)
  m2 <- dndsR:::.codeml_parse(fixture("model2.mlc"), tips8)
  expect_equal(m2$omega, c(0.39658, 1.88978)); expect_identical(m2$np, 16L)
})

test_that("codeml trees carry #1 marks and no internal names", {
  tree <- ape::unroot(ape::read.tree(text = "(((spA,spB),(spC,spD)),((spE,spF),(spG,spH)));"))
  marks <- logical(ape::Ntip(tree) + tree$Nnode)
  marks[match(c("spA", "spB"), tree$tip.label)] <- TRUE
  lines <- dndsR:::.paml_tree_lines(tree, marks)
  expect_identical(lines[1], " 8 1")
  expect_match(lines[2], "spA #1,spB #1", fixed = TRUE)
  expect_false(grepl("Node", lines[2]))
})

test_that("likelihood-ratio tests and their corrections", {
  lrt <- dndsR:::.paml_lrt(list(ok = TRUE, lnl = -100, np = 16L), list(ok = TRUE, lnl = -103, np = 15L))
  expect_equal(lrt$lrt, 6); expect_identical(lrt$df, 1L); expect_equal(lrt$p, stats::pchisq(6, 1, lower.tail = FALSE))
  neg <- dndsR:::.paml_lrt(list(ok = TRUE, lnl = -104, np = 16L), list(ok = TRUE, lnl = -103, np = 15L))
  expect_equal(neg$lrt, 0); expect_true(neg$negative)
  r <- data.frame(status = c("tested", "tested", "excluded", "tested"), p_two_vs_m0 = c(0.01, 0.04, NA, 0.03))
  # BH at m = 3 by hand: 0.03, 0.04, 0.04 (as in test-phylo-hyphy.R).
  expect_equal(dndsR:::.paml_adjust(r)$p_adj_two_vs_m0, c(0.03, 0.04, NA, 0.04))
})

test_that("pooled omega sums substitutions per site over genes", {
  # Gene 1: 10 nonsynonymous changes over 600 sites, 5 synonymous over 200; gene 2: 1 and 0.1.
  # (11/1200) / (5.1/400) = 0.719, whereas the mean of per-gene ratios is about 1.6.
  expect_equal(dndsR:::.phylo_pooled_omega(c(10, 1), c(600, 600), c(5, 0.1), c(200, 200)), (11 / 1200) / (5.1 / 400))
  expect_true(is.na(dndsR:::.phylo_pooled_omega(1, 600, 0, 200)))
})

test_that("the foreground file's set column, and RELAX's two sets after pruning", {
  tree <- dndsR:::.phylo_read_species_tree(system.file("extdata", "phylo_example", "species_tree.nwk", package = "dndsR"))
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("scope\ttaxa\tset", "clade\tspA,spB\ttest", "clade\tspC,spD\treference"), f)
  expect_error(dndsR:::.phylo_read_foreground(f, tree), "only phylo_relax")
  test <- dndsR:::.phylo_read_foreground(f, tree, allow_reference = TRUE)
  ref <- dndsR:::.phylo_read_foreground(f, tree, set = "reference", allow_reference = TRUE)
  m <- dndsR:::.phylo_map_relax(test, ref, tips8)
  expect_true(m$ok)
  expect_setequal(m$branches$branch[m$branches$set == "test"], c("spA", "spB", m$branches$branch[m$branches$taxa_side == "spA,spB"]))
  expect_identical(sum(m$branches$set == "reference"), 3L)
  expect_match(m$newick, "spA{TEST}", fixed = TRUE); expect_match(m$newick, "spC{REFERENCE}", fixed = TRUE)
  # Without reference rows, every other branch is the reference.
  m2 <- dndsR:::.phylo_map_relax(test, NULL, tips8)
  expect_identical(sum(m2$branches$set == "unclassified"), 0L)
  # A reference of spH alone is absent from an orthogroup without spH.
  writeLines(c("scope\ttaxa\tset", "clade\tspA,spB\ttest", "terminal\tspH\treference"), f)
  m3 <- dndsR:::.phylo_map_relax(dndsR:::.phylo_read_foreground(f, tree, allow_reference = TRUE),
                                 dndsR:::.phylo_read_foreground(f, tree, set = "reference", allow_reference = TRUE),
                                 setdiff(tips8, "spH"))
  expect_false(m3$ok); expect_identical(m3$reason, "reference_absent")
  expect_identical(dndsR:::.phylo_reason_basis("reference_absent"), "no_data")
  # Overlapping sets are refused.
  writeLines(c("scope\ttaxa\tset", "clade\tspA,spB\ttest", "terminal\tspA\treference"), f)
  m4 <- dndsR:::.phylo_map_relax(dndsR:::.phylo_read_foreground(f, tree, allow_reference = TRUE),
                                 dndsR:::.phylo_read_foreground(f, tree, set = "reference", allow_reference = TRUE), tips8)
  expect_false(m4$ok); expect_identical(m4$reason, "test_reference_overlap")
})

test_that("every branch is mapped and labelled for aBSREL without a foreground", {
  tree <- dndsR:::.phylo_read_species_tree(system.file("extdata", "phylo_example", "species_tree.nwk", package = "dndsR"))
  hyp <- dndsR:::.phylo_all_branches_hypothesis(tree)
  m <- dndsR:::.phylo_map_all_tested(hyp, setdiff(tips8, "spB"))
  expect_true(m$ok)
  expect_true(all(m$branches$foreground))
  # 7 taxa, unrooted: 11 branches, each labelled.
  expect_identical(nrow(m$branches), 11L)
  expect_identical(lengths(regmatches(m$newick, gregexpr("{FG}", m$newick, fixed = TRUE))), 11L)
  # The free-ratio mapping needs no foreground at all.
  expect_true(dndsR:::.phylo_map_foreground(hyp, tips8, require_foreground = FALSE)$ok)
})

test_that("aBSREL's per-lineage correction is one BH family per species-tree branch", {
  b <- data.frame(orthogroup = c("g1", "g2", "g3", "g1", "g2", "g4"),
                  species_tree_branches = c("spA", "spA", "spA", "spB", "spB", "spA"),
                  p_raw = c(0.01, 0.02, 0.5, 0.04, NA, 0.001))
  r <- data.frame(orthogroup = c("g1", "g2", "g3", "g4"), status = c("tested", "tested", "tested", "failed"))
  out <- dndsR:::.absrel_adjust(b, r)
  # spA over g1-g3 (g4 failed): 0.01, 0.02, 0.5 -> 0.03, 0.03, 0.5. spB: g1 only -> 0.04.
  expect_equal(out$p_bh_lineage, c(0.03, 0.03, 0.5, 0.04, NA, NA))
  expect_equal(out$lineage_family_size, c(3L, 3L, 3L, 1L, NA, NA))
})

test_that("branch identifiers follow the species tree, with ape's node numbers", {
  tree <- dndsR:::.phylo_read_species_tree(system.file("extdata", "phylo_example", "species_tree.nwk", package = "dndsR"))
  ids <- dndsR:::.phylo_branch_ids(tree, system.file("extdata", "phylo_example", "species_tree.nwk", package = "dndsR"))
  expect_identical(nrow(ids), 14L) # 8 tips, 7 internal nodes, rooted: 14 branches
  a <- ids[ids$branch_id == "spA", ]
  expect_identical(a$ape_node, match("spA", tree$tip.label))
  expect_identical(a$parent_id, "mrca(spA,spB)")
  expect_identical(sum(is.na(ids$parent_id)), 2L) # the two branches below the root
  expect_equal(a$species_tree_length, 0.10)
})

test_that("matched lineage pairs use only shared, reliable orthogroups", {
  long <- data.frame(orthogroup = c("g1", "g2", "g3", "g1", "g2", "g4"), branch_id = c("a", "a", "a", "b", "b", "b"),
                     NdN = c(2, 2, 2, 1, 1, 9), N = 100, SdS = 1, S = 50, omega = c(0.4, 0.4, 0.4, 0.2, 0.2, 9),
                     reliable = TRUE)
  p <- dndsR:::.phylo_lineage_pairs(long, c("a", "b"))
  expect_identical(p$n_shared, 2L)
  expect_identical(c(p$n_a_higher, p$n_b_higher), c(2L, 0L))
  expect_equal(p$median_log2_ratio, 1)
  expect_equal(p$pooled_omega_b, (2 / 200) / (2 / 100)) # g4 is not shared, so not counted
})

test_that("lineage coverage separates merged from absent branches", {
  # Three tested orthogroups. g1 has spA and spB; in g2, spB is missing, so spA joins the
  # clade's stem (merged); in g3, spA is missing (absent) and spB joins the stem.
  br <- data.frame(orthogroup = c("g1", "g1", "g1", "g2", "g3"),
                   species_tree_branches = c("spA", "spB", "ab", "spA+ab", "spB+ab"),
                   in_species_tree = TRUE, t_free = 0.1, N_free = 600, S_free = 200, dN_free = 0.01, dS_free = 0.05,
                   NdN_free = 6, SdS_free = 10, omega_free = 0.2, flags = "", reliable_free = TRUE)
  res <- data.frame(orthogroup = c("g1", "g2", "g3"), status = "tested", qc_flags = "")
  ids <- data.frame(branch_id = c("spA", "spB", "ab"))
  out <- dndsR:::.phylo_summarize_omega(list(results = res, branches = br), ids, bootstrap = 0, seed = 1)$lineages
  out <- out[match(ids$branch_id, out$branch_id), ]
  expect_identical(out$n_estimated, c(1L, 1L, 1L))
  expect_identical(out$n_merged, c(1L, 1L, 2L))
  expect_identical(out$n_absent, c(1L, 1L, 0L))
  expect_identical(out$n_merged_or_absent, out$n_merged + out$n_absent)
})

test_that("known calibration problems are flagged on tested rows of the outputs that have them", {
  res <- data.frame(orthogroup = c("g1", "g2"), status = c("tested", "excluded"), p_value = c(0.01, NA))
  out <- dndsR:::.calibration_column(res, "relax")
  expect_identical(out$calibration, c("p_value:anticonservative_in_simulation", NA))
  expect_match(names(dndsR:::.calibration_notes("relax", res))[1], "p_value")
  expect_warning(dndsR:::.calibration_warn("relax", res), "Calibration")
  # A branch-model run without the free-ratio model has no p_free_vs_m0, so no flag for it.
  bo <- data.frame(orthogroup = "g1", status = "tested", p_two_vs_m0 = 0.2, p_free_vs_m0 = NA_real_)
  expect_false("calibration" %in% names(dndsR:::.calibration_column(bo, "branch_omega")))
  expect_null(dndsR:::.calibration_notes("branch_omega", bo))
  # Both branch-model tests flagged when both were computed, in a fixed order.
  bo2 <- data.frame(orthogroup = "g1", status = "tested", p_free_vs_m0 = 0.2, p_fg_omega_gt1 = 0.4)
  expect_identical(dndsR:::.calibration_column(bo2, "branch_omega")$calibration,
                   "p_free_vs_m0:not_established_in_simulation;p_fg_omega_gt1:not_established_in_simulation")
  # BUSTED's outputs are unchanged.
  expect_false("calibration" %in% names(dndsR:::.calibration_column(res, "busted")))
})
