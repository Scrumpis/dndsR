# Foreground resolution, pruning, and identifiability. Expected answers are worked out by
# hand from the topology below; none comes from dndsR output.
#   ((((spA,spB),(spC,spD)),(spE,spF)),(spG,spH));
skip_if_not_installed("ape")

species <- "((((spA,spB),(spC,spD)),(spE,spF)),(spG,spH));"
all_taxa <- paste0("sp", LETTERS[1:8])

hypothesis <- function(rows, tree = species) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  writeLines(tree, file.path(dir, "tree.nwk"))
  writeLines(c("scope\ttaxa", paste(names(rows), rows, sep = "\t")), file.path(dir, "fg.tsv"))
  t <- .phylo_read_species_tree(file.path(dir, "tree.nwk"))
  .phylo_read_foreground(file.path(dir, "fg.tsv"), t)
}
fg_sets <- function(map) sort(map$branches$taxa_side[map$branches$foreground])

test_that("a clade hypothesis labels the clade's terminal branches and its stem", {
  m <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), all_taxa)
  expect_true(m$ok)
  expect_identical(fg_sets(m), c("spA", "spA,spB", "spB"))
  expect_identical(m$newick, "(((spA{FG},spB{FG})Node3{FG},(spC,spD)Node4)Node2,(spE,spF)Node5,(spG,spH)Node6);")
})

test_that("stem and terminal hypotheses label only their branches", {
  expect_identical(fg_sets(.phylo_map_foreground(hypothesis(c(stem = "spA,spB")), all_taxa)), "spA,spB")
  expect_identical(fg_sets(.phylo_map_foreground(hypothesis(c(terminal = "spA,spE")), all_taxa)), c("spA", "spE"))
  two <- .phylo_map_foreground(hypothesis(c(stem = "spA,spB", terminal = "spE")), all_taxa)
  expect_identical(fg_sets(two), c("spA,spB", "spE"))
})

test_that("a missing clade member makes a stem hypothesis unidentifiable, but not a clade hypothesis", {
  present <- setdiff(all_taxa, "spB")
  stem <- .phylo_map_foreground(hypothesis(c(stem = "spA,spB")), present)
  expect_false(stem$ok)
  expect_identical(stem$reason, "foreground_unidentifiable")
  expect_match(stem$detail, "mrca(spA,spB)+spA", fixed = TRUE)
  clade <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), present)
  expect_true(clade$ok)
  expect_identical(fg_sets(clade), "spA")
  expect_identical(clade$branches$species_tree_branches[clade$branches$foreground], "mrca(spA,spB)+spA")
  expect_identical(clade$info$fg_absent_branches, 1L)
})

test_that("a missing sister lineage: the stem is unidentifiable; a clade is excluded by policy", {
  present <- setdiff(all_taxa, c("spC", "spD"))
  # Stem scope: the only requested branch joins background, so nothing can be separated.
  stem <- .phylo_map_foreground(hypothesis(c(stem = "spA,spB")), present)
  expect_identical(stem$reason, "foreground_unidentifiable")
  # Clade scope: spA and spB stay separable, but testing them alone would narrow the hypothesis.
  clade <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), present)
  expect_identical(clade$reason, "foreground_partly_unidentifiable")
  expect_identical(.phylo_reason_basis(c(stem$reason, clade$reason)), c("mathematical", "policy"))
})

test_that("a terminal branch whose sister is missing is unidentifiable", {
  m <- .phylo_map_foreground(hypothesis(c(terminal = "spA")), setdiff(all_taxa, "spB"))
  expect_false(m$ok)
  expect_identical(m$reason, "foreground_unidentifiable")
  expect_true(.phylo_map_foreground(hypothesis(c(terminal = "spA")), setdiff(all_taxa, "spH"))$ok)
})

test_that("the stem of a clade beside the root cannot be told apart from the other root branch", {
  m <- .phylo_map_foreground(hypothesis(c(stem = "spG,spH")), all_taxa)
  expect_false(m$ok)
  expect_identical(m$reason, "foreground_unidentifiable")
})

test_that("a hypothesis row with no taxa present, or a terminal row lacking a taxon, is incomplete", {
  m <- .phylo_map_foreground(hypothesis(c(terminal = "spA,spE")), setdiff(all_taxa, "spE"))
  expect_identical(m$reason, "foreground_incomplete")
  expect_match(m$detail, "lacks: spE")
  two <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB", terminal = "spH")), setdiff(all_taxa, "spH"))
  expect_identical(two$reason, "foreground_incomplete")
  expect_identical(.phylo_reason_basis("foreground_incomplete"), "policy")
  expect_true(.phylo_map_foreground(hypothesis(c(clade = "spA,spB", terminal = "spH")), all_taxa)$ok)
})

test_that("missing taxa elsewhere leave the hypothesis intact", {
  m <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), setdiff(all_taxa, "spH"))
  expect_true(m$ok)
  expect_identical(fg_sets(m), c("spA", "spA,spB", "spB"))
  expect_false(grepl("spH", m$newick))
})

test_that("a foreground with no taxa present is reported as absent", {
  m <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), setdiff(all_taxa, c("spA", "spB")))
  expect_false(m$ok)
  expect_identical(m$reason, "foreground_absent")
})

test_that("invalid hypotheses stop with a clear message", {
  expect_error(hypothesis(c(stem = "spA,spC")), "do not form a clade.*spB, spD")
  expect_error(hypothesis(c(stem = "spA")), "at least two taxa")
  expect_error(hypothesis(c(clade = "spA,spZ")), "not in the species tree: spZ")
  expect_error(hypothesis(c(crown = "spA,spB")), "Unknown foreground scope")
  unrooted <- "(((spA,spB),(spC,spD)),(spE,spF),(spG,spH));"
  expect_error(hypothesis(c(stem = "spA,spB"), unrooted), "rooted species tree")
  expect_true(.phylo_map_foreground(hypothesis(c(terminal = "spA"), unrooted), all_taxa)$ok)
})

test_that("taxon names HyPhy would rename are rejected", {
  dir <- withr::local_tempdir()
  writeLines("((A.thaliana,B),(C,D));", file.path(dir, "t.nwk"))
  expect_error(.phylo_read_species_tree(file.path(dir, "t.nwk")), "letters, digits, and underscores")
})

test_that("split keys do not depend on rooting or order", {
  p <- c("spC", "spA", "spB", "spD")
  expect_identical(.phylo_split_key(c("spA", "spB"), p), .phylo_split_key(c("spC", "spD"), p))
  expect_identical(.phylo_split_key(c("spB", "spA"), p), "spC,spD")
  expect_true(is.na(.phylo_split_key(p, p)))
  expect_true(is.na(.phylo_split_key("spZ", p)))
})

test_that("a concordant gene tree gives the same labels as the species tree", {
  gt <- ape::read.tree(text = "((spG,spH),((spE,spF),((spB,spA),(spD,spC))));")
  m <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), all_taxa, gt)
  expect_true(m$ok)
  expect_identical(fg_sets(m), c("spA", "spA,spB", "spB"))
  expect_true(all(m$branches$topology_source == "gene_tree"))
})

test_that("a gene tree without the foreground clade cannot carry the hypothesis, and is not relabelled", {
  gt <- ape::read.tree(text = "((((spA,spC),(spB,spD)),(spE,spF)),(spG,spH));")
  m <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), all_taxa, gt)
  expect_false(m$ok)
  expect_identical(m$reason, "foreground_not_in_gene_tree")
  expect_identical(.phylo_reason_basis(m$reason), "mathematical")
  # In report mode the species tree is used, and the discordance is reported.
  r <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), all_taxa, gt, "report")
  expect_true(r$ok)
  expect_true(all(r$branches$topology_source == "species_tree"))
  expect_identical(fg_sets(r), c("spA", "spA,spB", "spB"))
  # Species splits {A,B} and {C,D} are missing from the gene tree; {A,B} is foreground.
  expect_identical(r$info$gt_missing_species_splits, 2L)
  expect_true(r$info$gt_discordant_near_foreground)
})

test_that("discordance away from the foreground is reported as such", {
  # Swapping spE and spG changes splits {E,F} and {G,H}, far from the stem of spA and spB.
  gt <- ape::read.tree(text = "((((spA,spB),(spC,spD)),(spG,spF)),(spE,spH));")
  m <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), all_taxa, gt)
  expect_true(m$ok)
  expect_identical(m$info$gt_missing_species_splits, 2L)
  expect_identical(m$info$gt_extra_splits, 2L)
  expect_false(m$info$gt_discordant_near_foreground)
  expect_identical(fg_sets(m), c("spA", "spA,spB", "spB"))
})

test_that("inside a clade-scope foreground, discordant gene-tree branches are foreground", {
  gt <- ape::read.tree(text = "((((spA,spC),(spB,spD)),(spE,spF)),(spG,spH));")
  m <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB,spC,spD")), all_taxa, gt)
  expect_true(m$ok)
  # The clade's stem splits 4 taxa from 4, so it is shown by the side without spA.
  expect_identical(fg_sets(m), sort(c("spA", "spB", "spC", "spD", "spA,spC", "spB,spD", "spE,spF,spG,spH")))
  # Under a stem scope, only the stem is labelled.
  s <- .phylo_map_foreground(hypothesis(c(stem = "spA,spB,spC,spD")), all_taxa, gt)
  expect_identical(fg_sets(s), "spE,spF,spG,spH")
})

test_that("gene tree tips may be taxa, sequence names, gene IDs, or OrthoFinder names, and must match", {
  dir <- withr::local_tempdir()
  seqs <- data.frame(seq_name = paste0(c("spA", "spB", "spC", "spD"), "|g", 1:4), taxon = c("spA", "spB", "spC", "spD"),
                     gene_id = paste0("g", 1:4))
  read <- function(text, s = seqs, aliases = NULL) {
    writeLines(text, file.path(dir, "g.nwk"))
    .phylo_read_gene_tree(file.path(dir, "g.nwk"), s, aliases)$tip.label
  }
  expect_identical(read("((spA|g1,spB|g2),(spC|g3,spD|g4));"), c("spA", "spB", "spC", "spD"))
  expect_identical(read("((g1,g2),(g3,g4));"), c("spA", "spB", "spC", "spD"))
  expect_identical(read("((spA_g1,spB_g2),(spC_g3,spD_g4));"), c("spA", "spB", "spC", "spD"))
  expect_identical(read("((A_prot_g1,spB_g2),(spC,spD));", aliases = c(A_prot = "spA")), c("spA", "spB", "spC", "spD"))
  expect_error(read("((spA,spB),(spC,spD));", seqs[1:3, ]), "gene_tree_taxa_mismatch")
  expect_error(read("((spA,spB),(spC,spZ));"), "gene_tree_taxa_mismatch")
})

test_that("every reason has a basis", {
  reasons <- c("foreground_unidentifiable", "foreground_not_in_gene_tree", "foreground_absent", "foreground_partly_unidentifiable",
               "foreground_incomplete", "multi_copy", "too_few_taxa", "cds_invalid", "internal_stop", "gene_not_found",
               "hyphy_error", "backtranslation_mismatch", "prepared_tree_mismatch", "label_mismatch")
  expect_false(anyNA(.phylo_reason_basis(reasons)))
  expect_identical(.phylo_reason_basis(c("multi_copy", "hyphy_error", "foreground_absent")), c("policy", "software", "no_data"))
})
