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

test_that("a missing clade member makes a stem hypothesis untestable, but not a clade hypothesis", {
  present <- setdiff(all_taxa, "spB")
  stem <- .phylo_map_foreground(hypothesis(c(stem = "spA,spB")), present)
  expect_false(stem$ok)
  expect_identical(stem$reason, "foreground_not_identifiable")
  expect_match(stem$detail, "mrca(spA,spB)+spA", fixed = TRUE)
  clade <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), present)
  expect_true(clade$ok)
  expect_identical(fg_sets(clade), "spA")
  expect_identical(clade$branches$species_tree_branches[clade$branches$foreground], "mrca(spA,spB)+spA")
})

test_that("a missing sister lineage makes the stem untestable under either scope", {
  present <- setdiff(all_taxa, c("spC", "spD"))
  for (scope in c("stem", "clade")) {
    m <- .phylo_map_foreground(hypothesis(stats::setNames("spA,spB", scope)), present)
    expect_false(m$ok)
    expect_identical(m$reason, "foreground_not_identifiable")
  }
})

test_that("a terminal branch whose sister is missing is untestable", {
  m <- .phylo_map_foreground(hypothesis(c(terminal = "spA")), setdiff(all_taxa, "spB"))
  expect_false(m$ok)
  expect_identical(m$reason, "foreground_not_identifiable")
  expect_true(.phylo_map_foreground(hypothesis(c(terminal = "spA")), setdiff(all_taxa, "spH"))$ok)
})

test_that("the stem of a clade beside the root cannot be told apart from the other root branch", {
  m <- .phylo_map_foreground(hypothesis(c(stem = "spG,spH")), all_taxa)
  expect_false(m$ok)
  expect_identical(m$reason, "foreground_not_identifiable")
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

test_that("a gene tree without the foreground clade is discordant, not relabelled", {
  gt <- ape::read.tree(text = "((((spA,spC),(spB,spD)),(spE,spF)),(spG,spH));")
  m <- .phylo_map_foreground(hypothesis(c(clade = "spA,spB")), all_taxa, gt)
  expect_false(m$ok)
  expect_identical(m$reason, "foreground_discordant")
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

test_that("gene tree tips may be sequence names, and must match the alignment", {
  dir <- withr::local_tempdir()
  writeLines("((spA|g1,spB|g2),(spC|g3,spD|g4));", file.path(dir, "g.nwk"))
  seqs <- c("spA|g1" = "spA", "spB|g2" = "spB", "spC|g3" = "spC", "spD|g4" = "spD")
  expect_setequal(.phylo_read_gene_tree(file.path(dir, "g.nwk"), seqs)$tip.label, c("spA", "spB", "spC", "spD"))
  expect_error(.phylo_read_gene_tree(file.path(dir, "g.nwk"), seqs[1:3]), "gene_tree_taxa_mismatch")
})
