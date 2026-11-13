# phylo_prepare on the synthetic genome example (inst/extdata/phylo_genomes), whose design,
# coding sequences, and one true alignment were made independently of dndsR. Needs MAFFT,
# not HyPhy.
skip_on_os("windows")
skip_if_not_installed("ape")
mafft <- tryCatch(.mafft_path(), error = function(e) NULL)
skip_if(is.null(mafft), "MAFFT not found")

ex <- withr::local_tempdir(.local_envir = teardown_env())
phylo_example(file.path(ex, "g"))
ex <- file.path(ex, "g")
prep <- function(..., run_name = "prepared", orthogroups = "orthogroup_table.tsv") {
  withr::local_dir(ex)
  suppressMessages(phylo_prepare(orthogroups, "species_tree.nwk", "foreground.tsv", run_name = run_name, threads = test_cores(), ...))
}
design <- utils::read.delim(file.path(ex, "truth", "design.tsv"), colClasses = "character")
base <- prep(genome_file = "genomes.tsv")
aln_md5 <- function(run) tools::md5sum(sort(list.files(file.path(ex, run, "alignments"), full.names = TRUE)))

test_that("every orthogroup is prepared or excluded as designed, with its basis", {
  got <- ifelse(base$status == "prepared", "prepared", paste0(base$status, ":", base$reason))
  expect_identical(got, design$expected_prepare[match(base$orthogroup, design$orthogroup)])
  ex_rows <- base$status == "excluded"
  expect_identical(base$exclusion_basis[ex_rows], design$expected_basis[match(base$orthogroup[ex_rows], design$orthogroup)])
  expect_match(base$detail[base$orthogroup == "OG07"], "spD \\(2\\)")
  expect_match(base$detail[base$orthogroup == "OG08"], "in-frame stop codon")
  expect_match(base$detail[base$orthogroup == "OG13"], "length not a multiple of 3")
})

test_that("the member table keeps every gene, including both paralogs and unused genes", {
  m <- utils::read.delim(file.path(ex, "prepared", "prepared_members.tsv"))
  table_in <- utils::read.delim(file.path(ex, "orthogroup_table.tsv"))
  expect_identical(nrow(m), nrow(table_in))
  expect_identical(m$copies[m$orthogroup == "OG07" & m$taxon == "spD"], c(2L, 2L))
  expect_false(any(m$used[m$orthogroup == "OG07"]))
  expect_true(all(m$used[m$orthogroup %in% base$orthogroup[base$status == "prepared"]]))
})

test_that("alignments hold exactly the designed coding sequences, and the indel alignment is the true one", {
  fa <- function(p) {
    l <- readLines(p)
    stats::setNames(l[c(FALSE, TRUE)], sub("^>", "", l[c(TRUE, FALSE)]))
  }
  truth <- fa(file.path(ex, "truth", "designed_cds.fasta"))
  for (og in base$orthogroup[base$status == "prepared"]) {
    a <- fa(file.path(ex, "prepared", "alignments", paste0(og, ".fasta")))
    ids <- sub("^[^|]*\\|", "", names(a))
    expect_identical(unname(gsub("-", "", a)), unname(substr(truth[ids], 1L, nchar(truth[ids]) - 3L)), label = og)
    expect_true(all(nchar(a) %% 3L == 0L))
  }
  got <- fa(file.path(ex, "prepared", "alignments", "OG10.fasta"))
  names(got) <- sub("\\|.*", "", names(got))
  true10 <- fa(file.path(ex, "truth", "OG10_true_codon_alignment.fasta"))
  expect_identical(unname(substring(got[names(true10)], 4L)), unname(true10))
})

test_that("preparation is deterministic", {
  again <- prep(genome_file = "genomes.tsv", run_name = "again")
  expect_identical(unname(aln_md5("again")), unname(aln_md5("prepared")))
  expect_identical(again$tree_md5, base$tree_md5)
})

test_that("OrthoFinder's Orthogroups.tsv and N0.tsv give the same alignments as the long table", {
  of <- prep(genome_file = "genomes.tsv", run_name = "of", orthogroups = "Orthogroups.tsv")
  expect_identical(of$status, base$status)
  expect_identical(unname(aln_md5("of")), unname(aln_md5("prepared")))
  wide <- utils::read.delim(file.path(ex, "Orthogroups.tsv"), check.names = FALSE, colClasses = "character")
  n0 <- data.frame(HOG = wide$Orthogroup, OG = wide$Orthogroup, "Gene Tree Parent Clade" = "n0", wide[, -1], check.names = FALSE)
  utils::write.table(n0, file.path(ex, "N0.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  h <- prep(genome_file = "genomes.tsv", run_name = "n0", orthogroups = "N0.tsv")
  expect_identical(unname(aln_md5("n0")), unname(aln_md5("prepared")))
  aliases <- stats::setNames(paste0("sp", LETTERS[1:6]), paste0("sp", LETTERS[1:6], "_proteins"))
  expect_identical(attr(.phylo_read_orthogroups(file.path(ex, "N0.tsv"), "auto", aliases), "format"), "orthofinder")
})

test_that("CDS FASTA input gives the same alignments; genes without a QC record become not found", {
  cds <- data.frame(taxon = paste0("sp", LETTERS[1:6]),
                    cds = file.path("prepared", "cds", paste0("sp", LETTERS[1:6]), paste0("sp", LETTERS[1:6], "_CDS.fasta")))
  utils::write.table(cds, file.path(ex, "cds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  c2 <- prep(cds_file = "cds.tsv", run_name = "cds")
  expect_identical(unname(aln_md5("cds")), unname(aln_md5("prepared")))
  expect_identical(c2$reason[c2$orthogroup %in% c("OG08", "OG13")], c("gene_not_found", "gene_not_found"))
})

test_that("terminal stops can exclude instead of being removed", {
  x <- prep(genome_file = "genomes.tsv", run_name = "stops", terminal_stops = "exclude")
  expect_true(all(x$status == "excluded"))
  expect_true(all(x$reason[!x$reason %in% c("multi_copy", "cds_invalid")] == "terminal_stop"))
})

test_that("input errors stop with a clear message", {
  withr::local_dir(ex)
  expect_error(phylo_prepare("orthogroup_table.tsv", "species_tree.nwk", "foreground.tsv"), "exactly one of genome_file or cds_file")
  bad <- utils::read.delim("genomes.tsv"); bad$taxon[1] <- "spZ"
  utils::write.table(bad, "bad_genomes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(phylo_prepare("orthogroup_table.tsv", "species_tree.nwk", "foreground.tsv", genome_file = "bad_genomes.tsv"), "not in the species tree: spZ")
  g <- utils::read.delim("genomes.tsv"); g$orthofinder_name <- NULL
  utils::write.table(g, "plain_genomes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(suppressMessages(phylo_prepare("Orthogroups.tsv", "species_tree.nwk", "foreground.tsv", genome_file = "plain_genomes.tsv", run_name = "x")),
               "OrthoFinder columns not matched to a taxon: spA_proteins")
  expect_error(phylo_prepare("orthogroup_table.tsv", "species_tree.nwk", "foreground.tsv", genome_file = "genomes.tsv", mafft = "/nonexistent/mafft",
                             run_name = "y"), "MAFFT not found at /nonexistent/mafft")
})

fake_mafft <- function(body) {
  p <- tempfile("mafft-", tmpdir = ex)
  writeLines(c("#!/bin/sh", 'if [ "$1" = --version ]; then echo "v7.526 (2024/Apr/26)" >&2; exit 0; fi', 'for a in "$@"; do in=$a; done', body), p)
  Sys.chmod(p, "0755")
  p
}

test_that("a MAFFT failure, or an alignment that alters a protein, fails the orthogroup", {
  expect_warning(crash <- prep(genome_file = "genomes.tsv", run_name = "crash", mafft = fake_mafft("exit 1")), "failed in preparation")
  expect_true(all(crash$reason[crash$status != "excluded"] == "alignment_failed"))
  expect_true(all(crash$status[crash$reason == "alignment_failed"] == "failed"))
  # Replaces the first residue of every protein with W: the codons no longer match.
  expect_warning(edit <- prep(genome_file = "genomes.tsv", run_name = "edit", mafft = fake_mafft("sed '/^>/!s/^./W/' \"$in\"")), "failed in preparation")
  expect_true(all(edit$reason[edit$status == "failed"] == "backtranslation_mismatch"))
  expect_identical(sum(edit$status == "failed"), sum(base$status == "prepared"))
  expect_identical(.phylo_reason_basis("backtranslation_mismatch"), "software")
})

test_that("a labelled tree changed after preparation fails the orthogroup before HyPhy runs", {
  copy <- withr::local_tempdir()
  file.copy(file.path(ex, "prepared"), copy, recursive = TRUE)
  t2 <- file.path(copy, "prepared")
  writeLines("((((spA{FG},spB)Node3,spC)Node2,(spD,spE)Node4),spF);", file.path(t2, "trees", "OG01.nwk"))
  fake <- file.path(ex, "hyphy-fake")
  writeLines(c("#!/bin/sh", 'if [ "$1" = --version ]; then echo "HYPHY 2.5.103(MP)"; exit 0; fi', "exit 1"), fake)
  Sys.chmod(fake, "0755")
  r <- suppressWarnings(suppressMessages(phylo_busted(prepared_dir = t2, output_dir = ex, run_name = "tampered", hyphy = fake)))
  expect_identical(r$reason[r$orthogroup == "OG01"], "prepared_tree_mismatch")
  expect_identical(r$reason[r$orthogroup == "OG02"], "hyphy_error")
  expect_identical(r$stage[r$orthogroup == "OG07"], "prepare")
})
