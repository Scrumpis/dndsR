# Codon alignment QC. Each case is designed by hand, with expected values worked out below.
taxa <- paste0("sp", LETTERS[1:5])
aln <- function(...) {
  s <- c(...)
  data.frame(seq_name = names(s), taxon = sub("\\|.*", "", names(s)), gene_id = NA_character_,
             sequence = unname(s), stringsAsFactors = FALSE)
}
check <- function(a, ...) .phylo_check_alignment(a, taxa, ...)

valid <- aln(spA = "ATGAAACCCTAA", spB = "ATGAAGCCC---", spC = "ATGAAACCG---", spD = "ATGAAACCC---")

test_that("a valid alignment passes, with its terminal stop masked and recorded", {
  r <- check(valid)
  expect_true(r$ok)
  expect_identical(r$aln$sequence[1], "ATGAAACCC---")
  expect_identical(r$metrics$terminal_stops_masked, 1L)
  expect_identical(r$flags, "terminal_stops_masked:1")
  # Pairwise p-distances over 9 sites: AB 1, AC 1, AD 0, BC 2, BD 1, CD 1 differences.
  expect_equal(r$sequences$mean_p_distance, round(c(2, 4, 4, 2) / 27, 6))
  expect_equal(r$metrics$median_p_distance, round(3 / 27, 6))
  expect_equal(r$metrics$max_p_distance, round(2 / 9, 6))
  expect_equal(r$sequences$gap_fraction, rep(0.25, 4))
  expect_equal(r$metrics$gappy_column_fraction, 0.25)
})

test_that("terminal stops can exclude instead of mask", {
  r <- check(valid, terminal_stops = "exclude")
  expect_false(r$ok)
  expect_identical(r$reason, "terminal_stop")
})

test_that("each hard failure excludes the orthogroup with its reason", {
  reason <- function(a, ...) check(a, ...)$reason
  expect_identical(reason(aln(spA = "ATGTGACCC", spB = "ATGAAACCC", spC = "ATGAAACCC", spD = "ATGAAACCC")), "internal_stop")
  expect_identical(reason(aln(spA = "ATGAAA", "spA|x" = "ATGAAA", spB = "ATGAAA", spC = "ATGAAA")), "multi_copy")
  expect_identical(reason(aln(spA = "ATGAAA", spZ = "ATGAAA", spB = "ATGAAA", spC = "ATGAAA")), "unknown_taxon")
  expect_identical(reason(aln(spA = "ATGAAA", spB = "ATGAA", spC = "ATGAAA", spD = "ATGAAA")), "not_aligned")
  expect_identical(reason(aln(spA = "ATGAAAC", spB = "ATGAAAC", spC = "ATGAAAC", spD = "ATGAAAC")), "length_not_multiple_of_3")
  expect_identical(reason(aln(spA = "ATGAXA", spB = "ATGAAA", spC = "ATGAAA", spD = "ATGAAA")), "invalid_character")
  expect_identical(reason(aln(spA = "AUGAAA", spB = "ATGAAA", spC = "ATGAAA", spD = "ATGAAA")), "invalid_character")
  expect_identical(reason(aln(spA = "A-GAAA", spB = "ATGAAA", spC = "ATGAAA", spD = "ATGAAA")), "split_codon_gap")
  expect_identical(reason(aln(spA = "------", spB = "ATGAAA", spC = "ATGAAA", spD = "ATGAAA")), "empty_sequence")
  expect_identical(reason(aln(spA = "NNN---", spB = "ATGAAA", spC = "ATGAAA", spD = "ATGAAA")), "empty_sequence")
  expect_identical(reason(aln(spA = "ATGAAA", spB = "ATGAAA", spC = "ATGAAA")), "too_few_taxa")
  expect_true(check(aln(spA = "ATGAAA", spB = "ATGAAA", spC = "ATGAAA"), min_taxa = 3L)$ok)
  expect_identical(reason(NULL), "empty_alignment")
})

test_that("ambiguous codons are reported, not changed", {
  a <- aln(spA = "ATGAAACCC", spB = "ATGNNNCCC", spC = "ATGAAACCC", spD = "ATGAAARCC")
  r <- check(a)
  expect_true(r$ok)
  expect_identical(r$aln$sequence, a$sequence)
  # spB: 1 of 3 determined codons is ambiguous; spD: 1 of 3.
  expect_equal(r$sequences$ambiguous_codon_fraction, round(c(0, 1 / 3, 0, 1 / 3), 6))
  expect_setequal(r$flags, c("ambiguous_codons:spB", "ambiguous_codons:spD"))
})

test_that("a divergent sequence and a dissimilar window are flagged, not masked", {
  # Four identical Ala sequences and one Gln sequence that differs at every site:
  # mean p is 0.25 for the four and 1 for spE; the median 0.25 times 3 is below 1.
  ala <- strrep("GCT", 10)
  a <- aln(spA = ala, spB = ala, spC = ala, spD = ala, spE = strrep("CAG", 10))
  r <- check(a)
  expect_true(r$ok)
  expect_identical(r$aln$sequence, a$sequence)
  expect_equal(r$sequences$mean_p_distance, c(0.25, 0.25, 0.25, 0.25, 1))
  expect_setequal(r$flags, c("divergent_sequence:spE", "high_divergence", "dissimilar_window:spE"))
  expect_identical(r$sequences$dissimilar_windows, c(0L, 0L, 0L, 0L, 1L))
})

test_that("FASTA headers give taxon and gene, and sequences are joined and upper-cased", {
  f <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">spA|gene1 a description", "atgaaa", "ccc", "", ">spB", "ATGAAACCC"), f)
  a <- .phylo_read_alignment(f)
  expect_identical(a$taxon, c("spA", "spB"))
  expect_identical(a$gene_id, c("gene1", NA))
  expect_identical(a$sequence, c("ATGAAACCC", "ATGAAACCC"))
  writeLines("not a fasta", f)
  expect_null(.phylo_read_alignment(f))
})
