test_that("a shared gene ID gets annotations from its own genome", {
  dir <- withr::local_tempdir()
  gff <- function(path, chrom, start, end, ipr) {
    writeLines(paste(chrom, "x", "mRNA", start, end, ".", "+", ".", paste0("ID=g1;Dbxref=InterPro:", ipr), sep = "\t"), path)
  }
  gff(file.path(dir, "q.gff"), "chr1", 1, 100, "IPR000001")
  gff(file.path(dir, "s.gff"), "chr9", 500, 900, "IPR000002")
  dnds <- file.path(dir, "d.tsv")
  writeLines(c("query_id\tsubject_id\tdN\tdS\tdNdS", "g1\tg1\t0.1\t0.2\t0.5"), dnds)
  got <- suppressMessages(append_annotations(dnds_file = dnds, query_gff = file.path(dir, "q.gff"), subject_gff = file.path(dir, "s.gff")))
  expect_equal(got$q_ipr, "IPR000001")
  expect_equal(got$s_ipr, "IPR000002")
  expect_equal(got$q_gff_seqname, "chr1")
  expect_equal(got$s_gff_seqname, "chr9")
  expect_equal(got$s_gff_start, 500)
})

test_that("feature names are not mistaken for functional terms", {
  dir <- withr::local_tempdir()
  q <- file.path(dir, "q.gff")
  writeLines(paste("chr1", "x", "mRNA", 1, 100, ".", "+", ".", "ID=AK12345.1;Note=K54321 and xPF00001", sep = "\t"), q)
  dnds <- file.path(dir, "d.tsv")
  writeLines(c("query_id\tsubject_id\tdN\tdS\tdNdS", "AK12345.1\tAK12345.1\t0.1\t0.2\t0.5"), dnds)
  got <- suppressMessages(append_annotations(dnds_file = dnds, query_gff = q, subject_gff = q))
  expect_equal(got$q_kegg, "K54321")
  expect_false("q_pfam" %in% names(got) && !is.na(got$q_pfam))
})

test_that("headerless BED files keep their first region", {
  bed <- withr::local_tempfile(fileext = ".bed")
  writeLines(c("chr1 100 200 R1", "chr2 100 200 R2"), bed)
  got <- .read_regions_bed(bed)
  expect_equal(nrow(got), 2L)
  expect_equal(got$seqname, c("chr1", "chr2"))
  expect_equal(got$start, c(101, 101)) # BED start is 0-based
  expect_equal(got$end, c(200, 200))

  writeLines("chr1 100 200", bed)
  expect_equal(nrow(.read_regions_bed(bed)), 1L)

  writeLines(c("seqname start end region_name", "chr1 100 200 R1"), bed)
  expect_equal(nrow(.read_regions_bed(bed)), 1L)
})

test_that("filter expressions fail loudly and NA rows are removed", {
  d <- data.frame(dNdS = c(0.5, 0.7, 20), x = c(1, NA, 1))
  got <- .filter_dnds(d, filter_expr = "x == 1")
  expect_equal(got$dNdS, 0.5)
  expect_false(anyNA(got$dNdS))
  expect_error(.filter_dnds(d, filter_expr = "no_such_column == 1"), "could not be evaluated")
  expect_error(.filter_dnds(d, filter_expr = "dNdS"), "TRUE or FALSE")
  expect_equal(nrow(.filter_dnds(d, filter_expr = "TRUE")), 2L)
})

test_that("non-finite ratios are excluded by the standard filter", {
  got <- .filter_dnds(data.frame(dNdS = c(0.5, Inf, NA, 12, 9.9)))
  expect_equal(got$dNdS, c(0.5, 9.9))
})

test_that("significant descendants claim genes through intermediate terms", {
  # grandparent -> parent -> child; only child is significant.
  claims <- .propagate_claims(
    order_terms = c("child", "parent", "grandparent", "lonely"),
    sig_terms = "child",
    term2rows = list(child = 1:10, parent = 1:20, grandparent = 1:30, lonely = 40:45),
    parents_map = list(child = "parent", parent = "grandparent"),
    nodes = c("child", "parent", "grandparent", "lonely")
  )
  expect_equal(claims$child, integer(0))
  expect_equal(claims$parent, 1:10)
  expect_equal(claims$grandparent, 1:10)
  expect_equal(claims$lonely, integer(0))
})

test_that("confidence intervals cope with one value", {
  got <- .calc_mean_ci(3.7, ci_method = "bootstrap", n_boot = 10)
  expect_equal(got[1], 3.7)
  expect_true(all(is.na(got[5:6])))
  x <- c(1, 2, 3, 4, 5)
  normal <- .calc_mean_ci(x)
  expect_equal(normal[1], 3)
  expect_equal(normal[6] - normal[5], 2 * 1.96 * sd(x) / sqrt(5))
})
