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

test_that("split_comparisons separates subgenomes and ignores FASTA descriptions", {
  skip_if_not_installed("Biostrings")
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  writeLines(c(">Chr1A a description", "ACGTACGTACGT", ">Chr1B", "TTTTACGTACGT"), "q.fasta")
  writeLines(c(">Chr1A", "ACGTACGTACGT", ">Chr1B", "TTTTACGTACGT"), "s.fasta")
  gff <- paste(c("Chr1A", "Chr1B"), "x", "gene", 1, 12, ".", "+", ".", c("ID=g1", "ID=g2"), sep = "\t")
  writeLines(c("##gff-version 3", gff), "q.gff3")
  writeLines(c("##gff-version 3", gff), "s.gff3")
  writeLines("c1 q.fasta q.gff3 s.fasta s.gff3", "c.txt")
  suppressMessages(split_comparisons("c.txt"))
  expect_true(all(file.exists(c("q_A.fasta", "q_B.fasta", "s_A.gff3", "s_B.gff3", "c_split.txt"))))
  expect_equal(sub("^>(\\S+).*", "\\1", readLines("q_A.fasta")[1]), "Chr1A")
  expect_false(any(grepl("Chr1B", readLines("q_A.gff3"))))
  expect_equal(length(readLines("c_split.txt")), 2L)
})

test_that("pairwise state contrasts match a directly fitted logistic model", {
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  # Group A: 30 positive of 100. Group B: 10 positive of 100. Others are purifying.
  table <- data.frame(
    group = rep(c("A", "B"), each = 100),
    dNdS = c(rep(2, 30), rep(0.1, 70), rep(2, 10), rep(0.1, 90))
  )
  utils::write.table(table, "t.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  suppressMessages(dnds_state_contrast(dnds_table_file = "t.tsv", output_dir = ".", make_plots = FALSE))
  files <- list.files(".", "state.*\\.tsv$", recursive = TRUE, full.names = TRUE)
  expect_gt(length(files), 0L)
  got <- utils::read.delim(files[1])
  row <- got[got$state == "Positive", ]
  expect_equal(nrow(row), 1L)
  expect_equal(row$odds_ratio, (30 / 70) / (10 / 90), tolerance = 1e-5)
})

test_that("RefSeq-style IDs connect sequences to their annotations", {
  skip_if_not_installed("txdbmaker")
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  example <- system.file("extdata", "example", package = "dndsR")
  file.copy(file.path(example, "spA.fasta"), ".")
  gff <- readLines(file.path(example, "spA.gff3"))
  # Rewrite to the prefixes used by NCBI: gene-, rna-, and cds-.
  gff <- sub("\tgene\t(.*)ID=", "\tgene\t\\1ID=gene-", gff)
  gff <- sub("\tmRNA\t(.*)ID=([^;]*);Parent=", "\tmRNA\t\\1ID=rna-\\2;Parent=gene-", gff)
  gff <- sub("\tCDS\t(.*)ID=([^;]*);Parent=", "\tCDS\t\\1ID=cds-\\2;Parent=rna-", gff)
  gff <- gff[!grepl("\texon\t", gff)]
  writeLines(gff, "ncbi.gff3")
  suppressMessages(extract_cds(comparison_name = "c", query_fasta = "spA.fasta", query_gff = "ncbi.gff3", output_dir = "."))
  ids <- sub(">", "", grep(">", readLines("c/spA_CDS.fasta"), value = TRUE))
  expect_true("spA_g01" %in% ids) # the prefix is gone
  writeLines(c("query_id\tsubject_id\tdN\tdS\tdNdS", paste(ids, ids, 0.1, 0.2, 0.5, sep = "\t")), "d.tsv")
  got <- suppressMessages(append_annotations(dnds_file = "d.tsv", query_gff = "ncbi.gff3", subject_gff = "ncbi.gff3"))
  expect_false(anyNA(got$q_gff_start))
  expect_equal(got$q_ipr[got$query_id == "spA_g01"], "IPR000007")
  expect_error(
    suppressMessages(append_annotations(dnds_file = "d.tsv", query_gff = file.path(example, "spB.gff3"), subject_gff = "ncbi.gff3")),
    "None of the"
  )
})

test_that("identifier prefixes are removed only at the start", {
  expect_equal(.normalize_gff_id(c("gene-AT1", "rna-NM_1.2", "cds-NP_1", "mygene-1", "g1", NA)),
               c("AT1", "NM_1.2", "NP_1", "mygene-1", "g1", NA))
})
