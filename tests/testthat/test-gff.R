gff <- function(lines) {
  path <- tempfile(fileext = ".gff3")
  writeLines(lines, path, sep = "\n")
  path
}
feature <- function(attrs, type = "gene", id = "chr1") paste(id, "src", type, 100, 200, ".", "+", ".", attrs, sep = "\t")

test_that("only a # at the start of a line is a comment", {
  g <- .read_gff_features(gff(c(
    "##gff-version 3",
    "# a comment GO:9999999",
    feature("ID=a;Note=serine#threonine;Ontology_term=GO:0000001"),
    "###",
    feature("ID=b;Note=#first")
  )))
  expect_equal(nrow(g), 2L)
  expect_equal(g$attributes[1], "ID=a;Note=serine#threonine;Ontology_term=GO:0000001")
  expect_equal(g$attributes[2], "ID=b;Note=#first")
  expect_type(g$start, "integer")
})

test_that("sequence after ##FASTA is not read as features", {
  g <- .read_gff_features(gff(c(feature("ID=a"), "##FASTA", ">chr1", "ACGT")))
  expect_equal(nrow(g), 1L)
})

test_that("line endings, missing fields, and whitespace-separated files are handled", {
  g <- .read_gff_features(gff(paste0(feature("ID=a"), "\r")))
  expect_equal(g$attributes, "ID=a")
  g <- .read_gff_features(gff(c(feature("."), "chr1\tsrc\tgene\t1\t2\t.\t+\t.")))
  expect_true(all(is.na(g$attributes)))
  g <- .read_gff_features(gff("chr1 src gene 100 200 . + . ID=a;Note=two words#here"))
  expect_equal(g$attributes, "ID=a;Note=two words#here")
  expect_equal(g$end, 200L)
})

test_that("a file with no features is an error", {
  expect_error(.read_gff_features(gff("##gff-version 3")), "no feature lines")
  expect_error(.read_gff_features(gff("not a gff")), "malformed")
})

test_that("GFF3 percent-encoding is decoded once, and last for the percent sign", {
  expect_equal(.decode_gff3("a%2Cb%3Bc%3Dd%26e%23f%3Ag"), "a,b;c=d&e#f:g")
  expect_equal(.decode_gff3("a%2cb"), "a,b")
  expect_equal(.decode_gff3("100%252C"), "100%2C")
  expect_equal(.decode_gff3("no codes"), "no codes")
})

test_that("append_annotations keeps terms after a # and adds none from comments", {
  dir <- withr::local_tempdir()
  q <- file.path(dir, "q.gff3"); s <- file.path(dir, "s.gff3")
  writeLines(c(
    "##gff-version 3",
    "# GO:0009999 appears only in this comment",
    feature("ID=gene-A;Note=serine#threonine;Ontology_term=GO:0000001;interpro=IPR000001"),
    feature("ID=gene-B;Note=encoded%23text;Ontology_term=GO%3A0000002%2CGO%3A0000003"),
    feature("ID=gene-C;Ontology_term=GO:0000004"),
    feature("ID=GO:0000005;Note=the ID is not a term"),
    "##FASTA", ">chr1", "ACGTGO:0000006"
  ), q)
  writeLines(feature("ID=gene-S1;Ontology_term=GO:0000007"), s)
  table <- file.path(dir, "t_dnds.tsv")
  utils::write.table(data.frame(query_id = c("A", "B", "C", "GO:0000005"), subject_id = "S1", dN = 0.1, dS = 0.2, dNdS = 0.5),
                     table, sep = "\t", quote = FALSE, row.names = FALSE)
  out <- file.path(dir, "out.tsv")
  suppressMessages(append_annotations(dnds_file = table, query_gff = q, subject_gff = s, output_file = out, threads = 1))
  got <- utils::read.delim(out, quote = "", comment.char = "")
  expect_equal(got$q_go, c("GO:0000001", "GO:0000002;GO:0000003", "GO:0000004", NA))
  expect_equal(got$q_ipr, c("IPR000001", NA, NA, NA))
  expect_equal(unique(got$s_go), "GO:0000007")
  expect_false(any(grepl("GO:0009999|GO:0000006", unlist(got))))
})
