# Runs the bundled example end to end. Expected values do not come from dndsR:
# truth/ holds the CDS as designed before the genomes were laid out, and an independent
# Nei-Gojobori calculation (tools/dev/validation/ng86.py) on those sequences.
skip_if_not_installed("orthologr")
skip_if_not_installed("txdbmaker")
skip_if(Sys.which("diamond") == "", "DIAMOND is not on PATH")

example <- system.file("extdata", "example", package = "dndsR")
work <- withr::local_tempdir(.local_envir = teardown_env())
file.copy(list.files(example, full.names = TRUE), work, recursive = TRUE)
old <- setwd(work)
withr::defer(setwd(old), teardown_env())

truth <- function(species) {
  Biostrings::readDNAStringSet(file.path("truth", paste0(species, "_longest_isoform_cds.fasta")))
}
reference <- utils::read.delim("truth/ng86_reference.tsv", na.strings = c("NA", ""))
suppressMessages(extract_cds(comparison_file = "comparison_file.txt", output_dir = "."))
suppressMessages(suppressWarnings(calculate_dnds(comparison_file = "comparison_file.txt", output_dir = ".", comp_cores = 2)))
result <- utils::read.delim("example_AvB/example_AvB_dnds.tsv")
result$pair <- sub("^spA_", "", result$query_id)
value <- function(pair, column) result[[column]][result$pair == pair]
expected <- function(pair, column) reference[[column]][reference$pair == pair]

test_that("the user-facing example check passes", {
  checks <- suppressMessages(verify_example("."))
  expect_true(all(checks$passed))
  expect_gte(nrow(checks), 13L)
})

test_that("the user-facing example check catches wrong results", {
  dir.create("broken")
  file.copy("truth", "broken", recursive = TRUE)
  dir.create("broken/example_AvB")
  file.copy(list.files("example_AvB", full.names = TRUE), "broken/example_AvB")
  wrong <- utils::read.delim("example_AvB/example_AvB_dnds.tsv")
  wrong$dNdS[wrong$query_id == "spA_g28"] <- 1.615
  utils::write.table(wrong, "broken/example_AvB/example_AvB_dnds.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(suppressMessages(verify_example("broken")), "example checks failed")
})

test_that("the example can be copied, once", {
  folder <- file.path(tempdir(), "copied-example")
  suppressMessages(dndsr_example(folder))
  expect_true(file.exists(file.path(folder, "comparison_file.txt")))
  expect_error(suppressMessages(dndsr_example(folder)), "already")
  expect_error(verify_example(tempdir()), "No results")
})

test_that("every written CDS matches the designed sequence exactly", {
  for (species in c("spA", "spB")) {
    got <- Biostrings::readDNAStringSet(file.path("example_AvB", paste0(species, "_CDS.fasta")))
    want <- truth(species)
    expect_true(all(names(got) %in% names(want)))
    # g30's longest valid isoform is deliberately 45 nt longer than the reference isoform.
    same <- !grepl("g30$", names(got))
    expect_identical(as.character(got[same]), as.character(want[names(got)[same]]))
    expect_equal(Biostrings::width(got[!same]), Biostrings::width(want[names(got)[!same]]) + 45L)
  }
})

test_that("multi-exon minus-strand genes keep their exon order", {
  got <- Biostrings::readDNAStringSet("example_AvB/spA_CDS.fasta")[["spA_g28"]]
  expect_identical(as.character(got), as.character(truth("spA")[["spA_g28"]]))
  expect_false(grepl("*", substr(as.character(Biostrings::translate(got)), 1, 301), fixed = TRUE))
})

test_that("the QC tables record exclusions and isoform choices", {
  a <- utils::read.delim("example_AvB/spA_CDS_qc.tsv")
  b <- utils::read.delim("example_AvB/spB_CDS_qc.tsv")
  expect_equal(a$id[!a$written], "spA_g32")
  expect_setequal(b$id[!b$written], c("spB_g31", "spB_g32"))
  expect_equal(a$transcript[a$id == "spA_g30"], "spA_g30.2") # the extended isoform is longest
  expect_equal(a$n_isoforms[a$id == "spA_g30"], 2L)
  expect_equal(a$transcript[a$id == "spA_g36"], "spA_g36.1") # the longer one has an in-frame stop
})

test_that("pairs with excluded or missing sequences give no result", {
  expect_false(any(c("g31", "g32", "g33", "g34") %in% result$pair))
  expect_equal(nrow(result), 32L)
})

test_that("zero-change cases behave exactly", {
  expect_equal(value("g01", "dN"), 0)
  expect_equal(value("g01", "dS"), 0)
  expect_true(is.na(value("g01", "dNdS")))
  expect_equal(value("g02", "dN"), 0)
  expect_gt(value("g02", "dS"), 0)
  expect_equal(value("g02", "dNdS"), 0)
  # g03 has nonsynonymous changes only. The Comeron method assigns site classes rather than
  # classifying each codon change, so its dS can be small but nonzero (see estimator validation).
  expect_gt(value("g03", "dN"), 0)
  expect_lt(value("g03", "dS"), value("g03", "dN"))
})

test_that("Comeron estimates agree with independent NG86 values within method differences", {
  # g30 is excluded: its longest isoform is extended, so it differs from the reference CDS.
  pairs <- setdiff(result$pair[result$pair %in% reference$pair], c("g01", "g02", "g03", "g30"))
  for (p in pairs) {
    expect_lt(abs(value(p, "dN") - expected(p, "dN")), 0.01 + 0.15 * expected(p, "dN"))
    expect_lt(abs(value(p, "dS") - expected(p, "dS")), 0.015 + 0.25 * expected(p, "dS"))
  }
})

test_that("designed selection classes are recovered", {
  designed <- utils::read.delim("truth/designed_substitutions.tsv")
  ratios <- function(class) vapply(designed$pair[designed$class == class], value, numeric(1), column = "dNdS")
  expect_true(all(ratios("conserved") < 0.2))
  expect_true(all(ratios("neutral_like") > 0.7 & ratios("neutral_like") < 1.3))
  expect_true(all(ratios("excess_nonsynonymous") > 2))
  expect_true(all(ratios("structure") < 0.2)) # includes g28, which was 1.6 before the exon-order fix
  expect_true(all(ratios("multi_isoform") < 0.2))
})

test_that("a provenance record identifies the inputs", {
  prov <- jsonlite::read_json("example_AvB/example_AvB_dnds_provenance.json")
  expect_equal(prov$query_md5, unname(tools::md5sum("example_AvB/spA_CDS.fasta")))
  expect_equal(prov$parameters$dnds_est.method, "Comeron")
  expect_null(.provenance_mismatch(
    "example_AvB/example_AvB_dnds_provenance.json", "example_AvB/spA_CDS.fasta", "example_AvB/spB_CDS.fasta"
  ))
  expect_match(
    .provenance_mismatch("example_AvB/example_AvB_dnds_provenance.json", "example_AvB/spB_CDS.fasta", "example_AvB/spB_CDS.fasta"),
    "changed"
  )
  expect_match(.provenance_mismatch("does-not-exist.json", "a", "b"), "older dndsR")
})

test_that("existing outputs without a QC table are not reused", {
  dir.create("stale")
  writeLines(c(">spA_g01", "ACGT"), "stale/spA_CDS.fasta")
  writeLines(c(">spB_g01", "ACGT"), "stale/spB_CDS.fasta")
  expect_message(
    extract_cds(comparison_name = "stale", query_fasta = "spA.fasta", query_gff = "spA.gff3", output_dir = "."),
    "Regenerating"
  )
  expect_gt(file.info("stale/spA_CDS.fasta")$size, 1000)
})

test_that("proteins can be exported and match the CDS", {
  res <- suppressMessages(extract_cds(
    comparison_name = "prot", query_fasta = "spA.fasta", query_gff = "spA.gff3",
    output_dir = ".", export_proteins = TRUE
  ))
  aa <- Biostrings::readAAStringSet("prot/spA_AA.fasta")
  cds <- Biostrings::readDNAStringSet("prot/spA_CDS.fasta")
  expect_identical(names(aa), names(cds))
  expect_false(any(grepl("*", as.character(aa), fixed = TRUE)))
})

test_that("group_by = 'tx' writes each isoform by transcript ID", {
  suppressMessages(extract_cds(
    comparison_name = "bytx", query_fasta = "spA.fasta", query_gff = "spA.gff3", output_dir = ".", group_by = "tx"
  ))
  got <- Biostrings::readDNAStringSet("bytx/spA_CDS.fasta")
  expect_true(all(c("spA_g29.1", "spA_g29.2", "spA_g30.1", "spA_g30.2", "spA_g36.1") %in% names(got)))
  expect_false("spA_g36.2" %in% names(got))
  expect_identical(as.character(got[["spA_g28.1"]]), as.character(truth("spA")[["spA_g28"]]))
})

test_that("a missing sequence file is an error, not a silent skip", {
  expect_error(
    calculate_dnds(comparison_name = "nowhere", query_seq = "no_q.fasta", subject_seq = "no_s.fasta", output_dir = "."),
    "not found"
  )
})

test_that("downstream steps run and match independent calculations", {
  suppressMessages(append_annotations(comparison_file = "comparison_file.txt", output_dir = "."))
  annot <- utils::read.delim("example_AvB/example_AvB_dnds_annot.tsv")
  expect_equal(nrow(annot), nrow(result))
  expect_equal(annot$q_kegg[annot$query_id == "spA_g20"], "K00001")
  expect_equal(annot$q_gff_seqname[annot$query_id == "spA_g01"], "Chr1")

  # Expected Fisher p-value from the hypergeometric distribution, using base R on this table.
  usable <- annot[!is.na(annot$dNdS) & annot$dNdS < 10, ]
  positive <- usable$dNdS > 1
  has_term <- usable$q_kegg == "K00001"
  a <- sum(positive & has_term)
  expected_p <- sum(stats::dhyper(a:min(sum(has_term), sum(positive)), sum(has_term), sum(!has_term), sum(positive)))

  suppressMessages(term_enrichment(comparison_file = "comparison_file.txt", output_dir = ".", make_plots = FALSE))
  kegg <- utils::read.delim("example_AvB/example_AvB_q_KEGG_enrichment.tsv")
  row <- kegg[kegg$term == "K00001", ]
  expect_equal(row$pos_count, a)
  expect_equal(row$total_count, sum(has_term))
  expect_equal(row$p_value, expected_p, tolerance = 1e-10)
  expect_equal(row$enrichment, (a / sum(positive)) / ((sum(has_term) - a) / sum(!positive)), tolerance = 1e-10)

  suppressMessages(ipr_enrichment(comparison_file = "comparison_file.txt", output_dir = ".", max_prop = 1, make_plots = FALSE))
  domain <- utils::read.delim("example_AvB/example_AvB_q_IPR_Domain_enrichment.tsv")
  expect_true("IPR000001" %in% domain$IPR)
  expect_gt(domain$odds_ratio[domain$IPR == "IPR000001"], 1)
  expect_false("IPR000003" %in% domain$IPR) # a Family term, tested only in its own stratum
})

test_that("region summaries count the right genes", {
  annot <- utils::read.delim("example_AvB/example_AvB_dnds_annot.tsv")
  suppressMessages(dnds_summary(comparison_file = "comparison_file.txt", output_dir = ".", regions_bed = "regions.bed", make_plots = FALSE))
  got <- utils::read.delim("example_AvB/example_AvB_dnds_summary.tsv")
  ok <- !is.na(annot$dNdS) & annot$dNdS < 10
  # BED 0-6000 is 1-6000 in GFF coordinates; a gene counts if it overlaps.
  inside <- ok & annot$q_gff_seqname == "Chr1" & annot$q_gff_start <= 6000 & annot$q_gff_end >= 1
  row <- got[got$scope == "regional" & got$side == "query" & got$region_status == "region", ]
  expect_equal(row$n, sum(inside))
  expect_equal(row$mean_dnds, mean(annot$dNdS[inside]), tolerance = 1e-6)
  everything <- got[got$scope == "global", ]
  expect_equal(everything$n, sum(ok))
  expect_equal(everything$n_removed_na, sum(is.na(annot$dNdS)))
})

test_that("contrasts have the right direction and pairing", {
  suppressMessages(suppressWarnings({
    extract_cds(comparison_file = "two_comparisons.txt", output_dir = ".")
    calculate_dnds(comparison_file = "two_comparisons.txt", output_dir = ".", comp_cores = 2)
    append_annotations(comparison_file = "two_comparisons.txt", output_dir = ".")
    dnds_contrast(comparison_file = "two_comparisons.txt", contrast_file = "contrast_file.txt", output_dir = ".", make_plots = FALSE)
  }))
  got <- utils::read.delim("contrast/AvB_vs_AvC_global_q_vs_q_dnds_contrast.tsv")
  a <- utils::read.delim("example_AvB/example_AvB_dnds_annot.tsv")
  b <- utils::read.delim("example_AvC/example_AvC_dnds_annot.tsv")
  shared <- merge(a, b, by = "query_id")
  shared <- shared[is.finite(shared$dNdS.x) & is.finite(shared$dNdS.y) & shared$dNdS.x < 10 & shared$dNdS.y < 10, ]
  delta <- shared$dNdS.x - shared$dNdS.y
  expect_equal(got$n, nrow(shared))
  expect_equal(got$mean_delta, mean(delta), tolerance = 1e-8)
  expect_gt(got$mean_delta, 0) # genes with excess nonsynonymous changes against B only
  expect_equal(got$wilcox_p_value, stats::wilcox.test(delta)$p.value, tolerance = 1e-8)
  expect_equal(got$n_pos_A, sum(shared$dNdS.x > 1))
  expect_equal(got$n_pos_B, sum(shared$dNdS.y > 1))
})

test_that("ideograms are written", {
  skip_if_not_installed("RIdeogram")
  suppressMessages(dnds_ideogram(comparison_file = "comparison_file.txt", output_dir = ".", sides = "query"))
  expect_true(length(list.files("example_AvB", "ideogram\\.svg$")) > 0L)
})
