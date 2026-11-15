# Tests for the McDonald-Kreitman machinery.
#
# Counting expectations here are derived independently from the implementations
# (codon tables and simple hand enumeration), not from running the code twice.

library(testthat)

skip_if_not_installed("orthologr")
skip_if_not_installed("txdbmaker")
skip_if_not_installed("Biostrings")

mk <- function(name) get(paste0(".mk_", name), envir = asNamespace("dndsR"))

# ---------- codon classifier ----------

# Expected values below are enumerated independently from the standard genetic code
# (all 3 alternatives per changed position, all shortest paths), not from dndsR.

test_that("codon degeneracy classes follow the genetic code", {
  deg <- mk("syn_alternatives")
  expect_equal(deg("GGT", 3L), 3L)   # Gly GGN: 4-fold
  expect_equal(deg("GGT", 1L), 0L)   # Gly: 0-fold
  expect_equal(deg("GAA", 3L), 1L)   # Glu GAR: 2-fold
  expect_equal(deg("ATT", 3L), 2L)   # Ile ATH: 3-fold
  expect_equal(deg("ATG", 1L), 0L)   # Met: 0-fold
})

test_that("codon path classifier averages 4-fold, 0-fold, and mid-fold steps", {
  cc <- mk("classify_codon")
  expect_equal(cc("GCT", "GCC"), c(1, 0, 0))          # Ala GCN, 4-fold
  expect_equal(cc("GCT", "GTT"), c(0, 1, 0))          # Ala -> Val, 0-fold
  expect_equal(cc("GAA", "GAG"), c(0, 0, 1))          # Glu GAR, 2-fold
  expect_equal(cc("ATT", "ATC"), c(0, 0, 1))          # Ile ATH, 3-fold
  # Two differences: both paths give one 0-fold and one mid-fold step.
  expect_equal(cc("CTT", "TTA"), c(0.5, 0.5, 1))
  # Three differences: no 4-fold step is possible, and one step is mid-fold.
  expect_equal(cc("AAA", "TTT"), c(0, 2, 1))
  expect_equal(cc("GGT", "GGT"), c(0, 0, 0))
})

test_that("divergence walk excludes SNPs and stops, and detects frame shifts", {
  ref <- "ATGGCTGCTTAA"
  oth <- "ATGGCTGCATAA"  # codon 3 GCT -> GCA (syn)
  w <- mk("pair_walk")(ref, oth, integer(0))
  expect_equal(unname(w[c("dn", "ds", "n_used", "frame_ok")]), c(0, 1, 3, 1))
  # Excluding codon 3 (CDS positions 7:9) removes the synonymous difference.
  w <- mk("pair_walk")(ref, oth, 7:9)
  expect_equal(unname(w[c("dn", "ds", "n_used", "n_excl")]), c(0, 0, 2, 1))
  # A reference gap inside a codon breaks the frame for everything after it.
  w <- mk("pair_walk")("ATG-CTGCTTAA", "ATGGCTGCTTAA", integer(0))
  expect_equal(unname(w["frame_ok"]), 0)
  # Every differing position lands in exactly one class: 4-fold, 0-fold, or mid-fold.
  # Glu GAA -> Phe TTT, three differences; enumerated independently: 0.25, 2, 0.75.
  w <- mk("pair_walk")("ATGGAA", "ATGTTT", integer(0))
  expect_equal(unname(w[c("ds", "dn", "mid")]), c(0.25, 2, 0.75))
  expect_equal(sum(unname(w[c("ds", "dn", "mid")])), 3)
})

# ---------- SNP parsing ----------

write_vcf <- function(path, lines) {
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=10000>",
    "##contig=<ID=chr2,length=2000>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    lines
  ), path)
}

test_that("VCF reader keeps only bi-allelic SNPs, dedups, and reads header lengths", {
  vcf <- file.path(tempdir(), "t.vcf")
  write_vcf(vcf, c(
    "chr1\t100\t.\tA\tG\t.\t.\t.",        # keep
    "chr1\t100\t.\tA\tG\t.\t.\t.",        # duplicate -> dedup
    "chr1\t200\t.\tC\tT,G\t.\t.\t.",      # multi-allelic -> drop
    "chr1\t300\t.\tA\tAT\t.\t.\t.",       # indel -> drop
    "chr1\t400\t.\t<DEL>\tA\t.\t.\t.",    # symbolic REF -> drop
    "chr1\t500\t.\tN\tA\t.\t.\t.",        # non-ACGT REF -> drop
    "chr1\t600\t.\tA\tA\t.\t.\t.",        # REF == ALT -> drop
    "chr2\t100\t.\tG\tC\t.\t.\t."         # keep (other contig)
  ))
  sn <- mk("read_snps")(vcf)
  expect_equal(paste(sn$snps$chrom, sn$snps$pos), c("chr1 100", "chr2 100"))
  expect_equal(sn$meta$format, "vcf")
  expect_equal(sn$meta$n_duplicates, 1L)
  expect_equal(sn$meta$n_dropped, 6L)  # 5 rejected records + 1 duplicate
  expect_equal(sn$meta$contig_len[["chr1"]], 10000)
})

test_that("BCF is rejected with a clear error", {
  f <- file.path(tempdir(), "x.bcf")
  file.create(f)
  expect_error(mk("read_snps")(f), "BCF")
})

test_that("BED input keeps only 1-bp intervals and is 0-based to 1-based", {
  bed <- file.path(tempdir(), "t.bed")
  writeLines(c("track name=x", "chr1\t99\t100\tname1", "chr1\t100\t102\t2bp"), bed)
  sn <- mk("read_snps")(bed)
  expect_equal(sn$snps$pos, 100L)
  expect_equal(sn$meta$format, "bed")
  expect_equal(sn$meta$n_dropped, 1L)
})

# ---------- REF-allele audit and Poly counting need a mapped gene ----------
# Build a small synthetic genome on the fly and check the whole map -> audit -> count path.

# Built once per test file: extract_cds and txdbmaker cost seconds each time.
mini_cache <- new.env(parent = emptyenv())
mini_genome <- function() {
  if (is.null(mini_cache$g)) {
    dir <- tempfile("mk_mini_")
    dir.create(dir)
    mini_cache$g <- build_mini_genome(dir)
  }
  mini_cache$g
}

build_mini_genome <- function(dir) {
  genome <- paste0(rep("A", 1000), collapse = "")
  cds1 <- "ATGGGTGGCGGAGGGGCTGCTGCTGCTGCTGCTTAA"   # 11 codons + stop
  stopifnot(nchar(cds1) == 36)
  genome <- paste0(substr(genome, 1, 99), cds1, substr(genome, 136, 1000))
  fasta <- file.path(dir, "mini.fasta")
  writeLines(c(">chrM", genome), fasta)
  gff <- file.path(dir, "mini.gff3")
  writeLines(c(
    "##gff-version 3",
    "chrM\ttem\tgene\t100\t135\t.\t+\t.\tID=gA;Name=gA",
    "chrM\ttem\tmRNA\t100\t135\t.\t+\t.\tID=gA.1;Parent=gA",
    "chrM\ttem\tCDS\t100\t135\t.\t+\t0\tID=gA.1.cds;Parent=gA.1"
  ), gff)
  cds_fa <- file.path(dir, "mini_CDS.fasta")
  suppressMessages(dndsR::extract_cds(comparison_name = "mini", query_fasta = fasta,
                                      query_gff = gff, output_dir = dir))
  cds_fa <- file.path(dir, "mini", "mini_CDS.fasta")
  list(fasta = fasta, gff = gff, cds_fa = cds_fa, genome = genome, cds = cds1)
}

test_that("CDS map reconstructs the CDS FASTA exactly", {
  skip_if_not_installed("txdbmaker")
  g <- mini_genome()
  map <- mk("cds_map")(g$gff, g$fasta, g$cds_fa)
  pos <- map$pos[["gA"]]
  s <- g$genome
  got <- paste(vapply(pos, function(i) substr(s, i, i), character(1)), collapse = "")
  expect_identical(got, g$cds)
})

test_that("REF audit flags wrong-reference SNPs; clean SNPs pass", {
  skip_if_not_installed("txdbmaker")
  g <- mini_genome()
  map <- mk("cds_map")(g$gff, g$fasta, g$cds_fa)
  pos <- map$pos[["gA"]]
  good <- data.frame(chrom = "chrM", pos = pos[c(4, 13)],
                     ref = substr(g$genome, pos[c(4, 13)], pos[c(4, 13)]),
                     alt = "A", stringsAsFactors = FALSE)
  bad <- good; bad$ref <- "T"   # wrong on purpose (unless ref really is T)
  bad$ref <- ifelse(bad$ref == "A", "G", "A")
  mapped_good <- mk("map_snps")(good, map)
  mapped_bad <- mk("map_snps")(bad, map)
  expect_equal(unname(mapped_good$audit["n_mismatch"]), 0L)
  expect_equal(unname(mapped_bad$audit["n_mismatch"]), 2L)
})

test_that("SNP mapping and polymorphism counting on the mini genome", {
  skip_if_not_installed("txdbmaker")
  g <- mini_genome()
  map <- mk("cds_map")(g$gff, g$fasta, g$cds_fa)
  # CDS: ATG GGT GGC GGA GGG CTG CTG CTG CTG CTT AA(T)
  # codon 2 (GGN, 4-fold): site at CDS pos 6 -> ChrM 105  -> 4-fold degenerate -> Ps
  # codon 1: position 102 (T, 2nd base of ATG): nonsyn change -> Pn
  snps <- data.frame(
    chrom = "chrM",
    pos = c(map$pos[["gA"]][6], map$pos[["gA"]][2]),
    ref = c(g$genome[map$pos[["gA"]][6]], substr(g$genome, map$pos[["gA"]][2], map$pos[["gA"]][2])),
    alt = c("C", "C"),
    stringsAsFactors = FALSE
  )
  hits <- mk("map_snps")(snps, map)$hits
  poly <- mk("poly_counts")(snps, map, hits)
  expect_equal(poly$gA[["ps"]], 1L)
  expect_equal(poly$gA[["pn"]], 1L)
  expect_equal(poly$gA[["mid"]], 0L)
})

test_that("polymorphisms at 2- and 3-fold sites are excluded from Pn and Ps", {
  # A hand-built map keeps this independent of genome extraction: codons are
  # ATG (Met), GAA (Glu, 2-fold at position 3), ATT (Ile, 3-fold), TAA (stop).
  map <- list(seq = c(g1 = "ATGGAAATTTAAGGT"))
  # CDS indices: 4 = GAA pos 1 (0-fold), 6 = GAA pos 3 (2-fold),
  #              9 = ATT pos 3 (3-fold), 10 = TAA pos 1 (stop), 15 = GGT pos 3 (4-fold)
  hits <- data.frame(gene = rep("g1", 5), snp_i = 1:5,
                     cds_idx = c(4L, 6L, 9L, 10L, 15L))
  snps <- data.frame(chrom = "c", pos = 1:5, ref = NA_character_, alt = NA_character_,
                     stringsAsFactors = FALSE)
  poly <- mk("poly_counts")(snps, map, hits)$g1
  expect_equal(poly[["pn"]], 1L)      # GAA position 1 is 0-fold
  expect_equal(poly[["ps"]], 1L)      # GGT position 3 is 4-fold
  expect_equal(poly[["mid"]], 2L)     # GAA position 3 and ATT position 3
  expect_equal(poly[["dropped"]], 1L) # the stop codon
})

test_that("snp_side = subject counts the subject's SNPs and excludes its codons", {
  # Two genes with identical CDS so the same VCF position is valid on either side.
  cds1 <- "ATGGGTGGCGGAGGGGCTGCTGCTGCTGCTGCTTAA"
  g <- paste0(rep("A", 1000), collapse = "")
  g <- paste0(substr(g, 1, 99), cds1, substr(g, 136, 199), cds1, substr(g, 236, 1000))
  dir <- withr::local_tempdir()
  fasta <- file.path(dir, "two.fasta"); writeLines(c(">chrM", g), fasta)
  gff <- file.path(dir, "two.gff3")
  writeLines(c(
    "##gff-version 3",
    "chrM\ttem\tgene\t100\t135\t.\t+\t.\tID=gA;Name=gA",
    "chrM\ttem\tmRNA\t100\t135\t.\t+\t.\tID=gA.1;Parent=gA",
    "chrM\ttem\tCDS\t100\t135\t.\t+\t0\tID=gA.1.cds;Parent=gA.1",
    "chrM\ttem\tgene\t200\t235\t.\t+\t.\tID=gB;Name=gB",
    "chrM\ttem\tmRNA\t200\t235\t.\t+\t.\tID=gB.1;Parent=gB",
    "chrM\ttem\tCDS\t200\t235\t.\t+\t0\tID=gB.1.cds;Parent=gB.1"
  ), gff)
  suppressMessages(dndsR::extract_cds(comparison_name = "two", query_fasta = fasta,
                                      query_gff = gff, output_dir = dir))
  cds_fa <- file.path(dir, "two", "two_CDS.fasta")
  map <- mk("cds_map")(gff, fasta, cds_fa)
  posB <- map$pos[["gB"]]
  vcf <- file.path(dir, "subject.vcf")
  writeLines(c("##fileformat=VCFv4.2",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
               sprintf("chrM\t%d\t.\tT\tC\t.\t.\t.", posB[6])), vcf)  # GGT/GGC, synonymous
  aln <- file.path(dir, "pair_pair_pal2nal.aln")
  writeLines(c(">gA", cds1, ">gB", sub("GGT", "GGC", cds1)), aln)
  dnds_res <- data.frame(query_id = "gA", subject_id = "gB", dN = 0, dS = 1, dNdS = 0)
  out <- mk("run")(comp_dir = dir, comparison_name = "side_subject", dnds_res = dnds_res,
                   aln_files = aln, snp_file = vcf, snp_side = "subject",
                   query_fasta = fasta, subject_fasta = fasta,
                   query_gff = gff, subject_gff = gff,
                   query_cds_fasta = cds_fa, subject_cds_fasta = cds_fa)
  expect_equal(nrow(out), 1L)
  expect_equal(out$reference_side, "subject")
  expect_equal(out$ps, 1L)          # the subject's synonymous SNP
  expect_equal(out$pn, 0L)
  # The only aligned difference is that same SNP, in a codon excluded from divergence.
  expect_equal(out$ds, 0)
  expect_equal(out$codons_excluded_polymorphic, 1L)
})

# ---------- SNP-side validation ----------

test_that("snp_file without snp_side (and vice versa) is an error", {
  expect_error(mk("validate_snp_args")("x.vcf", NULL, "cmp"), "snp_side")
  expect_error(mk("validate_snp_args")(NULL, "query", "cmp"), "snp_file")
  expect_error(mk("validate_snp_args")("x.vcf", "bogus", "cmp"), "query")
  expect_null(mk("validate_snp_args")(NULL, NULL, "cmp"))
  expect_equal(mk("validate_snp_args")("x.vcf", "Query", "cmp")$side, "query")
})

# ---------- end to end through calculate_dnds ----------

test_that("calculate_dnds computes the MK table alongside dN/dS", {
  skip_if(fast_tests(), "slow: runs orthologr")
  skip_if(Sys.which("diamond") == "", "DIAMOND is not on PATH")
  example <- system.file("extdata", "example", package = "dndsR")
  work <- withr::local_tempdir()
  file.copy(list.files(example, full.names = TRUE), work, recursive = TRUE)
  old <- setwd(work); withr::defer(setwd(old))
  suppressMessages(dndsR::extract_cds(comparison_file = "comparison_file.txt", output_dir = "."))

  # Build a VCF on the query genome with one synonymous + one nonsynonymous SNP in spA_g04
  cds_fa <- file.path("example_AvB", "spA_CDS.fasta")
  map <- mk("cds_map")("spA.gff3", "spA.fasta", cds_fa)
  genome <- as.character(Biostrings::readDNAStringSet("spA.fasta")[["Chr1"]])
  pos_ref <- function(cds_i) {
    p <- map$pos[["spA_g04"]][cds_i]
    list(pos = p, ref = substr(genome, p, p))
  }
  cds <- Biostrings::readDNAStringSet(cds_fa)
  g04 <- strsplit(as.character(cds[["spA_g04"]]), "")[[1]]
  code <- Biostrings::GENETIC_CODE
  syn <- NULL
  for (i in seq(3, length(g04) - 2, by = 3)) {
    c3 <- paste(g04[(i - 2):i], collapse = "")
    alt <- c("A", "C", "G", "T")[c("A", "C", "G", "T") != g04[i]][1]
    if (code[[paste0(substr(c3, 1, 2), alt)]] == code[[c3]] &&
        grepl("^(GG|CC|CG|AC|GC|TC|CT|GT)[ACGT]$", c3)) { syn <- list(cds_i = i, alt = alt); break }
  }
  stopifnot(!is.null(syn))
  prs <- pos_ref(syn$cds_i)
  vcf <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    sprintf("Chr1\t%d\t.\t%s\t%s\t.\t.\t.", prs$pos, prs$ref, syn$alt)
  )
  writeLines(vcf, "snps.vcf")
  writeLines(paste("example_AvB", "spA.fasta", "spA.gff3", "spB.fasta", "spB.gff3", "snps.vcf", "query"),
             "comparison_snp.txt")
  suppressWarnings(suppressMessages(
    dndsR::calculate_dnds(comparison_file = "comparison_snp.txt", output_dir = ".", comp_cores = test_cores())
  ))
  mk_tbl <- utils::read.delim("example_AvB/example_AvB_mk.tsv")
  g04 <- mk_tbl[mk_tbl$query_id == "spA_g04", ]
  expect_equal(g04$pn, 0L)
  expect_equal(g04$ps, 1L)   # the chosen codon is 4-fold degenerate
  # The designed pair has 48 differing codon positions. Each is counted in exactly
  # one class, so the three columns must add up to 48 whatever the split is.
  expect_equal(g04$dn + g04$ds + g04$div_steps_midfold, 48)
  expect_gt(g04$dn, 0)
  expect_gt(g04$ds, 0)
  expect_gt(g04$div_steps_midfold, 0)   # 2-/3-fold sites are excluded, not dropped
  expect_equal(g04$snps_midfold, 0L)
  expect_equal(g04$codons_excluded_polymorphic, 1L)
  pooled <- utils::read.delim("example_AvB/example_AvB_mk_pooled.tsv")
  expect_equal(pooled$ps, 1L)
  expect_equal(pooled$pn, 0L)
  expect_gt(pooled$ds, 0)
})

# ---------- BED-driven polymorphism classification ----------

test_that("BED sites are classified by codon degeneracy", {
  skip_if_not_installed("txdbmaker")
  g <- mini_genome()
  map <- mk("cds_map")(g$gff, g$fasta, g$cds_fa)
  bed_snps <- data.frame(
    chrom = c("chrM", "chrM"),
    pos = c(map$pos[["gA"]][6],   # 3rd pos of GGT (4-fold degenerate) -> Ps
            map$pos[["gA"]][7]),  # 1st pos of GGC (0-fold)             -> Pn
    ref = NA_character_,
    alt = NA_character_,
    stringsAsFactors = FALSE
  )
  hits <- mk("map_snps")(bed_snps, map)$hits
  poly <- mk("poly_counts")(bed_snps, map, hits)
  expect_equal(poly$gA[["ps"]], 1L)
  expect_equal(poly$gA[["pn"]], 1L)
  expect_equal(poly$gA[["mid"]], 0L)
  # VCF and BED now share one rule, so the same two sites give the same counts.
  vcf_snps <- data.frame(chrom = bed_snps$chrom, pos = bed_snps$pos,
                         ref = substr(g$genome, bed_snps$pos, bed_snps$pos),
                         alt = c("C", "C"), stringsAsFactors = FALSE)
  vpoly <- mk("poly_counts")(vcf_snps, map, mk("map_snps")(vcf_snps, map)$hits)$gA
  expect_equal(vpoly[["pn"]], poly$gA[["pn"]])
  expect_equal(vpoly[["ps"]], poly$gA[["ps"]])
  expect_equal(vpoly[["mid"]], poly$gA[["mid"]])
})

# ---------- .mk_run: validation and 1-to-1 enforcement (no orthologr needed) ----------

fake_mk_inputs <- function(dir, g) {
  # two trivial "pal2nal" alignments for pair gA/gB with one synonymous difference
  cds <- g$cds
  alt_cds <- sub("GGT", "GGC", cds)               # synonymous change in codon 2
  aln <- file.path(dir, "pair1_pair_pal2nal.aln")
  writeLines(c(">gA", cds, ">gB", alt_cds), aln)
  list(aln = aln, cds = cds)
}

mk_run_args <- function(g, snp_file, comp_dir) {
  dnds_res <- data.frame(query_id = "gA", subject_id = "gB", dN = 0, dS = 1, dNdS = 0)
  fi <- fake_mk_inputs(dirname(snp_file), g)
  list(comp_dir = comp_dir, comparison_name = "mini", dnds_res = dnds_res,
       aln_files = fi$aln, snp_file = snp_file, snp_side = "query",
       query_fasta = g$fasta, subject_fasta = g$fasta,
       query_gff = g$gff, subject_gff = g$gff,
       query_cds_fasta = g$cds_fa, subject_cds_fasta = g$cds_fa)
}

test_that("MK stops when the SNP contigs do not match the reference genome", {
  skip_if_not_installed("txdbmaker")
  dir <- withr::local_tempdir()   # writable output dir for this test
  g <- mini_genome()
  vcf <- file.path(dir, "wrong_contigs.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chrZZ\t105\t.\tT\tC\t.\t.\t."
  ), vcf)
  expect_error(do.call(mk("run"), mk_run_args(g, vcf, dir)), "No contig")
})

test_that("MK stops when a VCF contig length disagrees with the FASTA index", {
  skip_if_not_installed("txdbmaker")
  dir <- withr::local_tempdir()   # writable output dir for this test
  g <- mini_genome()
  vcf <- file.path(dir, "bad_len.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chrM,length=999>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chrM\t105\t.\tT\tC\t.\t.\t."
  ), vcf)
  expect_error(do.call(mk("run"), mk_run_args(g, vcf, dir)), "Contig length mismatch")
})

test_that("MK stops when most REF alleles disagree with the reference genome", {
  skip_if_not_installed("txdbmaker")
  dir <- withr::local_tempdir()   # writable output dir for this test
  g <- mini_genome()
  map <- mk("cds_map")(g$gff, g$fasta, g$cds_fa)
  pos <- map$pos[["gA"]]
  # every REF allele is deliberately wrong (never the true base); ALT differs from REF
  wrong_ref <- function(p) { b <- substr(g$genome, p, p); ifelse(b == "A", "C", "A") }
  wrong_alt <- function(r) ifelse(r == "C", "G", "C")
  rows <- vapply(pos[1:12], function(p) { r <- wrong_ref(p); sprintf("chrM\t%d\t.\t%s\t%s\t.\t.\t.", p, r, wrong_alt(r)) }, character(1))
  # make the first CDS position's ALT valid biallelic and refs definitely wrong
  vcf <- file.path(dir, "wrong_ref.vcf")
  writeLines(c("##fileformat=VCFv4.2",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", rows), vcf)
  expect_error(do.call(mk("run"), mk_run_args(g, vcf, dir)), "do not match")
})

test_that("MK drops non-1-to-1 pairs and writes both tables", {
  skip_if_not_installed("txdbmaker")
  dir <- withr::local_tempdir()   # writable output dir for this test
  g <- mini_genome()
  vcf <- file.path(dir, "empty.vcf")
  writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), vcf)
  dnds_res <- data.frame(
    query_id = c("gA", "gA", "gA"), subject_id = c("gB", "gC", "gD"),
    dN = c(0, 0, 0), dS = c(1, 1, 1), dNdS = c(0, 0, 0))
  args <- mk_run_args(g, vcf, dir)
  args$dnds_res <- dnds_res
  # alignments for gB and gC as stand-ins for the two "orthologs" of gA
  cds <- g$cds
  alt <- sub("GGT", "GGC", cds)
  a1 <- file.path(dir, "a1_pal2nal.aln"); a2 <- file.path(dir, "a2_pal2nal.aln")
  writeLines(c(">gA", cds, ">gB", alt), a1)
  writeLines(c(">gA", cds, ">gC", alt), a2)
  args$aln_files <- c(a1, a2)
  out <- do.call(mk("run"), args)
  expect_equal(nrow(out), 0L)  # all three copies of gA are dropped as non-1-to-1
  expect_true(file.exists(file.path(dir, "mini_mk.tsv")))
  expect_true(file.exists(file.path(dir, "mini_mk_pooled.tsv")))
})

test_that("fisher/neutrality statistics behave at the margins", {
  stats_fn <- mk("stats")
  expect_true(is.na(stats_fn(0, 0, 5, 2)["alpha"]))          # ps == 0
  expect_true(is.na(stats_fn(1, 1, 0, 2)["neutrality_index"])) # dn == 0
  expect_equal(unname(stats_fn(1, 1, 4, 1)["alpha"]), 0.75)  # 1 - (1*1)/(4*1)
  expect_true(is.na(stats_fn(2, 8, 0, 0)["dos"]))            # dn + ds == 0 -> undefined
})
