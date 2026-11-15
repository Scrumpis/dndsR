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
# (every shortest path, each step compared with the amino acid of the codon it starts
# from), not from dndsR.

test_that("a polymorphism class is the amino-acid effect of its allele", {
  cls <- mk("aa")
  expect_identical(cls("GAA"), "E")
  expect_identical(cls("GAG"), "E")   # synonymous, 2-fold position
  expect_identical(cls("TAA"), "*")
})

test_that("codon pairs are classed by the amino-acid change of each step", {
  cc <- mk("classify_codon")
  expect_equal(cc("GCT", "GCC"), c(1, 0))   # Ala -> Ala, synonymous
  expect_equal(cc("GCT", "GTT"), c(0, 1))   # Ala -> Val
  expect_equal(cc("GAA", "GAG"), c(1, 0))   # Glu -> Glu at a 2-fold position: still Ps
  expect_equal(cc("ATT", "ATC"), c(1, 0))   # Ile -> Ile at a 3-fold position: still Ps
  # Leu CTT -> Leu TTA, two differences. Path C->T first goes through Phe, so both
  # steps are nonsynonymous; path T->A first goes through CTA, Leu, so both are
  # synonymous. Averaged: one of each.
  expect_equal(cc("CTT", "TTA"), c(1, 1))
  # Three differences; enumerated independently: 0.25 synonymous, 2.75 nonsynonymous.
  expect_equal(cc("AAA", "TTT"), c(0.25, 2.75))
  # Direction does not matter: both codons of a pair see the same amino acids.
  expect_equal(cc("AGA", "GGA"), cc("GGA", "AGA"))
  expect_equal(cc("AGA", "GGA"), c(0, 1))   # Arg -> Gly
  expect_equal(cc("GGT", "GGT"), c(0, 0))
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
  # Every differing position lands in exactly one class: ds + dn equals the number of
  # differing positions. Glu GAA -> Phe TTT, three differences: 0.25 + 2.75.
  w <- mk("pair_walk")("ATGGAA", "ATGTTT", integer(0))
  expect_equal(unname(w[c("ds", "dn")]), c(0.25, 2.75))
  expect_equal(sum(unname(w[c("ds", "dn")])), 3)
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

test_that("BED input reads the alternative allele in column 4", {
  bed <- file.path(tempdir(), "t.bed")
  writeLines(c("track name=x", "chr1\t99\t100\tG", "chr1\t100\t102\tT"), bed)
  sn <- mk("read_snps")(bed)
  expect_equal(sn$snps$pos, 100L)
  expect_equal(sn$snps$alt, "G")
  expect_equal(sn$meta$format, "bed")
  expect_equal(sn$meta$n_dropped, 1L)   # the 2-bp interval
})

test_that("a BED file without an allele is refused, because Pn/Ps need it", {
  bed <- file.path(tempdir(), "three.bed")
  writeLines(c("chr1\t99\t100", "chr1\t199\t200"), bed)
  expect_error(mk("read_snps")(bed), "fourth column with the alternative allele")
  # A name that is not a single base is not accepted as an allele either.
  writeLines(c("chr1\t99\t100\trs12345"), bed)
  expect_error(mk("read_snps")(bed), "fourth column with the alternative allele")
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
  # CDS pos 6 is the third base of GGT: T -> C gives GGC, still Gly -> Ps.
  # CDS pos 2 is the second base of ATG: T -> C gives ACG, Thr -> Pn.
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
  expect_equal(poly$gA[["dropped"]], 0L)
})

test_that("polymorphisms are classed by their allele, including 2- and 3-fold sites", {
  # A hand-built map keeps this independent of genome extraction: codons are
  # ATG (Met), GAA (Glu, 2-fold at position 3), ATT (Ile, 3-fold), TAA (stop),
  # GGT (Gly: GGC, GGA and GGG are all Gly).
  map <- list(seq = c(g1 = "ATGGAAATTTAAGGT"), strand = c(g1 = "+"))
  hits <- data.frame(gene = rep("g1", 7), snp_i = 1:7,
                     cds_idx = c(4L, 6L, 9L, 10L, 13L, 15L, 17L))
  snps <- data.frame(chrom = "c", pos = 1:7, ref = NA_character_,
                     alt = c("A",   # GAA pos 1 -> AAA, Lys: nonsynonymous
                             "G",   # GAA pos 3 -> GAG, Glu: synonymous at a 2-fold site
                             "C",   # ATT pos 3 -> ATC, Ile: synonymous at a 3-fold site
                             "C",   # TAA pos 1, stop codon
                             "A",   # GGT pos 1 -> AGT, Ser: nonsynonymous
                             "C",   # GGT pos 3 -> GGC, Gly: synonymous
                             NA),   # no allele at all
                     stringsAsFactors = FALSE)
  poly <- mk("poly_counts")(snps, map, hits)$g1
  expect_equal(poly[["pn"]], 2L)      # GAA pos 1, GGT pos 1
  expect_equal(poly[["ps"]], 3L)      # GAA pos 3, ATT pos 3, GGT pos 3
  expect_equal(poly[["dropped"]], 2L) # the stop codon and the site with no allele
})

test_that("an allele equal to the reference base is dropped, not counted", {
  map <- list(seq = c(g1 = "ATGGGT"), strand = c(g1 = "+"))
  hits <- data.frame(gene = "g1", snp_i = 1L, cds_idx = 6L)
  snps <- data.frame(chrom = "c", pos = 1L, ref = NA_character_, alt = "T",
                     stringsAsFactors = FALSE)   # T is already the base there
  poly <- mk("poly_counts")(snps, map, hits)$g1
  expect_equal(poly[["pn"]], 0L)
  expect_equal(poly[["ps"]], 0L)
  expect_equal(poly[["dropped"]], 1L)
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

# ---------- headerless comparison files (regression test for .check_mk_columns) ----------

test_that("a headerless comparison file with a non-MK sixth column stops before any work", {
  dir <- withr::local_tempdir()
  cmp <- file.path(dir, "comparison.txt")
  # Five proper columns, then a sixth that is not a SNP file. Read headerless, column 6
  # becomes snp_file, and there is no snp_side: the run must stop at once, not after
  # the first comparison's dN/dS.
  writeLines(c("cmp1 q.fa q.gff3 s.fa s.gff3 not_a_snp_file.txt"), cmp)
  err <- tryCatch(
    dndsR::calculate_dnds(comparison_file = cmp, output_dir = dir, aligner = "blast"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "snp_side is missing")
  expect_match(err, "column 6 is read as snp_file")
  expect_match(err, "Use 5 columns, or 7 with snp_file and snp_side")
  # Nothing was attempted: no output directory for a comparison was created.
  expect_false(file.exists(file.path(dir, "cmp1")))
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
  expect_equal(g04$ps, 1L)   # the chosen SNP is synonymous
  # The designed pair has 48 differing codon positions. Each is counted in exactly
  # one class, so the two columns must add up to 48 whatever the split is.
  expect_equal(g04$dn + g04$ds, 48)
  expect_gt(g04$dn, 0)
  expect_gt(g04$ds, 0)
  expect_equal(g04$codons_excluded_polymorphic, 1L)
  pooled <- utils::read.delim("example_AvB/example_AvB_mk_pooled.tsv")
  expect_equal(pooled$ps, 1L)
  expect_equal(pooled$pn, 0L)
  expect_gt(pooled$ds, 0)
})

# ---------- BED-driven polymorphism classification ----------

test_that("BED and VCF give the same counts for the same sites and alleles", {
  skip_if_not_installed("txdbmaker")
  g <- mini_genome()
  map <- mk("cds_map")(g$gff, g$fasta, g$cds_fa)
  # CDS pos 6 (GGT -> GGC) and pos 9 (GGC -> GGT): both synonymous, both Ps.
  pos <- map$pos[["gA"]][c(6, 9)]
  bed <- data.frame(chrom = c("chrM", "chrM"), pos = pos, ref = NA_character_,
                    alt = c("C", "T"), stringsAsFactors = FALSE)
  bedpoly <- mk("poly_counts")(bed, map, mk("map_snps")(bed, map)$hits)$gA
  expect_equal(bedpoly[["ps"]], 2L)
  expect_equal(bedpoly[["pn"]], 0L)
  vcf <- data.frame(chrom = c("chrM", "chrM"), pos = pos,
                    ref = substr(g$genome, pos, pos), alt = c("C", "T"),
                    stringsAsFactors = FALSE)
  vpoly <- mk("poly_counts")(vcf, map, mk("map_snps")(vcf, map)$hits)$gA
  expect_equal(vpoly[["ps"]], bedpoly[["ps"]])
  expect_equal(vpoly[["pn"]], bedpoly[["pn"]])
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
