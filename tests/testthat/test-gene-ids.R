# Gene and transcript names written by extract_cds are the GFF3 IDs (without a type
# prefix such as gene- or gene:), whatever the Name, gene_id, or transcript_id
# attributes hold, so append_annotations, which looks genes up by ID, finds every one.
# Before, Bioconductor named genes by Name when every Name was unique, so annotations
# were empty for such files. Expected names are written by hand from the GFF3 below.
skip_if_not_installed("txdbmaker")
skip_if_not_installed("Biostrings")

# A 1,200-base sequence holding two genes: g1 on + with two CDS exons (30-89, 120-269),
# g2 on - with one CDS (700-999). Codons are drawn without stops; each CDS starts with ATG
# and ends with TAA (written on the coding strand).
make_genome <- function(dir) {
  set.seed(7)
  sense <- setdiff(apply(expand.grid(c("A", "C", "G", "T"), c("A", "C", "G", "T"), c("A", "C", "G", "T")), 1,
                         paste0, collapse = ""), c("TAA", "TAG", "TGA"))
  orf <- function(n_codons) paste0("ATG", paste(sample(sense, n_codons - 2, TRUE), collapse = ""), "TAA")
  revcomp <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
  s <- strsplit(paste(sample(c("A", "C", "G", "T"), 1200, TRUE), collapse = ""), "")[[1]]
  g1 <- orf(70)                                # 210 nt: 60 in exon 1, 150 in exon 2
  s[30:89] <- strsplit(substr(g1, 1, 60), "")[[1]]
  s[120:269] <- strsplit(substr(g1, 61, 210), "")[[1]]
  g2 <- orf(100)                               # 300 nt on the minus strand
  s[700:999] <- strsplit(revcomp(g2), "")[[1]]
  fa <- file.path(dir, "genome.fasta")
  writeLines(c(">chr1", paste(s, collapse = "")), fa)
  list(fasta = fa, g1 = g1, g2 = g2)
}

# One GFF3 per attribute style. {G1} etc. are replaced by the gene names of the style.
styles <- list(
  # RefSeq-like, every Name unique and different from the ID: named by Name before.
  refseq_unique_names = list(
    lines = c(
      "chr1\tRefSeq\tgene\t30\t269\t.\t+\t.\tID=gene-LOC1;Name=ABC1;gene_biotype=protein_coding;Note=GO:0003677",
      "chr1\tRefSeq\tmRNA\t30\t269\t.\t+\t.\tID=rna-XM_1.1;Parent=gene-LOC1;Name=XM_1.1;transcript_id=XM_1.1",
      "chr1\tRefSeq\tCDS\t30\t89\t.\t+\t0\tID=cds-XP_1.1;Parent=rna-XM_1.1",
      "chr1\tRefSeq\tCDS\t120\t269\t.\t+\t0\tID=cds-XP_1.1;Parent=rna-XM_1.1",
      "chr1\tRefSeq\tgene\t700\t999\t.\t-\t.\tID=gene-LOC2;Name=XYZ2;gene_biotype=protein_coding;Note=GO:0005524",
      "chr1\tRefSeq\tmRNA\t700\t999\t.\t-\t.\tID=rna-XM_2.1;Parent=gene-LOC2;Name=XM_2.1;transcript_id=XM_2.1",
      "chr1\tRefSeq\tCDS\t700\t999\t.\t-\t0\tID=cds-XP_2.1;Parent=rna-XM_2.1"),
    genes = c("LOC1", "LOC2"), tx = c("XM_1.1", "XM_2.1")),
  # RefSeq-like with one Name shared by both genes: named by ID before and after.
  refseq_shared_name = list(
    lines = c(
      "chr1\tRefSeq\tgene\t30\t269\t.\t+\t.\tID=gene-LOC1;Name=ABC1;Note=GO:0003677",
      "chr1\tRefSeq\tmRNA\t30\t269\t.\t+\t.\tID=rna-XM_1.1;Parent=gene-LOC1;Name=XM_1.1",
      "chr1\tRefSeq\tCDS\t30\t89\t.\t+\t0\tID=cds-XP_1.1;Parent=rna-XM_1.1",
      "chr1\tRefSeq\tCDS\t120\t269\t.\t+\t0\tID=cds-XP_1.1;Parent=rna-XM_1.1",
      "chr1\tRefSeq\tgene\t700\t999\t.\t-\t.\tID=gene-LOC2;Name=ABC1;Note=GO:0005524",
      "chr1\tRefSeq\tmRNA\t700\t999\t.\t-\t.\tID=rna-XM_2.1;Parent=gene-LOC2;Name=XM_2.1",
      "chr1\tRefSeq\tCDS\t700\t999\t.\t-\t0\tID=cds-XP_2.1;Parent=rna-XM_2.1"),
    genes = c("LOC1", "LOC2"), tx = c("XM_1.1", "XM_2.1")),
  # Ensembl-like: type: prefixes, with gene_id and transcript_id attributes.
  ensembl = list(
    lines = c(
      "chr1\tsrc\tgene\t30\t269\t.\t+\t.\tID=gene:AT1G01010;Name=NAC001;gene_id=AT1G01010;Note=GO:0003677",
      "chr1\tsrc\tmRNA\t30\t269\t.\t+\t.\tID=transcript:AT1G01010.1;Parent=gene:AT1G01010;transcript_id=AT1G01010.1",
      "chr1\tsrc\tCDS\t30\t89\t.\t+\t0\tID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1",
      "chr1\tsrc\tCDS\t120\t269\t.\t+\t0\tID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1",
      "chr1\tsrc\tgene\t700\t999\t.\t-\t.\tID=gene:AT1G01020;Name=ARV1;gene_id=AT1G01020;Note=GO:0005524",
      "chr1\tsrc\tmRNA\t700\t999\t.\t-\t.\tID=transcript:AT1G01020.1;Parent=gene:AT1G01020;transcript_id=AT1G01020.1",
      "chr1\tsrc\tCDS\t700\t999\t.\t-\t0\tID=CDS:AT1G01020.1;Parent=transcript:AT1G01020.1;protein_id=AT1G01020.1"),
    genes = c("AT1G01010", "AT1G01020"), tx = c("AT1G01010.1", "AT1G01020.1")),
  # Gene-predictor style (BRAKER, Helixer): IDs only.
  ids_only = list(
    lines = c(
      "chr1\tpred\tgene\t30\t269\t.\t+\t.\tID=g1;Note=GO:0003677",
      "chr1\tpred\tmRNA\t30\t269\t.\t+\t.\tID=g1.t1;Parent=g1",
      "chr1\tpred\tCDS\t30\t89\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1",
      "chr1\tpred\tCDS\t120\t269\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1",
      "chr1\tpred\tgene\t700\t999\t.\t-\t.\tID=g2;Note=GO:0005524",
      "chr1\tpred\tmRNA\t700\t999\t.\t-\t.\tID=g2.t1;Parent=g2",
      "chr1\tpred\tCDS\t700\t999\t.\t-\t0\tID=g2.t1.cds;Parent=g2.t1"),
    genes = c("g1", "g2"), tx = c("g1.t1", "g2.t1"))
)

run_style <- function(style, group_by) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  g <- make_genome(dir)
  gff <- file.path(dir, "genome.gff3")
  writeLines(c("##gff-version 3", style$lines), gff)
  suppressMessages(suppressWarnings(
    extract_cds(comparison_name = "c", query_fasta = g$fasta, subject_fasta = g$fasta,
                query_gff = gff, subject_gff = gff, output_dir = dir, group_by = group_by)))
  cds <- Biostrings::readDNAStringSet(file.path(dir, "c", "genome_CDS.fasta"))
  list(dir = dir, gff = gff, genome = g, cds = cds)
}

for (nm in names(styles)) {
  test_that(paste0("genes and transcripts are named by their GFF3 ID: ", nm), {
    st <- styles[[nm]]
    r <- run_style(st, "gene")
    expect_setequal(names(r$cds), st$genes)
    # The sequences are the designed CDS, so naming did not change what was extracted.
    expect_identical(as.character(r$cds[[st$genes[1]]]), r$genome$g1)
    expect_identical(as.character(r$cds[[st$genes[2]]]), r$genome$g2)
    rt <- run_style(st, "tx")
    expect_setequal(names(rt$cds), st$tx)
  })

  test_that(paste0("append_annotations finds every extracted gene: ", nm), {
    st <- styles[[nm]]
    r <- run_style(st, "gene")
    # The names extract_cds wrote, as a dN/dS table would carry them.
    tab <- data.frame(query_id = names(r$cds), subject_id = names(r$cds), dN = 0.01, dS = 0.1, dNdS = 0.1)
    f <- file.path(r$dir, "dnds.tsv")
    utils::write.table(tab, f, sep = "\t", quote = FALSE, row.names = FALSE)
    msgs <- character(0)
    out <- withCallingHandlers(
      append_annotations(dnds_file = f, query_gff = r$gff, subject_gff = r$gff),
      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
    expect_false(any(grepl("not found in its GFF", msgs)))
    expect_false(anyNA(out$q_gff_start))
    expect_identical(out$q_gff_start[match(st$genes, out$query_id)], c(30, 700))
    expect_match(out$q_attributes[match(st$genes[1], out$query_id)], "GO:0003677")
  })
}

# ---- Names DIAMOND and orthologr cannot use (docs/validation/known-issues.md) ----

test_that("unusable sequence names are recognised, and usable ones are not", {
  b <- dndsR:::.unsafe_seq_names(c("g 1", "g\t2", "gnl|JCVI|m.1", "ref|NM_1|", "gi|1|ref|g", "sp|P1|g",
                                   "lcl|g1", "x|g1", "AT1G01010", "g;1", "gene:A", NA))
  expect_setequal(b$whitespace, c("g 1", "g\t2"))
  expect_setequal(b$ncbi_tag, c("gnl|JCVI|m.1", "ref|NM_1|", "gi|1|ref|g", "sp|P1|g"))
  expect_length(dndsR:::.unsafe_seq_names("g 1", whitespace = FALSE)$whitespace, 0)
})

pipe_gff <- function(gene1, gene2, tx1, tx2) c(
  sprintf("chr1\tp\tgene\t30\t269\t.\t+\t.\tID=%s", gene1),
  sprintf("chr1\tp\tmRNA\t30\t269\t.\t+\t.\tID=%s;Parent=%s", tx1, gene1),
  sprintf("chr1\tp\tCDS\t30\t89\t.\t+\t0\tParent=%s", tx1), sprintf("chr1\tp\tCDS\t120\t269\t.\t+\t0\tParent=%s", tx1),
  sprintf("chr1\tp\tgene\t700\t999\t.\t-\t.\tID=%s", gene2),
  sprintf("chr1\tp\tmRNA\t700\t999\t.\t-\t.\tID=%s;Parent=%s", tx2, gene2),
  sprintf("chr1\tp\tCDS\t700\t999\t.\t-\t0\tParent=%s", tx2))

test_that("extract_cds stops on gene IDs with whitespace, naming them", {
  st <- list(lines = pipe_gff("gene A1", "gene A2", "t A1", "t A2"))
  expect_error(run_style(st, "gene"), "contain whitespace, e.g. 'gene A1', 'gene A2'")
})

test_that("NCBI-tagged transcript IDs stop transcript mode, and gene mode still works", {
  st <- list(lines = pipe_gff("gene-A1", "gene-A2", "rna-gnl|JCVI|mRNA.A1.1", "rna-gnl|JCVI|mRNA.A2.1"))
  expect_error(run_style(st, "tx"), "NCBI database tag.*'gnl\\|JCVI\\|mRNA.A1.1'.*group_by = \"gene\"")
  expect_setequal(names(run_style(st, "gene")$cds), c("A1", "A2"))
})

test_that("calculate_dnds stops on NCBI-tagged names before any alignment", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "q.fasta")
  writeLines(c(">ref|g1| a description", "ATGAAATAA", ">ref|g2|", "ATGCCCTAA"), f)
  expect_error(calculate_dnds(comparison_name = "c", query_seq = f, subject_seq = f, output_dir = dir, aligner = "blast"),
               "start with an NCBI database tag.*'ref\\|g1\\|'")
})
