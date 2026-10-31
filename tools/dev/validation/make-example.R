# Generates the bundled synthetic example: two small genomes with a designed set of
# substitutions and gene structures. Base R only, and independent of dndsR code.
#   Rscript --vanilla tools/dev/validation/make-example.R inst/extdata/example
# The sequences are synthetic and carry no biological meaning.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "inst/extdata/example"
dir.create(file.path(out, "truth"), recursive = TRUE, showWarnings = FALSE)
RNGversion("4.0.0")
set.seed(20261031)

bases <- c("T", "C", "A", "G")
codons <- apply(expand.grid(bases, bases, bases, stringsAsFactors = FALSE)[, 3:1], 1, paste0, collapse = "")
aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
code <- setNames(aa, codons)
sense <- codons[code != "*"]

random_cds <- function(n) c("ATG", sample(sense, n, TRUE), "TAA")
revcomp <- function(s) paste(rev(chartr("ACGT", "TGCA", strsplit(s, "")[[1]])), collapse = "")
random_dna <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")

# All single-nucleotide neighbours of a codon that remain sense codons.
neighbours <- function(codon) {
  res <- character()
  for (p in 1:3) for (b in setdiff(bases, substr(codon, p, p))) {
    n <- codon
    substr(n, p, p) <- b
    if (code[[n]] != "*") res <- c(res, n)
  }
  res
}

# Applies exact counts of one-step synonymous and nonsynonymous changes, one per codon.
mutate <- function(cds, n_syn, n_non) {
  idx <- sample(2:(length(cds) - 1))
  for (kind in c("syn", "non")) {
    need <- if (kind == "syn") n_syn else n_non
    while (need > 0) {
      i <- idx[1]
      idx <- idx[-1]
      opts <- neighbours(cds[i])
      keep <- if (kind == "syn") code[opts] == code[[cds[i]]] else code[opts] != code[[cds[i]]]
      opts <- opts[keep]
      if (!length(opts)) next
      cds[i] <- opts[sample.int(length(opts), 1)]
      need <- need - 1
    }
  }
  cds
}

# Designed classes: name, synonymous changes, nonsynonymous changes, term annotation.
# Terms are real identifiers assigned arbitrarily, so that lookups work; they mean nothing here.
n_codons <- 300
spec <- rbind(
  data.frame(id = "g01", class = "identical", syn = 0, non = 0),
  data.frame(id = "g02", class = "synonymous_only", syn = 30, non = 0),
  data.frame(id = "g03", class = "nonsynonymous_only", syn = 0, non = 30),
  data.frame(id = sprintf("g%02d", 4:15), class = "conserved", syn = 40, non = 8),
  data.frame(id = sprintf("g%02d", 16:19), class = "neutral_like", syn = 12, non = 36),
  data.frame(id = sprintf("g%02d", 20:25), class = "excess_nonsynonymous", syn = 4, non = 60),
  data.frame(id = sprintf("g%02d", 26:28), class = "structure", syn = 30, non = 10),
  data.frame(id = c("g29", "g30", "g35", "g36"), class = "multi_isoform", syn = 30, non = 10)
)
spec$term <- ifelse(spec$class == "excess_nonsynonymous", "InterPro:IPR000001",
                    ifelse(spec$class == "conserved", "InterPro:IPR000003", "InterPro:IPR000007"))
spec$go <- ifelse(spec$class == "excess_nonsynonymous", "GO:0006952", "GO:0008150")
spec$kegg <- ifelse(spec$class == "excess_nonsynonymous", "KEGG:K00001", "KEGG:K00002")

genes <- list()
truth <- list()
for (i in seq_len(nrow(spec))) {
  a <- random_cds(n_codons)
  b <- mutate(a, spec$syn[i], spec$non[i])
  # A second subject for contrasts: conserved at every gene, whatever its class against B.
  genes[[spec$id[i]]] <- list(a = a, b = b, c = mutate(a, 40, 8))
  truth[[spec$id[i]]] <- data.frame(pair = spec$id[i], class = spec$class[i], syn_changes = spec$syn[i],
                                    nonsyn_changes = spec$non[i], note = "")
}

# Problem genes that a correct pipeline should exclude or flag.
stop_a <- random_cds(n_codons)
stop_b <- mutate(stop_a, 20, 5)
stop_b[150] <- "TGA"
genes$g31 <- list(a = stop_a, b = stop_b, c = mutate(stop_a, 40, 8))
truth$g31 <- data.frame(pair = "g31", class = "internal_stop_subject", syn_changes = 20, nonsyn_changes = 5,
                        note = "subject has an in-frame stop at codon 150")
frame_a <- random_cds(n_codons)
frame_b <- mutate(frame_a, 20, 5)
genes$g32 <- list(a = c(frame_a, "A"), b = c(frame_b, "A"), c = c(mutate(frame_a, 40, 8), "A"))
truth$g32 <- data.frame(pair = "g32", class = "length_not_multiple_of_3", syn_changes = 20, nonsyn_changes = 5,
                        note = "both CDS have one trailing extra base")
only_a <- random_cds(n_codons)
only_b <- random_cds(n_codons)

structure_of <- function(id) {
  switch(id,
    g26 = list(strand = "-", cuts = NULL, alt = NULL),
    g27 = list(strand = "+", cuts = c(301L, 601L), alt = NULL),
    g28 = list(strand = "-", cuts = c(301L, 601L), alt = NULL),
    g29 = list(strand = "+", cuts = c(300L, 600L), alt = "skip_exon2"),
    g30 = list(strand = "+", cuts = c(300L, 600L), alt = "extend_exon2"),
    g35 = list(strand = "+", cuts = c(300L, 600L), alt = "shorten_exon2"),
    g36 = list(strand = "+", cuts = c(300L, 600L), alt = "extend_exon2_with_stop"),
    list(strand = if (as.integer(sub("g", "", id)) %% 2L == 0L) "+" else "-", cuts = NULL, alt = NULL)
  )
}

# Lays one species out on two chromosomes and writes FASTA, GFF3, and truth CDS.
write_species <- function(prefix, side, ids, extra) {
  chroms <- list(Chr1 = list(seq = random_dna(300), len = 300L), Chr2 = list(seq = random_dna(300), len = 300L))
  gff <- "##gff-version 3"
  cds_out <- character()
  parts <- list(Chr1 = list(chroms$Chr1$seq), Chr2 = list(chroms$Chr2$seq))
  lens <- c(Chr1 = 300L, Chr2 = 300L)

  for (k in seq_along(ids)) {
    id <- ids[k]
    codon_vec <- if (id %in% names(extra)) extra[[id]] else genes[[id]][[side]]
    st <- structure_of(id)
    full <- paste(codon_vec, collapse = "")
    cuts <- if (is.null(st$cuts)) integer() else st$cuts # in nucleotides; g27 and g28 split codons
    bounds <- c(0L, cuts, nchar(full))
    exons <- lapply(seq_len(length(bounds) - 1L), function(j) substr(full, bounds[j] + 1L, bounds[j + 1L]))
    # A 45-nt in-frame extension directly upstream of exon 2, used by alternative isoforms.
    ext <- switch(if (is.null(st$alt)) "none" else st$alt,
      extend_exon2 = ext_coding,
      extend_exon2_with_stop = ext_stop,
      ""
    )
    introns <- if (length(exons) > 1L) c(ext, rep("", length(exons) - 2L)) else character()
    gene_plus <- exons[[1]]
    exon_pos <- list(c(1L, nchar(exons[[1]])))
    for (j in seq_along(introns)) {
      intron <- paste0("GT", random_dna(150L), "AG", introns[j])
      start <- nchar(gene_plus) + nchar(intron) + 1L
      gene_plus <- paste0(gene_plus, intron, exons[[j + 1L]])
      exon_pos[[j + 1L]] <- c(start, start + nchar(exons[[j + 1L]]) - 1L)
    }
    glen <- nchar(gene_plus)
    chrom <- if (k <= ceiling(length(ids) / 2)) "Chr1" else "Chr2"
    offset <- lens[[chrom]]
    genomic <- if (st$strand == "+") gene_plus else revcomp(gene_plus)
    spacer <- random_dna(sample(200:400, 1))
    parts[[chrom]] <- c(parts[[chrom]], genomic, spacer)
    lens[[chrom]] <- offset + glen + nchar(spacer)

    # Converts gene-relative coordinates to chromosome coordinates on either strand.
    to_chr <- function(p) {
      if (st$strand == "+") c(offset + p[1], offset + p[2]) else c(offset + glen - p[2] + 1L, offset + glen - p[1] + 1L)
    }
    gid <- paste0(prefix, "_", id)
    s <- spec[spec$id == id, ]
    attrs <- paste0("ID=", gid, ";Name=", gid)
    if (nrow(s)) attrs <- paste0(attrs, ";Dbxref=", s$term, ",", s$kegg, ";Ontology_term=", s$go)
    g <- to_chr(c(1L, glen))
    line <- function(type, r, at, phase = ".") {
      paste(chrom, "synthetic", type, r[1], r[2], ".", st$strand, phase, at, sep = "\t")
    }
    gff <- c(gff, line("gene", g, attrs))

    write_tx <- function(tx, pos) {
      tid <- paste0(gid, ".", tx)
      gff <<- c(gff, line("mRNA", g, paste0("ID=", tid, ";Parent=", gid)))
      emit <- if (st$strand == "+") seq_along(pos) else rev(seq_along(pos))
      done <- 0L
      phases <- integer(length(pos))
      for (j in seq_along(pos)) {
        phases[j] <- (3L - done %% 3L) %% 3L
        done <- done + pos[[j]][2] - pos[[j]][1] + 1L
      }
      for (j in emit) {
        r <- to_chr(pos[[j]])
        gff <<- c(gff, line("exon", r, paste0("ID=", tid, ".exon", j, ";Parent=", tid)))
        gff <<- c(gff, line("CDS", r, paste0("ID=", tid, ".cds", ";Parent=", tid), phases[j]))
      }
    }
    write_tx(1L, exon_pos)
    if (identical(st$alt, "skip_exon2")) write_tx(2L, exon_pos[c(1L, 3L)])
    if (!is.null(st$alt) && startsWith(st$alt, "extend_exon2")) {
      p2 <- exon_pos
      p2[[2]][1] <- p2[[2]][1] - 45L
      write_tx(2L, p2)
    }
    if (identical(st$alt, "shorten_exon2")) {
      p2 <- exon_pos
      p2[[2]][2] <- p2[[2]][2] - 45L
      write_tx(2L, p2)
    }
    cds_out <- c(cds_out, paste0(">", gid), full)
  }

  fasta <- unlist(lapply(names(parts), function(ch) {
    s <- paste(unlist(parts[[ch]]), collapse = "")
    n <- nchar(s)
    starts <- seq(1L, n, by = 70L)
    c(paste0(">", ch), substring(s, starts, pmin(starts + 69L, n)))
  }))
  writeLines(fasta, file.path(out, paste0(prefix, ".fasta")))
  writeLines(gff, file.path(out, paste0(prefix, ".gff3")))
  writeLines(cds_out, file.path(out, "truth", paste0(prefix, "_longest_isoform_cds.fasta")))
}

ext_coding <- paste(sample(sense, 15, TRUE), collapse = "")
ext_stop <- paste0(paste(sample(sense, 7, TRUE), collapse = ""), "TAA", paste(sample(sense, 7, TRUE), collapse = ""))
ids_a <- c(names(genes), "g33")
ids_b <- c(names(genes), "g34")
write_species("spA", "a", ids_a, list(g33 = only_a))
write_species("spB", "b", ids_b, list(g34 = only_b))
write_species("spC", "c", names(genes), list())

truth_tbl <- do.call(rbind, truth)
rownames(truth_tbl) <- NULL
write.table(truth_tbl, file.path(out, "truth", "designed_substitutions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines("example_AvB spA.fasta spA.gff3 spB.fasta spB.gff3", file.path(out, "comparison_file.txt"))
writeLines(c("example_AvB spA.fasta spA.gff3 spB.fasta spB.gff3", "example_AvC spA.fasta spA.gff3 spC.fasta spC.gff3"),
           file.path(out, "two_comparisons.txt"))
writeLines("AvB_vs_AvC example_AvB query example_AvC query", file.path(out, "contrast_file.txt"))
writeLines("Chr1 0 6000 first_6kb", file.path(out, "regions.bed"))
cat("Wrote example to", out, "\n")
