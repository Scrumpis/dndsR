# Generates a synthetic six-genome example for the phylogenetic workflow: genome FASTA and
# GFF3 files, orthogroup tables in dndsR's long format and OrthoFinder's format, a species
# tree, a foreground, gene trees, and the designed truth. Base R only, independent of dndsR.
# The sequences are synthetic and carry no biological meaning.
#   Rscript --vanilla tools/dev/validation/make-phylo-genomes.R inst/extdata/phylo_genomes
#
# Coding sequences are simulated on the species tree under a BUSTED-form site-class model
# (see make-phylo-example.R). Each orthogroup is designed to test one behavior; see
# truth/design.tsv.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "inst/extdata/phylo_genomes"
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "codon-model.R")) # sense, code, rate_matrix
RNGversion("4.0.0")
set.seed(20261114)
for (d in c("", "gene_trees", "truth")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

taxa <- paste0("sp", LETTERS[1:6])
# ((((spA,spB)nAB,spC)nABC,(spD,spE)nDE)nAE,spF)root;
edges <- data.frame(
  parent = c("root", "root", "nAE", "nAE", "nABC", "nABC", "nAB", "nAB", "nDE", "nDE"),
  child  = c("nAE", "spF", "nABC", "nDE", "nAB", "spC", "spA", "spB", "spD", "spE"),
  length = c(0.04, 0.18, 0.05, 0.07, 0.06, 0.12, 0.10, 0.12, 0.10, 0.09),
  stringsAsFactors = FALSE
)
fg_children <- c("nAB", "spA", "spB")
writeLines("((((spA,spB),spC),(spD,spE)),spF);", file.path(out, "species_tree.nwk"))
writeLines(c("scope\ttaxa", "clade\tspA,spB"), file.path(out, "foreground.tsv"))

pi <- rep(1 / length(sense), length(sense))
kappa <- 2
neutral_scale <- { q <- rate_matrix(kappa, 1, pi); sum(sweep(q, 1, pi, "*")) }
eig <- list()
pmat <- function(omega, t) {
  k <- format(omega)
  if (is.null(eig[[k]])) {
    q <- rate_matrix(kappa, omega, pi) / neutral_scale
    diag(q) <- -rowSums(q)
    e <- eigen(q)
    eig[[k]] <<- list(v = e$vectors, vi = solve(e$vectors), d = e$values)
  }
  e <- eig[[k]]
  p <- Re(e$v %*% diag(exp(e$d * t)) %*% e$vi)
  p[p < 0] <- 0
  p / rowSums(p)
}
# Returns a codon matrix: one row per taxon, one column per site. (`n` belongs to codon-model.R.)
simulate <- function(k, classes) {
  cls <- sample.int(nrow(classes), k, TRUE, prob = classes$prop)
  st <- list(root = sample.int(length(sense), k, TRUE))
  for (i in seq_len(nrow(edges))) {
    e <- edges[i, ]
    fg <- e$child %in% fg_children
    to <- integer(k)
    for (cl in seq_len(nrow(classes))) {
      p <- pmat(if (fg) classes$omega_fg[cl] else classes$omega_bg[cl], e$length)
      for (s in which(cls == cl)) to[s] <- sample.int(length(sense), 1L, prob = p[st[[e$parent]][s], ])
    }
    st[[e$child]] <- to
  }
  t(vapply(taxa, function(tx) sense[st[[tx]]], character(k)))
}
null_classes <- data.frame(omega_bg = c(0.05, 0.5, 1), omega_fg = c(0.05, 0.5, 1), prop = c(0.55, 0.25, 0.20))
episodic_classes <- data.frame(omega_bg = c(0.05, 0.5, 1), omega_fg = c(0.05, 0.5, 8), prop = c(0.55, 0.25, 0.20))
stops <- c("TAA", "TAG", "TGA")

genes <- list()      # one entry per gene: taxon, og, cds (with ATG and stop)
design <- list()
true_aln <- list()
add_og <- function(og, cds_by_taxon, note, prepare, busted, basis = "") {
  for (tx in names(cds_by_taxon)) for (s in cds_by_taxon[[tx]]) genes[[length(genes) + 1L]] <<- list(taxon = tx, og = og, cds = s)
  design[[length(design) + 1L]] <<- data.frame(orthogroup = og, note = note, expected_prepare = prepare,
                                               expected_basis = basis, expected_busted = busted)
}
as_cds <- function(codon_row) paste0("ATG", paste(codon_row, collapse = ""), sample(stops, 1L))
og_from <- function(m, keep = taxa) setNames(lapply(keep, function(tx) as_cds(m[tx, ])), keep)

len <- 200L
for (i in 1:3) add_og(sprintf("OG%02d", i), og_from(simulate(len, episodic_classes)), "episodic on the foreground", "prepared", "tested")
for (i in 4:6) add_og(sprintf("OG%02d", i), og_from(simulate(len, null_classes)), "no foreground difference; gene trees for 04-06", "prepared", "tested")

m <- simulate(len, null_classes)
cds <- og_from(m)
dup <- m["spD", ]
change <- sample.int(len, 10L)
dup[change] <- sample(sense, 10L, TRUE)
cds$spD <- c(cds$spD, as_cds(dup))
add_og("OG07", cds, "two spD copies (a paralog)", "excluded:multi_copy", "excluded:multi_copy", "policy")

m <- simulate(len, null_classes)
cds <- og_from(m)
substr(cds$spC, 3L * 120L + 1L, 3L * 120L + 3L) <- "TGA" # codon 121 counting the ATG
add_og("OG08", cds, "internal TGA in spC's gene model", "excluded:cds_invalid", "excluded:cds_invalid", "data_quality")

add_og("OG09", og_from(simulate(len, null_classes), setdiff(taxa, "spB")), "no spB gene; clade hypothesis still testable",
       "prepared", "tested")

# Codon indels with a known alignment: 2 codons inserted in spA after site 50; sites 100-102
# deleted on the branch to spD and spE.
m <- simulate(len, null_classes)
ins <- sample(sense, 2L, TRUE)
aligned <- lapply(taxa, function(tx) {
  row <- m[tx, ]
  if (tx %in% c("spD", "spE")) row[100:102] <- "---"
  c(row[1:50], if (tx == "spA") ins else c("---", "---"), row[51:len])
})
names(aligned) <- taxa
true_aln$OG10 <- aligned
cds <- setNames(lapply(taxa, function(tx) { r <- aligned[[tx]]; as_cds(r[r != "---"]) }), taxa)
add_og("OG10", cds, "codon insertion in spA and deletion in spD and spE; true alignment in truth/", "prepared", "tested")

add_og("OG11", og_from(simulate(len, null_classes), setdiff(taxa, c("spA", "spB"))), "no spA or spB gene",
       "excluded:foreground_absent", "excluded:foreground_absent", "no_data")
add_og("OG12", og_from(simulate(len, null_classes), setdiff(taxa, "spC")),
       "no spC gene: the clade's stem joins the background branch above it; spA and spB stay separable",
       "excluded:foreground_partly_unidentifiable", "excluded:foreground_partly_unidentifiable", "policy")
add_og("OG14", og_from(simulate(len, null_classes), setdiff(taxa, c("spB", "spC"))),
       "no spB or spC gene: spA's branch joins the stem and the background branch above; nothing separable",
       "excluded:foreground_unidentifiable", "excluded:foreground_unidentifiable", "mathematical")
m <- simulate(len, null_classes)
cds <- og_from(m)
cds$spE <- paste0(substr(cds$spE, 1L, 299L), substr(cds$spE, 301L, nchar(cds$spE))) # 1-nt deletion
add_og("OG13", cds, "1-nt deletion in spE's gene model (frameshift)", "excluded:cds_invalid", "excluded:cds_invalid", "data_quality")

# Gene IDs that do not reveal orthogroups, so the tables are needed.
gene_tab <- do.call(rbind, lapply(genes, function(g) data.frame(taxon = g$taxon, orthogroup = g$og, cds = g$cds)))
gene_tab <- gene_tab[sample.int(nrow(gene_tab)), ]
gene_tab$gene_id <- ave(gene_tab$taxon, gene_tab$taxon, FUN = function(x) sprintf("%s%04d", x, seq_along(x) * 7L))

revcomp <- function(s) paste(rev(strsplit(chartr("ACGT", "TGCA", s), "")[[1]]), collapse = "")
random_dna <- function(k) paste(sample(c("A", "C", "G", "T"), k, TRUE), collapse = "")
for (tx in taxa) {
  g <- gene_tab[gene_tab$taxon == tx, ]
  g <- g[order(g$gene_id), ]
  chrom <- random_dna(300L)
  gff <- "##gff-version 3"
  for (i in seq_len(nrow(g))) {
    id <- g$gene_id[i]; cds <- g$cds[i]
    minus <- i %% 2L == 0L
    two_exons <- i %% 3L == 0L
    cut <- 151L # first exon length when split
    intron <- paste0("GT", random_dna(96L), "AG")
    plus <- if (two_exons) paste0(substr(cds, 1L, cut), intron, substr(cds, cut + 1L, nchar(cds))) else cds
    start <- nchar(chrom) + 1L
    chrom <- paste0(chrom, if (minus) revcomp(plus) else plus, random_dna(250L))
    end <- start + nchar(plus) - 1L
    strand <- if (minus) "-" else "+"
    line <- function(type, a, b, attrs, phase = ".") paste("Chr1", "synthetic", type, a, b, ".", strand, phase, attrs, sep = "\t")
    gff <- c(gff, line("gene", start, end, paste0("ID=", id)), line("mRNA", start, end, paste0("ID=", id, ".1;Parent=", id)))
    if (!two_exons) {
      parts <- list(c(start, end), c(start, end))
      exons <- list(c(start, end))
    } else if (!minus) {
      exons <- list(c(start, start + cut - 1L), c(start + cut + 100L, end))
    } else {
      # On the minus strand the first coding exon is at the right end.
      exons <- list(c(end - cut + 1L, end), c(start, end - cut - 100L))
    }
    done <- 0L
    for (j in seq_along(exons)) {
      r <- exons[[j]]
      phase <- (3L - done %% 3L) %% 3L
      gff <- c(gff, line("exon", r[1], r[2], sprintf("ID=%s.1.exon%d;Parent=%s.1", id, j, id)),
               line("CDS", r[1], r[2], sprintf("ID=%s.1.cds;Parent=%s.1", id, id), phase))
      done <- done + (r[2] - r[1] + 1L)
    }
  }
  writeLines(c(paste0(">Chr1"), substring(chrom, seq(1L, nchar(chrom), 80L), pmin(seq(80L, nchar(chrom) + 79L, 80L), nchar(chrom)))),
             file.path(out, paste0(tx, ".fasta")))
  writeLines(gff, file.path(out, paste0(tx, ".gff3")))
}

writeLines(c("taxon\tfasta\tgff\torthofinder_name", sprintf("%s\t%s.fasta\t%s.gff3\t%s_proteins", taxa, taxa, taxa, taxa)),
           file.path(out, "genomes.tsv"))
long <- gene_tab[order(gene_tab$orthogroup, gene_tab$taxon, gene_tab$gene_id), c("orthogroup", "taxon", "gene_id")]
write.table(long, file.path(out, "orthogroup_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# The same orthogroups as OrthoFinder writes them: one column per proteome, genes joined by ", ".
ogs <- sort(unique(long$orthogroup))
wide <- data.frame(Orthogroup = ogs, check.names = FALSE)
for (tx in taxa) {
  wide[[paste0(tx, "_proteins")]] <- vapply(ogs, function(o) paste(long$gene_id[long$orthogroup == o & long$taxon == tx], collapse = ", "), character(1))
}
write.table(wide, file.path(out, "Orthogroups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Gene trees, tips named as OrthoFinder names them: <proteome>_<gene id>.
tip <- function(og, tx) paste0(tx, "_proteins_", long$gene_id[long$orthogroup == og & long$taxon == tx])
gt <- function(og, pattern) {
  s <- pattern
  for (tx in taxa) s <- gsub(paste0("\\b", tx, "\\b"), tip(og, tx), s, perl = TRUE)
  writeLines(s, file.path(out, "gene_trees", paste0(og, ".nwk")))
}
gt("OG04", "((((spA,spB),spC),(spD,spF)),spE);")   # discordant away from the foreground
gt("OG05", "(((spA,spB),(spC,spD)),spE,spF);")     # discordant beside the foreground stem
gt("OG06", "((((spA,spC),spB),(spD,spE)),spF);")   # no spA+spB clade: the foreground is absent

gene_tree_design <- data.frame(
  orthogroup = c("OG04", "OG05", "OG06"),
  expected_gene_tree = c("tested; discordant, not near the foreground", "tested; discordant near the foreground",
                         "excluded:foreground_not_in_gene_tree"),
  missing_species_splits = c(1L, 2L, 1L)
)
write.table(gene_tree_design, file.path(out, "truth", "gene_tree_design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, design), file.path(out, "truth", "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(unlist(lapply(taxa, function(tx) c(paste0(">", tx), paste(true_aln$OG10[[tx]], collapse = "")))),
           file.path(out, "truth", "OG10_true_codon_alignment.fasta"))
truth_cds <- gene_tab[order(gene_tab$taxon, gene_tab$gene_id), ]
writeLines(as.vector(rbind(paste0(">", truth_cds$gene_id), truth_cds$cds)), file.path(out, "truth", "designed_cds.fasta"))
write.table(edges, file.path(out, "truth", "tree_edges.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
