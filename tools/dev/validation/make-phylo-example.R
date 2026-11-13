# Generates the bundled synthetic phylogenetic example: codon alignments simulated on an
# 8-taxon tree under a branch-site model, plus designed failure cases. Base R only, and
# independent of dndsR code. The sequences are synthetic and carry no biological meaning.
#   Rscript --vanilla tools/dev/validation/make-phylo-example.R inst/extdata/phylo_example
#
# Site classes follow the BUSTED generative form: each class has one omega on background
# branches and one on foreground branches. In "episodic" orthogroups, the third class has
# omega 6 on the foreground; in "null" orthogroups, foreground and background are equal.
# Branch lengths are expected substitutions per codon under the neutral (omega = 1) rate.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "inst/extdata/phylo_example"
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "codon-model.R")) # codons, sense, code, rate_matrix
RNGversion("4.0.0")
set.seed(20261113)
dir.create(file.path(out, "alignments"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "truth"), recursive = TRUE, showWarnings = FALSE)

# Tree as an edge list: parent, child, length. Root is "root".
# ((((spA,spB)nAB,(spC,spD)nCD)nABCD,(spE,spF)nEF)nAF,(spG,spH)nGH)root;
edges <- data.frame(
  parent = c("root", "root", "nAF", "nAF", "nABCD", "nABCD", "nAB", "nAB", "nCD", "nCD", "nEF", "nEF", "nGH", "nGH"),
  child  = c("nAF", "nGH", "nABCD", "nEF", "nAB", "nCD", "spA", "spB", "spC", "spD", "spE", "spF", "spG", "spH"),
  length = c(0.05, 0.10, 0.04, 0.08, 0.06, 0.05, 0.10, 0.12, 0.09, 0.11, 0.10, 0.09, 0.08, 0.10),
  stringsAsFactors = FALSE
)
taxa <- paste0("sp", LETTERS[1:8])
# Foreground hypothesis of the example: the clade of spA and spB, with its stem.
fg_children <- c("nAB", "spA", "spB")
newick <- "((((spA:0.10,spB:0.12):0.06,(spC:0.09,spD:0.11):0.05):0.04,(spE:0.10,spF:0.09):0.08):0.05,(spG:0.08,spH:0.10):0.10);"
writeLines(newick, file.path(out, "species_tree.nwk"))
writeLines(c("scope\ttaxa", "clade\tspA,spB"), file.path(out, "foreground.tsv"))

pi <- rep(1 / length(sense), length(sense))
kappa <- 2
neutral_scale <- {
  q <- rate_matrix(kappa, 1, pi); sum(sweep(q, 1, pi, "*"))
}
eig_cache <- list()
transition_probs <- function(omega, t) {
  key <- format(omega)
  if (is.null(eig_cache[[key]])) {
    q <- rate_matrix(kappa, omega, pi) / neutral_scale
    diag(q) <- -rowSums(q)
    e <- eigen(q)
    eig_cache[[key]] <<- list(v = e$vectors, vi = solve(e$vectors), d = e$values)
  }
  e <- eig_cache[[key]]
  p <- Re(e$v %*% diag(exp(e$d * t)) %*% e$vi)
  p[p < 0] <- 0
  p / rowSums(p)
}

simulate_og <- function(n_codons, classes) {
  cls <- sample.int(nrow(classes), n_codons, TRUE, prob = classes$prop)
  states <- list(root = sample.int(length(sense), n_codons, TRUE, prob = pi))
  for (i in seq_len(nrow(edges))) {
    e <- edges[i, ]
    fg <- e$child %in% fg_children
    from <- states[[e$parent]]
    to <- integer(n_codons)
    for (k in seq_len(nrow(classes))) {
      omega <- if (fg) classes$omega_fg[k] else classes$omega_bg[k]
      p <- transition_probs(omega, e$length)
      for (s in which(cls == k)) to[s] <- sample.int(length(sense), 1L, prob = p[from[s], ])
    }
    states[[e$child]] <- to
  }
  setNames(lapply(taxa, function(tx) paste(sense[states[[tx]]], collapse = "")), taxa)
}

write_fasta <- function(seqs, path, gene_suffix) {
  lines <- unlist(lapply(names(seqs), function(tx) c(paste0(">", tx, "|", tx, "_", gene_suffix), seqs[[tx]])))
  writeLines(lines, path)
}

null_classes <- data.frame(omega_bg = c(0.05, 0.5, 1), omega_fg = c(0.05, 0.5, 1), prop = c(0.55, 0.25, 0.20))
episodic_classes <- data.frame(omega_bg = c(0.05, 0.5, 1), omega_fg = c(0.05, 0.5, 6), prop = c(0.55, 0.25, 0.20))
n_codons <- 300L

design <- list()
add <- function(og, model, note, expect) {
  design[[length(design) + 1L]] <<- data.frame(orthogroup = og, model = model, note = note, expected_status = expect)
}

# 1. Episodic selection on the foreground, all taxa, with a terminal stop codon in every sequence.
s <- simulate_og(n_codons, episodic_classes)
s <- lapply(s, function(x) paste0(x, "TAA"))
write_fasta(s, file.path(out, "alignments", "OG0001.fasta"), "OG0001")
add("OG0001", "episodic", "all taxa; terminal TAA in every sequence", "tested")

# 2. No difference between foreground and background, all taxa.
s <- simulate_og(n_codons, null_classes)
write_fasta(s, file.path(out, "alignments", "OG0002.fasta"), "OG0002")
add("OG0002", "null", "all taxa", "tested")

# 3. Episodic, spH missing: the tree is pruned, and the foreground clade is intact.
s <- simulate_og(n_codons, episodic_classes)
s$spH <- NULL
write_fasta(s, file.path(out, "alignments", "OG0003.fasta"), "OG0003")
add("OG0003", "episodic", "spH missing; foreground clade intact", "tested")

# 4. Null, spB missing: the stem of {spA, spB} merges with spA's terminal branch. Under the
#    clade hypothesis every merged branch is foreground, so it stays identifiable.
s <- simulate_og(n_codons, null_classes)
s$spB <- NULL
write_fasta(s, file.path(out, "alignments", "OG0004.fasta"), "OG0004")
add("OG0004", "null", "spB missing; clade hypothesis still identifiable, stem would not be", "tested")

# 5. Both foreground taxa missing: the hypothesis cannot be tested.
s <- simulate_og(n_codons, null_classes)
s$spA <- NULL; s$spB <- NULL
write_fasta(s, file.path(out, "alignments", "OG0005.fasta"), "OG0005")
add("OG0005", "null", "spA and spB missing", "excluded:foreground_absent")

# 6. An internal stop codon in spC: excluded, not repaired.
s <- simulate_og(n_codons, null_classes)
substr(s$spC, 151, 153) <- "TGA"
write_fasta(s, file.path(out, "alignments", "OG0006.fasta"), "OG0006")
add("OG0006", "null", "internal TGA in spC at codon 51", "excluded:internal_stop")

# 7. Two sequences from spD: not single-copy.
s <- simulate_og(n_codons, null_classes)
seqs <- c(s, list(spD = s$spD))
lines <- unlist(lapply(seq_along(seqs), function(i) c(paste0(">", names(seqs)[i], "|", names(seqs)[i], "_OG0007_", i), seqs[[i]])))
writeLines(lines, file.path(out, "alignments", "OG0007.fasta"))
add("OG0007", "null", "two spD sequences", "excluded:multi_copy")

write.table(do.call(rbind, design), file.path(out, "truth", "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(edges, file.path(out, "truth", "tree_edges.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rbind(cbind(model = "null", null_classes), cbind(model = "episodic", episodic_classes)),
            file.path(out, "truth", "site_classes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
