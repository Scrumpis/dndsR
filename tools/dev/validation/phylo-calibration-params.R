# Derives simulation parameters from the real-genome benchmark, by the rules fixed in
# docs/dev/phylo-calibration-design.md ("Parameters from the real-genome benchmark").
#   Rscript --vanilla phylo-calibration-params.R <p3a results dir> <params.json> [c]
# Reads busted/busted_results.tsv, prepared.tar.gz (alignments), and busted_hyphy.tar.gz
# (HyPhy JSON). c (default 1) converts HyPhy branch lengths to the simulator's units; it is
# set afterwards by phylo-calibration-scale.sb, which looks only at tree lengths.
a <- commandArgs(TRUE); p3a <- a[1]; out <- a[2]; c_factor <- if (length(a) > 2) as.numeric(a[3]) else 1
tmp <- tempfile("p3a-"); dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
res <- utils::read.delim(file.path(p3a, "busted", "busted_results.tsv"), na.strings = "NA")
full <- res[res$status == "tested" & res$n_taxa == 5, ]
q <- round(stats::quantile(full$n_codons, c(0.25, 0.5, 0.75)) / 10) * 10

# F3x4 from all prepared alignments, pooled.
utils::untar(file.path(p3a, "prepared.tar.gz"), exdir = tmp, extras = "--wildcards", files = "prepared/alignments/*")
counts <- matrix(0, 3, 4, dimnames = list(NULL, c("T", "C", "A", "G")))
for (f in list.files(file.path(tmp, "prepared", "alignments"), full.names = TRUE)) {
  l <- readLines(f); s <- l[!startsWith(l, ">")]
  for (x in s) {
    ch <- strsplit(x, "")[[1]]
    for (p in 1:3) {
      b <- ch[seq(p, length(ch), 3)]
      counts[p, ] <- counts[p, ] + table(factor(b, levels = colnames(counts)))
    }
  }
}
nuc <- counts / rowSums(counts)

# Median unconstrained branch length of each branch, over five-taxon tested orthogroups.
utils::untar(file.path(p3a, "busted_hyphy.tar.gz"), exdir = tmp)
lens <- do.call(rbind, lapply(full$orthogroup, function(og) {
  j <- jsonlite::read_json(file.path(tmp, "busted", "hyphy", og, "busted.json"))
  b <- j[["branch attributes"]][["0"]]
  tree <- readLines(file.path(tmp, "busted", "hyphy", og, "tree.nwk"))
  # Name internal branches by the taxa they separate, so names do not depend on HyPhy's.
  t <- ape::read.tree(text = gsub("\\{FG\\}", "", tree))
  desc <- c(as.list(t$tip.label), lapply(ape::prop.part(t), function(p) t$tip.label[p]))
  nm <- c(t$tip.label, t$node.label)
  key <- vapply(seq_along(nm), function(i) paste(sort(desc[[i]]), collapse = ","), character(1))
  vals <- vapply(nm, function(n) if (is.null(b[[n]]$unconstrained)) NA_real_ else as.numeric(b[[n]]$unconstrained), numeric(1))
  stats::setNames(vals, key)[c("At", "Al", "Al,At", "Cr", "Es", "Br", "Br,Es")]
}))
med <- apply(lens, 2, stats::median, na.rm = TRUE)
# The species tree's two root branches (lineage I and II stems) are one branch unrooted; split it equally.
r <- med[["Br,Es"]] / 2
tree <- sprintf("(((At:%g,Al:%g):%g,Cr:%g):%g,(Es:%g,Br:%g):%g);", med[["At"]], med[["Al"]], med[["Al,At"]], med[["Cr"]], r, med[["Es"]], med[["Br"]], r)
# Discordant gene tree for M5: Al and Cr swapped, keeping branch lengths.
discordant <- sprintf("(((At:%g,Cr:%g):%g,Al:%g):%g,(Es:%g,Br:%g):%g);", med[["At"]], med[["At"]], med[["Al,At"]], med[["Al"]] + med[["Al,At"]], r, med[["Es"]], med[["Br"]], r)
params <- list(
  source = "P3A real-genome benchmark; rules in docs/dev/phylo-calibration-design.md",
  orthogroups_used = nrow(full), L25 = q[[1]], L50 = q[[2]], L75 = q[[3]],
  nuc_f3x4 = lapply(1:3, function(i) unname(nuc[i, ])),
  kappa = 2, c = c_factor, tree = tree, discordant_tree = discordant,
  species_topology = "(((At,Al),Cr),(Es,Br));", foreground = c("At", "Al"),
  branch_medians_hyphy = as.list(med),
  empirical_tree_length_median = stats::median(full$tree_length)
)
jsonlite::write_json(params, out, auto_unbox = TRUE, digits = NA, pretty = TRUE)
cat(jsonlite::toJSON(params[c("orthogroups_used", "L25", "L50", "L75", "tree", "empirical_tree_length_median")], auto_unbox = TRUE, pretty = TRUE), "\n")
