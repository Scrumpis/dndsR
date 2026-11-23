# Checks one simulated data set (phylo-lineage-sim2.R output) against its cell's design, with the
# expected values derived here from the design documents' numbers, not from the simulator's code:
#   Rscript --vanilla sim-check.R <sim dir> <cell>
# Prints one line per check and exits nonzero if any fails. Checks: every unselected gene's branch
# omega means equal the design's; selected genes are raised (classes 0.05/0.3/10, mean 1.12) on exactly
# the intended branches; selection rates; L-A1's loss rule; T16m's four absent taxa; no gaps where the
# cell has no indels; RELAX cells' test-branch omegas (gamma classes from numerical integration).
a <- commandArgs(TRUE); d <- a[1]; cell <- a[2]; tol <- 1e-9
genes <- utils::read.delim(file.path(d, "genes.tsv"), stringsAsFactors = FALSE, colClasses = c(dropped = "character"))
gbt <- utils::read.delim(file.path(d, "gene_branch_truth.tsv"), stringsAsFactors = FALSE)
bt <- utils::read.delim(file.path(d, "branch_truth.tsv"), stringsAsFactors = FALSE)
fg <- bt$branch_id[bt$foreground]
fails <- 0L
ok <- function(what, cond, detail = "") { cat(sprintf("%s %s: %s %s\n", if (isTRUE(cond)) "ok  " else "FAIL", cell, what, detail)); if (!isTRUE(cond)) fails <<- fails + 1L }
b0 <- sum(c(0.05, 0.3, 1) * c(0.6, 0.3, 0.1))       # 0.22
sel_mean <- sum(c(0.05, 0.3, 10) * c(0.6, 0.3, 0.1)) # 1.12
# Gamma(shape 0.5, given mean), 8 equal-probability classes, conditional means by numerical integration
gclass <- function(mean, shape = 0.5, k = 8) { rate <- shape / mean; q <- c(0, stats::qgamma(seq_len(k - 1) / k, shape, rate), Inf)
  vapply(seq_len(k), function(i) stats::integrate(function(x) x * stats::dgamma(x, shape, rate), q[i], q[i + 1])$value * k, 0) }
base <- switch(cell, LS0 = 0.2, LS2 = 0.2, LG1 = 0.2, LA1 = 0.2, RG0 = 0.3, RG05 = 0.3, RG2 = 0.3, b0)
eff <- stats::setNames(bt$branch_effect, bt$branch_id)
if (cell %in% c("LG1", "LA1")) {
  ok("branch effects vary (sd 0.5 on the log scale)", stats::sd(log(eff)) > 0.1, sprintf("sd %.2f", stats::sd(log(eff))))
} else ok("no branch effects", all(abs(eff - 1) < tol))
m <- stats::setNames(genes$gene_effect, genes$gene)
if (cell %in% c("LG1", "LA1")) { ok("gene effects vary", stats::sd(log(m)) > 0.1) } else ok("no gene effects", all(abs(m - 1) < tol))
relax_k <- switch(cell, RG05 = 0.5, RG2 = 2, 1)
fg_mean_relax <- mean(gclass(0.3)^relax_k)
expect <- function(g, b) {
  v <- base * m[[g]] * eff[[b]]
  if (cell %in% c("RG05", "RG2") && b %in% fg) v <- fg_mean_relax
  v
}
gbt$expected <- mapply(expect, gbt$gene, gbt$branch_id)
sel <- genes$gene[genes$selected]
raised <- gbt[abs(gbt$omega_mean - sel_mean) < 1e-6, ]
uns <- gbt[!gbt$gene %in% sel, ]
ok("unselected genes: every branch at the design's omega mean", all(abs(uns$omega_mean - uns$expected) < 1e-6),
   sprintf("(%d of %d branch values differ)", sum(abs(uns$omega_mean - uns$expected) >= 1e-6), nrow(uns)))
ok("no unselected gene raised to 1.12", !any(raised$gene %in% genes$gene[!genes$selected]))
sets <- split(raised$branch_id, raised$gene)
if (cell %in% c("BMIX", "BMIXh", "BE3hS")) {
  ok("selected genes raised on exactly the foreground (t1, t2, their stem)", length(sel) == 0 || all(vapply(sel, function(g) setequal(sets[[g]], fg), TRUE)))
  if (cell == "BE3hS") ok("every gene selected", all(genes$selected))
  else ok("selection rate consistent with 10% (binomial p > 0.001)", stats::binom.test(length(sel), nrow(genes), 0.1)$p.value > 0.001, sprintf("%d of %d", length(sel), nrow(genes)))
} else if (cell == "BMIXv") {
  ok("each selected gene raised on 1 or 2 branches", all(vapply(sel, function(g) length(sets[[g]]) %in% 1:2, TRUE)))
  ok("not the fixed foreground set (partial-matching regression)", !any(vapply(sel, function(g) setequal(sets[[g]], fg), TRUE)))
  ok("selected genes' other branches at the design's mean", all(abs(gbt$omega_mean[gbt$gene %in% sel & abs(gbt$omega_mean - sel_mean) > 1e-6] - b0) < 1e-6))
  ok("selection rate consistent with 10%", stats::binom.test(length(sel), nrow(genes), 0.1)$p.value > 0.001, sprintf("%d of %d", length(sel), nrow(genes)))
} else {
  ok("no gene selected", length(sel) == 0)
}
if (cell == "LA1") {
  cut <- exp(-0.7^2 / 2 + 0.7 * stats::qnorm(2 / 3))
  lost <- grepl("(^|,)t1(,|$)", genes$dropped)
  ok("t1 absent only above the gene-effect cutoff", all(genes$gene_effect[lost] > cut), sprintf("%d lost", sum(lost)))
}
tree_name <- genes$tree[1]
if (tree_name == "T16m") {
  ok("four taxa absent per gene", all(lengths(strsplit(genes$dropped, ",")) == 4L))
} else if (cell != "LA1") ok("no taxon absent", all(!nzchar(genes$dropped) | is.na(genes$dropped)))
fa <- list.files(file.path(d, "alignments"), full.names = TRUE)
if (length(fa)) {
  gaps <- vapply(fa, function(f) any(grepl("-", readLines(f)[-grep("^>", readLines(f))])), TRUE)
  if (cell %in% c("BE3", "BE3h", "BE3hS", "BMIXh")) ok("indel cell: some alignments have gaps", any(gaps))
  else if (tree_name != "T16m") ok("no gaps where the cell has no indels", !any(gaps))
}
if (fails) { cat(sprintf("%s: %d check(s) failed\n", cell, fails)); quit(status = 1) }
