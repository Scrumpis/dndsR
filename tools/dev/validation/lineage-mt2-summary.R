# Multiple-testing stage 2 summary (docs/dev/lineage-multiple-testing.md, "Stage 2"; criteria fixed there).
#   Rscript --vanilla lineage-mt2-summary.R <out prefix> <label>=<unit dir>:<run> ...
# Each data set is a unit folder (with sim/genes.tsv and sim/gene_branch_truth.tsv) and an aBSREL run
# name in it; labels starting with "realigned" are reported but not judged.
# Procedures on the gene-by-branch family (F-all), from p-values already computed (no rerun):
#   P-BH  one BH over all branch p_raw;  P-BY  Benjamini-Yekutieli over the same;
#   P-BB  Benjamini-Bogomolov: BH at q over genes on gene p_value (R of m selected), then BH within each
#         selected gene over its branch p_raw at q R / m;  P-lineage-pooled  p_bh_lineage <= q pooled.
# Also F-gene (gene p_adj) and F-lineage (p_bh_lineage within each family). Truth: a branch test is
# true when any species-tree branch it covers has omega_mean > 1 for that gene.
a <- commandArgs(TRUE); prefix <- a[1]; sets <- a[-1]; qv <- 0.05; B <- 2000L
fdp_ci <- function(gene, disc, false) {
  ids <- unique(gene); f <- factor(gene, levels = ids); d <- tapply(disc, f, sum); fa <- tapply(disc & false, f, sum)
  if (sum(d) == 0) return(c(NA, NA)); set.seed(1L)
  bs <- vapply(seq_len(B), function(i) { s <- sample(seq_along(ids), replace = TRUE); D <- sum(d[s]); if (D == 0) NA_real_ else sum(fa[s]) / D }, 0)
  unname(stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE))
}
judge <- function(lo, hi, D) if (D == 0) "no discoveries" else if (hi <= 0.075) "controls" else if (lo > 0.05) "fails" else "inconclusive"
out <- list(); fam <- list()
for (s in sets) {
  label <- sub("=.*$", "", s); x <- sub("^[^=]*=", "", s); ud <- sub(":[^:]*$", "", x); run <- sub("^.*:", "", x)
  r <- utils::read.delim(file.path(ud, run, paste0(run, "_results.tsv")), stringsAsFactors = FALSE)
  b <- utils::read.delim(file.path(ud, run, paste0(run, "_branches.tsv")), stringsAsFactors = FALSE)
  t <- r[r$status == "tested", ]; b <- b[b$orthogroup %in% t$orthogroup & !is.na(b$p_raw), ]
  genes <- utils::read.delim(file.path(ud, "sim", "genes.tsv"), stringsAsFactors = FALSE)
  gbt <- utils::read.delim(file.path(ud, "sim", "gene_branch_truth.tsv"), stringsAsFactors = FALSE)
  tk <- paste(gbt$gene, gbt$branch_id)[gbt$omega_mean > 1]
  b$true <- vapply(seq_len(nrow(b)), function(i) any(paste(b$orthogroup[i], strsplit(b$species_tree_branches[i], "+", fixed = TRUE)[[1]]) %in% tk), TRUE)
  t$true <- t$orthogroup %in% genes$gene[genes$selected]
  m <- nrow(t)
  sel_genes <- t$orthogroup[!is.na(t$p_adj) & t$p_adj <= qv]; R <- length(sel_genes)
  bb <- rep(FALSE, nrow(b))
  for (g in sel_genes) { i <- which(b$orthogroup == g); bb[i] <- stats::p.adjust(b$p_raw[i], "BH") <= qv * R / m }
  within <- if (R) vapply(sel_genes, function(g) { i <- b$orthogroup == g; D <- sum(bb[i]); if (D == 0) 0 else sum(bb[i] & !b$true[i]) / D }, 0) else numeric()
  procs <- list(`P-BH` = stats::p.adjust(b$p_raw, "BH") <= qv, `P-BY` = stats::p.adjust(b$p_raw, "BY") <= qv, `P-BB` = bb,
                `P-lineage-pooled` = !is.na(b$p_bh_lineage) & b$p_bh_lineage <= qv)
  for (p in names(procs)) {
    d <- procs[[p]]; D <- sum(d); F <- sum(d & !b$true); ci <- fdp_ci(b$orthogroup, d, !b$true)
    out[[length(out) + 1L]] <- data.frame(data = label, family = "F-all", procedure = p, tests = nrow(b), true_tests = sum(b$true), discoveries = D, false = F,
      fdp = if (D) F / D else NA, fdp_lo = ci[1], fdp_hi = ci[2], judgement = judge(ci[1], ci[2], D), power = sum(d & b$true) / sum(b$true),
      bb_selected_genes = if (p == "P-BB") R else NA, bb_mean_within_gene_fdp = if (p == "P-BB" && R) mean(within) else NA)
  }
  d <- !is.na(t$p_adj) & t$p_adj <= qv; ci <- fdp_ci(t$orthogroup, d, !t$true)
  out[[length(out) + 1L]] <- data.frame(data = label, family = "F-gene", procedure = "gene p_adj", tests = m, true_tests = sum(t$true), discoveries = sum(d), false = sum(d & !t$true),
    fdp = if (sum(d)) sum(d & !t$true) / sum(d) else NA, fdp_lo = ci[1], fdp_hi = ci[2], judgement = judge(ci[1], ci[2], sum(d)), power = sum(d & t$true) / sum(t$true),
    bb_selected_genes = NA, bb_mean_within_gene_fdp = NA)
  for (f in unique(b$species_tree_branches)) {
    i <- b$species_tree_branches == f; d <- !is.na(b$p_bh_lineage[i]) & b$p_bh_lineage[i] <= qv
    fam[[length(fam) + 1L]] <- data.frame(data = label, lineage = f, tests = sum(i), true_tests = sum(b$true[i]), discoveries = sum(d), false = sum(d & !b$true[i]))
  }
}
res <- do.call(rbind, out); fam <- do.call(rbind, fam)
utils::write.table(res, paste0(prefix, "_fdr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(fam, paste0(prefix, "_families.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, digits = 3, row.names = FALSE)
# F-lineage: families without any true test that rejected anything; families with true tests: FDP
for (dl in unique(fam$data)) { x <- fam[fam$data == dl, ]; n0 <- x$true_tests == 0
  cat(sprintf("%s F-lineage: lineages without selection rejecting anything %d of %d; lineages with selection: %d discoveries, %d false\n",
              dl, sum(x$discoveries[n0] > 0), sum(n0), sum(x$discoveries[!n0]), sum(x$false[!n0]))) }
# Proposal rule: among P-BH, P-BY, P-BB, those that control on every error-free data set; highest mean power.
ef <- res[res$family == "F-all" & !grepl("^realigned", res$data) & res$procedure %in% c("P-BH", "P-BY", "P-BB"), ]
ok <- tapply(ef$judgement == "controls", ef$procedure, all); pw <- tapply(ef$power, ef$procedure, mean)
cat("\nControls F-all on every error-free data set:\n"); print(data.frame(procedure = names(ok), controls_everywhere = ok, mean_power = pw[names(ok)]), row.names = FALSE)
cand <- names(ok)[ok]
# Ties: within 1 point of the best mean power, the simplest procedure (P-BH, then P-BB, then P-BY) is chosen.
pick <- if (length(cand)) { best <- max(pw[cand]); near <- cand[pw[cand] >= best - 0.01]; intersect(c("P-BH", "P-BB", "P-BY"), near)[1] } else NULL
cat("Proposal rule:", if (length(cand)) pick else "none controls everywhere: nothing proposed", "\n")
