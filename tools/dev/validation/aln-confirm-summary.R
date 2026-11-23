# Alignment confirmation summary (docs/dev/alignment-confirmation.md; measures and criteria fixed there).
#   Rscript --vanilla aln-confirm-summary.R <root> <out prefix>
# Units named C-N* (null), C-S* (every gene selected), C-M* (mixture); arms A0 (phylo_prepare, MAFFT),
# A1 (PRANK codon), T (true alignment, diagnostic). Writes <prefix>_{accuracy,size,power,fdr,diff,cost,
# repro,criteria}.tsv and prints them. Truth per branch: omega_mean > 1 in sim/gene_branch_truth.tsv.
a <- commandArgs(TRUE); root <- a[1]; prefix <- a[2]
arms <- c("A0", "A1", "T"); B <- 2000L
cp <- function(x, n) if (n == 0) c(NA, NA) else stats::binom.test(x, n)$conf.int[1:2]
category <- function(x, n) { ci <- cp(x, n); if (n == 0) return("no data")
  if (ci[2] <= 0.05) "calibrated" else if (ci[1] >= 0.025 && ci[2] <= 0.075) "within tolerance" else if (ci[1] > 0.05) "anticonservative" else "inconclusive" }
read_fasta <- function(path) {
  l <- readLines(path); h <- grep("^>", l)
  seqs <- vapply(seq_along(h), function(i) paste(l[seq.int(h[i] + 1L, length.out = (if (i < length(h)) h[i + 1L] else length(l) + 1L) - h[i] - 1L)], collapse = ""), "")
  stats::setNames(toupper(seqs), sub("\\|.*$", "", sub("^>", "", l[h])))
}
untar_arm <- function(u, arm) { d <- file.path(tempdir(), u); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  utils::untar(file.path(root, "units", u, paste0(arm, ".tar.gz")), exdir = d); file.path(d, arm) }
codon_idx <- function(s) { tri <- substring(s, seq(1, nchar(s), 3), seq(3, nchar(s), 3)); res <- !grepl("^-+$", tri)
  idx <- rep(NA_integer_, length(tri)); idx[res] <- seq_len(sum(res)); idx }
aln_m <- function(seqs) do.call(rbind, lapply(seqs, codon_idx))
pairs_of <- function(m) { out <- character(); tx <- rownames(m)
  for (i in seq_len(nrow(m) - 1L)) for (j in (i + 1L):nrow(m)) { k <- !is.na(m[i, ]) & !is.na(m[j, ]); if (any(k)) out <- c(out, paste(tx[i], tx[j], m[i, k], m[j, k])) }
  out }
columns_of <- function(m) { m <- m[order(rownames(m)), , drop = FALSE]; keep <- colSums(!is.na(m)) >= 2L; apply(m[, keep, drop = FALSE], 2, paste, collapse = ",") }
boot <- function(ids, stat) { set.seed(1L); vapply(seq_len(B), function(i) stat(sample(ids, replace = TRUE)), 0) }
q <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))
rd <- function(u, run, what) utils::read.delim(file.path(root, "units", u, run, sprintf("%s_%s.tsv", run, what)), stringsAsFactors = FALSE)

units <- basename(dirname(Sys.glob(file.path(root, "units", "*", "DONE"))))
acc <- list(); size <- list(); power <- list(); fdr <- list(); diffs <- list(); cost <- list(); repro <- list()
for (u in units) {
  ud <- file.path(root, "units", u); kind <- substr(u, 3, 3)
  dirs <- stats::setNames(lapply(arms, function(x) untar_arm(u, x)), arms); inputs <- untar_arm(u, "unaligned")
  genes_all <- utils::read.delim(file.path(ud, "sim", "genes.tsv"), stringsAsFactors = FALSE)
  kept <- sub("\\.fasta$", "", list.files(inputs, pattern = "\\.fasta$"))
  genes <- genes_all[genes_all$gene %in% kept, ]
  gbt <- utils::read.delim(file.path(ud, "sim", "gene_branch_truth.tsv"), stringsAsFactors = FALSE)
  truekey <- paste(gbt$gene, gbt$branch_id)[gbt$omega_mean > 1]
  # Accuracy and validity
  for (g in genes$gene) {
    raw <- read_fasta(file.path(inputs, paste0(g, ".fasta"))); m <- list(); pr <- list()
    for (arm in arms) {
      f <- file.path(dirs[[arm]], paste0(g, ".fasta"))
      if (!file.exists(f)) { acc[[length(acc) + 1L]] <- data.frame(unit = u, gene = g, arm = arm, present = FALSE, valid = FALSE, length_codons = NA,
        ungapped_columns = NA, sp_score = NA, tc_score = NA); next }
      s <- read_fasta(f)[names(raw)]; m[[arm]] <- aln_m(s); pr[[arm]] <- pairs_of(m[[arm]])
      acc[[length(acc) + 1L]] <- data.frame(unit = u, gene = g, arm = arm, present = TRUE, valid = all(gsub("-", "", s) == raw),
        length_codons = ncol(m[[arm]]), ungapped_columns = sum(colSums(is.na(m[[arm]])) == 0L),
        sp_score = if (arm == "T") 1 else NA, tc_score = if (arm == "T") 1 else NA)
    }
    for (arm in c("A0", "A1")) if (!is.null(pr[[arm]]) && !is.null(pr[["T"]])) {
      k <- length(acc) - (length(arms) - match(arm, arms))
      acc[[k]]$sp_score <- mean(pr[["T"]] %in% pr[[arm]]); acc[[k]]$tc_score <- mean(columns_of(m[["T"]]) %in% columns_of(m[[arm]]))
    }
  }
  per <- list()
  for (arm in arms) {
    r <- rd(u, paste0("absrel_", arm), "results"); b <- rd(u, paste0("absrel_", arm), "branches")
    t <- r[r$status == "tested", ]; b <- b[b$orthogroup %in% t$orthogroup & !is.na(b$p_raw), ]
    b$true <- paste(b$orthogroup, b$species_tree_branches) %in% truekey
    t$true <- t$orthogroup %in% genes$gene[genes$selected]
    per[[arm]] <- list(b = b, t = t, n_input = nrow(genes), n_tested = nrow(t))
    if (kind == "N") {
      for (x in list(list("p_raw per branch", b$p_raw), list("gene p_value", t$p_value), list("gene p_adj (rejections)", t$p_adj), list("p_bh_lineage (rejections)", b$p_bh_lineage))) {
        p <- x[[2]][!is.na(x[[2]])]; k <- sum(p <= 0.05); n <- length(p); ci <- cp(k, n)
        size[[length(size) + 1L]] <- data.frame(unit = u, arm = arm, test = x[[1]], n = n, positives = k, rate = k / n, lo = ci[1], hi = ci[2], category = category(k, n))
      }
    }
    if (kind == "S") {
      fgdet <- tapply(b$p_raw[b$true] <= 0.05, factor(b$orthogroup[b$true], levels = genes$gene), any); fgdet[is.na(fgdet)] <- FALSE
      gp <- (t$p_value <= 0.05)[match(genes$gene, t$orthogroup)]; gp[is.na(gp)] <- FALSE
      for (x in list(list("foreground p_raw (any fg branch)", fgdet), list("gene p_value", gp))) {
        k <- sum(x[[2]]); n <- length(x[[2]]); ci <- cp(k, n)
        power[[length(power) + 1L]] <- data.frame(unit = u, arm = arm, test = x[[1]], n = n, detected = k, power = k / n, lo = ci[1], hi = ci[2])
      }
    }
    if (kind == "M") {
      fdp_row <- function(output, gene, p, truth) {
        ok <- !is.na(p); gene <- gene[ok]; p <- p[ok]; truth <- truth[ok]; disc <- p <= 0.05; D <- sum(disc); F <- sum(disc & !truth)
        ids <- unique(gene); f <- factor(gene, levels = ids); dd <- tapply(disc, f, sum); ff <- tapply(disc & !truth, f, sum)
        bs <- if (D) boot(seq_along(ids), function(s) { x <- sum(dd[s]); if (x == 0) NA_real_ else sum(ff[s]) / x }) else NA
        lo <- if (D) q(bs, 0.025) else NA; hi <- if (D) q(bs, 0.975) else NA
        data.frame(unit = u, arm = arm, output = output, discoveries = D, false = F, fdp = if (D) F / D else NA, fdp_lo = lo, fdp_hi = hi,
          judgement = if (!D) "no discoveries" else if (hi <= 0.075) "controls" else if (lo > 0.05) "fails" else "inconclusive",
          power = if (sum(truth)) sum(disc & truth) / sum(truth) else NA)
      }
      fdr[[length(fdr) + 1L]] <- fdp_row("gene p_adj", t$orthogroup, t$p_adj, t$true)
      fdr[[length(fdr) + 1L]] <- fdp_row("p_bh_lineage, pooled over lineages", b$orthogroup, b$p_bh_lineage, b$true)
      nullfam <- tapply(b$true, b$species_tree_branches, function(x) !any(x)); rej <- tapply(b$p_bh_lineage <= 0.05, b$species_tree_branches, any, na.rm = TRUE)
      fdr[[length(fdr) + 1L]] <- data.frame(unit = u, arm = arm, output = sprintf("lineages without selection rejecting anything: %d of %d", sum(rej[nullfam], na.rm = TRUE), sum(nullfam)),
        discoveries = NA, false = NA, fdp = NA, fdp_lo = NA, fdp_hi = NA, judgement = NA, power = NA)
    }
  }
  # Bootstrap differences A1 - A0 (and T - A0), paired over genes
  ids <- genes$gene
  if (kind == "N") {
    rate_b <- function(b, s) { pos <- tapply(b$p_raw <= 0.05, factor(b$orthogroup, levels = ids), sum); n <- tapply(b$p_raw, factor(b$orthogroup, levels = ids), length)
      sum(pos[s], na.rm = TRUE) / sum(n[s], na.rm = TRUE) }
    rate_g <- function(t, s) { x <- (t$p_value <= 0.05)[match(s, t$orthogroup)]; mean(x, na.rm = TRUE) }
    for (arm in c("A1", "T")) {
      b0 <- per$A0$b; b1 <- per[[arm]]$b; t0 <- per$A0$t; t1 <- per[[arm]]$t
      bs <- boot(ids, function(s) rate_b(b1, s) - rate_b(b0, s))
      diffs[[length(diffs) + 1L]] <- data.frame(unit = u, measure = "p_raw false-positive rate", arm = arm, diff = rate_b(b1, ids) - rate_b(b0, ids), lo = q(bs, 0.025), hi = q(bs, 0.975))
      bs <- boot(ids, function(s) rate_g(t1, s) - rate_g(t0, s))
      diffs[[length(diffs) + 1L]] <- data.frame(unit = u, measure = "gene p_value false-positive rate", arm = arm, diff = rate_g(t1, ids) - rate_g(t0, ids), lo = q(bs, 0.025), hi = q(bs, 0.975))
    }
  }
  if (kind == "S") {
    det <- function(b) { x <- tapply(b$p_raw[b$true] <= 0.05, factor(b$orthogroup[b$true], levels = ids), any); x[is.na(x)] <- FALSE; x }
    gdet <- function(t) { x <- (t$p_value <= 0.05)[match(ids, t$orthogroup)]; x[is.na(x)] <- FALSE; stats::setNames(x, ids) }
    for (arm in c("A1", "T")) {
      d0 <- det(per$A0$b); d1 <- det(per[[arm]]$b); bs <- boot(ids, function(s) mean(d1[s]) - mean(d0[s]))
      diffs[[length(diffs) + 1L]] <- data.frame(unit = u, measure = "foreground p_raw power", arm = arm, diff = mean(d1) - mean(d0), lo = q(bs, 0.025), hi = q(bs, 0.975))
    }
    g1 <- gdet(per$A1$t); gT <- gdet(per$T$t); bs <- boot(ids, function(s) mean(g1[s]) - mean(gT[s]))
    diffs[[length(diffs) + 1L]] <- data.frame(unit = u, measure = "gene p_value power, A1 against T", arm = "A1", diff = mean(g1) - mean(gT), lo = q(bs, 0.025), hi = q(bs, 0.975))
  }
  # Cost
  ac <- utils::read.delim(file.path(ud, "align_cost.tsv"), stringsAsFactors = FALSE)
  as <- utils::read.delim(file.path(ud, "absrel_seconds.tsv"), stringsAsFactors = FALSE)
  a1 <- ac[ac$arm == "A1", ]; a0 <- ac[ac$arm == "A0", ]; abs0 <- as[as$arm == "A0", ]
  cost[[u]] <- data.frame(unit = u, genes = nrow(genes),
    prank_core_s_per_gene = mean(a1$seconds), prank_max_core_s = max(a1$seconds), prank_max_rss_gb = max(a1$max_rss_kb) / 1024^2,
    mafft_core_s_per_gene_upper = a0$seconds * 8 / nrow(genes_all), mafft_rss_gb = a0$max_rss_kb / 1024^2,
    absrel_core_s_per_gene = abs0$seconds * abs0$cpus / abs0$genes)
  cost[[u]]$prank_over_absrel <- cost[[u]]$prank_core_s_per_gene / cost[[u]]$absrel_core_s_per_gene
  rp <- utils::read.delim(file.path(ud, "reproducibility.tsv"), stringsAsFactors = FALSE)
  repro[[u]] <- data.frame(unit = u, arm = c("A0", "A1"), genes = c(sum(rp$arm == "A0"), sum(rp$arm == "A1")),
    identical = c(sum(rp$identical[rp$arm == "A0"]), sum(rp$identical[rp$arm == "A1"])))
}
acc <- do.call(rbind, acc)
w <- function(x, what) if (length(x)) { x <- if (is.data.frame(x)) x else do.call(rbind, x); utils::write.table(x, paste0(prefix, "_", what, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE); x }
acc_sum <- stats::aggregate(cbind(valid, sp_score, tc_score, length_codons, ungapped_columns) ~ unit + arm, acc[acc$present, ], mean)
acc_sum$present <- stats::aggregate(present ~ unit + arm, acc, mean)$present
w(acc, "accuracy_genes"); w(acc_sum, "accuracy"); size <- w(size, "size"); power <- w(power, "power"); fdr <- w(fdr, "fdr"); diffs <- w(diffs, "diff")
cost <- w(cost, "cost"); repro <- w(repro, "repro")

# Criteria (alignment-confirmation.md)
crit <- list(); add <- function(c, scope, ok, detail) crit[[length(crit) + 1L]] <<- data.frame(criterion = c, scope = scope, pass = ok, detail = detail)
for (u in grep("^C-N", names(cost), value = TRUE)) {
  d <- diffs[diffs$unit == u & diffs$arm == "A1", ]; s <- size[size$unit == u & size$arm == "A1", ]
  x <- d[d$measure == "p_raw false-positive rate", ]; y <- s[s$test == "p_raw per branch", ]
  add("1 branch Type I", u, x$hi < 0 && y$hi <= 0.075, sprintf("diff %.4f (%.4f to %.4f); A1 upper %.4f", x$diff, x$lo, x$hi, y$hi))
  x <- d[d$measure == "gene p_value false-positive rate", ]; y <- s[s$test == "gene p_value", ]
  add("2 gene Type I", u, x$hi < 0 && y$lo <= 0.05, sprintf("diff %.4f (%.4f to %.4f); A1 lower %.4f", x$diff, x$lo, x$hi, y$lo))
}
for (u in grep("^C-M", names(cost), value = TRUE)) for (o in c("gene p_adj", "p_bh_lineage, pooled over lineages")) {
  f1 <- fdr[fdr$unit == u & fdr$arm == "A1" & fdr$output == o, ]; f0 <- fdr[fdr$unit == u & fdr$arm == "A0" & fdr$output == o, ]
  ok <- (is.na(f1$fdp_lo) || f1$fdp_lo <= 0.05) && (is.na(f1$fdp) || is.na(f0$fdp) || f1$fdp < f0$fdp)
  add("3 false discoveries", paste(u, o), ok, sprintf("A1 FDP %.3f (%.3f-%.3f), A0 %.3f", f1$fdp, f1$fdp_lo, f1$fdp_hi, f0$fdp))
}
for (u in grep("^C-S", names(cost), value = TRUE)) {
  x <- diffs[diffs$unit == u & diffs$arm == "A1" & diffs$measure == "foreground p_raw power", ]
  add("4 power", u, x$diff > -0.05 && x$lo > -0.10, sprintf("diff %.3f (%.3f to %.3f)", x$diff, x$lo, x$hi))
}
for (u in names(cost)) {
  a <- acc[acc$unit == u & acc$arm == "A1", ]; rp <- repro[repro$unit == u, ]
  add("5 completeness and reproducibility", u, mean(a$present & a$valid) >= 0.99 && all(rp$identical == rp$genes),
      sprintf("A1 valid %.3f; identical reruns A0 %d/%d, A1 %d/%d", mean(a$present & a$valid), rp$identical[1], rp$genes[1], rp$identical[2], rp$genes[2]))
}
for (u in grep("16", names(cost), value = TRUE)) {
  x <- cost[cost$unit == u, ]
  add("6 cost (16 taxa)", u, x$prank_over_absrel <= 1 && x$prank_max_rss_gb <= 4, sprintf("PRANK %.1f core-s/gene against aBSREL %.1f (ratio %.2f); peak %.2f GB", x$prank_core_s_per_gene, x$absrel_core_s_per_gene, x$prank_over_absrel, x$prank_max_rss_gb))
}
crit <- w(crit, "criteria")
for (x in list(acc_sum, size, power, fdr, diffs, cost, repro, crit)) if (!is.null(x)) print(x, digits = 3, row.names = FALSE)
