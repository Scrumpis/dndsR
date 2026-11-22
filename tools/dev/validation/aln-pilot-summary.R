# Alignment-method pilot summary (docs/dev/alignment-pilot.md; measures and criteria fixed there).
#   Rscript --vanilla aln-pilot-summary.R <pilot root> <out prefix>
# Writes <prefix>_accuracy.tsv (per gene and arm), _size.tsv, _power.tsv, _diff.tsv (bootstrap
# differences against A0), _cost.tsv, _uncertainty.tsv, and prints a short report.
a <- commandArgs(TRUE); root <- a[1]; prefix <- a[2]
arms <- c("A0", "A1", "A2", "T")
cp <- function(x, n) if (n == 0) c(NA, NA) else stats::binom.test(x, n)$conf.int[1:2]
category <- function(x, n) {
  ci <- cp(x, n); if (n == 0) return("no data")
  if (ci[2] <= 0.05) "calibrated" else if (ci[1] >= 0.025 && ci[2] <= 0.075) "within tolerance" else if (ci[1] > 0.05) "anticonservative" else "inconclusive"
}
read_fasta <- function(path) {
  l <- readLines(path); h <- grep("^>", l)
  seqs <- vapply(seq_along(h), function(i) paste(l[seq.int(h[i] + 1L, length.out = (if (i < length(h)) h[i + 1L] else length(l) + 1L) - h[i] - 1L)], collapse = ""), "")
  stats::setNames(toupper(seqs), sub("\\|.*$", "", sub("^>", "", l[h])))
}
untar_arm <- function(set, arm) {
  d <- file.path(tempdir(), set, arm); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  utils::untar(file.path(root, "units", set, paste0(arm, ".tar.gz")), exdir = dirname(d)); d
}
# Codon index of each aligned codon column per taxon (NA = gap); a triplet with any base counts.
codon_index <- function(s) {
  tri <- substring(s, seq(1, nchar(s), 3), seq(3, nchar(s), 3))
  res <- !grepl("^-+$", tri); partial <- sum(grepl("-", tri) & res)
  idx <- rep(NA_integer_, length(tri)); idx[res] <- seq_len(sum(res)); list(idx = idx, partial = partial)
}
aln_index <- function(seqs) {
  ci <- lapply(seqs, codon_index)
  list(m = do.call(rbind, lapply(ci, `[[`, "idx")), partial = sum(vapply(ci, `[[`, 0L, "partial")))
}
pairs_of <- function(m) {
  out <- character(); tx <- rownames(m)
  for (i in seq_len(nrow(m) - 1L)) for (j in (i + 1L):nrow(m)) {
    k <- !is.na(m[i, ]) & !is.na(m[j, ])
    if (any(k)) out <- c(out, paste(tx[i], tx[j], m[i, k], m[j, k]))
  }
  out
}
columns_of <- function(m) {
  m <- m[order(rownames(m)), , drop = FALSE]
  keep <- colSums(!is.na(m)) >= 2L
  apply(m[, keep, drop = FALSE], 2, paste, collapse = ",")
}
sets <- intersect(c("null", "selected"), list.files(file.path(root, "units")))
acc <- list(); size <- list(); power <- list(); unc <- list(); per_gene <- list(); cost <- list()
for (set in sets) {
  u <- file.path(root, "units", set)
  genes <- utils::read.delim(file.path(u, "sim", "genes.tsv"), stringsAsFactors = FALSE)
  dirs <- stats::setNames(lapply(arms, function(x) untar_arm(set, x)), arms)
  inputs <- untar_arm(set, "unaligned")
  idx <- list(); pr <- list()
  for (g in genes$gene) {
    raw <- read_fasta(file.path(inputs, paste0(g, ".fasta")))
    for (arm in arms) {
      f <- file.path(dirs[[arm]], paste0(g, ".fasta"))
      if (!file.exists(f)) { acc[[length(acc) + 1L]] <- data.frame(set = set, gene = g, arm = arm, present = FALSE); next }
      s <- read_fasta(f); s <- s[names(raw)]
      valid <- all(gsub("-", "", s) == raw)
      ai <- aln_index(s); idx[[arm]][[g]] <- ai$m; pr[[arm]][[g]] <- pairs_of(ai$m)
      m <- ai$m; nt <- nrow(m)
      acc[[length(acc) + 1L]] <- data.frame(set = set, gene = g, arm = arm, present = TRUE, valid = valid, partial_codons = ai$partial,
        length_codons = ncol(m), ungapped_columns = sum(colSums(is.na(m)) == 0L), columns_half_or_less_gapped = sum(colSums(is.na(m)) <= nt / 2))
    }
    tp <- pr[["T"]][[g]]; tc <- columns_of(idx[["T"]][[g]])
    for (arm in setdiff(arms, "T")) {
      if (is.null(pr[[arm]][[g]])) next
      k <- which(vapply(acc, function(x) x$set == set && x$gene == g && x$arm == arm, TRUE))
      acc[[k]]$sp_score <- mean(tp %in% pr[[arm]][[g]])
      acc[[k]]$tc_score <- mean(tc %in% columns_of(idx[[arm]][[g]]))
    }
    a0 <- pr[["A0"]][[g]]
    unc[[length(unc) + 1L]] <- data.frame(set = set, gene = g,
      dis_A0_A1 = if (length(a0) && !is.null(pr[["A1"]][[g]])) mean(!(a0 %in% pr[["A1"]][[g]])) else NA_real_,
      dis_A0_A2 = if (length(a0) && !is.null(pr[["A2"]][[g]])) mean(!(a0 %in% pr[["A2"]][[g]])) else NA_real_,
      A0_exact = if (length(a0)) setequal(a0, tp) else NA)
  }
  for (arm in arms) {
    rd <- function(what, run) { f <- file.path(u, run, sprintf("%s_%s.tsv", run, what)); utils::read.delim(f, stringsAsFactors = FALSE) }
    r <- rd("results", paste0("absrel_", arm)); b <- rd("branches", paste0("absrel_", arm))
    t <- r[r$status == "tested", ]; b <- b[b$orthogroup %in% t$orthogroup & !is.na(b$p_raw), ]
    bo <- rd("results", paste0("branch_omega_", arm)); bo <- bo[bo$status == "tested", ]
    per_gene[[set]][[arm]] <- list(branches = b, genes = t, bo = bo)
    if (set == "null") {
      for (x in list(list("absrel p_raw per branch", b$p_raw), list("absrel p_bh_lineage per branch", b$p_bh_lineage),
                     list("absrel p_value per gene", t$p_value), list("absrel p_adj per gene", t$p_adj),
                     list("p_two_vs_m0", bo$p_two_vs_m0), list("p_fg_omega_gt1", bo$p_fg_omega_gt1))) {
        p <- x[[2]][!is.na(x[[2]])]; k <- sum(p <= 0.05); n <- length(p); ci <- cp(k, n)
        size[[length(size) + 1L]] <- data.frame(arm = arm, test = x[[1]], n = n, positives = k, rate = k / n, lo = ci[1], hi = ci[2], category = category(k, n))
      }
    } else {
      bt <- utils::read.delim(file.path(u, "sim", "branch_truth.tsv"), stringsAsFactors = FALSE)
      isfg <- b$species_tree_branches %in% bt$branch_id[bt$foreground]   # not b$foreground: TRUE for every tested branch
      fgdet <- tapply(b$p_raw[isfg] <= 0.05, b$orthogroup[isfg], any)
      for (x in list(list("absrel foreground p_raw (any fg branch)", fgdet[genes$gene]), list("absrel p_value per gene", (t$p_value <= 0.05)[match(genes$gene, t$orthogroup)]),
                     list("p_two_vs_m0", (bo$p_two_vs_m0 <= 0.05)[match(genes$gene, bo$orthogroup)]), list("p_fg_omega_gt1", (bo$p_fg_omega_gt1 <= 0.05)[match(genes$gene, bo$orthogroup)]))) {
        d <- x[[2]]; d[is.na(d)] <- FALSE; k <- sum(d); n <- length(d); ci <- cp(k, n)
        power[[length(power) + 1L]] <- data.frame(arm = arm, test = x[[1]], n = n, detected = k, power = k / n, lo = ci[1], hi = ci[2])
      }
    }
  }
  sec <- utils::read.delim(file.path(u, "align_seconds.tsv"), stringsAsFactors = FALSE)
  cpus <- 8L
  cost[[set]] <- data.frame(set = set, arm = c("A0", "A1", "A2"),
    core_seconds_per_gene = c(sec$seconds[sec$arm == "A0"] * cpus / nrow(genes), mean(sec$seconds[sec$arm == "A1"]), mean(sec$seconds[sec$arm == "A2"])),
    note = c("wall time of phylo_prepare on 8 threads x 8 / genes (upper bound: includes start-up)", "one core per gene, 8 in parallel", "one core per gene, 8 in parallel"))
}
acc <- do.call(rbind, lapply(acc, function(x) { for (c in c("valid", "partial_codons", "length_codons", "ungapped_columns", "columns_half_or_less_gapped", "sp_score", "tc_score")) if (is.null(x[[c]])) x[[c]] <- NA; x }))
w <- function(x, what) utils::write.table(x, paste0(prefix, "_", what, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
w(acc, "accuracy"); w(do.call(rbind, cost), "cost")
if (length(size)) w(do.call(rbind, size), "size")
if (length(power)) w(do.call(rbind, power), "power")

# Paired gene-bootstrap differences against A0 (2,000 resamples, seed 1).
boot <- function(gene_ids, stat, B = 2000L) { set.seed(1L); vapply(seq_len(B), function(i) stat(sample(gene_ids, replace = TRUE)), 0) }
diffs <- list()
if (!is.null(per_gene$null)) {
  g <- unique(per_gene$null$A0$branches$orthogroup)
  rate <- function(b, ids) { tab <- table(factor(b$orthogroup, levels = g)); pos <- tapply(b$p_raw <= 0.05, factor(b$orthogroup, levels = g), sum); sum(pos[ids], na.rm = TRUE) / sum(tab[ids]) }
  for (arm in c("A1", "A2", "T")) {
    b0 <- per_gene$null$A0$branches; b1 <- per_gene$null[[arm]]$branches
    est <- rate(b1, g) - rate(b0, g); bs <- boot(g, function(ids) rate(b1, ids) - rate(b0, ids))
    diffs[[length(diffs) + 1L]] <- data.frame(set = "null", measure = "absrel p_raw false-positive rate", arm = arm, diff = est, lo = unname(stats::quantile(bs, 0.025)), hi = unname(stats::quantile(bs, 0.975)))
  }
}
if (!is.null(per_gene$selected)) {
  g <- unique(per_gene$selected$A0$genes$orthogroup)
  bt <- utils::read.delim(file.path(root, "units", "selected", "sim", "branch_truth.tsv"), stringsAsFactors = FALSE)
  det <- function(b) { isfg <- b$species_tree_branches %in% bt$branch_id[bt$foreground]; x <- tapply(b$p_raw[isfg] <= 0.05, factor(b$orthogroup[isfg], levels = g), any); x[is.na(x)] <- FALSE; x }
  d0 <- det(per_gene$selected$A0$branches)
  for (arm in c("A1", "A2", "T")) {
    d1 <- det(per_gene$selected[[arm]]$branches)
    bs <- boot(g, function(ids) mean(d1[ids]) - mean(d0[ids]))
    diffs[[length(diffs) + 1L]] <- data.frame(set = "selected", measure = "absrel foreground power", arm = arm, diff = mean(d1) - mean(d0), lo = unname(stats::quantile(bs, 0.025)), hi = unname(stats::quantile(bs, 0.975)))
  }
}
if (length(diffs)) w(do.call(rbind, diffs), "diff")

# P4: A0-A1 and A0-A2 disagreement against (a) an inexact A0 alignment and (b) a realignment-induced positive.
unc <- do.call(rbind, unc)
if (!is.null(per_gene$null)) {
  key <- function(b) paste(b$orthogroup, b$species_tree_branches)
  b0 <- per_gene$null$A0$branches; bt <- per_gene$null$T$branches
  pt <- bt$p_raw[match(key(b0), key(bt))]
  induced <- unique(b0$orthogroup[b0$p_raw <= 0.05 & !is.na(pt) & pt > 0.05])
  unc$induced_positive <- ifelse(unc$set == "null", unc$gene %in% induced, NA)
}
w(unc, "uncertainty")
rep_rule <- function(flag, target) {
  ok <- !is.na(flag) & !is.na(target)
  sprintf("sensitivity %d/%d, flags on unaffected %d/%d", sum(flag[ok] & target[ok]), sum(target[ok]), sum(flag[ok] & !target[ok]), sum(!target[ok]))
}
for (d in c("dis_A0_A1", "dis_A0_A2")) for (set in sets) {
  x <- unc[unc$set == set, ]
  cat(sprintf("%s %s vs inexact A0: %s\n", set, d, rep_rule(x[[d]] > 0, !x$A0_exact)))
  if (set == "null") cat(sprintf("%s %s vs induced positive: %s\n", set, d, rep_rule(x[[d]] > 0, x$induced_positive)))
}
print(stats::aggregate(cbind(sp_score, tc_score, length_codons, ungapped_columns, columns_half_or_less_gapped, partial_codons) ~ set + arm, acc[acc$present, ], mean))
cat("invalid alignments (do not reproduce the CDS):", sum(!acc$valid, na.rm = TRUE), "; missing:", sum(!acc$present), "\n")
if (length(size)) print(do.call(rbind, size), digits = 3)
if (length(power)) print(do.call(rbind, power), digits = 3)
if (length(diffs)) print(do.call(rbind, diffs), digits = 3)
print(do.call(rbind, cost), digits = 3)
