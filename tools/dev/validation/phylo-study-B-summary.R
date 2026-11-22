# Study B summary (docs/dev/lineage-robustness-study.md): size under each error, false discoveries in
# the mixture, masking, and detection, from phylo-lineage-study2.sb units.
#   Rscript --vanilla phylo-study-B-summary.R <results root> <out prefix>
# Size definitions as the first lineage study (phylo-lineage-calibration-summary.R): aBSREL p_raw and
# p_bh_lineage per tested branch, aBSREL p_value per gene (Holm minimum); branch models p_two_vs_m0
# and p_fg_omega_gt1 per gene; BUSTED p_value per gene. At 0.05, Clopper-Pearson intervals.
a <- commandArgs(TRUE); root <- a[1]; prefix <- a[2]
cp <- function(x, n) if (n == 0) c(NA, NA) else stats::binom.test(x, n)$conf.int[1:2]
category <- function(x, n) {
  ci <- cp(x, n); if (n == 0) return("no data")
  if (ci[2] <= 0.05) "calibrated" else if (ci[1] >= 0.025 && ci[2] <= 0.075) "within tolerance" else if (ci[1] > 0.05) "anticonservative" else if (ci[1] > 0.075 || ci[2] < 0.025) "outside tolerance" else "inconclusive"
}
rd <- function(u, run, what) { f <- file.path(root, "units", u, run, sprintf("%s_%s.tsv", run, what)); if (file.exists(f)) utils::read.delim(f, stringsAsFactors = FALSE) else NULL }
units <- basename(dirname(Sys.glob(file.path(root, "units", "*", "DONE"))))
size <- list(); detect <- list(); mixture <- list()
tests <- function(u, run) {
  m <- sub("_mask$", "", run); out <- list()
  r <- rd(u, run, "results"); if (is.null(r)) return(out)
  t <- r[r$status == "tested", ]
  if (m == "absrel") {
    b <- rd(u, run, "branches"); b <- b[b$orthogroup %in% t$orthogroup & !is.na(b$p_raw), ]
  # Truth per branch comes from the simulator (sim/branch_truth.tsv), not from the branches table's
  # `foreground` column, which is TRUE for every tested branch when aBSREL tests all branches.
    bt <- utils::read.delim(file.path(root, "units", u, "sim", "branch_truth.tsv"), stringsAsFactors = FALSE)
    fg <- b$species_tree_branches %in% bt$branch_id[bt$foreground]
    out[["absrel p_raw per branch"]] <- list(p = b$p_raw, og = b$orthogroup, fg = fg)
    out[["absrel p_bh_lineage per branch"]] <- list(p = b$p_bh_lineage, og = b$orthogroup, fg = fg)
    out[["absrel p_value per gene"]] <- list(p = t$p_value, og = t$orthogroup)
  } else if (m == "branch_omega") {
    out[["p_two_vs_m0"]] <- list(p = t$p_two_vs_m0, og = t$orthogroup)
    out[["p_fg_omega_gt1"]] <- list(p = t$p_fg_omega_gt1, og = t$orthogroup)
  } else if (m == "busted") out[["busted p_value"]] <- list(p = t$p_value, og = t$orthogroup)
  out
}
for (u in units) {
  genes <- utils::read.delim(file.path(root, "units", u, "sim", "genes.tsv"), stringsAsFactors = FALSE)
  sel <- genes$gene[genes$selected]
  for (run in basename(list.dirs(file.path(root, "units", u), recursive = FALSE))) {
    if (run == "sim") next
    for (nm in names(tt <- tests(u, run))) {
      x <- tt[[nm]]; ok <- !is.na(x$p)
      if (!grepl("^B-MIX", u)) {
        k <- sum(x$p[ok] <= 0.05); n <- sum(ok); ci <- cp(k, n)
        size[[length(size) + 1L]] <- data.frame(unit = u, run = run, masked = grepl("_mask$", run), test = nm, n = n, positives = k,
                                                rate = k / n, lo = ci[1], hi = ci[2], category = category(k, n))
      } else { # B2: discoveries after the documented correction; truth = selected gene (and, per branch, a foreground branch)
        p <- x$p[ok]; og <- x$og[ok]
        q <- if (grepl("p_bh_lineage", nm)) p else stats::p.adjust(p, "BH")
        disc <- q <= 0.05
        true <- og %in% sel & (if (is.null(x$fg)) TRUE else x$fg[ok] %in% TRUE)
        nd <- sum(disc); nf <- sum(disc & !true); ci <- cp(nf, nd)
        mixture[[length(mixture) + 1L]] <- data.frame(unit = u, run = run, test = nm, correction = if (grepl("p_bh_lineage", nm)) "BH per lineage (dndsR)" else "BH across genes",
          discoveries = nd, false = nf, fdp = if (nd) nf / nd else NA, lo = ci[1], hi = ci[2],
          power = if (any(true)) mean(disc[true]) else NA, true_tests = sum(true))
      }
    }
  }
  # B4: detection signals per gene, from the unmasked inputs and fits
  prof <- utils::read.delim(file.path(root, "units", u, "sim", "windows_profile.tsv"), stringsAsFactors = FALSE)
  prof$gene <- sub("\\.fasta$", "", prof$alignment)
  win <- tapply(prof$windows, prof$gene, function(v) any(v > 0))
  b <- rd(u, "branch_omega", "branches"); r <- rd(u, "branch_omega", "results")
  if (!is.null(b)) {
    term <- b[b$terminal %in% TRUE, ]
    ratio3 <- tapply(seq_len(nrow(term)), term$orthogroup, function(i) { d <- term$dS_m0[i]; any(d > 3 * stats::median(d), na.rm = TRUE) })
    hds <- tapply(b$flags, b$orthogroup, function(f) any(grepl("high_ds", f)))
    fsc <- tapply(b$flags, b$orthogroup, function(f) any(grepl("few_synonymous", f)))
    qc <- if (!is.null(r)) stats::setNames(nzchar(ifelse(is.na(r$qc_flags), "", r$qc_flags)), r$orthogroup) else NULL
    g <- sort(unique(b$orthogroup))
    detect[[u]] <- data.frame(unit = u, gene = g, dissimilar_window = as.vector(win[g]) %in% TRUE, high_ds = as.vector(hds[g]) %in% TRUE,
      few_synonymous = as.vector(fsc[g]) %in% TRUE, tip_ds_3x_median = as.vector(ratio3[g]) %in% TRUE, qc_flag = if (is.null(qc)) NA else as.vector(qc[g]) %in% TRUE)
  }
}
sz <- do.call(rbind, size)
utils::write.table(sz, paste0(prefix, "_size.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (length(mixture)) utils::write.table(do.call(rbind, mixture), paste0(prefix, "_mixture.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (length(detect)) {
  dt <- do.call(rbind, detect)
  utils::write.table(dt, paste0(prefix, "_detection_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  sig <- c("dissimilar_window", "high_ds", "few_synonymous", "tip_ds_3x_median", "qc_flag")
  base <- dt[grepl("^B-E0", dt$unit), ]
  det <- do.call(rbind, lapply(split(dt, dt$unit), function(x) data.frame(unit = x$unit[1], signal = sig,
    share_flagged = vapply(sig, function(s) mean(x[[s]], na.rm = TRUE), numeric(1)),
    share_flagged_E0 = vapply(sig, function(s) mean(base[[s]][base$unit == sub("^B-E[0-9]+-", "B-E0-", x$unit[1])], na.rm = TRUE), numeric(1)))))
  utils::write.table(det, paste0(prefix, "_detection.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
print(sz, row.names = FALSE)
if (length(mixture)) print(do.call(rbind, mixture), row.names = FALSE)
