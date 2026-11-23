# Study B, "upstream limitation or dndsR defect" (docs/dev/lineage-robustness-study.md, fixed
# before running): for a method and error cell outside tolerance, the n genes with the smallest
# p-values are rerun directly in HyPhy or codeml, without dndsR code, from the alignment and tree
# dndsR wrote. Agreement (LRT within 1e-3 or a relative 1e-4) means the failure belongs to the
# model or engine on that input; disagreement is a dndsR defect.
#   Rscript --vanilla phylo-study-B-direct.R <unit dir> <run> <absrel|branch_omega|busted> <n> <out dir> <hyphy> <codeml>
# BUSTED (added for the remaining arms, 2026-11-23): alignment.fas and tree.nwk as dndsR wrote them
# (<run>/hyphy.tar.gz), HyPhy called with the settings in the run's provenance (C1: --srv Yes,
# --rates 3, --syn-rates 3, --multiple-hits Double, --error-sink No, --grid-size 250,
# --starting-points 3, CPU=1, RANDOM_SEED=1); the LRT is compared with the results' `lrt`.
# aBSREL: alignment.fas and tree.nwk as dndsR wrote them (kept in <run>/hyphy.tar.gz), HyPhy
# called as documented (CPU=1, RANDOM_SEED=1, --branches FG, --srv No, --multiple-hits None).
# Branch models: codeml's inputs are not kept (keep_runs = FALSE), so they are written here by
# hand: PHYLIP from the same alignment aBSREL received in the matching run (absrel or
# absrel_mask), the species tree unrooted with the foreground clade (t1, t2 and its stem) marked
# #1, and a control file from PAML's documentation with the settings in the provenance
# (F3x4, kappa 2, omega starts 0.4, 0.1, 2; best of the starts).
a <- commandArgs(TRUE)
unit <- a[1]; run <- a[2]; method <- a[3]; n <- as.integer(a[4]); out <- a[5]; hyphy <- a[6]; codeml <- a[7]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
agree <- function(x, y) abs(x - y) <= 1e-3 | abs(x - y) <= 1e-4 * abs(y)
untar1 <- function(tarball, member, dest) {
  st <- system2("tar", c("-xzf", tarball, "-C", dest, "--strip-components=2", member))
  if (st != 0) stop("cannot extract ", member)
}
absrel_run <- if (method %in% c("absrel", "busted")) run else sub("^branch_omega", "absrel", run)
tarball <- file.path(unit, absrel_run, "hyphy.tar.gz")

if (method == "busted") {
  r <- utils::read.delim(file.path(unit, run, paste0(run, "_results.tsv")), stringsAsFactors = FALSE)
  r <- r[r$status == "tested", ]; og <- head(r$orthogroup[order(r$p_value)], n)
  res <- list()
  for (g in og) {
    d <- file.path(out, g); dir.create(d, showWarnings = FALSE)
    untar1(tarball, file.path("hyphy", g, "alignment.fas"), d); untar1(tarball, file.path("hyphy", g, "tree.nwk"), d)
    old <- setwd(d)
    system2(hyphy, c("CPU=1", "ENV=RANDOM_SEED=1", "busted", "--alignment", "alignment.fas", "--tree", "tree.nwk",
                     "--branches", "FG", "--code", "Universal", "--srv", "Yes", "--rates", "3", "--syn-rates", "3",
                     "--multiple-hits", "Double", "--error-sink", "No", "--grid-size", "250", "--starting-points", "3",
                     "--output", "direct.json"), stdout = "direct.log", stderr = "direct.log")
    j <- jsonlite::fromJSON("direct.json", simplifyVector = FALSE); setwd(old)
    res[[g]] <- data.frame(unit = basename(unit), run = run, method = method, orthogroup = g, branch = "fg",
                           lrt_dndsr = r$lrt[r$orthogroup == g], lrt_direct = as.numeric(j[["test results"]][["LRT"]]))
  }
} else if (method == "absrel") {
  r <- utils::read.delim(file.path(unit, run, paste0(run, "_results.tsv")), stringsAsFactors = FALSE)
  b <- utils::read.delim(file.path(unit, run, paste0(run, "_branches.tsv")), stringsAsFactors = FALSE)
  r <- r[r$status == "tested", ]; og <- head(r$orthogroup[order(r$p_value)], n)
  res <- list()
  for (g in og) {
    d <- file.path(out, g); dir.create(d, showWarnings = FALSE)
    untar1(tarball, file.path("hyphy", g, "alignment.fas"), d); untar1(tarball, file.path("hyphy", g, "tree.nwk"), d)
    old <- setwd(d)
    system2(hyphy, c("CPU=1", "ENV=RANDOM_SEED=1", "absrel", "--alignment", "alignment.fas", "--tree", "tree.nwk",
                     "--branches", "FG", "--code", "Universal", "--srv", "No", "--multiple-hits", "None",
                     "--output", "direct.json"), stdout = "direct.log", stderr = "direct.log")
    j <- jsonlite::fromJSON("direct.json", simplifyVector = FALSE); setwd(old)
    ba <- j[["branch attributes"]][["0"]]
    lrt <- vapply(ba, function(x) if (is.null(x$LRT)) NA_real_ else as.numeric(x$LRT), numeric(1))
    bg <- b[b$orthogroup == g & !is.na(b$lrt), c("branch", "lrt")]
    m <- merge(bg, data.frame(branch = names(lrt), lrt_direct = lrt), by = "branch")
    res[[g]] <- data.frame(unit = basename(unit), run = run, method = method, orthogroup = g, branch = m$branch,
                           lrt_dndsr = m$lrt, lrt_direct = m$lrt_direct)
  }
} else {
  r <- utils::read.delim(file.path(unit, run, paste0(run, "_results.tsv")), stringsAsFactors = FALSE)
  r <- r[r$status == "tested" & is.finite(r$p_two_vs_m0), ]; og <- head(r$orthogroup[order(r$p_two_vs_m0)], n)
  # Species tree unrooted at its root; foreground (t1, t2) and stem. Written by hand per tree:
  # T8b (((t1,t2),(t3,t4)),((t5,t6),(t7,t8))); T16 (added 2026-11-23) the balanced 16-taxon tree.
  ntaxa <- length(utils::read.delim(file.path(unit, "sim", "branch_truth.tsv"))$branch_id) / 2 + 1
  if (ntaxa == 8) {
    tree0 <- "(((t1,t2),(t3,t4)),(t5,t6),(t7,t8));"
    tree2 <- "(((t1 #1,t2 #1) #1,(t3,t4)),(t5,t6),(t7,t8));"
  } else if (ntaxa == 16) {
    tree0 <- "((((t1,t2),(t3,t4)),((t5,t6),(t7,t8))),((t9,t10),(t11,t12)),((t13,t14),(t15,t16)));"
    tree2 <- "((((t1 #1,t2 #1) #1,(t3,t4)),((t5,t6),(t7,t8))),((t9,t10),(t11,t12)),((t13,t14),(t15,t16)));"
  } else stop("no hand-written tree for ", ntaxa, " taxa")
  ctl <- function(model, omega, fix) c("seqfile = alignment.phy", "treefile = tree.nwk", "outfile = out.mlc", "noisy = 0",
    "verbose = 0", "runmode = 0", "seqtype = 1", "CodonFreq = 2", "clock = 0", "model = " %s% model, "NSsites = 0",
    "icode = 0", "fix_kappa = 0", "kappa = 2", "fix_omega = " %s% fix, "omega = " %s% omega, "fix_alpha = 1",
    "alpha = 0", "ncatG = 1", "getSE = 0", "RateAncestor = 0", "cleandata = 0", "method = 0")
  `%s%` <- function(x, y) paste0(x, y)
  fit <- function(d, tree, model, omega, fix) {
    f <- file.path(d, sprintf("m%d_w%s_f%d", model, omega, fix)); dir.create(f, showWarnings = FALSE)
    file.copy(file.path(d, "alignment.phy"), f, overwrite = TRUE); writeLines(c(sprintf(" %d 1", ntaxa), tree), file.path(f, "tree.nwk"))
    writeLines(ctl(model, omega, fix), file.path(f, "codeml.ctl"))
    old <- setwd(f); system2(codeml, "codeml.ctl", stdout = "log", stderr = "log", stdin = "/dev/null"); setwd(old)
    x <- grep("^lnL\\(ntime:", readLines(file.path(f, "out.mlc")), value = TRUE)
    as.numeric(sub("^.*\\):\\s*(-?[0-9.]+).*$", "\\1", x[length(x)]))
  }
  res <- list()
  for (g in og) {
    d <- file.path(out, g); dir.create(d, showWarnings = FALSE)
    untar1(tarball, file.path("hyphy", g, "alignment.fas"), d)
    fa <- readLines(file.path(d, "alignment.fas")); h <- grep("^>", fa)
    seqs <- vapply(seq_along(h), function(i) paste(fa[(h[i] + 1):(if (i < length(h)) h[i + 1] - 1 else length(fa))], collapse = ""), "")
    writeLines(c(sprintf(" %d %d", length(h), nchar(seqs[1])), paste0(sub("^>", "", fa[h]), "  ", seqs)), file.path(d, "alignment.phy"))
    l0 <- max(vapply(c(0.4, 0.1, 2), function(w) fit(d, tree0, 0, w, 0), 0))
    l2 <- max(vapply(c(0.4, 0.1, 2), function(w) fit(d, tree2, 2, w, 0), 0))
    res[[g]] <- data.frame(unit = basename(unit), run = run, method = method, orthogroup = g, branch = "fg",
                           lrt_dndsr = r$lrt_two_vs_m0[r$orthogroup == g], lrt_direct = 2 * (l2 - l0))
  }
}
res <- do.call(rbind, res); res$agree <- agree(res$lrt_direct, res$lrt_dndsr)
utils::write.table(res, file.path(out, "direct_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s %s %s: %d of %d statistics agree (genes %d); max abs diff %.3g\n", basename(unit), run, method,
            sum(res$agree), nrow(res), length(unique(res$orthogroup)), max(abs(res$lrt_direct - res$lrt_dndsr))))
