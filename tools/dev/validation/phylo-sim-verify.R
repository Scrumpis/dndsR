# Verifies the calibration simulator (phylo-sim.R) before any pilot, as required by
# docs/dev/phylo-calibration-design.md. Needs ape and PAML's codeml.
#   Rscript --vanilla phylo-sim-verify.R <codeml> <work dir>
a <- commandArgs(TRUE); codeml <- a[1]; work <- a[2]
dir.create(work, recursive = TRUE, showWarnings = FALSE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "phylo-sim.R"))
RNGversion("4.0.0"); set.seed(20261114)
nuc <- rbind(c(0.25, 0.20, 0.30, 0.25), c(0.30, 0.22, 0.28, 0.20), c(0.33, 0.17, 0.25, 0.25))
pi <- f3x4(nuc)
ok <- TRUE
report <- function(name, pass, detail) { cat(sprintf("[%s] %s: %s\n", if (pass) "PASS" else "FAIL", name, detail)); ok <<- ok && pass }

# 1. Matrices
worst <- 0
for (om in c(0.1, 1, 3)) for (r in c(0.3, 2.2)) for (d in c(0, 0.2)) {
  p <- p_matrix(2, om, r, d, pi, 0.3)
  worst <- max(worst, max(abs(rowSums(p) - 1)), max(abs(pi %*% p - pi)), max(abs(pi * p - t(pi * p))))
}
report("transition matrices are stochastic, stationary, reversible", worst < 1e-10, sprintf("largest error %.2g", worst))

# 2. codeml recovers omega, kappa, t on long pairwise alignments (no site variation):
# the mean of 10 replicates within 5% of the truth (amended; see the design's change log).
pair <- ape::read.tree(text = "(A:0.15,B:0.15);")
one <- list(rate = 1, weight = 1)
fit_pair <- function(aln, d) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  writeLines(c(sprintf("2 %d", 3L * ncol(aln)), paste0("A  ", paste(aln["A", ], collapse = "")), paste0("B  ", paste(aln["B", ], collapse = ""))), file.path(d, "seq.phy"))
  writeLines(c("seqfile = seq.phy", "outfile = out.txt", "runmode = -2", "seqtype = 1", "CodonFreq = 2",
               "model = 0", "NSsites = 0", "icode = 0", "fix_kappa = 0", "kappa = 2", "fix_omega = 0", "omega = 0.4",
               "cleandata = 0", "verbose = 0"), file.path(d, "codeml.ctl"))
  old <- setwd(d); on.exit(setwd(old)); system2(codeml, "codeml.ctl", stdout = FALSE, stderr = FALSE)
  out <- readLines("out.txt")
  # The line after the pairwise log-likelihood holds t, kappa, and omega.
  as.numeric(strsplit(trimws(out[max(grep("^lnL =", out)) + 1L]), " +")[[1]])
}
for (om in c(0.1, 0.5, 1, 2)) {
  est <- t(vapply(1:10, function(r) fit_pair(sim_gene(pair, 20000L, pi, 2, om, om, 1, one), file.path(work, sprintf("w%g_%d", om, r))), numeric(3)))
  m <- colMeans(est)
  true_t <- 0.3 * sum(-pi * diag(rate_q(2, om, 1, 0, pi)))
  pass <- abs(m[3] / om - 1) < 0.05 && abs(m[2] / 2 - 1) < 0.05 && abs(m[1] / true_t - 1) < 0.05
  report(sprintf("codeml recovers omega %g (mean of 10)", om), pass,
         sprintf("omega %.4f (sd %.4f), kappa %.3f, t %.4f (true %.4f)", m[3], stats::sd(est[, 3]), m[2], m[1], true_t))
}

# 3. Double changes: share of two-position differences on a short branch.
short <- ape::read.tree(text = "(A:0.0025,B:0.0025);")
aln <- sim_gene(short, 100000L, pi, 2, 0.5, 0.5, 1, one, delta = 0.2)
diff_pos <- vapply(seq_len(ncol(aln)), function(k) sum(strsplit(aln["A", k], "")[[1]] != strsplit(aln["B", k], "")[[1]]), integer(1))
n_diff <- sum(diff_pos > 0); n_two <- sum(diff_pos == 2)
q <- rate_q(2, 0.5, 1, 0.2, pi); off <- sweep(q, 1, pi, "*"); diag(off) <- 0
expected <- sum(off[ndiff == 2]) / sum(off)
ci <- stats::binom.test(n_two, n_diff)$conf.int
report("double-change share matches the rate matrix", expected >= ci[1] && expected <= ci[2],
       sprintf("%d of %d differences (%.4f, 95%% CI %.4f-%.4f); predicted %.4f", n_two, n_diff, n_two / n_diff, ci[1], ci[2], expected))

# 4. Indel bookkeeping.
par <- list(L50 = 300, tree = "(((At:0.1,Al:0.1):0.05,Cr:0.15):0.05,(Es:0.2,Br:0.2):0.05);", c = 1, kappa = 2,
            nuc_f3x4 = as.list(as.data.frame(t(nuc))), foreground = c("At", "Al"),
            species_topology = "(((At,Al),Cr),(Es,Br));", discordant_tree = "(((At:0.1,Cr:0.1):0.05,Al:0.15):0.05,(Es:0.2,Br:0.2):0.05);")
par$nuc_f3x4 <- lapply(1:3, function(i) nuc[i, ])
jsonlite::write_json(par, file.path(work, "params.json"), auto_unbox = TRUE, digits = NA)
run_cell(file.path(work, "params.json"), "M4_indel", 20L, 7L, file.path(work, "m4"))
cds <- read.delim(file.path(work, "m4", "orthogroup_table.tsv"))
good <- TRUE; gaps <- 0L
for (g in unique(cds$orthogroup)) {
  l <- readLines(file.path(work, "m4", "alignments", paste0(g, ".true.fasta")))
  s <- stats::setNames(l[c(FALSE, TRUE)], sub(">", "", l[c(TRUE, FALSE)]))
  gaps <- gaps + sum(grepl("---", s))
  for (tx in names(s)) {
    f <- readLines(file.path(work, "m4", "cds", paste0(tx, ".fasta")))
    want <- f[which(f == paste0(">", tx, "_", g)) + 1L]
    good <- good && identical(gsub("-", "", s[[tx]]), want) && nchar(want) %% 3L == 0L
  }
}
report("indel alignments, without gaps, equal the simulated CDS", good && gaps > 0, sprintf("%d sequences with gaps across 20 genes", gaps))
if (!ok) quit(status = 1)
cat("Simulator verified.\n")
