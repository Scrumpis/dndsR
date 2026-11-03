# Benchmarks dN/dS estimators against sequences with known truth, one grid cell per Slurm task.
#   Rscript estimator-benchmark.R <dir with yn00, codeml, KaKs_Calculator> <CDS FASTA for codon usage> <task> <out dir> [reps] [boot reps]
# Estimators: Comeron (the dndsR default, via orthologr), and, independently of dndsR, PAML yn00
# (YN00 and NG86), PAML codeml (pairwise ML), KaKs_Calculator (MLWL and LPB).
# Regimes: uniform, one omega per gene; sites, per-codon omega drawn from a gamma with shape 0.5
# (four equal-probability categories) and mean omega, so most sites are constrained.
# Truth is the model's dN and dS, from mutational opportunity, so their ratio is omega. For `sites`
# it is averaged over categories. Also records a percentile bootstrap interval, resampling codons,
# for the Comeron ratio. That is an investigation of coverage, not a dndsR feature.
a <- commandArgs(TRUE)
tools_dir <- a[1]; usage_fasta <- a[2]; task <- as.integer(a[3]); out <- a[4]
reps <- if (length(a) > 4L) as.integer(a[5]) else 200L
boots <- if (length(a) > 5L) as.integer(a[6]) else 100L
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(self), "codon-model.R"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
grid <- expand.grid(codons = c(100L, 200L, 400L, 800L, 1600L),
                    ds = c(0.005, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 3),
                    omega = c(0.1, 0.3, 0.6, 1, 1.5, 3),
                    regime = c("uniform", "sites"), stringsAsFactors = FALSE)
if (task > nrow(grid)) quit(save = "no")
x <- grid[task, ]
RNGversion("4.0.0"); set.seed(20261103L + task)
kappa <- 2
seqs <- as.character(Biostrings::readDNAStringSet(usage_fasta))
seqs <- seqs[nchar(seqs) >= 30L & nchar(seqs) %% 3L == 0L]
cod <- unlist(lapply(seqs, function(s) substring(s, seq(1, nchar(s) - 2, 3), seq(3, nchar(s), 3))))
counts <- table(factor(toupper(cod), levels = sense))
pi <- as.numeric(counts) / sum(counts)

# The tree scale t is shared by all categories, so a site's synonymous rate does not depend on
# its omega. Rates here are not rescaled per category, which keeps dS the same across categories.
omegas <- if (x$regime == "sites") {
  m <- stats::qgamma((2 * (1:4) - 1) / 8, shape = 0.5, rate = 0.5); x$omega * m / mean(m)
} else x$omega
raw <- function(w) { q <- rate_matrix(kappa, w, pi); q }
base <- flux(raw(1), pi, one_step) # scale so that t is substitutions per codon at omega = 1
syn_sites <- flux(raw(1), pi, one_step & syn) / base
non_sites <- 1 - syn_sites
per_t <- function(w) {
  q <- raw(w)
  c(dN = flux(q, pi, one_step & !syn) / base / (3 * non_sites), dS = flux(q, pi, one_step & syn) / base / (3 * syn_sites))
}
rates <- sapply(omegas, per_t)
t <- x$ds / mean(rates["dS", ])
true_dN <- t * mean(rates["dN", ]); true_dS <- t * mean(rates["dS", ])
tp <- lapply(omegas, function(w) {
  q <- raw(w) / base; diag(q) <- -rowSums(q)
  e <- eigen(q); p <- Re(e$vectors %*% (exp(e$values * t) * solve(e$vectors))); p[p < 0] <- 0; p
})

wd <- tempfile("bench"); dir.create(wd); setwd(wd)
comeron <- function(s, q) {
  writeLines(c(">q", s, ">s", q), "c.fa")
  r <- tryCatch(orthologr::substitutionrate("c.fa", est.method = "Comeron", quiet = TRUE), error = function(e) NULL)
  if (is.null(r)) c(NA, NA) else c(r$dN, r$dS)
}
yn00 <- function(s, q) {
  writeLines(c(sprintf(" 2 %d", nchar(s)), "q", s, "s", q), "p.phy")
  writeLines(c("seqfile = p.phy", "outfile = yn.out", "verbose = 0", "icode = 0", "weighting = 0", "commonf3x4 = 0"), "yn00.ctl")
  system2(file.path(tools_dir, "yn00"), "yn00.ctl", stdout = FALSE, stderr = FALSE)
  res <- readLines("yn.out")
  v <- as.numeric(strsplit(trimws(grep("^\\s*2\\s+1\\s+\\d", res, value = TRUE)[1]), "\\s+")[[1]])
  ng <- grep("\\(", res, value = TRUE); ng <- ng[grepl("^s\\s", ng)][1]
  g <- as.numeric(regmatches(ng, gregexpr("-?[0-9]+\\.[0-9]+", ng))[[1]])
  c(v[8], v[11], g[2], g[3])
}
codeml <- function(s, q) {
  writeLines(c(sprintf(" 2 %d", nchar(s)), "q", s, "s", q), "p.phy")
  writeLines(c("seqfile = p.phy", "outfile = ml.out", "noisy = 0", "verbose = 0", "runmode = -2", "seqtype = 1",
               "CodonFreq = 2", "model = 0", "NSsites = 0", "icode = 0", "fix_kappa = 0", "kappa = 2", "fix_omega = 0", "omega = 0.4", "cleandata = 1"), "codeml.ctl")
  system2(file.path(tools_dir, "codeml"), "codeml.ctl", stdout = FALSE, stderr = FALSE, input = "\n", timeout = 120)
  if (!file.exists("ml.out")) return(c(NA, NA))
  res <- grep("dN/dS", readLines("ml.out"), value = TRUE)
  file.remove("ml.out")
  if (!length(res)) return(c(NA, NA))
  dn <- as.numeric(sub(".*dN *= *([0-9.]+).*", "\\1", res[length(res)]))
  ds <- as.numeric(sub(".*dS *= *([0-9.]+).*", "\\1", res[length(res)]))
  c(dn, ds)
}
kaks <- function(s, q) {
  writeLines(c("q-s", s, q, ""), "p.axt")
  system2(file.path(tools_dir, "KaKs_Calculator"), c("-i", "p.axt", "-o", "p.kaks", "-m", "MLWL", "-m", "LPB"), stdout = FALSE, stderr = FALSE)
  r <- tryCatch(utils::read.delim("p.kaks"), error = function(e) NULL)
  if (is.null(r)) return(rep(NA_real_, 4))
  as.numeric(c(r$Ka[r$Method == "MLWL"], r$Ks[r$Method == "MLWL"], r$Ka[r$Method == "LPB"], r$Ks[r$Method == "LPB"]))
}
rows <- vector("list", reps)
for (r in seq_len(reps)) {
  anc <- sample.int(n, x$codons, TRUE, prob = pi)
  cat <- sample.int(length(omegas), x$codons, TRUE)
  der <- integer(x$codons)
  for (k in seq_along(omegas)) for (i in unique(anc[cat == k])) {
    w <- which(anc == i & cat == k); der[w] <- sample.int(n, length(w), TRUE, prob = tp[[k]][i, ])
  }
  s <- paste(sense[anc], collapse = ""); q <- paste(sense[der], collapse = "")
  ci <- c(NA, NA, NA)
  # A bootstrap per replicate is the costly part, so it runs on the first 50 replicates.
  if (boots > 0L && r <= 50L) {
    b <- vapply(seq_len(boots), function(i) {
      j <- sample.int(x$codons, x$codons, TRUE); e <- comeron(paste(sense[anc[j]], collapse = ""), paste(sense[der[j]], collapse = ""))
      if (is.finite(e[2]) && e[2] > 0) e[1] / e[2] else NA_real_
    }, 1)
    ci <- c(stats::quantile(b, c(0.025, 0.975), na.rm = TRUE, names = FALSE), mean(is.na(b)))
  }
  rows[[r]] <- data.frame(x, true_dN = true_dN, true_dS = true_dS, rep = r,
    t(setNames(c(comeron(s, q), yn00(s, q), codeml(s, q), kaks(s, q)),
               c("Comeron_dN", "Comeron_dS", "YN00_dN", "YN00_dS", "NG86_dN", "NG86_dS", "ML_dN", "ML_dS", "MLWL_dN", "MLWL_dS", "LPB_dN", "LPB_dS"))),
    boot_lo = ci[1], boot_hi = ci[2], boot_undefined = ci[3])
}
write.table(do.call(rbind, rows), file.path(out, sprintf("cell_%04d.tsv", task)), sep = "\t", quote = FALSE, row.names = FALSE)
