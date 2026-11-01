# Why does the Comeron estimator read low in simulation? Long simulated pairs remove most
# sampling noise, so differences between estimators and truth are systematic. Varies transition
# bias, omega, divergence, and codon usage, and compares two definitions of truth:
#   mutational sites: S and N weighted by mutation rates (Goldman-Yang, Yang-Nielsen).
#                     True dN/dS then equals omega exactly.
#   physical sites:   S and N counted without regard to transition bias (kappa = 1).
#                     True dN/dS can differ from omega when kappa is not 1.
#   Rscript comeron-bias.R <dir with yn00, KaKs_Calculator> <CDS FASTA for codon usage> <output dir> [codons]
args <- commandArgs(trailingOnly = TRUE)
tools_dir <- args[1]; usage_fasta <- args[2]; out <- args[3]
n_codons <- if (length(args) >= 4) as.integer(args[4]) else 50000L
dir.create(out, recursive = TRUE, showWarnings = FALSE)
RNGversion("4.0.0"); set.seed(20261101)

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "codon-model.R"))

tmp <- tempfile("bias"); dir.create(tmp)
estimate <- function(pair) {
  old <- setwd(tmp); on.exit(setwd(old))
  writeLines(c(">q", pair$a, ">s", pair$b), "pair.fa")
  r <- orthologr::substitutionrate(file.path(tmp, "pair.fa"), est.method = "Comeron", quiet = TRUE)
  writeLines(c(sprintf(" 2 %d", nchar(pair$a)), "q", pair$a, "s", pair$b), "pair.phy")
  writeLines(c("seqfile = pair.phy", "outfile = yn.out", "verbose = 0", "icode = 0", "weighting = 0", "commonf3x4 = 0"), "yn00.ctl")
  system2(file.path(tools_dir, "yn00"), "yn00.ctl", stdout = FALSE, stderr = FALSE)
  res <- readLines("yn.out")
  v <- as.numeric(strsplit(trimws(grep("^\\s*2\\s+1\\s+\\d", res, value = TRUE)[1]), "\\s+")[[1]])
  ng <- grep("^s\\s", grep("\\(", res, value = TRUE), value = TRUE)[1]
  w <- as.numeric(regmatches(ng, gregexpr("-?[0-9]+\\.[0-9]+", ng))[[1]])
  writeLines(c("q-s", pair$a, pair$b, ""), "pair.axt")
  system2(file.path(tools_dir, "KaKs_Calculator"), c("-i", "pair.axt", "-o", "pair.kaks", "-m", "LWL", "-m", "MLWL", "-m", "LPB"), stdout = FALSE, stderr = FALSE)
  k <- utils::read.delim("pair.kaks")
  kv <- as.numeric(c(rbind(k$Ka, k$Ks)))
  names(kv) <- c(rbind(paste0(k$Method, "_dN"), paste0(k$Method, "_dS")))
  c(Comeron_dN = r$dN, Comeron_dS = r$dS, YN00_dN = v[8], YN00_dS = v[11], NG86_dN = w[2], NG86_dS = w[3], kv)
}

# Codon usage from a real CDS FASTA (public genome).
seqs <- as.character(Biostrings::readDNAStringSet(usage_fasta))
seqs <- seqs[nchar(seqs) >= 30L & nchar(seqs) %% 3L == 0L]
all_codons <- unlist(lapply(seqs, function(x) substring(x, seq(1, nchar(x) - 2, 3), seq(3, nchar(x), 3))))
counts <- table(factor(toupper(all_codons), levels = sense))
usage <- list(uniform = rep(1 / n, n), empirical = as.numeric(counts) / sum(counts))

grid <- expand.grid(usage = names(usage), kappa = c(1, 2, 4, 8), omega = c(0.1, 0.5, 1, 2), t = c(0.05, 0.2, 0.6, 1.5), stringsAsFactors = FALSE)
rows <- list()
for (g in seq_len(nrow(grid))) {
  x <- grid[g, ]; pi <- usage[[x$usage]]
  tr <- truth(x$kappa, x$omega, x$t, pi)
  est <- estimate(simulate_pair(x$kappa, x$omega, x$t, pi))
  rows[[g]] <- data.frame(x, t(tr), t(est))
}
d <- do.call(rbind, rows)
utils::write.table(d, file.path(out, "comeron_bias.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

d$true_ratio_mut <- d$dN_mut / d$dS_mut
d$true_ratio_phys <- d$dN_phys / d$dS_phys
for (m in c("Comeron", "YN00", "NG86", "LWL", "MLWL", "LPB")) d[[paste0(m, "_ratio")]] <- d[[paste0(m, "_dN")]] / d[[paste0(m, "_dS")]]
cols <- c("usage", "kappa", "omega", "t", "true_ratio_mut", "true_ratio_phys", "Comeron_ratio", "YN00_ratio", "NG86_ratio", "MLWL_ratio", "LPB_ratio")
print(format(d[d$omega %in% c(0.5, 1) & d$t %in% c(0.2, 0.6), cols], digits = 3), row.names = FALSE)
cat("\nMean of estimated ratio / omega, by kappa, codon usage, and method (t <= 0.6)\n")
s <- d[d$t <= 0.6, ]
print(stats::aggregate(cbind(Comeron = Comeron_ratio / omega, YN00 = YN00_ratio / omega, NG86 = NG86_ratio / omega, MLWL = MLWL_ratio / omega, LPB = LPB_ratio / omega) ~ usage + kappa, s, mean), digits = 3)
cat("\nMean of Comeron ratio / physical-site true ratio, by kappa and usage (t <= 0.6)\n")
print(stats::aggregate(cbind(Comeron = Comeron_ratio / true_ratio_phys, NG86 = NG86_ratio / true_ratio_phys) ~ usage + kappa, s, mean), digits = 3)
cat("\nMean relative error in dS and dN against mutational truth, Comeron and YN00 (t <= 0.6)\n")
print(stats::aggregate(cbind(Comeron_dS = Comeron_dS / dS_mut, Comeron_dN = Comeron_dN / dN_mut, YN00_dS = YN00_dS / dS_mut, YN00_dN = YN00_dN / dN_mut) ~ usage + kappa, s, mean), digits = 3)
