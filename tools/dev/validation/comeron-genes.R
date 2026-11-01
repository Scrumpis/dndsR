# Per-gene behavior of dN/dS at realistic lengths and realistic codon usage: how often is a
# gene called above 1 when it evolves neutrally, and how often when omega is above 1?
#   Rscript comeron-genes.R <dir with yn00> <CDS FASTA for codon usage> <output dir> [replicates]
args <- commandArgs(trailingOnly = TRUE)
tools_dir <- args[1]; usage_fasta <- args[2]; out <- args[3]
reps <- if (length(args) >= 4) as.integer(args[4]) else 100L
dir.create(out, recursive = TRUE, showWarnings = FALSE)
RNGversion("4.0.0"); set.seed(20261102)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "codon-model.R"))
seqs <- as.character(Biostrings::readDNAStringSet(usage_fasta))
seqs <- seqs[nchar(seqs) >= 30L & nchar(seqs) %% 3L == 0L]
cod <- unlist(lapply(seqs, function(x) substring(x, seq(1, nchar(x) - 2, 3), seq(3, nchar(x), 3))))
counts <- table(factor(toupper(cod), levels = sense))
pi <- as.numeric(counts) / sum(counts)

tmp <- tempfile("genes"); dir.create(tmp)
one <- function(a, b) {
  old <- setwd(tmp); on.exit(setwd(old))
  writeLines(c(">q", a, ">s", b), "p.fa")
  r <- orthologr::substitutionrate(file.path(tmp, "p.fa"), est.method = "Comeron", quiet = TRUE)
  writeLines(c(sprintf(" 2 %d", nchar(a)), "q", a, "s", b), "p.phy")
  writeLines(c("seqfile = p.phy", "outfile = yn.out", "verbose = 0", "icode = 0", "weighting = 0", "commonf3x4 = 0"), "yn00.ctl")
  system2(file.path(tools_dir, "yn00"), "yn00.ctl", stdout = FALSE, stderr = FALSE)
  res <- readLines("yn.out")
  v <- as.numeric(strsplit(trimws(grep("^\\s*2\\s+1\\s+\\d", res, value = TRUE)[1]), "\\s+")[[1]])
  c(c_dN = r$dN, c_dS = r$dS, y_dN = v[8], y_dS = v[11])
}
grid <- expand.grid(codons = c(150L, 300L, 600L, 1500L), t = c(0.05, 0.2, 0.6, 1.5), omega = c(0.2, 0.5, 1, 1.5, 2, 4), kappa = c(2, 4))
rows <- list()
for (g in seq_len(nrow(grid))) {
  x <- grid[g, ]
  q <- rate_matrix(x$kappa, x$omega, pi); q <- q / sum(sweep(q, 1, pi, "*")); diag(q) <- -rowSums(q)
  e <- eigen(q); p <- Re(e$vectors %*% diag(exp(e$values * x$t)) %*% solve(e$vectors)); p[p < 0] <- 0
  est <- matrix(NA_real_, reps, 4)
  for (r in seq_len(reps)) {
    a <- sample.int(n, x$codons, TRUE, prob = pi)
    b <- vapply(a, function(i) sample.int(n, 1L, prob = p[i, ]), 1L)
    est[r, ] <- one(paste(sense[a], collapse = ""), paste(sense[b], collapse = ""))
  }
  cr <- ifelse(est[, 2] > 0, est[, 1] / est[, 2], NA); yr <- ifelse(est[, 4] > 0, est[, 3] / est[, 4], NA)
  tr <- truth(x$kappa, x$omega, x$t, pi)
  rows[[g]] <- data.frame(x, true_dS = tr[["dS_mut"]], undefined_ratio = mean(is.na(cr)),
    comeron_median = stats::median(cr, na.rm = TRUE), comeron_iqr = stats::IQR(cr, na.rm = TRUE),
    comeron_above1 = mean(cr > 1, na.rm = TRUE), yn00_median = stats::median(yr, na.rm = TRUE), yn00_above1 = mean(yr > 1, na.rm = TRUE),
    median_ds = stats::median(est[, 2], na.rm = TRUE))
}
d <- do.call(rbind, rows)
utils::write.table(d, file.path(out, "comeron_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(format(d[d$kappa == 2 & d$omega %in% c(1, 2), ], digits = 3), row.names = FALSE)
