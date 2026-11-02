# Compares orthologr's Comeron estimator, as used by dndsR, with independent tools on
# (a) the designed example and (b) sequences simulated under a codon model with known truth.
#   Rscript compare-estimators.R <dir with yn00 and KaKs_Calculator> <example dir> <output dir>
# Run in the dndsR image, which provides orthologr. The simulation uses base R only.
args <- commandArgs(trailingOnly = TRUE)
tools_dir <- args[1]; example <- args[2]; out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
RNGversion("4.0.0")
set.seed(4242)

bases <- c("T", "C", "A", "G")
codons <- apply(expand.grid(bases, bases, bases, stringsAsFactors = FALSE)[, 3:1], 1, paste0, collapse = "")
aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
code <- setNames(aa, codons)
sense <- codons[code != "*"]
n <- length(sense)

# Goldman-Yang style codon model with equal codon frequencies.
transition <- function(x, y) paste(sort(c(x, y)), collapse = "") %in% c("AG", "CT")
rate_matrix <- function(kappa, omega) {
  q <- matrix(0, n, n, dimnames = list(sense, sense))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    a <- strsplit(sense[i], "")[[1]]
    b <- strsplit(sense[j], "")[[1]]
    d <- which(a != b)
    if (length(d) != 1L) next
    r <- if (transition(a[d], b[d])) kappa else 1
    q[i, j] <- r * if (code[[sense[i]]] == code[[sense[j]]]) 1 else omega
  }
  q
}
syn_rate <- function(q) sum(q[outer(code[sense], code[sense], "==")]) / n

# True dN and dS for branch length t, in substitutions per codon. S and N come from the
# mutational opportunity at omega = 1, so dN/dS equals omega exactly.
truth <- function(kappa, omega, t) {
  q <- rate_matrix(kappa, omega)
  q1 <- rate_matrix(kappa, 1)
  scale_q <- sum(q) / n
  syn <- syn_rate(q) / scale_q
  non <- 1 - syn
  p_s <- syn_rate(q1) / (sum(q1) / n)
  c(dN = t * non / (3 * (1 - p_s)), dS = t * syn / (3 * p_s))
}

simulate_pair <- function(kappa, omega, t, n_codons) {
  q <- rate_matrix(kappa, omega)
  q <- q / (sum(q) / n)
  diag(q) <- -rowSums(q)
  e <- eigen(q)
  p <- Re(e$vectors %*% diag(exp(e$values * t)) %*% solve(e$vectors))
  p[p < 0] <- 0
  a <- sample.int(n, n_codons, TRUE)
  b <- vapply(a, function(i) sample.int(n, 1, prob = p[i, ]), integer(1))
  list(a = paste(sense[a], collapse = ""), b = paste(sense[b], collapse = ""))
}

# --- the three estimators on one aligned pair ------------------------------------------
comeron <- function(a, b, dir) {
  f <- file.path(dir, "pair.fasta")
  writeLines(c(">q", a, ">s", b), f)
  res <- orthologr::substitutionrate(f, est.method = "Comeron", quiet = TRUE)
  c(dN = res$dN, dS = res$dS)
}
yn00 <- function(a, b, dir) {
  old <- setwd(dir); on.exit(setwd(old))
  writeLines(c(sprintf(" 2 %d", nchar(a)), "q", a, "s", b), "pair.phy")
  writeLines(c("seqfile = pair.phy", "outfile = yn.out", "verbose = 0", "icode = 0", "weighting = 0", "commonf3x4 = 0"), "yn00.ctl")
  system2(file.path(tools_dir, "yn00"), "yn00.ctl", stdout = FALSE, stderr = FALSE)
  res <- readLines("yn.out")
  # (B) Yang & Nielsen (2000) table: S N t kappa omega dN +- SE dS +- SE
  row <- grep("^\\s*2\\s+1\\s+\\d", res, value = TRUE)[1]
  v <- as.numeric(strsplit(trimws(row), "\\s+")[[1]])
  yn <- c(dN = v[8], dS = v[11])
  # (A) Nei-Gojobori table: lower-triangle "w (dN dS)".
  ng_line <- grep("\\(", res, value = TRUE)
  ng_line <- ng_line[grepl("^s\\s", ng_line)][1]
  ng <- as.numeric(regmatches(ng_line, gregexpr("-?[0-9]+\\.[0-9]+", ng_line))[[1]])
  c(YN00_dN = yn[["dN"]], YN00_dS = yn[["dS"]], NG86_dN = ng[2], NG86_dS = ng[3])
}
kaks <- function(a, b, dir) {
  old <- setwd(dir); on.exit(setwd(old))
  writeLines(c("q-s", a, b, ""), "pair.axt")
  system2(file.path(tools_dir, "KaKs_Calculator"), c("-i", "pair.axt", "-o", "pair.kaks", "-m", "LWL", "-m", "MLWL", "-m", "LPB", "-m", "MYN"),
          stdout = FALSE, stderr = FALSE)
  r <- utils::read.delim("pair.kaks")
  out <- as.numeric(c(rbind(r$Ka, r$Ks)))
  names(out) <- c(rbind(paste0(r$Method, "_dN"), paste0(r$Method, "_dS")))
  out
}
all_estimates <- function(a, b, dir) {
  c(Comeron = setNames(comeron(a, b, dir), c("dN", "dS")), yn00(a, b, dir), kaks(a, b, dir))
}
strip_stop <- function(x) if (code[[substr(x, nchar(x) - 2L, nchar(x))]] == "*") substr(x, 1L, nchar(x) - 3L) else x

work <- tempfile("est"); dir.create(work)

# --- designed example --------------------------------------------------------------------
read_fa <- function(path) {
  l <- readLines(path); h <- grepl("^>", l)
  setNames(l[!h], sub("^>spA_|^>spB_", "", l[h]))
}
qa <- read_fa(file.path(example, "truth", "spA_longest_isoform_cds.fasta"))
qb <- read_fa(file.path(example, "truth", "spB_longest_isoform_cds.fasta"))
designed <- utils::read.delim(file.path(example, "truth", "designed_substitutions.tsv"))
keep <- designed$pair[!designed$class %in% c("internal_stop_subject", "length_not_multiple_of_3")]
rows <- lapply(keep, function(p) {
  data.frame(pair = p, class = designed$class[designed$pair == p],
             t(all_estimates(strip_stop(qa[[p]]), strip_stop(qb[[p]]), work)), check.names = FALSE)
})
ex <- do.call(rbind, rows)
utils::write.table(ex, file.path(out, "example_estimates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- simulation --------------------------------------------------------------------------
grid <- expand.grid(kappa = c(2, 4), omega = c(0.1, 0.5, 1, 2), t = c(0.1, 0.5, 1.5))
sims <- list()
for (g in seq_len(nrow(grid))) {
  cell <- grid[g, ]
  tr <- truth(cell$kappa, cell$omega, cell$t)
  for (rep in 1:20) {
    pair <- simulate_pair(cell$kappa, cell$omega, cell$t, 600L)
    est <- all_estimates(pair$a, pair$b, work)
    sims[[length(sims) + 1L]] <- data.frame(kappa = cell$kappa, omega = cell$omega, t = cell$t, rep = rep,
                                            true_dN = tr[["dN"]], true_dS = tr[["dS"]], t(est), check.names = FALSE)
  }
}
sim <- do.call(rbind, sims)
utils::write.table(sim, file.path(out, "simulation_estimates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

methods <- sub("_dN$", "", sub("\\.dN$", "", grep("dN$", setdiff(names(sim), "true_dN"), value = TRUE)))
dn_col <- function(m) if (m == "Comeron") "Comeron.dN" else paste0(m, "_dN")
ds_col <- function(m) if (m == "Comeron") "Comeron.dS" else paste0(m, "_dS")
summary_rows <- lapply(split(sim, list(sim$kappa, sim$omega, sim$t), drop = TRUE), function(d) {
  do.call(rbind, lapply(methods, function(m) {
    dn <- d[[dn_col(m)]]; ds <- d[[ds_col(m)]]
    data.frame(kappa = d$kappa[1], omega = d$omega[1], t = d$t[1], true_dN = d$true_dN[1], true_dS = d$true_dS[1], method = m,
               mean_dN = mean(dn, na.rm = TRUE), mean_dS = mean(ds, na.rm = TRUE),
               median_ratio = stats::median(dn / ds, na.rm = TRUE), frac_ratio_above_1 = mean(dn / ds > 1, na.rm = TRUE), usable = sum(is.finite(dn / ds)))
  }))
})
summ <- do.call(rbind, summary_rows)
rownames(summ) <- NULL
utils::write.table(summ, file.path(out, "simulation_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ[summ$method %in% c("Comeron", "NG86", "YN00", "MLWL", "LPB"), ], digits = 3)
cat("\n--- example: nonsynonymous-only and synonymous-only genes\n")
print(ex[ex$class %in% c("nonsynonymous_only", "synonymous_only", "identical"), ], digits = 3)
