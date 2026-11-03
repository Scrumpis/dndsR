# Simulates a whole genome's worth of sequence pairs with known omega, and estimates dN/dS from
# them with the real estimator, so that a gene's label comes from estimator noise, not a shuffle.
#   Rscript sim-genome.R <dnds_annot.tsv> <query CDS FASTA> <scenario> <seed> <out.tsv.gz> [kappa]
# Each gene keeps its own codons (as the ancestor), length, observed dS (as true dS), and terms.
# Scenarios choose omega, independently of any term unless a term is spiked:
#   long     omega drawn from observed dN/dS of long genes, whatever the gene's length
#   bylen    the same, scaled by the median dN/dS of the gene's length decile
#   spike_F_W  bylen, then fraction F of the genes in 5 chosen InterPro terms get omega W
# A chosen term's true status is in the output's comment line, so that enrichment can be scored.
# Limits: no indels, no recombination, one omega per gene, stationary codon usage.
a <- commandArgs(TRUE)
annot <- a[1]; cds <- a[2]; scenario <- a[3]; seed <- as.integer(a[4]); out <- a[5]
kappa <- if (length(a) > 5L) as.numeric(a[6]) else 2
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(self), "codon-model.R"))
RNGversion("4.0.0"); set.seed(seed)

d <- read.table(annot, sep = "\t", header = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE)
d$dNdS <- suppressWarnings(as.numeric(d$dNdS)); d$dS <- suppressWarnings(as.numeric(d$dS))
d <- d[is.finite(d$dNdS) & d$dNdS < 10 & is.finite(d$dS), ]
fa <- Biostrings::readDNAStringSet(cds)
names(fa) <- sub(" .*", "", names(fa))
d <- d[d$query_id %in% names(fa), ]
seqs <- toupper(as.character(fa[d$query_id]))
codon_ids <- function(x) {
  i <- match(substring(x, seq(1, nchar(x) - 2L, 3L), seq(3L, nchar(x), 3L)), sense)
  i[!is.na(i)] # drops the terminal stop and any ambiguous codon
}
ancestors <- lapply(seqs, codon_ids)
keep <- lengths(ancestors) >= 30L
d <- d[keep, ]; ancestors <- ancestors[keep]; ng <- nrow(d)
counts <- tabulate(unlist(ancestors), n); pi <- counts / sum(counts)

n_aa <- lengths(ancestors)
long <- d$dNdS[n_aa >= 600L & d$dS >= 0.02]
decile <- cut(rank(n_aa, ties.method = "first"), 10, labels = FALSE)
scale_bin <- tapply(d$dNdS, decile, stats::median) / stats::median(long)
omega <- sample(long, ng, TRUE)
spiked <- character(0); spiked_genes <- integer(0)
if (scenario != "long") omega <- omega * scale_bin[decile]
if (startsWith(scenario, "spike_")) {
  parts <- strsplit(scenario, "_", fixed = TRUE)[[1]]
  frac <- as.numeric(parts[2]); w <- as.numeric(parts[3])
  members <- split(rep(seq_len(ng), lengths(strsplit(ifelse(is.na(d$q_ipr), "", d$q_ipr), ";", fixed = TRUE))),
                   unlist(strsplit(ifelse(is.na(d$q_ipr), "", d$q_ipr), ";", fixed = TRUE)))
  size <- lengths(members)
  pool <- names(size)[size >= 20L & size <= 300L]
  spiked <- sample(pool, 5L)
  spiked_genes <- unique(unlist(lapply(members[spiked], function(g) sample(g, ceiling(frac * length(g))))))
  omega[spiked_genes] <- w
}
omega <- pmax(omega, 0.001)

# Cached by rounded omega: the eigen decomposition is the slow part.
grid <- exp(seq(log(0.001), log(20), length.out = 60L))
bin <- findInterval(omega, grid, all.inside = TRUE)
models <- vector("list", length(grid))
model <- function(k) {
  if (!is.null(models[[k]])) return(models[[k]])
  q <- rate_matrix(kappa, grid[k], pi); q <- q / sum(sweep(q, 1, pi, "*")); diag(q) <- -rowSums(q)
  e <- eigen(q)
  models[[k]] <<- list(v = e$vectors, l = e$values, inv = solve(e$vectors),
                       ds_per_t = truth(kappa, grid[k], 1, pi)[["dS_mut"]],
                       dn_per_t = truth(kappa, grid[k], 1, pi)[["dN_mut"]])
  models[[k]]
}
wd <- tempfile("sim"); dir.create(wd); setwd(wd)
est <- matrix(NA_real_, ng, 2L)
for (g in seq_len(ng)) {
  m <- model(bin[g])
  t <- pmax(d$dS[g], 0.002) / m$ds_per_t
  p <- Re(m$v %*% (exp(m$l * t) * m$inv)); p[p < 0] <- 0
  x <- ancestors[[g]]; y <- integer(length(x))
  for (i in unique(x)) { w <- which(x == i); y[w] <- sample.int(n, length(w), TRUE, prob = p[i, ]) }
  writeLines(c(">q", paste(sense[x], collapse = ""), ">s", paste(sense[y], collapse = "")), "p.fa")
  r <- tryCatch(orthologr::substitutionrate("p.fa", est.method = "Comeron", quiet = TRUE), error = function(e) NULL)
  if (!is.null(r)) est[g, ] <- c(r$dN, r$dS)
}
res <- data.frame(query_id = d$query_id, true_omega = signif(grid[bin], 4), true_dS = pmax(d$dS, 0.002),
                  n_codons = n_aa, dN = est[, 1], dS = est[, 2], spiked = seq_len(ng) %in% spiked_genes)
res$dNdS <- ifelse(res$dS > 0, res$dN / res$dS, NA_real_)
con <- gzfile(out, "w")
writeLines(paste0("# scenario=", scenario, " seed=", seed, " kappa=", kappa, " spiked_terms=", paste(spiked, collapse = ",")), con)
write.table(res, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
