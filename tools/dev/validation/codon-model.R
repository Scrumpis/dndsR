# Codon substitution model shared by the estimator validation scripts. Base R only.
bases <- c("T", "C", "A", "G")
codons <- apply(expand.grid(bases, bases, bases, stringsAsFactors = FALSE)[, 3:1], 1, paste0, collapse = "")
aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
code <- setNames(aa, codons)
sense <- codons[code != "*"]; n <- length(sense)
syn <- outer(code[sense], code[sense], "==")
one_step <- matrix(FALSE, n, n); transition <- matrix(FALSE, n, n)
for (i in seq_len(n)) for (j in seq_len(n)) {
  a <- strsplit(sense[i], "")[[1]]; b <- strsplit(sense[j], "")[[1]]; d <- which(a != b)
  if (length(d) == 1L) {
    one_step[i, j] <- TRUE
    transition[i, j] <- paste(sort(c(a[d], b[d])), collapse = "") %in% c("AG", "CT")
  }
}

rate_matrix <- function(kappa, omega, pi) {
  q <- matrix(0, n, n)
  q[one_step] <- 1
  q[transition] <- kappa
  q[one_step & !syn] <- q[one_step & !syn] * omega
  q <- sweep(q, 2, pi, "*")
  q
}
flux <- function(q, pi, keep) sum(sweep(q, 1, pi, "*")[keep])

truth <- function(kappa, omega, t, pi) {
  q <- rate_matrix(kappa, omega, pi)
  total <- flux(q, pi, one_step)
  syn_share <- flux(q, pi, one_step & syn) / total # substitutions that are synonymous
  # Proportion of mutational opportunity that is synonymous, at omega = 1.
  q1 <- rate_matrix(kappa, 1, pi)
  p_mut <- flux(q1, pi, one_step & syn) / flux(q1, pi, one_step)
  # Physical sites: the same at kappa = 1 and omega = 1.
  q0 <- rate_matrix(1, 1, pi)
  p_phys <- flux(q0, pi, one_step & syn) / flux(q0, pi, one_step)
  c(dN_mut = t * (1 - syn_share) / (3 * (1 - p_mut)), dS_mut = t * syn_share / (3 * p_mut),
    dN_phys = t * (1 - syn_share) / (3 * (1 - p_phys)), dS_phys = t * syn_share / (3 * p_phys))
}

simulate_pair <- function(kappa, omega, t, pi) {
  q <- rate_matrix(kappa, omega, pi)
  q <- q / sum(sweep(q, 1, pi, "*"))
  diag(q) <- -rowSums(q)
  e <- eigen(q)
  p <- Re(e$vectors %*% diag(exp(e$values * t)) %*% solve(e$vectors))
  p[p < 0] <- 0
  a <- sample.int(n, n_codons, TRUE, prob = pi)
  b <- integer(n_codons)
  for (i in unique(a)) {
    w <- which(a == i)
    b[w] <- sample.int(n, length(w), TRUE, prob = p[i, ])
  }
  list(a = paste(sense[a], collapse = ""), b = paste(sense[b], collapse = ""))
}

