# Checks each dnds_method of calculate_dnds against tools that share no code with dndsR, on
# the codon alignments dndsR itself used.
#   Rscript --vanilla dnds-methods.R <pairs.tsv> <tools bin> <output dir> [cores]
# pairs.tsv: comparison_name, query_seq, subject_seq (CDS FASTAs, as from extract_cds).
# <tools bin> holds KaKs_Calculator, yn00 (PAML), and python3 is used for ng86.py.
# Writes <output dir>/dnds_methods_pairs.tsv (one row per pair and method) and
# dnds_methods_summary.tsv.
#
# Checks:
# 1. Alignments: every method's run produces byte-identical codon alignments.
# 2. Reproduction: dndsR's values equal KaKs_Calculator run directly on those alignments,
#    wherever KaKs_Calculator gives a number. Where it gives NA, dndsR's 0 or NA is counted.
# 3. Independent: NG against PAML yn00 NG86 and ng86.py; YN against yn00 YN00; LWL against
#    yn00 LWL85; LPB against yn00 LPB93. Same codons as KaKs_Calculator uses (codons with a
#    gap, ambiguity, or stop in either sequence removed).
suppressPackageStartupMessages(library(dndsR))
a <- commandArgs(TRUE)
pairs_file <- a[1]; bin <- normalizePath(a[2]); out <- a[3]; cores <- if (length(a) > 3) as.integer(a[4]) else 2L
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
specs <- utils::read.delim(pairs_file, stringsAsFactors = FALSE)
methods <- dndsR:::.dnds_methods()
kaks_bin <- file.path(bin, "KaKs_Calculator")
work <- tempfile("dm-"); dir.create(work)
code <- Biostrings::GENETIC_CODE

read_pair <- dndsR:::.read_pair_fasta
clean <- function(p) {   # codons KaKs_Calculator compares, without stops
  s <- seq(1L, nchar(p[[1]]), 3L)
  cx <- substring(p[[1]], s, s + 2L); cy <- substring(p[[2]], s, s + 2L)
  ok <- grepl("^[ACGT]{3}$", cx) & grepl("^[ACGT]{3}$", cy)
  ok[ok] <- code[cx[ok]] != "*" & code[cy[ok]] != "*"
  c(paste(cx[ok], collapse = ""), paste(cy[ok], collapse = ""))
}
yn00 <- function(p) {
  d <- tempfile(tmpdir = work); dir.create(d); old <- setwd(d); on.exit(setwd(old))
  writeLines(c(sprintf(" 2 %d", nchar(p[1])), "q", p[1], "s", p[2]), "pair.phy")
  writeLines(c("seqfile = pair.phy", "outfile = yn.out", "verbose = 0", "icode = 0", "weighting = 0", "commonf3x4 = 0"), "yn00.ctl")
  system2(file.path(bin, "yn00"), "yn00.ctl", stdout = FALSE, stderr = FALSE)
  r <- readLines("yn.out")
  num <- function(x) as.numeric(regmatches(x, gregexpr("-?[0-9]+\\.[0-9]+", x))[[1]])
  ng <- num(grep("^s\\s", r[grep("Nei & Gojobori", r):length(r)], value = TRUE)[1])
  yn <- num(grep("^\\s*2\\s+1\\s+[0-9]", r, value = TRUE)[1])
  lwl <- num(grep("^LWL85:", r, value = TRUE)); lpb <- num(grep("^LPB93:", r, value = TRUE))
  c(NG_dN = ng[2], NG_dS = ng[3], YN_dN = yn[6], YN_dS = yn[8],
    LWL_dN = lwl[2], LWL_dS = lwl[1], LPB_dN = lpb[2], LPB_dS = lpb[1])
}

# Ratios of dndsR to the independent tool, over values above 0.005 (rounding matters below).
.ratios <- function(u, v) { ok <- is.finite(u) & is.finite(v) & v > 0.005; u[ok] / v[ok] }
.med_ratio <- function(u, v) { r <- .ratios(u, v); if (length(r)) stats::median(r) else NA_real_ }
.rng_ratio <- function(u, v) { r <- .ratios(u, v); if (length(r)) sprintf("%.3f-%.3f", min(r), max(r)) else NA_character_ }
rows <- list(); summ <- list()
for (i in seq_len(nrow(specs))) {
  cmp <- specs$comparison_name[i]
  res <- list(); aln <- list()
  for (m in methods) {
    od <- file.path(work, cmp, m)
    calculate_dnds(comparison_name = cmp, query_seq = specs$query_seq[i], subject_seq = specs$subject_seq[i],
                   output_dir = od, dnds_method = m, kaks_calculator = kaks_bin, comp_cores = cores, overwrite = TRUE)
    res[[m]] <- utils::read.delim(file.path(od, cmp, paste0(cmp, "_dnds.tsv")), stringsAsFactors = FALSE)
    files <- list.files(file.path(tempdir(), "_alignment", "codon_aln"), pattern = "_pal2nal\\.aln$", full.names = TRUE)
    aln[[m]] <- sort(unname(tools::md5sum(files)))
    if (m == "Comeron") {
      pairs <- Filter(Negate(is.null), lapply(files, read_pair))
      keys <- vapply(pairs, function(p) paste(names(p), collapse = "\t"), character(1))
    }
  }
  same_alignments <- all(vapply(aln, identical, logical(1), aln[["Comeron"]]))
  # Direct KaKs_Calculator, all methods at once, on the raw alignments.
  axt <- file.path(work, paste0(cmp, ".axt"))
  writeLines(unlist(lapply(seq_along(pairs), function(k) c(paste0("p", k), pairs[[k]][[1]], pairs[[k]][[2]], ""))), axt)
  system2(kaks_bin, c("-i", axt, "-o", paste0(axt, ".kaks"), unlist(lapply(dndsR:::.kaks_methods(), function(m) c("-m", m)))),
          stdout = FALSE, stderr = FALSE)
  direct <- utils::read.delim(paste0(axt, ".kaks"), check.names = FALSE, na.strings = "NA")
  # Independent tools on the cleaned alignments.
  cl <- lapply(pairs, clean)
  ind <- do.call(rbind, lapply(cl, function(p) if (nchar(p[1]) >= 3) yn00(p) else rep(NA_real_, 8)))
  fq <- file.path(work, "q.fa"); fs <- file.path(work, "s.fa")
  writeLines(unlist(lapply(seq_along(cl), function(k) c(paste0(">q_", k), cl[[k]][1]))), fq)
  writeLines(unlist(lapply(seq_along(cl), function(k) c(paste0(">s_", k), cl[[k]][2]))), fs)
  ng <- utils::read.delim(text = system2("python3", c(file.path(here, "ng86.py"), fq, fs), stdout = TRUE))
  ngi <- match(as.character(seq_along(cl)), sub("^.*_", "", ng[[1]]))
  for (m in setdiff(methods, "Comeron")) {
    r <- res[[m]]; k <- match(keys, paste(r$query_id, r$subject_id, sep = "\t"))
    d <- direct[direct$Method == m, ][match(paste0("p", seq_along(pairs)), direct$Sequence[direct$Method == m]), ]
    x <- data.frame(comparison = cmp, method = m, query_id = sub("\t.*", "", keys), subject_id = sub(".*\t", "", keys),
                    codons = nchar(vapply(cl, `[`, "", 1)) / 3, dN = r$dN[k], dS = r$dS[k], dNdS = r$dNdS[k],
                    direct_dN = d$Ka, direct_dS = d$Ks,
                    ind_dN = if (m %in% c("NG", "YN", "LWL", "LPB")) ind[, paste0(m, "_dN")] else NA_real_,
                    ind_dS = if (m %in% c("NG", "YN", "LWL", "LPB")) ind[, paste0(m, "_dS")] else NA_real_,
                    ng86py_dN = if (m == "NG") ng$dN[ngi] else NA_real_, ng86py_dS = if (m == "NG") ng$dS[ngi] else NA_real_)
    rows[[length(rows) + 1L]] <- x
    reld <- function(u, v) { ok <- is.finite(u) & is.finite(v); if (!any(ok)) NA_real_ else max(abs(u[ok] - v[ok]) / pmax(abs(v[ok]), 1e-3)) }
    absd <- function(u, v) { ok <- is.finite(u) & is.finite(v); if (!any(ok)) NA_real_ else max(abs(u[ok] - v[ok])) }
    summ[[length(summ) + 1L]] <- data.frame(
      comparison = cmp, method = m, pairs = nrow(x), same_alignments = same_alignments,
      direct_numeric = sum(is.finite(x$direct_dN)) + sum(is.finite(x$direct_dS)),
      direct_max_abs_diff = max(absd(x$dN, x$direct_dN), absd(x$dS, x$direct_dS)),
      direct_NA_to_zero = sum(is.na(x$direct_dN) & x$dN %in% 0) + sum(is.na(x$direct_dS) & x$dS %in% 0),
      direct_NA_kept = sum(is.na(x$direct_dN) & is.na(x$dN)) + sum(is.na(x$direct_dS) & is.na(x$dS)),
      ind_max_abs_diff_dN = absd(x$dN, x$ind_dN), ind_max_abs_diff_dS = absd(x$dS, x$ind_dS),
      # yn00 prints 4 decimals, and LWL-type estimates can be slightly negative.
      ind_zero_agree = sum(x$dN %in% 0 & abs(x$ind_dN) <= 1e-4, na.rm = TRUE) + sum(x$dS %in% 0 & abs(x$ind_dS) <= 1e-4, na.rm = TRUE),
      ind_zero_disagree = sum(x$dN %in% 0 & abs(x$ind_dN) > 1e-4, na.rm = TRUE) + sum(x$dS %in% 0 & abs(x$ind_dS) > 1e-4, na.rm = TRUE),
      ind_median_ratio_dN = .med_ratio(x$dN, x$ind_dN), ind_range_ratio_dN = .rng_ratio(x$dN, x$ind_dN),
      ind_median_ratio_dS = .med_ratio(x$dS, x$ind_dS), ind_range_ratio_dS = .rng_ratio(x$dS, x$ind_dS),
      ng86py_max_abs_diff = max(absd(x$dN, x$ng86py_dN), absd(x$dS, x$ng86py_dS)),
      spearman_ratio_ind = tryCatch(suppressWarnings(stats::cor(x$dNdS, ifelse(x$ind_dS > 0, x$ind_dN / x$ind_dS, NA),
                                                                method = "spearman", use = "complete.obs")), error = function(e) NA_real_)
    )
  }
}
x <- do.call(rbind, rows); s <- do.call(rbind, summ)
utils::write.table(x, file.path(out, "dnds_methods_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(s, file.path(out, "dnds_methods_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(s, digits = 3)
unlink(work, recursive = TRUE)
