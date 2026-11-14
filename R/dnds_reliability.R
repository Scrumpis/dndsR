# Per-pair reliability flags for pairwise dN/dS. They change no value: the dN/dS table is
# written exactly as before, and the flags go to a separate table. Criteria and evidence:
# docs/validation/estimators.md, "Reliability flags".

#' Flag thresholds for pairwise dN/dS
#'
#' \describe{
#'   \item{min_codons}{Codons compared (both sequences A, C, G, T, neither a stop) below this:
#'     \code{few_codons}. In simulation, undefined and wrong-side ratios are common at 100 codons.}
#'   \item{min_syn_differences}{NG86 synonymous differences below this: \code{few_synonymous_differences}.
#'     With fewer than 3 events, the Poisson relative error of dS alone exceeds 58 percent.}
#'   \item{ds_count_fraction}{A dS below this fraction of the observed synonymous p-distance
#'     (NG86 differences over sites), with at least one synonymous difference:
#'     \code{ds_below_observed}. The NG86 distance correction only raises dS above this p-distance,
#'     so a much lower dS is either a failed estimate or a method whose sites differ greatly
#'     from those of NG86. Check it with another method.}
#'   \item{extreme_ratio}{dN/dS above this: \code{extreme_ratio}.}
#'   \item{low_ds, high_ds, outlier_fold}{As in the calculate_dnds message: dS below 0.02
#'     (\code{low_ds}), above 0.8 (\code{high_ds}), above 3 times the table's median (\code{outlier_ds}).}
#' }
#' @keywords internal
.dnds_flag_limits <- list(min_codons = 100L, min_syn_differences = 3, ds_count_fraction = 0.5,
                          extreme_ratio = 10, low_ds = 0.02, high_ds = 0.8, outlier_fold = 3)

#' NG86 synonymous sites of each sense codon, changes to stop codons excluded
#' @keywords internal
.ng86_sites <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    code <- Biostrings::GENETIC_CODE
    sense <- names(code)[code != "*"]
    s <- vapply(sense, function(cd) {
      x <- strsplit(cd, "")[[1]]
      sum(vapply(1:3, function(p) {
        alt <- setdiff(c("A", "C", "G", "T"), x[p])
        mut <- vapply(alt, function(b) { y <- x; y[p] <- b; paste(y, collapse = "") }, "")
        aa <- code[mut]; ok <- aa != "*"
        if (any(ok)) mean(aa[ok] == code[[cd]]) else 0
      }, numeric(1)))
    }, numeric(1))
    cache <<- s
    s
  }
})

#' NG86 synonymous and nonsynonymous differences between two sense codons
#'
#' Averaged over the shortest pathways that avoid stop codons. Cached per pair.
#' @return \code{c(sd, nd)}.
#' @keywords internal
.ng86_codon_diffs <- local({
  cache <- new.env(parent = emptyenv())
  perms <- function(v) if (length(v) <= 1L) list(v) else
    unlist(lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))), recursive = FALSE)
  function(a, b) {
    key <- paste0(a, b)
    if (!is.null(hit <- cache[[key]])) return(hit)
    code <- Biostrings::GENETIC_CODE
    x <- strsplit(a, "")[[1]]; y <- strsplit(b, "")[[1]]
    pos <- which(x != y)
    out <- c(0, 0)
    if (length(pos)) {
      tot <- c(0, 0); n_paths <- 0
      for (o in perms(pos)) {
        cur <- x; aa <- code[[a]]; s <- 0; n <- 0; ok <- TRUE
        for (p in o) {
          cur[p] <- y[p]; nxt <- code[[paste(cur, collapse = "")]]
          if (nxt == "*") { ok <- FALSE; break }
          if (nxt == aa) s <- s + 1 else n <- n + 1
          aa <- nxt
        }
        if (ok) { tot <- tot + c(s, n); n_paths <- n_paths + 1 }
      }
      if (n_paths) out <- tot / n_paths
    }
    assign(key, out, envir = cache)
    out
  }
})

#' NG86 counts for an aligned codon pair
#'
#' Codons with a gap or any character other than A, C, G, T in either sequence, or a stop
#' in either, are not compared.
#' @return A list: \code{codons}, \code{syn_sites}, \code{syn_differences}, \code{nonsyn_differences}.
#' @keywords internal
.ng86_counts <- function(x, y) {
  starts <- seq(1L, nchar(x), by = 3L)
  cx <- substring(toupper(x), starts, starts + 2L); cy <- substring(toupper(y), starts, starts + 2L)
  sites <- .ng86_sites()
  ok <- cx %in% names(sites) & cy %in% names(sites)
  cx <- cx[ok]; cy <- cy[ok]
  d <- cx != cy
  diffs <- if (any(d)) rowSums(mapply(.ng86_codon_diffs, cx[d], cy[d])) else c(0, 0)
  list(codons = length(cx), syn_sites = sum((sites[cx] + sites[cy]) / 2),
       syn_differences = diffs[1], nonsyn_differences = diffs[2])
}

#' Reliability flags for every pair of a dN/dS table
#'
#' @param res The dN/dS table (query_id, subject_id, dS, dNdS).
#' @param aln_files The pairs' codon alignments (two sequences each).
#' @return A data frame: query_id, subject_id, codons_compared, syn_sites, syn_differences,
#'   nonsyn_differences, flags (semicolon-separated, empty when none).
#' @keywords internal
.dnds_pair_flags <- function(res, aln_files) {
  lim <- .dnds_flag_limits
  pairs <- Filter(Negate(is.null), lapply(aln_files, .read_pair_fasta))
  keys <- vapply(pairs, function(p) paste(names(p), collapse = "\r"), "")
  counts <- lapply(pairs, function(p) .ng86_counts(p[[1]], p[[2]]))
  i <- match(paste(res$query_id, res$subject_id, sep = "\r"), keys)
  get <- function(f) vapply(seq_along(i), function(k) if (is.na(i[k])) NA_real_ else as.numeric(counts[[i[k]]][[f]]), numeric(1))
  codons <- get("codons"); S <- get("syn_sites"); Sd <- get("syn_differences"); Nd <- get("nonsyn_differences")
  ds <- suppressWarnings(as.numeric(res$dS)); ratio <- suppressWarnings(as.numeric(res$dNdS))
  med <- stats::median(ds[is.finite(ds)])
  ps <- ifelse(is.finite(S) & S > 0, Sd / S, NA_real_)
  f <- list(
    no_alignment = is.na(i),
    few_codons = !is.na(codons) & codons < lim$min_codons,
    few_synonymous_differences = !is.na(Sd) & Sd < lim$min_syn_differences,
    ds_below_observed = is.finite(ds) & !is.na(Sd) & Sd >= 1 & is.finite(ps) & ds < lim$ds_count_fraction * ps,
    extreme_ratio = is.finite(ratio) & ratio > lim$extreme_ratio,
    low_ds = is.finite(ds) & ds < lim$low_ds,
    high_ds = is.finite(ds) & ds > lim$high_ds,
    outlier_ds = is.finite(ds) & is.finite(med) & ds > lim$outlier_fold * med
  )
  flags <- vapply(seq_len(nrow(res)), function(k) paste(names(f)[vapply(f, `[`, logical(1), k)], collapse = ";"), "")
  data.frame(query_id = res$query_id, subject_id = res$subject_id, codons_compared = codons,
             syn_sites = round(S, 4), syn_differences = round(Sd, 4), nonsyn_differences = round(Nd, 4),
             flags = flags, stringsAsFactors = FALSE)
}

#' Summarize the reliability flags in a message
#' @keywords internal
.report_pair_flags <- function(flags, label, path) {
  f <- unlist(strsplit(flags$flags[nzchar(flags$flags)], ";", fixed = TRUE))
  n <- function(x) sum(f == x)
  message(sprintf(paste0("[calculate_dnds] %s: per-pair reliability flags in %s. Values are unchanged. ",
                         "Fewer than %d codons compared: %d. Fewer than %g synonymous differences: %d. ",
                         "dS far below the observed synonymous differences: %d. dN/dS above %g: %d. No codon alignment: %d."),
                  label, basename(path), .dnds_flag_limits$min_codons, n("few_codons"), .dnds_flag_limits$min_syn_differences,
                  n("few_synonymous_differences"), n("ds_below_observed"), .dnds_flag_limits$extreme_ratio, n("extreme_ratio"),
                  n("no_alignment")))
  invisible(flags)
}
