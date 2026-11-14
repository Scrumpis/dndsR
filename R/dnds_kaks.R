# Pairwise dN/dS methods computed with KaKs_Calculator, on the codon alignments orthologr
# makes for the default method. See docs/validation/estimators.md ("Other methods").

#' dN/dS methods dndsR supports
#'
#' \code{"Comeron"} runs inside orthologr. The others run in KaKs_Calculator 2.0.
#' @keywords internal
.dnds_methods <- function() c("Comeron", .kaks_methods())

# KaKs_Calculator's GY and MLPB are left out: GY gives a near-zero dS and a large ratio
# where no synonymous change occurred, and MLPB has had no check against another tool.
#' @rdname dot-dnds_methods
#' @keywords internal
.kaks_methods <- function() c("NG", "LWL", "LPB", "MLWL", "YN", "MYN")

#' Find KaKs_Calculator and its version
#'
#' @param path Path to the executable, or \code{NULL} to search \code{PATH}.
#' @return A list with \code{bin} and \code{version}. Stops if it cannot run.
#' @keywords internal
.kaks_find <- function(path = NULL) {
  bin <- if (is.null(path)) Sys.which("KaKs_Calculator") else path
  if (!nzchar(bin) || !file.exists(bin)) {
    stop("KaKs_Calculator is needed for this dN/dS method but was not found",
         if (!is.null(path)) paste0(" at ", path), ". It is in the dndsR image; elsewhere install ",
         "kakscalculator2 (Bioconda) or give its path with kaks_calculator.", call. = FALSE)
  }
  out <- tryCatch(suppressWarnings(system2(bin, "-h", stdout = TRUE, stderr = TRUE)), error = function(e) character())
  v <- regmatches(out, regexpr("Version: *[^,]+", out))
  if (!any(grepl("KaKs_Calculator", out))) stop("Could not run KaKs_Calculator at ", bin, call. = FALSE)
  list(bin = unname(bin), version = if (length(v)) trimws(sub("Version:", "", v[1])) else NA_character_)
}

#' Read a two-sequence FASTA alignment
#' @keywords internal
.read_pair_fasta <- function(path) {
  l <- readLines(path, warn = FALSE)
  h <- grep("^>", l)
  if (length(h) != 2L) return(NULL)
  ends <- c(h[-1] - 1L, length(l))
  seqs <- vapply(1:2, function(i) toupper(paste(l[seq_len(ends[i] - h[i]) + h[i]], collapse = "")), character(1))
  names(seqs) <- trimws(sub("\\s.*$", "", sub("^>", "", l[h])))
  if (nchar(seqs[1]) != nchar(seqs[2]) || nchar(seqs[1]) %% 3L != 0L) return(NULL)
  seqs
}

#' Classify the mutational pathways between two sense codons
#'
#' Every order of the single-base changes between \code{a} and \code{b} is a pathway. Returns
#' \code{"same"}, \code{"syn"} (every step on every pathway keeps the amino acid),
#' \code{"nonsyn"} (every step changes it), \code{"mixed"}, or \code{"stop"} (a pathway passes
#' through a stop codon, which methods treat differently). Cached per codon pair.
#' @keywords internal
.codon_paths <- local({
  cache <- new.env(parent = emptyenv())
  perms <- function(v) if (length(v) <= 1L) list(v) else
    unlist(lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))), recursive = FALSE)
  function(a, b) {
    key <- paste0(a, b)
    if (!is.null(hit <- cache[[key]])) return(hit)
    code <- Biostrings::GENETIC_CODE
    x <- strsplit(a, "")[[1]]; y <- strsplit(b, "")[[1]]
    pos <- which(x != y)
    out <- if (!length(pos)) "same" else {
      steps <- unlist(lapply(perms(pos), function(o) {
        cur <- x; aa <- code[[a]]; r <- logical(length(o))
        for (j in seq_along(o)) {
          cur[o[j]] <- y[o[j]]; nxt <- code[[paste(cur, collapse = "")]]
          if (nxt == "*") return(NA)
          r[j] <- nxt == aa; aa <- nxt
        }
        r
      }))
      if (anyNA(steps)) "stop" else if (all(steps)) "syn" else if (!any(steps)) "nonsyn" else "mixed"
    }
    assign(key, out, envir = cache)
    out
  }
})

#' Whether a pair must have dN = 0 or dS = 0 under a KaKs_Calculator method
#'
#' Used where KaKs_Calculator writes "NA", which it writes both for a value of 0 and for an
#' undefined one. Codons with any character other than A, C, G, T in either sequence are
#' skipped, as KaKs_Calculator does, as are stop codons.
#'
#' For NG, YN, and MYN, which count differences along mutational pathways: \code{no_n} when
#' every step of every pathway between every compared codon pair keeps the amino acid, and
#' \code{no_s} when every step changes it. A pathway through a stop codon decides neither.
#' Identical proteins are not enough: TTA and CTC both encode Leu, but one pathway passes
#' through TTC (Phe).
#'
#' For LWL, LPB, and MLWL, which classify sites by degeneracy, a synonymous transversion at a
#' twofold site counts toward dN, and corrections can be negative, so only identical compared
#' codons give 0 (both values).
#' @param method A KaKs_Calculator method.
#' @return A list with logicals \code{no_n} and \code{no_s}.
#' @keywords internal
.pair_zero_changes <- function(x, y, method = "NG") {
  starts <- seq(1L, nchar(x), by = 3L)
  cx <- substring(x, starts, starts + 2L)
  cy <- substring(y, starts, starts + 2L)
  code <- Biostrings::GENETIC_CODE
  ok <- grepl("^[ACGT]{3}$", cx) & grepl("^[ACGT]{3}$", cy)
  cx <- cx[ok]; cy <- cy[ok]
  keep <- code[cx] != "*" & code[cy] != "*"
  cx <- cx[keep]; cy <- cy[keep]
  d <- cx != cy
  if (!any(d)) return(list(no_n = TRUE, no_s = TRUE))
  if (!method %in% c("NG", "YN", "MYN")) return(list(no_n = FALSE, no_s = FALSE))
  cls <- mapply(.codon_paths, cx[d], cy[d], USE.NAMES = FALSE)
  list(no_n = all(cls == "syn"), no_s = all(cls == "nonsyn"))
}

#' Estimate dN and dS with KaKs_Calculator for every codon alignment in a folder
#'
#' Each alignment holds one query and one subject sequence (orthologr's pal2nal output).
#' Pairs are named \code{p<k>} in the AXT file, so sequence names are never parsed.
#' KaKs_Calculator writes "NA" both for a value of 0 and for an undefined one; dndsR sets
#' dN or dS to 0 only when \code{.pair_zero_changes} shows it must be. KaKs_Calculator's YN
#' code writes 99 where it cannot estimate a distance (the cap PAML's yn00, from which it
#' derives, uses); values of 99 or more are set to NA, and their number is the attribute
#' \code{capped}. dN/dS is \code{dN / dS} when dS > 0, otherwise NA.
#'
#' @param aln_files Codon alignment files.
#' @param method A KaKs_Calculator method.
#' @param bin KaKs_Calculator executable.
#' @param work Folder for the AXT input and output.
#' @return A data frame: query_id, subject_id, dN, dS, dNdS.
#' @keywords internal
.kaks_estimate <- function(aln_files, method, bin, work) {
  empty <- data.frame(query_id = character(), subject_id = character(), dN = numeric(), dS = numeric(), dNdS = numeric())
  pairs <- Filter(Negate(is.null), lapply(aln_files, .read_pair_fasta))
  if (!length(pairs)) return(empty)
  ids <- paste0("p", seq_along(pairs))
  axt <- file.path(work, "pairs.axt"); out <- file.path(work, paste0("pairs_", method, ".kaks"))
  writeLines(unlist(lapply(seq_along(pairs), function(k) c(ids[k], pairs[[k]][[1]], pairs[[k]][[2]], ""))), axt)
  unlink(out)
  log <- suppressWarnings(system2(bin, c("-i", shQuote(axt), "-o", shQuote(out), "-m", method), stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  if (!is.null(status) && status != 0L || !file.exists(out)) {
    stop("KaKs_Calculator failed (method ", method, "):\n", paste(utils::tail(log, 10), collapse = "\n"), call. = FALSE)
  }
  k <- utils::read.delim(out, check.names = FALSE, na.strings = "NA", colClasses = "character")
  num <- function(v) suppressWarnings(as.numeric(v))
  ka <- num(k$Ka)[match(ids, k$Sequence)]
  ks <- num(k$Ks)[match(ids, k$Sequence)]
  na_ka <- is.na(ka); na_ks <- is.na(ks)
  capped <- sum(ka >= 99 | ks >= 99, na.rm = TRUE)
  ka[!na_ka & ka >= 99] <- NA; ks[!na_ks & ks >= 99] <- NA
  zero <- lapply(pairs, function(p) .pair_zero_changes(p[[1]], p[[2]], method))
  ka[na_ka & vapply(zero, `[[`, logical(1), "no_n")] <- 0
  ks[na_ks & vapply(zero, `[[`, logical(1), "no_s")] <- 0
  ka[!is.na(ka) & ka == 0] <- 0   # "-0"
  ks[!is.na(ks) & ks == 0] <- 0
  est <- data.frame(query_id = vapply(pairs, function(p) names(p)[1], character(1)),
                    subject_id = vapply(pairs, function(p) names(p)[2], character(1)),
                    dN = ka, dS = ks, dNdS = ifelse(!is.na(ks) & ks > 0, ka / ks, NA_real_))
  attr(est, "capped") <- capped
  est
}

#' Replace orthologr's Comeron values with another method's
#'
#' @param res orthologr's result table.
#' @param est \code{.kaks_estimate()} output.
#' @return \code{res} with dN, dS, and dNdS from \code{est}, matched by query and subject.
#'   Rows with no alignment get NA, and their number is the attribute \code{unmatched}.
#' @keywords internal
.dnds_replace_values <- function(res, est) {
  res <- as.data.frame(res)
  i <- match(paste(res$query_id, res$subject_id, sep = "\r"), paste(est$query_id, est$subject_id, sep = "\r"))
  for (col in c("dN", "dS", "dNdS")) res[[col]] <- est[[col]][i]
  attr(res, "unmatched") <- sum(is.na(i))
  res
}
