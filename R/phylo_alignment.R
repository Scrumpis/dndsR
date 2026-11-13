# Codon alignment input and QC for the phylogenetic module. Method-agnostic.
#
# Policy (docs/phylo-methods.md): nothing is repaired or masked, except that a stop codon
# ending a sequence is replaced by gaps when `terminal_stops = "mask"`, and recorded.
# Hard failures exclude the orthogroup with a reason. Heuristic flags are reported only.

# Heuristic flag thresholds. Reported, not validated, and never used to exclude data.
.phylo_qc_limits <- list(
  ambiguous_codon_fraction = 0.05, # per sequence, of non-gap codons
  gap_fraction = 0.5,              # per sequence, of codon columns
  gappy_column_fraction = 0.3,     # of codon columns with gaps in more than half the sequences
  divergent_fold = 3,              # mean p-distance of a sequence over the median of all sequences
  max_p_distance = 0.35,           # largest pairwise nucleotide p-distance
  window_codons = 10L,             # window for locally dissimilar sequence
  window_min_codons = 8L,          # non-gap codons a window needs to be judged
  window_mismatch = 0.7            # share of a window's amino acids unlike every other sequence
)

#' Read an alignment FASTA into named character sequences
#'
#' Sequence names are the first word of each header. A name `taxon|gene` gives the taxon
#' and gene; a bare name is the taxon.
#' @keywords internal
.phylo_read_alignment <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  starts <- grep("^>", lines)
  if (!length(starts) || starts[1] != 1L) return(NULL)
  ends <- c(starts[-1] - 1L, length(lines))
  seqs <- mapply(function(s, e) {
    if (e <= s) return("")
    paste(gsub("[[:space:]]", "", lines[(s + 1L):e]), collapse = "")
  }, starts, ends)
  ids <- sub("^>[[:space:]]*([^[:space:]]*).*", "\\1", lines[starts])
  data.frame(
    seq_name = ids,
    taxon = sub("\\|.*$", "", ids),
    gene_id = ifelse(grepl("|", ids, fixed = TRUE), sub("^[^|]*\\|", "", ids), NA_character_),
    sequence = toupper(unname(seqs)),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.phylo_codon_table <- function() {
  gc <- Biostrings::GENETIC_CODE
  stats::setNames(as.character(gc), names(gc))
}

#' Check one codon alignment
#'
#' @param aln Data frame from [.phylo_read_alignment()].
#' @param taxa Taxa in the species tree.
#' @param min_taxa Minimum number of sequences.
#' @param terminal_stops `"mask"` or `"exclude"`.
#' @return A list with `ok`, `reason`, `detail`, `aln` (with masked sequences), `metrics`
#'   (one-row data frame), `sequences` (per-sequence metrics), and `flags`.
#' @keywords internal
.phylo_check_alignment <- function(aln, taxa, min_taxa = 4L, terminal_stops = c("mask", "exclude")) {
  terminal_stops <- match.arg(terminal_stops)
  out <- function(ok, reason = NA_character_, detail = NA_character_, metrics = NULL, sequences = NULL, flags = character()) {
    list(ok = ok, reason = reason, detail = detail, aln = aln, metrics = metrics, sequences = sequences, flags = flags)
  }
  if (is.null(aln) || !nrow(aln)) return(out(FALSE, "empty_alignment", "no FASTA records"))
  unknown <- setdiff(aln$taxon, taxa)
  if (length(unknown)) return(out(FALSE, "unknown_taxon", paste("not in the species tree:", paste(unknown, collapse = ", "))))
  dup <- unique(aln$taxon[duplicated(aln$taxon)])
  if (length(dup)) return(out(FALSE, "multi_copy", paste("more than one sequence from:", paste(dup, collapse = ", "))))
  len <- nchar(aln$sequence)
  if (length(unique(len)) != 1L) return(out(FALSE, "not_aligned", sprintf("sequence lengths range from %d to %d", min(len), max(len))))
  if (len[1] %% 3L != 0L || len[1] == 0L) return(out(FALSE, "length_not_multiple_of_3", sprintf("alignment length %d", len[1])))
  bad <- grepl("[^ACGTRYKMSWBDHVN?-]", aln$sequence)
  if (any(bad)) return(out(FALSE, "invalid_character", paste("characters other than nucleotides, IUPAC codes, '-' or '?' in:", paste(aln$seq_name[bad], collapse = ", "))))

  n_codon <- len[1] %/% 3L
  starts <- seq.int(1L, by = 3L, length.out = n_codon)
  codons <- t(vapply(aln$sequence, function(s) substring(s, starts, starts + 2L), character(n_codon), USE.NAMES = FALSE))
  if (n_codon == 1L) codons <- matrix(codons, ncol = 1L)
  gap_n <- nchar(gsub("[^-?]", "", codons))
  dim(gap_n) <- dim(codons)
  split <- gap_n > 0L & gap_n < 3L & nchar(gsub("[^-]", "", codons)) > 0L
  if (any(split)) {
    where <- which(split, arr.ind = TRUE)[1, ]
    return(out(FALSE, "split_codon_gap", sprintf("%d codons contain one or two gap characters, first in %s at codon %d",
                                                  sum(split), aln$seq_name[where[1]], where[2])))
  }
  is_gap <- gap_n == 3L
  empty <- rowSums(!is_gap & codons != "NNN") == 0L
  if (any(empty)) return(out(FALSE, "empty_sequence", paste("no determined codons in:", paste(aln$seq_name[empty], collapse = ", "))))

  table <- .phylo_codon_table()
  aa <- matrix(unname(table[codons]), nrow(codons)) # NA for gaps and ambiguous codons
  is_stop <- !is.na(aa) & aa == "*"
  last <- apply(!is_gap, 1L, function(r) max(which(r)))
  terminal <- is_stop[cbind(seq_len(nrow(codons)), last)]
  internal <- is_stop
  internal[cbind(seq_len(nrow(codons)), last)] <- FALSE
  if (any(internal)) {
    where <- which(internal, arr.ind = TRUE)[1, ]
    return(out(FALSE, "internal_stop", sprintf("%d stop codons before the end of a sequence, first %s in %s at codon %d",
                                                sum(internal), codons[where[1], where[2]], aln$seq_name[where[1]], where[2])))
  }
  n_masked <- sum(terminal)
  if (n_masked && terminal_stops == "exclude") {
    return(out(FALSE, "terminal_stop", paste("sequences ending in a stop codon:", paste(aln$seq_name[terminal], collapse = ", "))))
  }
  for (i in which(terminal)) {
    codons[i, last[i]] <- "---"
    is_gap[i, last[i]] <- TRUE
    aa[i, last[i]] <- NA
    aln$sequence[i] <- paste(codons[i, ], collapse = "")
  }
  if (nrow(aln) < min_taxa) {
    return(out(FALSE, "too_few_taxa", sprintf("%d sequences, fewer than min_taxa = %d", nrow(aln), min_taxa)))
  }

  # Reported metrics and heuristic flags.
  lim <- .phylo_qc_limits
  determined <- !is_gap
  ambiguous <- determined & grepl("[^ACGT]", codons)
  gap_frac <- rowMeans(is_gap)
  ambig_frac <- ifelse(rowSums(determined) > 0, rowSums(ambiguous) / pmax(rowSums(determined), 1), 0)
  gappy_cols <- mean(colMeans(is_gap) > 0.5)

  nt <- do.call(rbind, strsplit(aln$sequence, "", fixed = TRUE))
  acgt <- nt %in% c("A", "C", "G", "T")
  dim(acgt) <- dim(nt)
  n <- nrow(nt)
  pdist <- matrix(NA_real_, n, n)
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    both <- acgt[i, ] & acgt[j, ]
    if (any(both)) pdist[i, j] <- pdist[j, i] <- mean(nt[i, both] != nt[j, both])
  }
  mean_p <- rowMeans(pdist, na.rm = TRUE)
  med_p <- stats::median(mean_p, na.rm = TRUE)
  max_p <- suppressWarnings(max(pdist, na.rm = TRUE))

  windows <- integer(n)
  w <- lim$window_codons
  if (n_codon >= w) {
    unlike <- matrix(FALSE, n, n_codon)
    for (i in seq_len(n)) {
      others <- aa[-i, , drop = FALSE]
      unlike[i, ] <- !is.na(aa[i, ]) & vapply(seq_len(n_codon), function(k) {
        o <- others[, k]
        o <- o[!is.na(o)]
        length(o) > 0L && !(aa[i, k] %in% o)
      }, logical(1))
    }
    has_aa <- !is.na(aa)
    for (i in seq_len(n)) {
      for (s in seq_len(n_codon - w + 1L)) {
        idx <- s:(s + w - 1L)
        k <- sum(has_aa[i, idx])
        if (k >= lim$window_min_codons && sum(unlike[i, idx]) / k >= lim$window_mismatch) windows[i] <- windows[i] + 1L
      }
    }
  }

  flags <- character()
  if (any(ambig_frac > lim$ambiguous_codon_fraction)) flags <- c(flags, paste0("ambiguous_codons:", aln$taxon[ambig_frac > lim$ambiguous_codon_fraction]))
  if (any(gap_frac > lim$gap_fraction)) flags <- c(flags, paste0("gappy_sequence:", aln$taxon[gap_frac > lim$gap_fraction]))
  if (gappy_cols > lim$gappy_column_fraction) flags <- c(flags, "gappy_columns")
  if (n >= 4L && is.finite(med_p) && med_p > 0) {
    div <- which(mean_p > lim$divergent_fold * med_p)
    if (length(div)) flags <- c(flags, paste0("divergent_sequence:", aln$taxon[div]))
  }
  if (is.finite(max_p) && max_p > lim$max_p_distance) flags <- c(flags, "high_divergence")
  if (any(windows > 0L)) flags <- c(flags, paste0("dissimilar_window:", aln$taxon[windows > 0L]))
  if (n_masked) flags <- c(flags, sprintf("terminal_stops_masked:%d", n_masked))

  sequences <- data.frame(
    seq_name = aln$seq_name, taxon = aln$taxon, gene_id = aln$gene_id,
    gap_fraction = round(gap_frac, 6), ambiguous_codon_fraction = round(ambig_frac, 6),
    mean_p_distance = round(mean_p, 6), dissimilar_windows = windows,
    terminal_stop_masked = terminal, stringsAsFactors = FALSE
  )
  metrics <- data.frame(
    n_taxa = n, n_codons = n_codon, terminal_stops_masked = n_masked,
    max_gap_fraction = round(max(gap_frac), 6), max_ambiguous_codon_fraction = round(max(ambig_frac), 6),
    gappy_column_fraction = round(gappy_cols, 6), median_p_distance = round(med_p, 6),
    max_p_distance = round(if (is.finite(max_p)) max_p else NA_real_, 6),
    dissimilar_windows = sum(windows)
  )
  out(TRUE, metrics = metrics, sequences = sequences, flags = flags)
}

#' Write an alignment for HyPhy, named by taxon
#' @keywords internal
.phylo_write_alignment <- function(aln, path) {
  writeLines(as.vector(rbind(paste0(">", aln$taxon), aln$sequence)), path)
}
