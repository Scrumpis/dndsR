# McDonald-Kreitman (MK) initial implementation.
#
# One MK table is written per comparison, beside the dN/dS outputs, when the comparison
# carries `snp_file` and `snp_side` columns:
#   <comp>_mk.tsv        one row per 1:1 pair (query_id, subject_id)
#   <comp>_mk_pooled.tsv one row summed over all pairs
#
# Method, in brief (see docs/mcdonald-kreitman.md for the full description):
#  * Sites are classified by their allele, as in the original McDonald-Kreitman
#    test: a polymorphism is Ps when its ALT leaves the amino acid unchanged and Pn
#    when it changes it; a fixed difference is Ds or Dn by the same rule. VCF gives
#    the ALT allele; a BED file gives it in column 4.
#  * Dn/Ds are counted from the same pal2nal codon alignments orthologr used, one
#    codon at a time, with multi-site codons averaged over all shortest mutational
#    paths (equal weight). Codons with gaps, non-ACGT bases, a stop codon on either
#    side, or a SNP in the reference species are excluded from the counts.
#  * Pn/Ps are counted from bi-allelic SNP positions mapped onto the reference
#    species' CDS. Indels, multi-allelic records, and symbolic alleles are ignored.
#    dndsR does no SNP filtering: filter calls (PASS, allele count, depth) beforehand.
#  * alpha and the neutrality index use the path-averaged Dn/Ds; Fisher's exact test
#    uses their rounded (integer) values, since exact tests need counts.

# ---------- codon utilities ----------

.mk_aa_lookup <- local({
  env <- new.env(parent = emptyenv())
  code <- Biostrings::GENETIC_CODE
  for (nm in names(code)) assign(nm, unname(code[[nm]]), envir = env)
  env
})

.mk_aa <- function(codon) {
  x <- .mk_aa_lookup[[codon]]
  if (is.null(x)) "X" else x
}

.mk_classify_cache <- new.env(parent = emptyenv())

.mk_permutations <- function(v) {
  if (length(v) <= 1L) return(list(v))
  out <- list()
  for (i in seq_along(v)) {
    for (rest in .mk_permutations(v[-i])) out[[length(out) + 1L]] <- c(v[i], rest)
  }
  out
}

#' Classify a codon pair as synonymous and nonsynonymous differences
#'
#' Allele-aware, as the original McDonald-Kreitman test is: a difference counts as
#' synonymous when it does not change the amino acid of the codon, nonsynonymous
#' when it does. Endpoints must be pure ACGT and neither may encode a stop codon;
#' callers filter those out. A codon pair that differs at more than one position is
#' split over every shortest mutational path, each step compared with the codon it
#' starts from; paths through a stop codon are discarded. Returns c(synonymous,
#' nonsynonymous) averaged over the surviving paths, summing to the number of
#' differing positions.
#' @keywords internal
.mk_classify_codon <- function(a, b) {
  if (a == b) return(c(0, 0))
  xa <- strsplit(a, "", fixed = TRUE)[[1]]
  xb <- strsplit(b, "", fixed = TRUE)[[1]]
  diff <- which(xa != xb)
  if (length(diff) == 1L) return(if (.mk_aa(a) == .mk_aa(b)) c(1, 0) else c(0, 1))
  key <- paste(a, b, sep = ":")
  cached <- .mk_classify_cache[[key]]
  if (!is.null(cached)) return(cached)
  syn <- 0; nonsyn <- 0; ok <- 0L
  for (p in .mk_permutations(diff)) {
    cur <- xa; aa_cur <- .mk_aa(a); s <- 0; n <- 0; ok_path <- TRUE
    for (pos in p) {
      cur[pos] <- xb[pos]
      aa_next <- .mk_aa(paste(cur, collapse = ""))
      if (aa_next == "*") { ok_path <- FALSE; break }
      if (aa_next == aa_cur) s <- s + 1 else n <- n + 1
      aa_cur <- aa_next
    }
    if (ok_path) { syn <- syn + s; nonsyn <- nonsyn + n; ok <- ok + 1L }
  }
  res <- if (ok == 0L) c(0, 0) else c(syn / ok, nonsyn / ok)
  assign(key, res, envir = .mk_classify_cache)
  res
}

# ---------- SNP input ----------

# Count header lines (lines starting with #, track, or browser) and capture the
# contig lengths embedded in a VCF header, if any. Also reports whether any
# non-header line exists, so callers can tolerate an empty data section.
.mk_header_info <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  n <- 0L
  contig_len <- NULL
  has_body <- FALSE
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (!length(ln)) break
    if (grepl("^##contig=", ln)) {
      m <- regmatches(ln, regexec("ID=([^,>]+).*?length=([0-9]+)", ln))[[1]]
      if (length(m) == 3L) contig_len[m[2]] <- as.numeric(m[3])
      n <- n + 1L; next
    }
    if (grepl("^#", ln) || grepl("^(track|browser)\\b", ln, ignore.case = TRUE)) { n <- n + 1L; next }
    has_body <- TRUE
    break
  }
  list(n_header = n, contig_len = contig_len, has_body = has_body)
}

# Empty SNP table carrying the format metadata.
.mk_empty_snps <- function(fmt, contig_len = NULL) {
  out <- data.frame(chrom = character(0), pos = integer(0),
                    ref = character(0), alt = character(0), stringsAsFactors = FALSE)
  list(snps = out,
       meta = list(format = fmt, n_input = 0L, n_kept = 0L, n_dropped = 0L,
                   n_duplicates = 0L, contig_len = contig_len))
}

#' Read a SNP position file (VCF or BED)
#'
#' `snp_file` may be a VCF (optionally bgzipped/gzipped) or a BED file of single-base
#' intervals. For VCF, only bi-allelic SNPs are kept: REF is a single ACGT base and there is
#' exactly one ALT allele, also a single ACGT base. Multi-allelic records, indels, and
#' symbolic alleles are dropped. For BED, only intervals of exactly one base are kept;
#' coordinates are converted from 0-based half-open to 1-based. Alleles are not available
#' for BED input, so sites are classified by degeneracy (see .mk_poly_counts).
#'
#' @return list with `pos` (1-based positions), `chrom`, `ref`, `alt` (single uppercase
#'   ALT base, `NA` for BED), plus counts in `meta`.
#' @keywords internal
.mk_read_snps <- function(path) {
  if (!file.exists(path)) stop("snp_file not found: ", path, call. = FALSE)
  lower <- tolower(path)
  if (grepl("\\.bcf$", lower)) {
    stop("BCF is not supported; write a (bgzipped) VCF instead: ", path, call. = FALSE)
  }
  is_vcf <- grepl("\\.vcf(\\.gz)?$", lower)
  is_bed <- grepl("\\.bed(\\.gz)?$", lower)
  if (!is_vcf && !is_bed) {
    stop("snp_file must be a .vcf, .vcf.gz, .bed, or .bed.gz file: ", path, call. = FALSE)
  }
  hd <- .mk_header_info(path)
  if (!isTRUE(hd$has_body)) {
    return(.mk_empty_snps(if (is_vcf) "vcf" else "bed", hd$contig_len))
  }

  if (is_vcf) {
    raw <- data.table::fread(path, skip = hd$n_header, header = FALSE, sep = "\t",
                             select = 1:5, showProgress = FALSE, data.table = FALSE,
                             colClasses = list(character = 1:5))
    if (ncol(raw) < 5L) stop("VCF has fewer than 5 columns: ", path, call. = FALSE)
    names(raw) <- c("chrom", "pos", "id", "ref", "alt")
    pos <- suppressWarnings(as.integer(raw$pos))
    bad_pos <- is.na(pos) | pos < 1L
    ref <- toupper(raw$ref)
    alt <- toupper(raw$alt)
    ref_ok <- grepl("^[ACGT]$", ref)
    # Only bi-allelic SNPs are used at all: REF is one ACGT base and there is exactly one
    # ALT allele, also one ACGT base, and REF != ALT. Multi-allelic records, indels,
    # symbolic alleles, and REF == ALT records are dropped and counted in the tally.
    biallelic <- grepl("^[ACGT]$", alt) & alt != ref
    drop_n <- sum(bad_pos | !ref_ok | !biallelic)
    keep <- !bad_pos & ref_ok & biallelic
    out <- data.frame(chrom = raw$chrom[keep], pos = pos[keep],
                      ref = ref[keep], alt = alt[keep], stringsAsFactors = FALSE)
    fmt <- "vcf"
    contig_len <- hd$contig_len
  } else {
    raw <- data.table::fread(path, skip = hd$n_header, header = FALSE, sep = "\t",
                             showProgress = FALSE, data.table = FALSE,
                             colClasses = "character")
    if (ncol(raw) < 3L) stop("BED file needs at least 3 columns: ", path, call. = FALSE)
    # Counting needs the allele: a BED file supplies it in column 4, which standard
    # BED calls the name. Positions without it cannot be classed as Pn or Ps.
    if (ncol(raw) < 4L) {
      stop("A BED SNP file needs a fourth column with the alternative allele, one of ",
           "A, C, G, T. Standard BED column 4 (the name) is read as that allele here, so ",
           "write chrom, start, end, alt.", call. = FALSE)
    }
    if (ncol(raw) > 4L) raw <- raw[, 1:4, drop = FALSE]
    names(raw) <- c("chrom", "start", "end", "alt")
    start <- suppressWarnings(as.integer(raw$start)); end <- suppressWarnings(as.integer(raw$end))
    alt <- toupper(trimws(raw$alt))
    if (nrow(raw) && !any(grepl("^[ACGT]$", alt))) {
      stop("A BED SNP file needs a fourth column with the alternative allele, one of ",
           "A, C, G, T; none of the ", nrow(raw), " rows has one. Standard BED column 4 ",
           "(the name) is read as that allele here, so write chrom, start, end, alt.",
           call. = FALSE)
    }
    ok <- !is.na(start) & !is.na(end) & start >= 0L & end - start == 1L &
      grepl("^[ACGT]$", alt)
    drop_n <- sum(!ok)
    out <- data.frame(chrom = raw$chrom[ok], pos = start[ok] + 1L,
                      ref = NA_character_, alt = alt[ok], stringsAsFactors = FALSE)
    fmt <- "bed"
    contig_len <- NULL
  }
  if (nrow(out)) {
    dup <- duplicated(paste(out$chrom, out$pos))
    if (any(dup)) out <- out[!dup, , drop = FALSE]
  } else {
    dup <- logical(0)
  }
  list(snps = out,
       meta = list(format = fmt, n_input = nrow(out) + drop_n + sum(dup),
                   n_kept = nrow(out), n_dropped = drop_n + sum(dup),
                   n_duplicates = sum(dup), contig_len = contig_len))
}

# ---------- CDS coordinate map ----------

# Build, for one genome, per-gene CDS coordinate maps identical to the CDS FASTAs
# written by extract_cds(): the same isoform choice, the same CDS order, gaps removed.
.mk_cds_map <- function(gff_path, genome_fasta, cds_fasta, by = "gene") {
  if (!file.exists(gff_path)) stop("GFF not found: ", gff_path, call. = FALSE)
  if (!file.exists(genome_fasta)) stop("Genome FASTA not found: ", genome_fasta, call. = FALSE)
  if (!file.exists(cds_fasta)) {
    stop("CDS FASTA not found: ", cds_fasta,
         "\nRun extract_cds() for this comparison before calculate_dnds().", call. = FALSE)
  }
  if (!file.exists(paste0(genome_fasta, ".fai"))) Rsamtools::indexFa(genome_fasta)
  fa <- Rsamtools::FaFile(genome_fasta)
  chroms <- GenomeInfoDb::seqnames(Rsamtools::scanFaIndex(fa))

  tg <- .tx_cds_groups(gff_path, by = by, verbose = FALSE)
  tg <- .keep_valid_groups(tg, chroms)
  all_cds <- GenomicFeatures::extractTranscriptSeqs(fa, tg$cds)
  q <- .cds_quality(all_cds)
  valid <- q$multiple_of_3 & !is.na(q$internal_stop) & !q$internal_stop
  chosen <- .select_longest_isoform(names(all_cds), tg$gene, Biostrings::width(all_cds), valid)
  cds <- all_cds[match(chosen$tx, names(all_cds))]
  names(cds) <- chosen$gene

  # Intersect with what extract_cds actually wrote; anything else cannot be mapped.
  written <- Biostrings::readDNAStringSet(cds_fasta)
  keep <- names(cds) %in% names(written)
  dropped <- names(cds)[!keep]
  if (length(dropped)) {
    message("  - MK: ", length(dropped), " CDS selected from the GFF are not in ",
            basename(cds_fasta), " and were dropped (e.g. ",
            paste(utils::head(dropped, 3), collapse = ", "), ").")
    cds <- cds[keep]
    chosen <- chosen[match(names(cds), chosen$gene), , drop = FALSE]
  }

  # Per gene: CDS-ordered integer vector of genomic positions, plus strand/chrom/parts
  pos_list <- vector("list", length(cds))
  strand_v <- character(length(cds))
  chrom_v <- character(length(cds))
  gr_parts <- vector("list", length(cds))
  names(pos_list) <- names(cds)
  for (i in seq_along(cds)) {
    gr <- tg$cds[[chosen$tx[i]]]
    strand_v[i] <- as.character(GenomicRanges::strand(gr))[1]
    chrom_v[i] <- as.character(GenomicRanges::seqnames(gr))[1]
    ord <- order(BiocGenerics::start(gr))
    if (strand_v[i] == "-") ord <- rev(ord)
    idx <- integer(sum(BiocGenerics::width(gr)))
    k <- 0L
    for (j in ord) {
      s <- BiocGenerics::start(gr)[j]; e <- BiocGenerics::end(gr)[j]
      seq_v <- if (strand_v[i] == "+") seq.int(s, e) else seq.int(e, s)
      idx[k + seq_len(length(seq_v))] <- seq_v
      k <- k + length(seq_v)
    }
    pos_list[[i]] <- idx
    gr_parts[[i]] <- gr
  }
  names(gr_parts) <- names(cds)
  cds_char <- as.character(cds)
  names(cds_char) <- names(cds)
  list(pos = pos_list, strand = stats::setNames(strand_v, names(cds)),
       chrom = stats::setNames(chrom_v, names(cds)),
       ranges = GenomicRanges::GRangesList(gr_parts),
       seq = cds_char)
}

# chrom lengths from the FASTA index (creates the .fai if absent, like extract_cds does)
.mk_fasta_lengths <- function(genome_fasta) {
  if (!file.exists(paste0(genome_fasta, ".fai"))) Rsamtools::indexFa(genome_fasta)
  fa <- Rsamtools::FaFile(genome_fasta)
  idx <- Rsamtools::scanFaIndex(fa)
  stats::setNames(as.integer(GenomeInfoDb::seqlengths(idx)), names(GenomeInfoDb::seqlengths(idx)))
}

# Validate the snp_file / snp_side pair from a comparison row.
# Returns list(file, side) or NULL when the comparison has no SNP columns. Stops on
# half-specified or invalid input.
.mk_validate_snp_args <- function(snp_file, snp_side, comparison) {
  val <- function(x) {
    if (is.null(x)) return(NA_character_)
    x <- trimws(as.character(x)[1])
    if (!nzchar(x) || toupper(x) %in% c("NA", "NA_CHARACTER_", "NULL", "NONE")) NA_character_ else x
  }
  f <- val(snp_file); s <- val(snp_side)
  if (is.na(f) && is.na(s)) return(NULL)
  if (is.na(f) && !is.na(s)) {
    stop("snp_side is set but snp_file is empty for comparison '", comparison, "'.", call. = FALSE)
  }
  if (is.na(s)) {
    stop("snp_file is set for comparison '", comparison, "' but snp_side is missing. ",
         "Use snp_side = query or subject to say which species' coordinates the SNP file uses.",
         call. = FALSE)
  }
  s <- tolower(s)
  if (!s %in% c("query", "subject")) {
    stop("snp_side must be 'query' or 'subject', not '", s, "' (comparison '", comparison, "').", call. = FALSE)
  }
  list(file = f, side = s)
}

# Map SNPs onto the CDS of one genome. Returns, per SNP hit: gene, snp_i, cds_idx.
# Also carries a per-gene set of CDS indices carrying a polymorphism (for divergence
# exclusion) and the audit counters for the REF check.
.mk_map_snps <- function(snps, map) {
  empty <- data.frame(gene = character(), snp_i = integer(), cds_idx = integer(),
                      stringsAsFactors = FALSE)
  audit <- c(n_checked = 0, n_mismatch = 0)
  if (!nrow(snps)) {
    return(list(hits = empty, excl = list(), audit = audit))
  }
  map_chroms <- unique(map$chrom)
  keep_snp <- snps$chrom %in% map_chroms
  snps <- snps[keep_snp, , drop = FALSE]
  if (!nrow(snps)) {
    return(list(hits = empty, excl = list(), audit = audit))
  }
  snp_gr <- GenomicRanges::GRanges(snps$chrom, IRanges::IRanges(snps$pos, width = 1L))
  hits <- GenomicRanges::findOverlaps(snp_gr, map$ranges, ignore.strand = TRUE)
  if (!length(hits)) {
    return(list(hits = empty, excl = list(), audit = audit))
  }
  q <- S4Vectors::queryHits(hits)
  s <- S4Vectors::subjectHits(hits)
  gene_ids <- names(map$ranges)
  hit_list <- lapply(split(seq_along(q), s), function(rows) {
    g <- s[rows[1]]
    gname <- gene_ids[g]
    p <- map$pos[[gname]]
    ord <- order(p)
    ps <- p[ord]
    qs <- snps$pos[q[rows]]
    k <- findInterval(qs, ps)
    ok <- k >= 1L & ps[pmax(k, 1L)] == qs
    data.frame(gene = rep(gname, sum(ok)),
               snp_i = q[rows[ok]],
               cds_idx = ord[k[ok]],
               stringsAsFactors = FALSE)
  })
  hits_df <- do.call(rbind, hit_list)
  rownames(hits_df) <- NULL

  # REF audit against the CDS sequence (complemented on minus strand)
  if (!is.na(snps$ref[1])) {
    # REF allele must match the reference genome. Check it against the CDS base
    # (complemented on minus-strand genes), which we already have in memory.
    mismatch <- 0L; checked <- 0L
    for (i in seq_len(nrow(hits_df))) {
      g <- hits_df$gene[i]; idx <- hits_df$cds_idx[i]
      base <- substr(map$seq[[g]], idx, idx)
      if (!grepl("^[ACGT]$", base)) next
      if (identical(map$strand[[g]], "-")) base <- chartr("ACGT", "TGCA", base)
      refb <- snps$ref[hits_df$snp_i[i]]
      checked <- checked + 1L
      if (!identical(base, refb)) mismatch <- mismatch + 1L
    }
    audit <- c(n_checked = checked, n_mismatch = mismatch)
  }
  excl <- tapply(hits_df$cds_idx, hits_df$gene, function(x) sort(unique(as.integer(x))), simplify = FALSE)
  list(hits = hits_df, excl = excl, audit = audit)
}

# Per-gene polymorphism counts from SNP hits, allele-aware: a site is Ps when its
# ALT allele leaves the amino acid unchanged and Pn when it changes it, the way the
# original McDonald-Kreitman test classifies a polymorphism. The ALT allele is on the
# genome's strand, so it is complemented for a minus-strand gene before the codon is
# rebuilt. Sites whose reference codon is not pure ACGT or encodes a stop are dropped,
# and so is an ALT equal to the reference base (not a polymorphism).
.mk_poly_counts <- function(snps, map, hits) {
  genes <- unique(hits$gene)
  per <- lapply(genes, function(g) c(pn = 0L, ps = 0L, dropped = 0L))
  names(per) <- genes
  hl <- split(hits, hits$gene)
  for (g in names(hl)) {
    hh <- hl[[g]]
    seqg <- map$seq[[g]]
    minus <- identical(map$strand[[g]], "-")
    for (r in seq_len(nrow(hh))) {
      idx <- hh$cds_idx[r]
      ci <- (idx - 1L) %/% 3L + 1L
      ph <- (idx - 1L) %% 3L
      codon <- substr(seqg, 3L * (ci - 1L) + 1L, 3L * (ci - 1L) + 3L)
      alt <- snps$alt[hh$snp_i[r]]
      if (minus) alt <- chartr("ACGT", "TGCA", alt)
      if (!grepl("^[ACGT]{3}$", codon) || .mk_aa(codon) == "*" || is.na(alt) ||
          alt == substr(codon, ph + 1L, ph + 1L)) {
        per[[g]]["dropped"] <- per[[g]]["dropped"] + 1L
        next
      }
      alt_codon <- paste0(substr(codon, 1L, ph), alt,
                          substr(codon, ph + 2L, 3L))
      if (.mk_aa(alt_codon) == .mk_aa(codon)) {
        per[[g]]["ps"] <- per[[g]]["ps"] + 1L
      } else {
        per[[g]]["pn"] <- per[[g]]["pn"] + 1L
      }
    }
  }
  per
}

# Walk one aligned codon pair; counts substitution classes for the reference side.
# excl: sorted integer vector of reference CDS indices carrying a polymorphism.
.mk_pair_walk <- function(ref_seq, other_seq, excl) {
  n <- nchar(ref_seq)
  dn <- 0; ds <- 0; n_diff <- 0L; n_excl <- 0L; n_stop <- 0L
  n_gap <- 0L; n_invalid <- 0L; n_used <- 0L
  frame_ok <- TRUE
  ci <- 0L # CDS index of the next reference base
  rexcl <- if (length(excl)) excl else integer(0)
  in_excl <- function(i) {
    if (!length(rexcl)) return(FALSE)
    k <- findInterval(i, rexcl)
    k >= 1L && rexcl[k] == i
  }
  if (n >= 3L) for (col in seq.int(1L, n - 2L, by = 3L)) {
    a <- substr(ref_seq, col, col + 2L)
    b <- substr(other_seq, col, col + 2L)
    a_nongap <- sum(strsplit(a, "", fixed = TRUE)[[1]] != "-")
    if (a_nongap > 0L && a_nongap < 3L) {
      # A partial codon on the reference side is a frameshift; downstream codon
      # boundaries no longer line up with CDS codons, so the pair cannot be counted.
      frame_ok <- FALSE
      break
    }
    if (a_nongap == 0L) { n_gap <- n_gap + 1L; next }
    if (grepl("[^ACGT]", a) || grepl("[^ACGT]", b)) {
      ci <- ci + 3L
      if (grepl("[^ACGT]", b)) n_gap <- n_gap + 1L else n_invalid <- n_invalid + 1L
      next
    }
    aa_a <- .mk_aa(a); aa_b <- .mk_aa(b)
    if (aa_a == "*" || aa_b == "*") { ci <- ci + 3L; n_stop <- n_stop + 1L; next }
    idx <- ci + 1:3
    if (any(vapply(idx, in_excl, logical(1)))) { ci <- ci + 3L; n_excl <- n_excl + 1L; next }
    ci <- ci + 3L
    n_used <- n_used + 1L
    if (a == b) next
    res <- .mk_classify_codon(a, b)
    ds <- ds + res[1]; dn <- dn + res[2]; n_diff <- n_diff + 1L
  }
  c(dn = dn, ds = ds, n_diff = n_diff, n_used = n_used, n_excl = n_excl,
    n_stop = n_stop, n_gap = n_gap, n_invalid = n_invalid, frame_ok = as.numeric(frame_ok))
}

# 2x2 Fisher test on rounded counts; NA when either margin is zero.
.mk_fisher <- function(pn, ps, dn, ds) {
  if (!is.finite(dn) || !is.finite(ds)) return(NA_real_)
  dn_i <- round(dn); ds_i <- round(ds)
  if (pn + ps == 0L || dn_i + ds_i == 0L) return(NA_real_)
  tryCatch(
    stats::fisher.test(matrix(c(pn, ps, dn_i, ds_i), nrow = 2L, byrow = TRUE))$p.value,
    error = function(e) NA_real_
  )
}

# alpha, neutrality index and direction of selection from a count tuple.
.mk_stats <- function(pn, ps, dn, ds) {
  if (!is.finite(dn) || !is.finite(ds)) {
    return(c(alpha = NA_real_, neutrality_index = NA_real_, dos = NA_real_, fisher_p = NA_real_))
  }
  alpha <- if (is.finite(dn) && dn > 0 && ps > 0) 1 - (ds * pn) / (dn * ps) else NA_real_
  ni    <- if (is.finite(dn) && dn > 0 && ps > 0) (pn * ds) / (ps * dn) else NA_real_
  dos   <- if ((pn + ps) > 0 && (dn + ds) > 0) pn / (pn + ps) - dn / (dn + ds) else NA_real_
  c(alpha = alpha, neutrality_index = ni, dos = dos, fisher_p = .mk_fisher(pn, ps, dn, ds))
}

.mk_bh <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = "BH")
  out
}

.mk_write_frame <- function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

# Run the whole MK step for one comparison. Called by calculate_dnds().
# Returns the per-gene data.frame invisibly.
.mk_run <- function(comp_dir, comparison_name, dnds_res, aln_files,
                    snp_file, snp_side,
                    query_fasta, subject_fasta, query_gff, subject_gff,
                    query_cds_fasta, subject_cds_fasta) {
  side <- snp_side
  other_side <- if (side == "query") "subject" else "query"
  gff   <- if (side == "query") query_gff   else subject_gff
  fasta <- if (side == "query") query_fasta else subject_fasta
  other_fasta <- if (side == "query") subject_fasta else query_fasta
  cds_fa <- if (side == "query") query_cds_fasta else subject_cds_fasta

  snp <- .mk_read_snps(snp_file)
  message(sprintf("[%s] MK: %d SNP site(s), %d dropped at read time (%s, %s).",
                  comparison_name, snp$meta$n_kept, snp$meta$n_dropped,
                  snp$meta$format, snp_file))

  map <- .mk_cds_map(gff, fasta, cds_fa)
  if (!nrow(snp$snps)) {
    snps <- snp$snps
    mapped <- list(hits = data.frame(gene = character(), snp_i = integer(), cds_idx = integer()),
                   excl = list(), audit = c(n_checked = 0, n_mismatch = 0))
    message("  - MK: no usable SNP sites; Pn and Ps will be zero.")
    return(.mk_finish(comparison_name, comp_dir, dnds_res, aln_files, side, snps, map, mapped))
  }

  fa_len <- .mk_fasta_lengths(fasta)
  snp_chroms <- unique(snp$snps$chrom)
  shared <- intersect(snp_chroms, names(fa_len))
  if (!length(shared)) {
    other_len <- tryCatch(.mk_fasta_lengths(other_fasta), error = function(e) NULL)
    hint <- if (length(other_len) && length(intersect(snp_chroms, names(other_len)))) {
      paste0(" The contigs match the ", other_side, " genome instead - check snp_side.")
    } else ""
    stop("No contig in the SNP file exists in the ", side, " genome (", fasta, ").", hint, call. = FALSE)
  }
  dropped_contigs <- setdiff(snp_chroms, names(fa_len))
  if (length(dropped_contigs)) {
    message(sprintf("  - MK: %d contig(s) in the SNP file are absent from the %s genome; their SNPs were ignored (e.g. %s).",
                    length(dropped_contigs), side, paste(utils::head(dropped_contigs, 3), collapse = ", ")))
    snps <- snp$snps[snp$snps$chrom %in% shared, , drop = FALSE]
  } else {
    snps <- snp$snps
  }
  cl <- snp$meta$contig_len
  if (!is.null(cl)) {
    bad <- shared[vapply(shared, function(cn) {
      L <- cl[[cn]]; !is.null(L) && !is.na(L) && !identical(as.numeric(L), as.numeric(fa_len[[cn]]))
    }, logical(1))]
    if (length(bad)) {
      stop("Contig length mismatch between the SNP file and the ", side, " genome (",
           paste(utils::head(bad, 5), collapse = ", "),
           "). The SNP file is probably from a different assembly.", call. = FALSE)
    }
  }

  mapped <- .mk_map_snps(snps, map)

  # Audit REF alleles against the CDS (catches wrong species / off-by-one coordinates).
  if (snp$meta$format == "vcf" && mapped$audit["n_checked"] > 0) {
    frac <- mapped$audit["n_mismatch"] / mapped$audit["n_checked"]
    if (mapped$audit["n_checked"] >= 10L && frac > 0.10) {
      stop(sprintf(paste0("MK: %d of %d REF alleles (%.0f%%) do not match the %s genome. ",
                          "The SNP file is probably not from this genome, or uses different coordinates."),
                   mapped$audit["n_mismatch"], mapped$audit["n_checked"], 100 * frac, side),
           call. = FALSE)
    }
    if (frac > 0.01) {
      message(sprintf("  - MK: WARNING: %d of %d REF alleles (%.1f%%) do not match the %s genome.",
                      mapped$audit["n_mismatch"], mapped$audit["n_checked"], 100 * frac, side))
    }
  }

  .mk_finish(comparison_name, comp_dir, dnds_res, aln_files, side, snps, map, mapped)
}

# Shared tail of .mk_run: divergence counts, statistics, and output.
.mk_finish <- function(comparison_name, comp_dir, dnds_res, aln_files, side, snps, map, mapped) {
  poly <- .mk_poly_counts(snps, map, mapped$hits)

  # ---- divergence, per aligned pair ----
  pairs <- Filter(Negate(is.null), lapply(aln_files, .read_pair_fasta))
  pair_key <- vapply(pairs, function(p) paste(names(p), collapse = "\r"), character(1))
  dd <- as.data.frame(dnds_res, stringsAsFactors = FALSE)
  dd$query_id <- as.character(dd$query_id)
  dd$subject_id <- as.character(dd$subject_id)
  dup_q <- duplicated(dd$query_id) | duplicated(dd$query_id, fromLast = TRUE)
  dup_s <- duplicated(dd$subject_id) | duplicated(dd$subject_id, fromLast = TRUE)
  if (any(dup_q | dup_s)) {
    message(sprintf("  - MK: %d pair(s) dropped because a gene is not 1-to-1 (appears in multiple pairs).",
                    sum(dup_q | dup_s)))
  }
  dd <- dd[!(dup_q | dup_s), , drop = FALSE]

  n_no_alignment <- 0L
  n_missing_map <- 0L
  n_frame <- 0L
  rows <- vector("list", nrow(dd))
  for (i in seq_len(nrow(dd))) {
    qid <- dd$query_id[i]; sid <- dd$subject_id[i]
    j <- match(paste(qid, sid, sep = "\r"), pair_key)
    if (is.na(j)) j <- match(paste(sid, qid, sep = "\r"), pair_key)
    if (is.na(j)) { n_no_alignment <- n_no_alignment + 1L; next }
    p <- pairs[[j]]
    if (side == "query") {
      ref_seq <- p[1]; oth_seq <- p[2]; ref_gene <- qid
    } else {
      ref_seq <- p[2]; oth_seq <- p[1]; ref_gene <- sid
    }
    if (!ref_gene %in% names(map$pos)) {
      n_missing_map <- n_missing_map + 1L
      next
    }
    wk <- .mk_pair_walk(ref_seq, oth_seq, mapped$excl[[ref_gene]])
    if (!isTRUE(wk[["frame_ok"]] == 1)) {
      n_frame <- n_frame + 1L
      next
    }
    pc <- poly[[ref_gene]]
    pn <- if (is.null(pc)) 0L else as.integer(pc[["pn"]])
    ps <- if (is.null(pc)) 0L else as.integer(pc[["ps"]])
    rows[[i]] <- c(list(query_id = qid, subject_id = sid, reference_side = side,
                        pn = pn, ps = ps,
                        dn = unname(wk["dn"]), ds = unname(wk["ds"]),
                        codons_used = unname(wk["n_used"]),
                        codons_excluded_polymorphic = unname(wk["n_excl"])),
                   as.list(.mk_stats(pn, ps, unname(wk["dn"]), unname(wk["ds"]))))
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (n_no_alignment || n_missing_map || n_frame) {
    message(sprintf(paste0("  - MK: %d pair(s) skipped: %d without a codon alignment, ",
                           "%d whose %s gene is not in the CDS map, %d with a frame-shifted alignment."),
                    n_no_alignment + n_missing_map + n_frame, n_no_alignment, n_missing_map, side, n_frame))
  }
  out <- if (length(rows)) do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE)) else
    data.frame(query_id = character(), subject_id = character(), reference_side = character(),
               pn = integer(), ps = integer(), dn = numeric(), ds = numeric(),
               codons_used = integer(), codons_excluded_polymorphic = integer(),
               alpha = numeric(), neutrality_index = numeric(), dos = numeric(),
               fisher_p = numeric(), stringsAsFactors = FALSE)
  if (nrow(out)) {
    out$dn_int <- round(out$dn); out$ds_int <- round(out$ds)
    out$fisher_q <- .mk_bh(out$fisher_p)
    out <- out[, c("query_id", "subject_id", "reference_side", "pn", "ps", "dn", "ds",
                   "dn_int", "ds_int",
                   "codons_used", "codons_excluded_polymorphic",
                   "alpha", "neutrality_index", "dos", "fisher_p", "fisher_q")]
  } else {
    out$dn_int <- integer(0); out$ds_int <- integer(0); out$fisher_q <- numeric(0)
    out <- out[, c("query_id", "subject_id", "reference_side", "pn", "ps", "dn", "ds",
                   "dn_int", "ds_int",
                   "codons_used", "codons_excluded_polymorphic",
                   "alpha", "neutrality_index", "dos", "fisher_p", "fisher_q")]
  }

  sum_pn <- sum(out$pn); sum_ps <- sum(out$ps); sum_dn <- sum(out$dn); sum_ds <- sum(out$ds)
  pooled_stats <- .mk_stats(sum_pn, sum_ps, sum_dn, sum_ds)
  pooled <- data.frame(
    comparison = comparison_name, reference_side = side, n_pairs = nrow(out),
    pn = sum_pn, ps = sum_ps, dn = sum_dn, ds = sum_ds,
    dn_int = round(sum_dn), ds_int = round(sum_ds),
    codons_used = sum(out$codons_used), codons_excluded_polymorphic = sum(out$codons_excluded_polymorphic),
    alpha = unname(pooled_stats["alpha"]), neutrality_index = unname(pooled_stats["neutrality_index"]),
    dos = unname(pooled_stats["dos"]), fisher_p = unname(pooled_stats["fisher_p"]),
    row.names = NULL, check.names = FALSE)

  f1 <- .mk_write_frame(out, file.path(comp_dir, paste0(comparison_name, "_mk.tsv")))
  f2 <- .mk_write_frame(pooled, file.path(comp_dir, paste0(comparison_name, "_mk_pooled.tsv")))
  message(sprintf(paste0("[%s] MK: %d pairs; Pn=%d Ps=%d Dn=%.2f Ds=%.2f (rounded %d/%d); ",
                         "alpha=%s, NI=%s; wrote %s"),
                  comparison_name, nrow(out), sum_pn, sum_ps, sum_dn, sum_ds, pooled$dn_int, pooled$ds_int,
                  if (is.na(pooled$alpha)) "NA" else sprintf("%.4f", pooled$alpha),
                  if (is.na(pooled$neutrality_index)) "NA" else sprintf("%.4f", pooled$neutrality_index),
                  basename(f1)))
  invisible(out)
}
