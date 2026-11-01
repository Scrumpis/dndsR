# Compares two dndsR versions on real genomes, gene by gene, with gffread as an independent
# source of CDS. Writes summary tables to <work>/summary and gene tables to <work>/genes.
#   Rscript compare-versions.R <work dir> <gffread binary>
# <work> has old/, cur/, inputs/, and comparison.txt, as made by run-version.sb.
# Results stay outside Git. Only aggregate tables should be copied into the repository.
suppressPackageStartupMessages({
  library(Biostrings); library(GenomicFeatures); library(txdbmaker); library(GenomicRanges)
})
args <- commandArgs(trailingOnly = TRUE)
work <- args[1]; gffread <- args[2]
dir.create(file.path(work, "summary"), showWarnings = FALSE)
dir.create(file.path(work, "genes"), showWarnings = FALSE)
comparisons <- read.table(file.path(work, "comparison.txt"), stringsAsFactors = FALSE,
                          col.names = c("name", "qfa", "qgff", "sfa", "sgff"))

stem <- function(path) tools::file_path_sans_ext(basename(path))

# Structure of every gene, and gffread's CDS for every transcript.
genome_table <- function(fa, gff) {
  txdb <- suppressMessages(suppressWarnings(makeTxDbFromGFF(file.path(work, gff), format = "gff3")))
  cds <- cdsBy(txdb, by = "tx")
  map <- suppressMessages(AnnotationDbi::select(txdb, keys = names(cds), keytype = "TXID", columns = c("TXNAME", "GENEID")))
  map <- map[match(names(cds), map$TXID), ]
  named <- dndsR:::.label_transcripts(map$TXID, map$TXNAME, map$GENEID)
  tx <- data.frame(
    tx = named$tx, gene = named$gene,
    strand = as.character(unlist(unique(strand(cds)))[seq_along(cds)]),
    n_cds_parts = lengths(cds),
    stringsAsFactors = FALSE
  )
  out <- tempfile(fileext = ".fa")
  # gffread needs a writable, indexed FASTA, so it gets a private link.
  link <- tempfile(fileext = ".fa")
  file.symlink(normalizePath(file.path(work, fa)), link)
  # gffread stops on a strand of "?", used for some trans-spliced organelle genes.
  clean_gff <- tempfile(fileext = ".gff")
  gff_lines <- readLines(file.path(work, gff))
  writeLines(gff_lines[!grepl("^[^#\\t]*\\t[^\\t]*\\t[^\\t]*\\t[^\\t]*\\t[^\\t]*\\t[^\\t]*\\t\\?\\t", gff_lines)], clean_gff)
  status <- system2(gffread, c("-x", out, "-g", link, clean_gff), stdout = FALSE, stderr = FALSE)
  stopifnot(status == 0L)
  g <- as.character(readDNAStringSet(out))
  names(g) <- sub("^rna-", "", sub(" .*", "", names(g)))
  list(tx = tx, gffread = g)
}

classify <- function(n_isoforms, strand, parts) {
  ifelse(n_isoforms > 1L, "multi_isoform",
         ifelse(parts == 1L, "single_exon",
                ifelse(strand == "-", "multi_exon_minus", "multi_exon_plus")))
}

summaries <- list(); structure_rows <- list(); changes <- list(); problems <- list(); shared <- list()

for (i in seq_len(nrow(comparisons))) {
  cmp <- comparisons[i, ]
  message("== ", cmp$name)
  for (side in c("q", "s")) {
    fa <- cmp[[paste0(side, "fa")]]; gff <- cmp[[paste0(side, "gff")]]
    g <- genome_table(fa, gff)
    old <- as.character(readDNAStringSet(file.path(work, "old", cmp$name, paste0(stem(fa), "_CDS.fasta"))))
    cur <- as.character(readDNAStringSet(file.path(work, "cur", cmp$name, paste0(stem(fa), "_CDS.fasta"))))
    qc <- read.delim(file.path(work, "cur", cmp$name, paste0(stem(fa), "_CDS_qc.tsv")))
    qc$strand <- g$tx$strand[match(qc$transcript, g$tx$tx)]
    qc$parts <- g$tx$n_cds_parts[match(qc$transcript, g$tx$tx)]
    qc$class <- classify(qc$n_isoforms, qc$strand, qc$parts)
    qc$gffread_cds <- toupper(g$gffread[qc$transcript])
    qc$cur_cds <- cur[qc$id]
    qc$old_cds <- old[qc$id]
    # gffread is correct for any isoform choice, so compare with the transcript dndsR chose.
    qc$cur_equals_gffread <- ifelse(qc$written, toupper(qc$cur_cds) == qc$gffread_cds, NA)
    qc$old_equals_gffread <- toupper(qc$old_cds) == qc$gffread_cds
    qc$old_equals_cur <- ifelse(qc$written, toupper(qc$old_cds) == toupper(qc$cur_cds), NA)
    # Did old assemble the same bases in another order? Same length and base composition is the test.
    keep <- c("id", "transcript", "n_isoforms", "length_nt", "multiple_of_3", "internal_stop", "written", "strand", "parts", "class",
              "cur_equals_gffread", "old_equals_gffread", "old_equals_cur")
    utils::write.table(qc[, keep], file.path(work, "genes", sprintf("%s_%s_cds.tsv", cmp$name, side)), sep = "\t", quote = FALSE, row.names = FALSE)
    by_class <- split(qc, qc$class)
    structure_rows[[length(structure_rows) + 1L]] <- do.call(rbind, lapply(names(by_class), function(k) {
      d <- by_class[[k]]
      data.frame(comparison = cmp$name, side = side, class = k, genes = nrow(d),
                 cur_written = sum(d$written),
                 cur_matches_gffread = sum(d$cur_equals_gffread, na.rm = TRUE),
                 cur_differs_from_gffread = sum(!d$cur_equals_gffread, na.rm = TRUE),
                 cur_not_compared = sum(is.na(d$cur_equals_gffread) & d$written),
                 old_matches_gffread = sum(d$old_equals_gffread, na.rm = TRUE),
                 old_differs_from_gffread = sum(!d$old_equals_gffread, na.rm = TRUE),
                 cur_internal_stop = sum(d$internal_stop %in% TRUE),
                 cur_length_not_3n = sum(!d$multiple_of_3),
                 old_with_internal_stop_or_changed = sum(!d$old_equals_cur, na.rm = TRUE))
    }))
  }

  # dN/dS, matched by query gene.
  read_dnds <- function(tag) {
    read.delim(file.path(work, tag, cmp$name, paste0(cmp$name, "_dnds.tsv")))[, c("query_id", "subject_id", "dN", "dS", "dNdS")]
  }
  o <- read_dnds("old"); n <- read_dnds("cur")
  m <- merge(o, n, by = "query_id", all = TRUE, suffixes = c("_old", "_cur"))
  qcds <- read.delim(file.path(work, "genes", sprintf("%s_q_cds.tsv", cmp$name)))
  scds <- read.delim(file.path(work, "genes", sprintf("%s_s_cds.tsv", cmp$name)))
  m$q_class <- qcds$class[match(m$query_id, qcds$id)]
  m$s_class <- scds$class[match(m$subject_id_cur, scds$id)]
  # A pair is affected if either version's CDS differed for either gene.
  q_changed <- !(qcds$old_equals_cur[match(m$query_id, qcds$id)] %in% TRUE)
  s_changed <- !(scds$old_equals_cur[match(ifelse(is.na(m$subject_id_cur), m$subject_id_old, m$subject_id_cur), scds$id)] %in% TRUE)
  m$cds_changed <- q_changed | s_changed
  m$status <- ifelse(is.na(m$subject_id_old), "only_cur",
                     ifelse(is.na(m$subject_id_cur), "only_old",
                            ifelse(m$subject_id_old == m$subject_id_cur, "same_partner", "different_partner")))
  m$delta <- m$dNdS_cur - m$dNdS_old
  m$old_positive <- m$dNdS_old > 1 & m$dNdS_old < 10
  m$cur_positive <- m$dNdS_cur > 1 & m$dNdS_cur < 10
  utils::write.table(m, file.path(work, "genes", sprintf("%s_pairs.tsv", cmp$name)), sep = "\t", quote = FALSE, row.names = FALSE)

  both <- m[m$status == "same_partner" & !is.na(m$dNdS_old) & !is.na(m$dNdS_cur), ]
  describe <- function(d, label) {
    data.frame(
      comparison = cmp$name, set = label, pairs = nrow(d),
      identical_value = sum(d$dNdS_old == d$dNdS_cur),
      median_abs_change = stats::median(abs(d$delta)),
      p95_abs_change = if (nrow(d)) as.numeric(stats::quantile(abs(d$delta), 0.95)) else NA_real_,
      max_abs_change = if (nrow(d)) max(abs(d$delta)) else NA_real_,
      old_above_1 = sum(d$old_positive, na.rm = TRUE),
      cur_above_1 = sum(d$cur_positive, na.rm = TRUE),
      lost_above_1 = sum(d$old_positive & !d$cur_positive, na.rm = TRUE),
      gained_above_1 = sum(!d$old_positive & d$cur_positive, na.rm = TRUE),
      changed_above_0.5x = sum(abs(log2((d$dNdS_cur + 1e-4) / (d$dNdS_old + 1e-4))) > 1)
    )
  }
  changes[[length(changes) + 1L]] <- rbind(
    describe(both, "all_comparable"),
    describe(both[both$cds_changed, ], "cds_changed"),
    describe(both[!both$cds_changed, ], "cds_unchanged"),
    describe(both[both$q_class == "multi_exon_minus" & both$cds_changed, ], "q_multi_exon_minus_changed"),
    describe(both[both$q_class == "multi_isoform" & both$cds_changed, ], "q_multi_isoform_changed"),
    describe(both[both$q_class %in% c("single_exon", "multi_exon_plus") & both$s_class %in% c("single_exon", "multi_exon_plus"), ], "both_genes_unaffected_structure")
  )

  summaries[[length(summaries) + 1L]] <- data.frame(
    comparison = cmp$name, old_pairs = nrow(o), cur_pairs = nrow(n),
    same_partner = sum(m$status == "same_partner"), different_partner = sum(m$status == "different_partner"),
    only_old = sum(m$status == "only_old"), only_cur = sum(m$status == "only_cur")
  )

  value_problems <- function(d, tag) {
    data.frame(
      comparison = cmp$name, version = tag, pairs = nrow(d),
      ratio_na = sum(is.na(d$dNdS)), ds_zero = sum(d$dS == 0, na.rm = TRUE), dn_zero = sum(d$dN == 0, na.rm = TRUE),
      both_zero = sum(d$dN == 0 & d$dS == 0, na.rm = TRUE), dn_na = sum(is.na(d$dN)), ds_na = sum(is.na(d$dS)),
      ratio_infinite = sum(is.infinite(d$dNdS)), ratio_above_1 = sum(d$dNdS > 1, na.rm = TRUE),
      ratio_at_least_10 = sum(d$dNdS >= 10, na.rm = TRUE), ds_above_2 = sum(d$dS > 2, na.rm = TRUE),
      ds_below_0.01 = sum(d$dS < 0.01 & d$dS > 0, na.rm = TRUE),
      median_ds = stats::median(d$dS, na.rm = TRUE), median_ratio = stats::median(d$dNdS, na.rm = TRUE)
    )
  }
  problems[[length(problems) + 1L]] <- rbind(value_problems(o, "old"), value_problems(n, "cur"))

  # Duplicate identifiers between the two genomes, and the annotation result.
  ann_old <- read.delim(file.path(work, "old", cmp$name, paste0(cmp$name, "_dnds_annot.tsv")), quote = "", comment.char = "")
  ann_cur <- read.delim(file.path(work, "cur", cmp$name, paste0(cmp$name, "_dnds_annot.tsv")), quote = "", comment.char = "")
  ids <- function(gff) {
    x <- readLines(file.path(work, gff))
    x <- x[!startsWith(x, "#")]
    unique(sub(".*?ID=([^;]+).*", "\\1", x[grepl("ID=", x)], perl = TRUE))
  }
  qi <- ids(cmp$qgff); si <- ids(cmp$sgff)
  shared[[length(shared) + 1L]] <- data.frame(
    comparison = cmp$name, query_ids = length(qi), subject_ids = length(si), shared_ids = length(intersect(qi, si)),
    pairs_with_same_id = sum(ann_cur$query_id == ann_cur$subject_id),
    annot_rows_old = nrow(ann_old), annot_rows_cur = nrow(ann_cur)
  )
}

write_out <- function(x, name) utils::write.table(do.call(rbind, x), file.path(work, "summary", name), sep = "\t", quote = FALSE, row.names = FALSE)
write_out(structure_rows, "cds_by_gene_structure.tsv")
write_out(summaries, "pairing.tsv")
write_out(changes, "dnds_changes.tsv")
write_out(problems, "value_problems.tsv")
write_out(shared, "shared_identifiers.tsv")
message("Done")
