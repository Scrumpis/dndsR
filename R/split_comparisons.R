#' Split genomes by subgenome / haplotype / custom and emit a new comparison_file
#'
#' @param comparison_file Path to whitespace-delimited file (tabs/spaces; header or not)
#'   with columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#' @param mode One of "subgenome", "haplotype", or "custom".
#' @param custom_regex Used only when mode="custom". A regex with EXACTLY ONE capture group
#'   that extracts the label from sequence names / GFF seqids (e.g. "_hap(\\d+)$", "(A|B|C)$").
#' @param case_insensitive Logical; match labels case-insensitively (default TRUE).
#' @return Invisibly, the path to the new comparison file (with "_split" inserted before the extension).
#' @export
split_comparisons <- function(comparison_file,
                                       mode = c("subgenome","haplotype","custom"),
                                       custom_regex = NULL,
                                       case_insensitive = TRUE) {
  mode <- match.arg(mode)

  # --- helpers ---------------------------------------------------------------

  .read_ws <- function(path, header_try = TRUE) {
    utils::read.table(path, header = header_try, sep = "", quote = "\"",
                      stringsAsFactors = FALSE, comment.char = "",
                      strip.white = TRUE, blank.lines.skip = TRUE, check.names = FALSE)
  }
  .read_comparisons <- function(x) {
    req <- c("comparison_name","query_fasta","query_gff","subject_fasta","subject_gff")
    df1 <- try(.read_ws(x, header_try = TRUE), silent = TRUE)
    if (!inherits(df1, "try-error") && all(req %in% names(df1))) return(df1[, req, drop = FALSE])
    df2 <- .read_ws(x, header_try = FALSE)
    stopifnot(ncol(df2) >= 5)
    names(df2)[1:5] <- req
    df2[, req, drop = FALSE]
  }
  .split_ext <- function(path) {
    # returns list(base_noext, ext_full) preserving composite extensions like .fa.gz or .gff3.gz
    bn <- basename(path)
    m <- regexec("^(.*)\\.(fa|fasta|fna|gff3?|gff)(\\.gz)?$", bn, ignore.case = TRUE)
    mm <- regmatches(bn, m)[[1]]
    if (length(mm) >= 3) {
      base <- mm[2]
      ext  <- paste0(".", mm[3], ifelse(length(mm) >= 4 && nzchar(mm[4]), mm[4], ""))
    } else {
      base <- tools::file_path_sans_ext(bn)
      ext  <- paste0(".", tools::file_ext(bn))
      if (ext == ".") ext <- ""
    }
    list(base = base, ext = ext)
  }
  .label_from <- switch(mode,
                        subgenome = function(x) {
                          # trailing letter only if immediately preceded by a digit (e.g., Chr1A, chr10b)
                          # excludes chr_cp, Chr_MT, Chr00
                          m <- regexec("(?<=\\d)([A-Za-z])$", x, perl = TRUE, ignore.case = case_insensitive)
                          z <- regmatches(x, m)
                          vapply(z, function(y) if (length(y) == 2) toupper(y[2]) else NA_character_, character(1))
                        },
                        haplotype = function(x) {
                          m <- regexec("_(\\d+)$", x, ignore.case = case_insensitive)
                          z <- regmatches(x, m)
                          vapply(z, function(y) if (length(y) == 2) y[2] else NA_character_, character(1))
                        },
                        custom = {
                          if (is.null(custom_regex) || !nzchar(custom_regex))
                            stop("For mode='custom', provide custom_regex with exactly one capture group, e.g. \"_hap(\\\\d+)$\".")
                          function(x) {
                            m <- regexec(custom_regex, x, perl = TRUE, ignore.case = case_insensitive)
                            z <- regmatches(x, m)
                            vapply(z, function(y) if (length(y) == 2) y[2] else NA_character_, character(1))
                          }
                        }
  )

  .split_fasta <- function(fa_path) {
    if (!file.exists(fa_path)) stop("FASTA not found: ", fa_path)
    x <- Biostrings::readDNAStringSet(fa_path)
    nms <- names(x)
    labs <- .label_from(nms)
    lab_set <- unique(labs[!is.na(labs) & nzchar(labs)])
    if (!length(lab_set)) return(list(labels = character(0), files = character(0)))
    se <- .split_ext(fa_path)
    out_files <- character(0)
    for (L in lab_set) {
      idx <- which(labs == L)
      out <- file.path(dirname(fa_path), sprintf("%s_%s%s", se$base, L, se$ext))
      Biostrings::writeXStringSet(x[idx], filepath = out, compress = grepl("\\.gz$", se$ext, ignore.case = TRUE))
      out_files <- c(out_files, out)
    }
    list(labels = lab_set, files = out_files)
  }

  .split_gff <- function(gff_path) {
    if (!file.exists(gff_path)) stop("GFF not found: ", gff_path)
    lines <- readLines(gff_path)
    is_comment <- startsWith(lines, "#")
    data_lines <- lines[!is_comment]
    if (!length(data_lines)) {
      return(list(labels = character(0), files = character(0)))
    }
    # seqid is the first tab-delimited field
    seqids <- sub("\t.*$", "", data_lines, perl = TRUE)
    labs <- .label_from(seqids)
    lab_set <- unique(labs[!is.na(labs) & nzchar(labs)])
    if (!length(lab_set)) return(list(labels = character(0), files = character(0)))
    se <- .split_ext(gff_path)
    out_files <- character(0)
    for (L in lab_set) {
      keep <- which(labs == L)
      out_lines <- c(lines[is_comment], data_lines[keep])  # keep all comments + only data lines with this label
      out <- file.path(dirname(gff_path), sprintf("%s_%s%s", se$base, L, se$ext))
      writeLines(out_lines, out, useBytes = TRUE)
      out_files <- c(out_files, out)
    }
    list(labels = lab_set, files = out_files)
  }

  .dq <- function(x) paste0("\"", gsub("\"", "\\\\\"", x), "\"")  # quote only paths

  # --- main ------------------------------------------------------------------

  df <- .read_comparisons(comparison_file)
  out_rows <- list()

  for (i in seq_len(nrow(df))) {
    comp  <- df$comparison_name[i]
    qfa   <- df$query_fasta[i]
    qgff  <- df$query_gff[i]
    sfa   <- df$subject_fasta[i]
    sgff  <- df$subject_gff[i]

    message("Splitting: ", comp)

    qfa_res  <- .split_fasta(qfa)
    qgff_res <- .split_gff(qgff)
    sfa_res  <- .split_fasta(sfa)
    sgff_res <- .split_gff(sgff)

    q_labels <- intersect(qfa_res$labels, qgff_res$labels)
    s_labels <- intersect(sfa_res$labels, sgff_res$labels)
    both <- intersect(q_labels, s_labels)
    if (!length(both)) {
      warning("No shared labels for ", comp, " (query labels: ", paste(q_labels, collapse = ","),
              "; subject labels: ", paste(s_labels, collapse = ","), "). Skipping.")
      next
    }

    map_outfile <- function(path, labels) {
      se <- .split_ext(path)
      stats::setNames(file.path(dirname(path), sprintf("%s_%s%s", se$base, labels, se$ext)), labels)
    }
    qfa_map  <- map_outfile(qfa, q_labels)
    qgff_map <- map_outfile(qgff, q_labels)
    sfa_map  <- map_outfile(sfa, s_labels)
    sgff_map <- map_outfile(sgff, s_labels)

    for (L in both) {
      out_rows[[length(out_rows) + 1]] <- data.frame(
        comparison_name = paste0(comp, "_", L),
        query_fasta     = qfa_map[[L]],
        query_gff       = qgff_map[[L]],
        subject_fasta   = sfa_map[[L]],
        subject_gff     = sgff_map[[L]],
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(out_rows)) stop("No split comparisons produced.")

  out_df <- do.call(rbind, out_rows)

  # quote ONLY the path columns; no header
  out_df$query_fasta   <- .dq(out_df$query_fasta)
  out_df$query_gff     <- .dq(out_df$query_gff)
  out_df$subject_fasta <- .dq(out_df$subject_fasta)
  out_df$subject_gff   <- .dq(out_df$subject_gff)

  # write new comparison file next to the original, with "_split" inserted before extension
  comp_dir <- dirname(comparison_file)
  comp_se  <- .split_ext(comparison_file)
  out_path <- file.path(comp_dir, sprintf("%s_split%s", comp_se$base, comp_se$ext))
  utils::write.table(out_df,
                     file = out_path,
                     sep = "\t",
                     quote = FALSE,     # don't auto-quote anything
                     row.names = FALSE,
                     col.names = FALSE) # headerless

  message("Wrote split comparison file: ", out_path)
  invisible(out_path)
}
