#' Genome-wide dN/dS ideograms (RIdeogram) for query/subject
#'
#' Reads <comp>/<comp>_merged.tsv and generates chromosome ideograms with
#' two heatmaps: inverted negative selection (dNdS < 1) overlaid, and
#' positive selection (dNdS >= 1) as the label heatmap. Writes
#' <comp>_{q|s}_ideogram.{svg,png} into each comparison directory.
#'
#' @param dnds_merged_file Path to a single <comp>_merged.tsv (single mode).
#' @param comparison_file  Whitespace-delimited file with columns:
#'   comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#'   (Header optional.) In batch mode, the merged file is taken as
#'   file.path(output_dir, comparison_name, paste0(comparison_name, "_merged.tsv")).
#' @param output_dir Root directory containing per-comparison folders in batch mode.
#'   Default: current working directory.
#' @param sides Character vector among c("query","subject"). Default both.
#' @param window_size Integer window size for tiling genes (default 300000).
#' @param max_dnds Drop rows with dNdS >= max_dnds or NA (default 10).
#' @param filter_expr Optional character with a logical expression evaluated
#'   in the merged table (e.g. "q_chr == s_chr"). Rows failing are removed.
#' @param make_png Logical; if TRUE (default) also write a PNG via convertSVG().
#' @param overwrite Logical; if FALSE (default) skip when outputs exist.
#' @param verbose Logical; print progress.
#'
#' @return Invisibly, character vector of output paths written.
#' @export
dnds_ideogram <- function(dnds_merged_file = NULL,
                          comparison_file  = NULL,
                          output_dir       = getwd(),
                          sides            = c("query","subject"),
                          window_size      = 300000,
                          max_dnds         = 10,
                          filter_expr      = NULL,
                          make_png         = TRUE,
                          overwrite        = FALSE,
                          verbose          = TRUE) {
  sides <- match.arg(sides, choices = c("query","subject"), several.ok = TRUE)

  req_cli <- function(...) if (verbose) message("[dndsR::dnds_ideogram] ", sprintf(...))
  die <- function(...) stop(sprintf(...), call. = FALSE)

  # ---- utilities ----
  .read_ws <- function(path, header_try = TRUE) {
    utils::read.table(path, header = header_try, sep = "", quote = "\"",
                      stringsAsFactors = FALSE, comment.char = "",
                      strip.white = TRUE, blank.lines.skip = TRUE, check.names = FALSE)
  }

  .read_comparisons <- function(x) {
    req <- c("comparison_name","query_fasta","query_gff","subject_fasta","subject_gff")
    if (is.data.frame(x)) { stopifnot(all(req %in% names(x))); return(x[, req, drop = FALSE]) }
    df1 <- try(.read_ws(x, header_try = TRUE), silent = TRUE)
    if (!inherits(df1, "try-error") && all(req %in% names(df1))) return(df1[, req, drop = FALSE])
    df2 <- .read_ws(x, header_try = FALSE)
    if (ncol(df2) < 5) stop("comparison_file must have 5 columns (or header) with: ", paste(req, collapse = ", "))
    names(df2)[1:5] <- req
    df2[, req, drop = FALSE]
  }

  .ensure_fai <- function(fasta) {
    fai <- paste0(fasta, ".fai")
    if (!file.exists(fai)) {
      code <- suppressWarnings(system2("samtools", c("faidx", shQuote(fasta)), stdout = TRUE, stderr = TRUE))
      if (!file.exists(fai)) die("Failed to create FAI for %s. samtools said:\n%s", fasta, paste(code, collapse = "\n"))
    }
    fai
  }

  .read_karyotype <- function(fai_path) {
    fai <- utils::read.table(fai_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    # Strip leading non-digits (Nick's convention), start at 0, then custom order
    karyotype <- data.frame(
      Chr   = sub("^[^1-9]*([1-9].*)", "\\1", fai[[1]]),
      Start = 0L,
      End   = as.integer(fai[[2]]),
      stringsAsFactors = FALSE
    )
    karyotype[order(substring(karyotype$Chr, 2, 2), substring(karyotype$Chr, 1, 1)), , drop = FALSE]
  }

  .order_chr_tbl <- function(tbl, chr_col = "Chr") {
    tbl[order(substring(tbl[[chr_col]], 2, 2), substring(tbl[[chr_col]], 1, 1), tbl$Start), , drop = FALSE]
  }

  .one_side <- function(merged_path, side, fasta, window_size, max_dnds, filter_expr, make_png, overwrite) {
    if (!requireNamespace("RIdeogram", quietly = TRUE)) die("Please install RIdeogram: install.packages('RIdeogram')")
    if (!requireNamespace("xml2", quietly = TRUE))      die("Please install xml2: install.packages('xml2')")

    side <- match.arg(side, c("query","subject"))
    data_type_sub <- substr(side, 1, 1) # 'q' or 's'
    comp_dir  <- dirname(normalizePath(merged_path))
    comp_base <- sub("_merged\\.tsv$", "", basename(merged_path))
    out_svg   <- file.path(comp_dir, sprintf("%s_%s_ideogram.svg", comp_base, side))
    out_png   <- file.path(comp_dir, sprintf("%s_%s_ideogram.png", comp_base, side))

    if (!overwrite && file.exists(out_svg) && (!make_png || file.exists(out_png))) {
      req_cli("Exists -> skipping: %s (%s)", comp_base, side); return(out_svg)
    }

    req_cli("Loading merged: %s", merged_path)
    d <- utils::read.table(merged_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    need <- c("dNdS", sprintf("%s_chr", data_type_sub), sprintf("%s_gene_start", data_type_sub))
    miss <- setdiff(need, names(d))
    if (length(miss)) die("Merged file missing columns: %s", paste(miss, collapse = ", "))

    # Filter
    keep <- !is.na(d$dNdS) & d$dNdS < max_dnds
    if (!is.null(filter_expr) && nzchar(filter_expr)) {
      ok <- try(eval(parse(text = filter_expr), envir = d, enclos = parent.frame()), silent = TRUE)
      if (inherits(ok, "try-error")) die("Bad filter_expr: %s", filter_expr)
      ok[is.na(ok)] <- FALSE
      keep <- keep & ok
    }
    d <- d[keep, , drop = FALSE]
    if (!nrow(d)) { req_cli("No rows after filtering -> skip %s", comp_base); return(character(0)) }

    # Transform + windowing
    chr_col   <- sprintf("%s_chr", data_type_sub)
    start_col <- sprintf("%s_gene_start", data_type_sub)
    ln        <- log(d$dNdS + 1)
    d$ln_dNdS      <- ln
    d$ln_neg_dNdS  <- ifelse(d$dNdS >= 1, NA_real_, ln)
    d$ln_pos_dNdS  <- ifelse(d$dNdS <  1, NA_real_, ln)

    win <- aggregate(
      list(total_raw_dNdS   = d$dNdS,
           total_ln_dNdS    = d$ln_dNdS,
           total_ln_neg_dNdS= d$ln_neg_dNdS,
           total_ln_pos_dNdS= d$ln_pos_dNdS,
           n_genes          = rep(1L, nrow(d)),
           neg_genes        = as.integer(d$dNdS < 1),
           pos_genes        = as.integer(d$dNdS >= 1)),
      by = list(Chr = d[[chr_col]],
                Start = (d[[start_col]] %/% window_size) * window_size),
      FUN = function(x) if (is.numeric(x)) sum(x, na.rm = TRUE) else sum(x)
    )
    win$End <- win$Start + window_size
    win$avg_ln_dNdS     <- ifelse(win$n_genes == 0, NA_real_, win$total_ln_dNdS     / pmax(win$n_genes, 1))
    win$avg_ln_neg_dNdS <- ifelse(win$neg_genes== 0, NA_real_, win$total_ln_neg_dNdS/ pmax(win$neg_genes,1))
    win$avg_ln_pos_dNdS <- ifelse(win$pos_genes== 0, NA_real_, win$total_ln_pos_dNdS/ pmax(win$pos_genes,1))

    # Karyotype from FASTA index
    fai <- .ensure_fai(fasta)
    karyotype <- .read_karyotype(fai)

    # Fix last window to chrom end
    last_idx <- tapply(seq_len(nrow(win)), win$Chr, tail, n = 1)
    last_idx <- unlist(last_idx, use.names = FALSE)
    if (length(last_idx)) {
      m <- match(win$Chr[last_idx], karyotype$Chr)
      win$End[last_idx] <- karyotype$End[m]
    }

    # Inputs for RIdeogram
    neg <- win[!is.nan(win$avg_ln_neg_dNdS), c("Chr","Start","End","avg_ln_neg_dNdS")]
    pos <- win[!is.nan(win$avg_ln_pos_dNdS), c("Chr","Start","End","avg_ln_pos_dNdS")]

    if (!nrow(neg) && !nrow(pos)) { req_cli("No heatmap data -> skip %s", comp_base); return(character(0)) }

    # Order + invert neg values for yellow->red scale
    if (nrow(neg)) {
      neg <- .order_chr_tbl(neg, "Chr"); names(neg)[4] <- "Value"; neg$Value <- -1 * neg$Value
    }
    if (nrow(pos)) {
      pos <- .order_chr_tbl(pos, "Chr"); names(pos)[4] <- "Value"
    }

    req_cli("Rendering ideogram: %s (%s)", comp_base, side)
    RIdeogram::ideogram(
      karyotype = .order_chr_tbl(karyotype, "Chr"),
      overlaid  = if (nrow(neg)) neg else NULL,
      label     = if (nrow(pos)) pos else NULL,
      label_type = "heatmap",
      colorset1 = c("#FFFFCC", "#e34a33"),  # Neg (inverted) yellow->red
      colorset2 = c("#FFFFCC", "#2c7fb8"),  # Pos yellow->blue
      Ly = 3
    )

    # Patch legend/text in SVG, rename & convert
    if (!file.exists("chromosome.svg")) die("RIdeogram did not write chromosome.svg")
    svg <- xml2::read_xml("chromosome.svg")
    tn  <- xml2::xml_find_all(svg, ".//text")
    xml2::xml_attr(tn, "font-size")  <- "14"
    xml2::xml_attr(tn, "font-weight")<- "bold"

    txt <- xml2::xml_text(tn)
    low  <- which(txt == "Low")
    high <- which(txt == "High")
    if (length(low))  xml2::xml_text(tn[low]) <- "Neut"
    if (length(high) >= 2) {
      xml2::xml_text(tn[high[1]]) <- "Neg"
      xml2::xml_text(tn[high[2]]) <- "Pos"
      parent <- xml2::xml_parent(tn[high[1]])
      hdr <- xml2::xml_add_child(parent, "text", "Selection Pressure")
      xml2::xml_set_attr(hdr, "x", "600")
      xml2::xml_set_attr(hdr, "y", "0")
      xml2::xml_set_attr(hdr, "font-size", "14")
      xml2::xml_set_attr(hdr, "font-weight", "bold")
      xml2::xml_set_attr(hdr, "font-family", "Arial")
      xml2::xml_set_attr(hdr, "text-anchor", "middle")
    }

    xml2::write_xml(svg, "chromosome.svg")
    file.rename("chromosome.svg", out_svg)

    if (make_png) {
      RIdeogram::convertSVG(out_svg, device = "png")
      if (file.exists("chromosome.png")) file.rename("chromosome.png", out_png)
    }

    if (verbose) message("[dndsR::dnds_ideogram] Wrote: ", out_svg,
                         if (make_png && file.exists(out_png)) paste0(" and ", out_png) else "")
    out_svg
  }

  # ---- batch vs single ----
  outs <- character(0)

  if (!is.null(comparison_file)) {
    df <- .read_comparisons(comparison_file)
    for (i in seq_len(nrow(df))) {
      comp <- df$comparison_name[i]
      comp_dir <- file.path(output_dir, comp)
      merged_path <- file.path(comp_dir, sprintf("%s_dnds.tsv", comp))
      if (!file.exists(merged_path)) { if (verbose) message("Missing dnds file: ", merged_path); next }
      for (sd in sides) {
        fasta <- if (sd == "query") df$query_fasta[i] else df$subject_fasta[i]
        if (!file.exists(fasta)) { if (verbose) message("Missing FASTA (", sd, "): ", fasta); next }
        outs <- c(outs, .one_side(merged_path, sd, fasta, window_size, max_dnds, filter_expr, make_png, overwrite))
      }
    }
    if (verbose) message("All ideograms complete.")
    return(invisible(stats::na.omit(outs)))
  }

  # single mode
  if (is.null(dnds_merged_file)) die("Provide either comparison_file (batch) OR dnds_merged_file (single).")
  if (!file.exists(dnds_merged_file)) die("dnds_merged_file not found: %s", dnds_merged_file)

  comp_dir  <- dirname(dnds_merged_file)
  comp_base <- sub("_merged\\.tsv$", "", basename(dnds_merged_file))
  # side-specific FASTA is required in single mode; choose via 'sides'
  if (length(sides) != 1L) die("Single mode requires exactly one 'side' in 'sides'.")
  stop("In single mode, please call via CLI shim (below) or provide a FASTA; ",
       "the exported function is designed for batch via comparison_file.")
}
