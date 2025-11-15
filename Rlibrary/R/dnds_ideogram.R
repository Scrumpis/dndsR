#' Genome-wide dN/dS ideograms (RIdeogram) for query/subject
#'
#' For each comparison, reads a dN/dS table (first existing among
#' \code{<comp>_dnds.tsv}, \code{<comp>_dnds_annot.tsv}, or \code{<comp>_merged.tsv}),
#' derives genomic coordinates for \code{query_id}/\code{subject_id} from the
#' provided GFFs (first attribute match per ID), writes a temporary
#' \code{<comp>/<comp>_ideogram_input.tsv} containing \code{q_chr},
#' \code{q_gene_start}, \code{s_chr}, \code{s_gene_start}, and then generates
#' chromosome ideograms per side with two heatmaps: inverted negative selection
#' (\code{dNdS < 1}) overlaid, and positive selection (\code{dNdS >= 1}) as the
#' label heatmap. Outputs are \code{<comp>_{query|subject}_ideogram.{svg,png}}.
#'
#' @param dnds_merged_file Path to a single prebuilt intermediate table that
#'   already contains \code{q_chr/q_gene_start} or \code{s_chr/s_gene_start}
#'   (single mode). Most users should prefer batch mode via \code{comparison_file}.
#' @param comparison_file Whitespace-delimited file with columns:
#'   \code{comparison_name}, \code{query_fasta}, \code{query_gff},
#'   \code{subject_fasta}, \code{subject_gff}. (Header optional.) In batch mode,
#'   the function will look for the dN/dS table inside
#'   \code{file.path(output_dir, comparison_name)} using the first existing among
#'   \code{<comp>_dnds.tsv}, \code{<comp>_dnds_annot.tsv}, \code{<comp>_merged.tsv},
#'   build the intermediate by joining coordinates from the GFFs, and plot.
#' @param output_dir Root directory containing per-comparison folders in batch
#'   mode. Default: current working directory.
#' @param sides Character vector among \code{c("query","subject")}. Default both.
#' @param window_size Integer window size (bp) for binning genes (default \code{300000}).
#' @param max_dnds Drop rows with \code{dNdS >= max_dnds} or \code{NA} (default \code{10}).
#' @param filter_expr Optional character with a logical expression evaluated in
#'   the augmented table (e.g. \code{"q_chr == s_chr"}). Rows failing are removed.
#'   (No chromosome mapping is performed; this is a raw filter.)
#' @param make_png Logical; if \code{TRUE} (default) also write a PNG via \code{convertSVG()}.
#' @param overwrite Logical; if \code{FALSE} (default) skip when outputs exist.
#' @param keep_intermediate Logical; keep the generated
#'   \code{<comp>_ideogram_input.tsv} (default \code{FALSE}).
#' @param chr_strip_leading_chr0 Logical; optionally strip leading "chr" and optional "0"
#'   (e.g., \code{chr01}→\code{1}). Case-insensitive by default. Default \code{FALSE}.
#' @param chr_strip_leading Character (regex) to strip from the start (e.g., \code{"Cf0"}).
#'   Default \code{NULL} (no stripping).
#' @param chr_strip_trailing Character (regex) to strip from the end (e.g., \code{"_random"}).
#'   Default \code{NULL} (no stripping).
#' @param chr_case_insensitive Logical; apply the above stripping case-insensitively
#'   (default \code{TRUE}).
#' @param restrict_gff_to_gene Logical; if \code{TRUE}, only use GFF rows with \code{type=="gene"}
#'   when indexing coordinates. Default \code{FALSE}.
#' @param verbose Logical; print progress.
#' @param ... Additional arguments passed to \code{RIdeogram::ideogram()}, such
#'   as layout or color options. Arguments that would override core mappings
#'   (\code{karyotype}, \code{overlaid}, \code{label}) are ignored.
#'
#' @details
#' Coordinates are derived by parsing the GFF attribute column and matching common
#' keys (\code{ID}, \code{gene_id}, \code{Parent}, \code{Name}, \code{gene},
#' \code{locus_tag}, \code{geneID}); the first occurrence of an ID is used,
#' taking seqid (col 1) and start (col 4). Header lines beginning with \code{'#'}
#' are ignored. Chromosome labels are kept \emph{raw} by default; optional
#' normalization flags can be used to strip site-specific prefixes/suffixes.
#'
#' @return Invisibly, a character vector of output file paths written.
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
                          keep_intermediate= FALSE,
                          # label normalization (opt-in; raw by default)
                          chr_strip_leading_chr0 = FALSE,
                          chr_strip_leading      = NULL,
                          chr_strip_trailing     = NULL,
                          chr_case_insensitive   = TRUE,
                          # GFF parsing
                          restrict_gff_to_gene   = FALSE,
                          verbose          = TRUE,
                          ...) {

  sides <- match.arg(sides, choices = c("query","subject"), several.ok = TRUE)
  req_cli <- function(...) if (verbose) message("[dndsR::dnds_ideogram] ", sprintf(...))
  die <- function(...) stop(sprintf(...), call. = FALSE)

  # extra arguments to pass on to RIdeogram::ideogram()
  ideogram_dots <- list(...)

  # ---------- utilities ----------
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
  # Keep RAW names (no digit stripping). RIdeogram accepts arbitrary labels.
  .read_karyotype <- function(fai_path) {
    fai <- utils::read.table(fai_path, sep = "\t", header = FALSE,
                             stringsAsFactors = FALSE, quote = "", comment.char = "")
    if (!nrow(fai)) stop("FAI appears empty: ", fai_path, call. = FALSE)
    data.frame(
      Chr   = as.character(fai[[1]]),
      Start = 0L,
      End   = as.integer(fai[[2]]),
      stringsAsFactors = FALSE
    )
  }
  .order_by_karyotype <- function(tbl, karyo, chr_col = "Chr") {
    if (!nrow(tbl)) return(tbl)
    lev <- karyo$Chr
    f <- factor(tbl[[chr_col]], levels = lev, ordered = TRUE)
    tbl <- tbl[order(f, tbl$Start), , drop = FALSE]
    # ensure Chr is character, not factor
    tbl[[chr_col]] <- as.character(tbl[[chr_col]])
    tbl
  }
  .normalize_chr <- function(x,
                             strip_chr0 = FALSE,
                             strip_lead = NULL,
                             strip_trail= NULL,
                             ci = TRUE) {
    if (!length(x)) return(x)
    out <- as.character(x)
    if (strip_chr0) {
      rx <- if (ci) "(?i)^chr0?" else "^chr0?"
      out <- sub(rx, "", out, perl = TRUE)
    }
    if (!is.null(strip_lead) && nzchar(strip_lead)) {
      rx <- if (ci) paste0("(?i)^(", strip_lead, ")") else paste0("^(", strip_lead, ")")
      out <- sub(rx, "", out, perl = TRUE)
    }
    if (!is.null(strip_trail) && nzchar(strip_trail)) {
      rx <- if (ci) paste0("(?i)(", strip_trail, ")$") else paste0("(", strip_trail, ")$")
      out <- sub(rx, "", out, perl = TRUE)
    }
    out
  }
  .gff_index <- function(gff_path, restrict_gene = FALSE) {
    if (!file.exists(gff_path)) die("GFF not found: %s", gff_path)
    g <- try(utils::read.table(gff_path, sep = "\t", quote = "", comment.char = "#",
                               header = FALSE, stringsAsFactors = FALSE, fill = TRUE,
                               col.names = paste0("V", 1:9)), silent = TRUE)
    if (inherits(g, "try-error")) die("Could not read GFF: %s", gff_path)
    if (ncol(g) < 9) die("GFF appears malformed (<9 columns): %s", gff_path)
    if (restrict_gene && "V3" %in% names(g) && any(g$V3 == "gene", na.rm = TRUE)) {
      g <- g[g$V3 == "gene", , drop = FALSE]
    }
    attrs <- g$V9
    seqid <- g$V1
    start <- suppressWarnings(as.integer(g$V4))
    split_attrs <- strsplit(attrs, ";", fixed = TRUE)
    grab_val <- function(parts, key) {
      kv <- parts[startsWith(parts, paste0(key, "="))]
      if (!length(kv)) return(NA_character_)
      sub(paste0("^", key, "="), "", kv[1])
    }
    keys <- c("ID","gene_id","Parent","Name","gene","locus_tag","geneID")
    collect_first <- function(parts) {
      for (k in keys) {
        v <- grab_val(parts, k)
        if (!is.na(v) && nzchar(v)) return(v)
      }
      NA_character_
    }
    ids <- vapply(split_attrs, collect_first, character(1))
    keep <- !is.na(ids) & nzchar(ids)
    df <- data.frame(id = ids[keep], chr = seqid[keep], start = start[keep], stringsAsFactors = FALSE)
    df <- df[!duplicated(df$id), , drop = FALSE]
    rownames(df) <- df$id
    df
  }
  .augment_with_coords <- function(d, id_map, id_col, out_prefix) {
    miss <- setdiff(id_col, names(d))
    if (length(miss)) die("Input missing required column(s): %s", paste(miss, collapse = ", "))
    ids <- d[[id_col]]
    m   <- id_map[ids, , drop = FALSE]
    d[[paste0(out_prefix, "_chr")]]        <- m$chr
    d[[paste0(out_prefix, "_gene_start")]] <- m$start
    d
  }
  .find_dnds_table <- function(comp_dir, comp) {
    candidates <- file.path(comp_dir, c(
      sprintf("%s_dnds.tsv", comp),
      sprintf("%s_dnds_annot.tsv", comp),
      sprintf("%s_merged.tsv", comp)
    ))
    ex <- file.exists(candidates)
    if (!any(ex)) return(NA_character_)
    candidates[which(ex)[1]]
  }
  # Crop output SVG to remove bottom white space
  crop_svg_by_bottom_margin <- function(svg, cut_px) {
    root <- xml2::xml_root(svg)

    vb <- xml2::xml_attr(root, "viewBox")
    if (is.na(vb) || !nzchar(vb)) {
      w <- suppressWarnings(as.numeric(sub("px$", "", xml2::xml_attr(root, "width"))))
      h <- suppressWarnings(as.numeric(sub("px$", "", xml2::xml_attr(root, "height"))))
      if (!is.finite(w) || !is.finite(h)) { w <- 1200; h <- 800 }
      xml2::xml_attr(root, "viewBox") <- sprintf("0 0 %g %g", w, h)
      vb <- xml2::xml_attr(root, "viewBox")
    }

    nums <- as.numeric(strsplit(vb, "[ ,]+")[[1]])
    minX <- nums[1]; minY <- nums[2]; vw <- nums[3]; vh <- nums[4]

    new_vh <- max(10, vh - cut_px)

    xml2::xml_attr(root, "viewBox") <- sprintf("%g %g %g %g", minX, minY, vw, new_vh)

    old_h <- suppressWarnings(as.numeric(sub("px$", "", xml2::xml_attr(root, "height"))))
    if (is.finite(old_h) && vh > 0) {
      xml2::xml_attr(root, "height") <- paste0(old_h * (new_vh / vh), "px")
    }

    bg <- xml2::xml_find_first(root, ".//rect[@id='svg-bg-rect']")
    if (!inherits(bg, "xml_missing")) {
      xml2::xml_set_attr(bg, "x", as.character(minX))
      xml2::xml_set_attr(bg, "y", as.character(minY))
      xml2::xml_set_attr(bg, "width",  as.character(vw))
      xml2::xml_set_attr(bg, "height", as.character(new_vh))
    }

    svg
  }

  # ---------- rendering worker (per side) ----------
  .one_side <- function(merged_path, side, fasta, window_size, max_dnds, filter_expr,
                        make_png, overwrite,
                        chr_strip_leading_chr0, chr_strip_leading, chr_strip_trailing, chr_case_insensitive,
                        ideogram_args = NULL) {

    if (!requireNamespace("RIdeogram", quietly = TRUE)) die("Please install RIdeogram: install.packages('RIdeogram')")
    if (!requireNamespace("xml2", quietly = TRUE))      die("Please install xml2: install.packages('xml2')")

    side <- match.arg(side, c("query","subject"))
    short <- substr(side, 1, 1) # 'q' or 's'
    comp_dir  <- dirname(normalizePath(merged_path))
    comp_base <- sub("\\.(tsv|txt)$", "", basename(merged_path))
    out_svg   <- file.path(comp_dir, sprintf("%s_%s_ideogram.svg", comp_base, side))
    out_png   <- file.path(comp_dir, sprintf("%s_%s_ideogram.png", comp_base, side))

    if (!overwrite && file.exists(out_svg) && (!make_png || file.exists(out_png))) {
      req_cli("Exists -> skipping: %s (%s)", comp_base, side); return(out_svg)
    }

    req_cli("Loading: %s", merged_path)
    d <- utils::read.table(merged_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")

    # Expect at minimum dNdS + id cols + our coord cols
    need0 <- c("dNdS","query_id","subject_id")
    miss0 <- setdiff(need0, names(d))
    if (length(miss0)) die("Table missing columns: %s", paste(miss0, collapse = ", "))

    # Early filter on dNdS
    keep <- !is.na(d$dNdS) & d$dNdS < max_dnds
    d <- d[keep, , drop = FALSE]
    if (!nrow(d)) { req_cli("No rows after dNdS filtering -> skip"); return(character(0)) }

    need <- c(sprintf("%s_chr", short), sprintf("%s_gene_start", short))
    if (!all(need %in% names(d))) die("Internal error: expected coord columns %s in intermediate.", paste(need, collapse = ", "))

    # Optional user filter (raw logical; no mapping implied)
    if (!is.null(filter_expr) && nzchar(filter_expr)) {
      ok <- try(eval(parse(text = filter_expr), envir = d, enclos = parent.frame()), silent = TRUE)
      if (inherits(ok, "try-error")) die("Bad filter_expr: %s", filter_expr)
      ok[is.na(ok)] <- FALSE
      d <- d[ok, , drop = FALSE]
    }
    if (!nrow(d)) { req_cli("No rows after filter_expr -> skip"); return(character(0)) }

    # Transform + windowing
    chr_col   <- sprintf("%s_chr", short)
    start_col <- sprintf("%s_gene_start", short)
    ln        <- log(d$dNdS + 1)
    d$ln_dNdS      <- ln
    d$ln_neg_dNdS  <- ifelse(d$dNdS >= 1, NA_real_, ln)
    d$ln_pos_dNdS  <- ifelse(d$dNdS <  1, NA_real_, ln)

    win <- aggregate(
      list(total_raw_dNdS    = d$dNdS,
           total_ln_dNdS     = d$ln_dNdS,
           total_ln_neg_dNdS = d$ln_neg_dNdS,
           total_ln_pos_dNdS = d$ln_pos_dNdS,
           n_genes           = rep(1L, nrow(d)),
           neg_genes         = as.integer(d$dNdS < 1),
           pos_genes         = as.integer(d$dNdS >= 1)),
      by = list(Chr = d[[chr_col]],
                Start = (d[[start_col]] %/% window_size) * window_size),
      FUN = function(x) if (is.numeric(x)) sum(x, na.rm = TRUE) else sum(x)
    )

    # add End and the averaged columns
    win$End <- win$Start + window_size
    win$avg_ln_dNdS      <- ifelse(win$n_genes  == 0, NA_real_, win$total_ln_dNdS      / pmax(win$n_genes,  1))
    win$avg_ln_neg_dNdS  <- ifelse(win$neg_genes== 0, NA_real_, win$total_ln_neg_dNdS  / pmax(win$neg_genes,1))
    win$avg_ln_pos_dNdS  <- ifelse(win$pos_genes== 0, NA_real_, win$total_ln_pos_dNdS  / pmax(win$pos_genes,1))
    if (!nrow(win)) { req_cli("No windows formed -> skip"); return(character(0)) }

    # Karyotype (RAW names), then optional normalization for both karyotype & win
    # fai <- .ensure_fai(fasta)
    # karyotype <- .read_karyotype(fai)
    karyotype <- .fasta_karyotype(fasta)

    # Optional normalization (opt-in; default no change)
    karyotype$Chr <- .normalize_chr(
      karyotype$Chr,
      strip_chr0 = chr_strip_leading_chr0,
      strip_lead = chr_strip_leading,
      strip_trail= chr_strip_trailing,
      ci = chr_case_insensitive
    )
    win$Chr <- .normalize_chr(
      win$Chr,
      strip_chr0 = chr_strip_leading_chr0,
      strip_lead = chr_strip_leading,
      strip_trail= chr_strip_trailing,
      ci = chr_case_insensitive
    )

    # Align to karyotype; clamp windows to chrom ends
    win <- win[win$Chr %in% karyotype$Chr, , drop = FALSE]
    if (!nrow(win)) { req_cli("No windows remain after matching chromosomes to FASTA -> skip"); return(character(0)) }
    chr_len <- setNames(karyotype$End, karyotype$Chr)
    win$End <- pmin(win$End, chr_len[win$Chr])
    win$Start[is.na(win$Start)] <- 0L
    win$End[is.na(win$End)]     <- win$Start[is.na(win$End)]
    win <- win[is.finite(win$Start) & is.finite(win$End) & win$End > win$Start, , drop = FALSE]
    if (!nrow(win)) { req_cli("All windows invalid after clamping -> skip"); return(character(0)) }

    # Build heatmaps
    neg <- subset(win, is.finite(avg_ln_neg_dNdS), select = c("Chr","Start","End","avg_ln_neg_dNdS"))
    pos <- subset(win, is.finite(avg_ln_pos_dNdS), select = c("Chr","Start","End","avg_ln_pos_dNdS"))
    if (!nrow(neg) && !nrow(pos)) { req_cli("No heatmap data -> skip"); return(character(0)) }

    if (nrow(neg)) { neg <- .order_by_karyotype(neg, karyotype, "Chr"); names(neg)[4] <- "Value"; neg$Value <- -1 * neg$Value }
    if (nrow(pos)) { pos <- .order_by_karyotype(pos, karyotype, "Chr"); names(pos)[4] <- "Value" }
    karyotype <- .order_by_karyotype(karyotype, karyotype, "Chr")

    req_cli(sprintf("Rendering ideogram (%s)", side))

    if (is.null(ideogram_args)) ideogram_args <- list()

    # prevent user from stomping on core mapping args
    protected <- c("karyotype","overlaid","label")
    bad <- intersect(names(ideogram_args), protected)
    if (length(bad)) {
      warning("[dndsR::dnds_ideogram] Ignoring arguments in ... that conflict with core RIdeogram params: ",
              paste(bad, collapse = ", "))
      ideogram_args[bad] <- NULL
    }

    base_args <- list(
      karyotype = karyotype,
      overlaid  = if (nrow(neg)) neg else NULL,
      label     = if (nrow(pos)) pos else NULL,
      label_type = "heatmap",
      colorset1 = c("#FFFFCC", "#e34a33"),  # inverted negative: yellow->red
      colorset2 = c("#FFFFCC", "#2c7fb8"),  # positive: yellow->blue
      Ly = 3
    )

    do.call(RIdeogram::ideogram, c(base_args, ideogram_args))

    if (!file.exists("chromosome.svg")) die("RIdeogram did not write chromosome.svg")
    svg <- xml2::read_xml("chromosome.svg")
    tn  <- xml2::xml_find_all(svg, ".//text")
    xml2::xml_attr(tn, "font-size")  <- "14"
    xml2::xml_attr(tn, "font-weight")<- "bold"
    txt <- xml2::xml_text(tn)
    low  <- which(txt == "Low"); high <- which(txt == "High")
    if (length(low))  xml2::xml_text(tn[low]) <- "Neut"
    if (length(high) >= 2) {
      xml2::xml_text(tn[high[1]]) <- "Neg"
      xml2::xml_text(tn[high[2]]) <- "Pos"
      parent <- xml2::xml_parent(tn[high[1]])
      hdr <- xml2::xml_add_child(parent, "text", "Selection Pressure")
      xml2::xml_set_attr(hdr, "x", "600"); xml2::xml_set_attr(hdr, "y", "0")
      xml2::xml_set_attr(hdr, "font-size", "14"); xml2::xml_set_attr(hdr, "font-weight", "bold")
      xml2::xml_set_attr(hdr, "font-family", "Arial"); xml2::xml_set_attr(hdr, "text-anchor", "middle")
    }

    # === Center chromosome labels around their current placement ===
    chr_names <- karyotype$Chr
    lab_nodes <- tn[xml2::xml_text(tn) %in% chr_names]
    if (length(lab_nodes)) {

      # helpers: parse translate(...) and compute absolute X
      parse_tx <- function(tr) {
        if (is.na(tr) || !nzchar(tr)) return(0)
        m <- regmatches(tr, regexpr("translate\\s*\\(([^)]*)\\)", tr))
        if (!length(m)) return(0)
        nums <- as.numeric(strsplit(sub(".*\\(([^)]*)\\).*", "\\1", m), "[ ,]+")[[1]])
        if (length(nums) >= 1 && is.finite(nums[1])) nums[1] else 0
      }
      abs_x <- function(node) {
        x <- suppressWarnings(as.numeric(xml2::xml_attr(node, "x"))); if (is.na(x)) x <- 0
        cur <- node
        repeat {
          x <- x + parse_tx(xml2::xml_attr(cur, "transform"))
          par <- xml2::xml_parent(cur)
          if (inherits(par, "xml_missing")) break
          cur <- par
        }
        x
      }

      # average glyph width in em for Arial-ish fonts; adjust if needed
      char_width_em <- 0.56

      for (i in seq_along(lab_nodes)) {
        node <- lab_nodes[i]
        txt   <- xml2::xml_text(node)
        nchar_txt <- nchar(txt)

        # current styling
        fs <- suppressWarnings(as.numeric(xml2::xml_attr(node, "font-size")))
        if (is.na(fs)) fs <- 14  # your code sets 14; keep in sync
        old_anchor <- xml2::xml_attr(node, "text-anchor")
        if (is.na(old_anchor) || !nzchar(old_anchor)) old_anchor <- "start"

        # we always want center anchoring going forward
        xml2::xml_attr(node, "text-anchor") <- "middle"
        xml2::xml_attr(node, "dx") <- NULL  # avoid compounding shifts

        # compute an adjustment so the label is centered *around its current placement*
        # If RIdeogram used left-edge (start), shift right by ~half text width.
        # If it used right-edge (end), shift left by ~half text width.
        # If it already used middle, no shift needed.
        half_px <- 0.5 * nchar_txt * fs * char_width_em
        dx <- switch(tolower(old_anchor),
                     "start"  =  +half_px,
                     "end"    =  -half_px,
                     "middle" =  0,
                     # unknown -> assume start
                     +half_px
        )

        dx <- dx - 4

        # Apply dx via a local translate that preserves Y and all ancestors.
        old_tr <- xml2::xml_attr(node, "transform")
        new_tr <- if (!is.na(old_tr) && nzchar(old_tr)) {
          sprintf("translate(%g,0) %s", dx, old_tr)
        } else {
          sprintf("translate(%g,0)", dx)
        }
        xml2::xml_attr(node, "transform") <- new_tr
      }
    }

    # === Add top padding by shifting the viewBox (no reparenting) ===
    y_pad <- 80  # px; increase if your header still touches the top

    root <- xml2::xml_root(svg)

    # Ensure viewBox exists; if not, synthesize from width/height
    vb <- xml2::xml_attr(root, "viewBox")
    if (is.na(vb) || !nzchar(vb)) {
      w <- suppressWarnings(as.numeric(sub("px$", "", xml2::xml_attr(root, "width"))))
      h <- suppressWarnings(as.numeric(sub("px$", "", xml2::xml_attr(root, "height"))))
      if (is.na(w) || is.na(h)) {
        # fallback to a reasonable default if dimensions are missing
        w <- 1200; h <- 800
      }
      xml2::xml_attr(root, "viewBox") <- sprintf("0 0 %g %g", w, h)
      vb <- xml2::xml_attr(root, "viewBox")
    }

    # Shift min-y up by y_pad and extend height by y_pad
    nums <- as.numeric(strsplit(vb, "[ ,]+")[[1]])
    if (length(nums) == 4) {
      nums[2] <- nums[2] - y_pad   # minY becomes negative -> adds top space
      nums[4] <- nums[4] + y_pad   # increase viewBox height
      xml2::xml_attr(root, "viewBox") <- paste(nums, collapse = " ")
    }

    # If explicit pixel height is set, bump it so nothing is clipped when exporting
    old_h <- suppressWarnings(as.numeric(sub("px$", "", xml2::xml_attr(root, "height"))))
    if (!is.na(old_h)) {
      xml2::xml_attr(root, "height") <- paste0(old_h + y_pad, "px")
    }

    # === Solid white background behind everything (robust) ===
    root <- xml2::xml_root(svg)

    # 1) Ensure viewBox exists, then parse it
    vb <- xml2::xml_attr(root, "viewBox")
    if (is.na(vb) || !nzchar(vb)) {
      w <- suppressWarnings(as.numeric(sub("px$", "", xml2::xml_attr(root, "width"))))
      h <- suppressWarnings(as.numeric(sub("px$", "", xml2::xml_attr(root, "height"))))
      if (is.na(w) || is.na(h)) { w <- 1200; h <- 800 }
      xml2::xml_attr(root, "viewBox") <- sprintf("0 0 %g %g", w, h)
      vb <- xml2::xml_attr(root, "viewBox")
    }
    nums <- as.numeric(strsplit(vb, "[ ,]+")[[1]])
    stopifnot(length(nums) == 4)
    minX <- nums[1]; minY <- nums[2]; vw <- nums[3]; vh <- nums[4]

    # 2) Remove any previous bg rect(s)
    xml2::xml_find_all(root, ".//rect[@id='svg-bg-rect']") |> xml2::xml_remove()

    # 3) Find the FIRST non-<defs> child directly under <svg>
    first_draw <- xml2::xml_find_first(root, "./*[not(self::defs)][1]")

    # 4) Build the background node as a tiny fragment and insert it
    bg_markup <- sprintf(
      "<rect id='svg-bg-rect' x='%g' y='%g' width='%g' height='%g' fill='#ffffff' fill-opacity='1' stroke='none' pointer-events='none'/>",
      minX, minY, vw, vh
    )
    bg_node <- xml2::read_xml(bg_markup)

    if (!inherits(first_draw, "xml_missing")) {
      # Insert BEFORE the first drawable child so it paints underneath everything else
      xml2::xml_add_sibling(first_draw, bg_node, .where = "before")
    } else {
      # No drawable children? Just append under <svg>; it will be the only child anyway
      xml2::xml_add_child(root, bg_node)
    }
    svg <- crop_svg_by_bottom_margin(svg, cut_px = 350)
    # --- save & return ---
    xml2::write_xml(svg, "chromosome.svg")
    file.rename("chromosome.svg", out_svg)

    if (make_png) {
      RIdeogram::convertSVG(out_svg, device = "png")
      if (file.exists("chromosome.png")) file.rename("chromosome.png", out_png)
    }

    if (verbose) message("[dndsR::dnds_ideogram] Wrote: ", out_svg,
                         if (make_png && file.exists(out_png)) paste0(" and ", out_png) else "")
    out_svg
  }  # <- this closes .one_side()

  # ---------- batch vs single ----------
  outs <- character(0)

  if (!is.null(comparison_file)) {
    df <- .read_comparisons(comparison_file)
    for (i in seq_len(nrow(df))) {
      comp <- df$comparison_name[i]
      comp_dir <- file.path(output_dir, comp)
      if (!dir.exists(comp_dir)) { if (verbose) message("Missing comparison dir: ", comp_dir); next }

      dnds_path <- .find_dnds_table(comp_dir, comp)
      if (is.na(dnds_path)) { if (verbose) message("Missing dN/dS table for ", comp); next }

      inter_path <- file.path(comp_dir, paste0(comp, "_ideogram_input.tsv"))
      need_build <- TRUE
      if (file.exists(inter_path)) need_build <- file.mtime(dnds_path) > file.mtime(inter_path)

      if (need_build) {
        req_cli("Indexing GFFs for %s", comp)
        q_map <- .gff_index(df$query_gff[i], restrict_gene = restrict_gff_to_gene)
        s_map <- .gff_index(df$subject_gff[i], restrict_gene = restrict_gff_to_gene)

        req_cli("Joining coords \u2192 intermediate: %s", basename(inter_path))
        D <- utils::read.table(dnds_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
        need0 <- c("dNdS","query_id","subject_id")
        miss0 <- setdiff(need0, names(D))
        if (length(miss0)) die("[%s] table missing columns: %s", comp, paste(miss0, collapse = ", "))

        D <- .augment_with_coords(D, q_map, "query_id",   "q")
        D <- .augment_with_coords(D, s_map, "subject_id", "s")

        # Optional normalization of chr labels in the intermediate (opt-in; raw by default)
        D$q_chr <- .normalize_chr(D$q_chr,
                                  strip_chr0 = chr_strip_leading_chr0,
                                  strip_lead = chr_strip_leading,
                                  strip_trail= chr_strip_trailing,
                                  ci = chr_case_insensitive)
        D$s_chr <- .normalize_chr(D$s_chr,
                                  strip_chr0 = chr_strip_leading_chr0,
                                  strip_lead = chr_strip_leading,
                                  strip_trail= chr_strip_trailing,
                                  ci = chr_case_insensitive)

        utils::write.table(D, inter_path, sep = "\t", row.names = FALSE, quote = FALSE)
      } else {
        req_cli("Using existing intermediate: %s", basename(inter_path))
      }

      for (sd in sides) {
        fasta <- if (sd == "query") df$query_fasta[i] else df$subject_fasta[i]
        if (!file.exists(fasta)) { if (verbose) message("Missing FASTA (", sd, "): ", fasta); next }
        outs <- c(outs, .one_side(
          inter_path, sd, fasta,
          window_size, max_dnds, filter_expr,
          make_png, overwrite,
          chr_strip_leading_chr0, chr_strip_leading, chr_strip_trailing, chr_case_insensitive,
          ideogram_args = ideogram_dots
        ))
      }

      if (!keep_intermediate && file.exists(inter_path)) unlink(inter_path)
    }
    if (verbose) message("All ideograms complete.")
    return(invisible(stats::na.omit(outs)))
  }

  # single mode (expects a prebuilt intermediate with q_/s_ coords; choose exactly one side)
  if (is.null(dnds_merged_file)) die("Provide either comparison_file (batch) OR dnds_merged_file (single).")
  if (!file.exists(dnds_merged_file)) die("dnds_merged_file not found: %s", dnds_merged_file)
  if (length(sides) != 1L) die("Single mode requires exactly one 'side' in 'sides'.")

  # Require caller to provide the correct side's FASTA via a wrapper/CLI; this
  # exported function is primarily intended for batch mode.
  stop("Single-mode expects prebuilt intermediate and a FASTA supplied by a wrapper. ",
       "Run via batch (comparison_file) or build your own wrapper that calls .one_side().")
}
