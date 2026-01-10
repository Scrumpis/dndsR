#' GO enrichment with topGO (DAG-aware) for query/subject columns
#'
#' Uses topGO algorithms (classic/elim/weight/weight01) to test enrichment of GO terms
#' among positively selected pairs (dNdS > pos_threshold) vs background, with an
#' annotation-aware universe (drop rows with no GO for that side) by default.
#'
#' @param dnds_annot_file Path to single <comp>_dnds_annot.tsv (single mode).
#' @param comparison_file Batch mode: whitespace-delimited file with
#'   columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#' @param output_dir Root directory that contains per-comparison folders (batch mode).
#' @param sides c("query","subject") which side(s) to analyze.
#' @param ontologies c("BP","MF","CC") which GO branches.
#' @param algorithm topGO algorithm: one of "weight01","elim","classic","weight".
#' @param statistic topGO statistic: usually "fisher" (default) or "ks".
#' @param pos_threshold Numeric; dNdS > pos_threshold defines "positive" (default 1).
#' @param max_dnds Drop rows with dNdS >= max_dnds or NA (default 10).
#' @param filter_expr Optional character filter evaluated in the data
#'   (e.g., "q_gff_seqname == s_gff_seqname").
#' @param drop_rows_without_go If TRUE (default), universe is genes with >=1 GO for that side.
#' @param node_min Minimum term size (topGO nodeSize) to keep (default 10).
#' @param node_max Maximum term size (filtered after test; default Inf).
#' @param p_adjust "none"|"BH" for multiple-testing correction across tested nodes (default "none").
#' @param make_plots If TRUE, writes a top-N bubble plot per result.
#' @param top_n Number of rows to plot.
#' @param exclude_gos Character vector of GO IDs to exclude globally
#'   from positives/background and from tested nodes (default NULL).
#' @param exclude_descendants If TRUE, also exclude descendants of each GO ID (default FALSE).
#' @param exclude_descendants_depth Integer depth limit for descendant exclusion;
#'   1 = direct children only, Inf = entire subtree (default Inf).
#' @param exclude_descendants_limit Hard cap on total excluded nodes per (side,ontology)
#'   run to avoid accidentally nuking huge swaths (default 5000).
#' @param include_definition If TRUE (default), append TERM definition from GO.db.
#' @param alpha Numeric FDR reference level for color scale (default 0.05).
#' @param verbose If TRUE, prints summary of GO filtering (default FALSE).
#' @param threads Integer; maximum number of parallel workers (forked processes; default 4). 
#'   In batch mode, workers are applied across (comparison, side) jobs. 
#'   On Windows, runs sequentially (no forking).
#' @param ... Additional arguments passed on to 'topGO::runTest()', e.g. `score`,
#'   `ties.method`, etc. Core arguments `object`, `algorithm`, and `statistic`
#'   are always set by dndsR and cannot be overridden.
#'
#' @return Invisibly: vector of output TSV paths (batch or single).
#' @export
go_enrichment <- function(
  dnds_annot_file = NULL,
  comparison_file = NULL,
  output_dir = getwd(),
  sides = c("query", "subject"),
  ontologies = c("BP", "MF", "CC"),
  algorithm = c("weight01", "elim", "classic", "weight"),
  statistic = c("fisher", "ks"),
  pos_threshold = 1,
  max_dnds = 10,
  filter_expr = NULL,
  drop_rows_without_go = TRUE,
  node_min = 10,
  node_max = Inf,
  p_adjust = c("none", "BH"),
  make_plots = TRUE,
  top_n = 20,
  exclude_gos = NULL,
  exclude_descendants = FALSE,
  exclude_descendants_depth = Inf,
  exclude_descendants_limit = 5000,
  include_definition = TRUE,
  alpha = 0.05,
  verbose = FALSE,
  threads = 4,
  ...
) {

  # ---- deps ----
  if (!requireNamespace("topGO", quietly = TRUE)) {
    stop("Package 'topGO' is required. Install via Bioconductor.", call. = FALSE)
  }
  if (!requireNamespace("GO.db", quietly = TRUE)) {
    stop("Package 'GO.db' is required. Install via Bioconductor.", call. = FALSE)
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required. Install via Bioconductor.", call. = FALSE)
  }

  # Ensure topGO is ATTACHED so its on-attach side effects (GOBPTerm/GOMFTerm/GOCCTerm)
  # exist for topGO internals that refer to them as plain symbols.
  if (!"package:topGO" %in% search()) {
    suppressMessages(base::attachNamespace("topGO"))
  }
  
  # sanity check (now safe to look in the attached package env)
  topgo_pkg <- as.environment("package:topGO")
  ok_terms <- exists("GOBPTerm", envir = topgo_pkg, inherits = FALSE) &&
              exists("GOMFTerm", envir = topgo_pkg, inherits = FALSE) &&
              exists("GOCCTerm", envir = topgo_pkg, inherits = FALSE)
  
  if (!ok_terms) {
    stop(
      "topGO did not create GOBPTerm/GOMFTerm/GOCCTerm on attach. ",
      "This usually indicates a broken/partial topGO installation.",
      call. = FALSE
    )
  }

  # extra args to pass to topGO::runTest()
  topgo_dots <- list(...)

  sides      <- match.arg(sides, choices = c("query", "subject"), several.ok = TRUE)
  ontologies <- match.arg(ontologies, choices = c("BP", "MF", "CC"), several.ok = TRUE)
  algorithm  <- match.arg(algorithm, choices = c("weight01", "elim", "classic", "weight"))
  statistic  <- match.arg(statistic, choices = c("fisher", "ks"))
  p_adjust   <- match.arg(p_adjust, choices = c("none", "BH"))
  threads <- suppressWarnings(as.integer(threads))
  if (is.na(threads) || threads < 1L) threads <- 1L
  message(sprintf("[go_enrichment] threads=%d (pid=%d)", threads, Sys.getpid()))
  if (isTRUE(make_plots) && threads > 8L && !is.null(comparison_file)) {
    warning(
      "[go_enrichment] make_plots=TRUE with threads>8 may be I/O-heavy (many SVG writes). ",
      "Consider lowering threads.",
      call. = FALSE
    )
  }

  # --- helpers ---
  .read_ws <- function(path, header_try = TRUE) {
    utils::read.table(
      path,
      header = header_try,
      sep = "",
      quote = "\"",
      stringsAsFactors = FALSE,
      comment.char = "",
      strip.white = TRUE,
      blank.lines.skip = TRUE,
      check.names = FALSE
    )
  }

  .read_comparisons <- function(x) {
    req <- c("comparison_name", "query_fasta", "query_gff", "subject_fasta", "subject_gff")
    if (is.data.frame(x)) {
      stopifnot(all(req %in% names(x)))
      return(x[, req, drop = FALSE])
    }
    df1 <- try(.read_ws(x, header_try = TRUE), silent = TRUE)
    if (!inherits(df1, "try-error") && all(req %in% names(df1))) {
      return(df1[, req, drop = FALSE])
    }
    df2 <- .read_ws(x, header_try = FALSE)
    stopifnot(ncol(df2) >= 5)
    names(df2)[1:5] <- req
    df2[, req, drop = FALSE]
  }

  .apply_filter <- function(d) {
    keep <- !is.na(d$dNdS) & d$dNdS < max_dnds

    if (!is.null(filter_expr) && nzchar(filter_expr)) {
      ok <- try(eval(parse(text = filter_expr), envir = d, enclos = parent.frame()), silent = TRUE)

      if (inherits(ok, "try-error")) {
        warning("[go_enrichment] filter_expr failed: ", filter_expr, call. = FALSE)
      } else {
        ok <- as.logical(ok)

        if (length(ok) == 1L) ok <- rep(ok, nrow(d))

        if (length(ok) != nrow(d)) {
          warning("[go_enrichment] filter_expr returned length ", length(ok),
                  " but expected ", nrow(d), "; ignoring filter.", call. = FALSE)
        } else {
          keep <- keep & ok
        }
      }
    }

    d[keep, , drop = FALSE]
  }

  # Named logical vector (GOID = TRUE) for fast membership checks
  .build_exclude_set <- function(excl) {
    if (is.null(excl) || !length(excl)) return(NULL)
    x <- unique(as.character(excl))
    stats::setNames(rep(TRUE, length(x)), x)
  }

  # Depth-/limit-aware descendant expansion using CHILDREN maps
  .expand_descendants <- function(seeds, ontology, depth = Inf, limit = Inf) {
    if (!length(seeds)) return(seeds)
    ontology <- as.character(ontology)
    if (length(ontology) != 1L) ontology <- ontology[[1L]]

    child_map <- switch(
      ontology,
      BP = GO.db::GOBPCHILDREN,
      MF = GO.db::GOMFCHILDREN,
      CC = GO.db::GOCCCHILDREN,
      stop("Unknown ontology: ", ontology, call. = FALSE)
    )

    seen <- unique(seeds)
    frontier <- unique(seeds)
    cur_depth <- 0L

    while (length(frontier) && cur_depth < depth) {
      kids <- unique(
        unlist(
          AnnotationDbi::mget(frontier, child_map, ifnotfound = NA),
          use.names = FALSE
        )
      )
      kids <- kids[!is.na(kids)]
      kids <- setdiff(kids, seen)
      if (!length(kids)) break
      seen <- c(seen, kids)
      if (is.finite(limit) && length(seen) > limit) {
        warning(
          sprintf(
            "Descendant exclusion capped at %d nodes for ontology %s.",
            limit, ontology
          ),
          call. = FALSE
        )
        seen <- unique(seen)[seq_len(limit)]
        break
      }
      frontier <- kids
      cur_depth <- cur_depth + 1L
    }
    unique(seen)
  }

  # Extract GO IDs from one row's GO string (semicolon-delimited), with exclusions applied
  .extract_go_ids <- function(s, exclude_set = NULL) {
    s <- as.character(s)
    if (is.na(s) || !nzchar(s)) return(character(0))
    v <- strsplit(s, ";", fixed = TRUE)[[1L]]
    ids <- unique(grep("^GO:\\d{7}$", v, value = TRUE))
    if (!length(ids)) return(character(0))

    if (!is.null(exclude_set)) {
      ex <- as.logical(exclude_set[ids])
      ex[is.na(ex)] <- FALSE
      ids <- ids[!ex]
    }
    ids
  }

  # ---- plotting helpers (match IPR style) ----
  .upper_padj <- function(d, alpha) {
    max_p <- suppressWarnings(max(d$p_adj, na.rm = TRUE))
    if (!is.finite(max_p)) max_p <- alpha
    eps <- max(1e-12, alpha * 1e-6)
    max(max_p, alpha + eps)
  }

  .padj_scale <- function(alpha) {
    # Clamp colors above alpha to alpha (so everything non-significant shares the same top color)
    oob_fun <- if (requireNamespace("scales", quietly = TRUE)) scales::squish else NULL
    ggplot2::scale_color_viridis_c(
      option = "viridis",
      direction = -1,
      limits = c(0, alpha),
      oob = oob_fun,
      name = "adj p"
    )
  }

  # ---------- font & plotting helpers (match ipr_enrichment) ----------
  .bundled_arial_path <- function() {
    ttf <- system.file("fonts", "ArialBold.ttf", package = "dndsR")
    if (!nzchar(ttf) || !file.exists(ttf)) return(NULL)
    ttf
  }

  .setup_showtext <- function() {
    ttf <- .bundled_arial_path()
    if (is.null(ttf)) return(NULL)

    if (!requireNamespace("showtext", quietly = TRUE)) return(NULL)
    if (!requireNamespace("sysfonts", quietly = TRUE)) return(NULL)

    ok <- try({
      sysfonts::font_add(family = "Arial Bold", regular = ttf)
      showtext::showtext_auto(TRUE)
      TRUE
    }, silent = TRUE)

    if (isTRUE(ok)) "Arial Bold" else NULL
  }

  .register_svglite_mapping <- function() {
    if (!requireNamespace("svglite", quietly = TRUE)) return(NULL)
    ttf <- .bundled_arial_path()
    if (is.null(ttf)) return(NULL)

    function(file, ...) svglite::svglite(
      file,
      user_fonts = list(`Arial Bold` = ttf),
      ...
    )
  }

  .pick_sans_family <- function() {
    fam <- .setup_showtext()
    if (!is.null(fam)) return(fam)
    "sans"
  }

  .svg_device <- function() {
    dev_map <- .register_svglite_mapping()
    if (!is.null(dev_map)) return(dev_map)
    if (requireNamespace("svglite", quietly = TRUE)) return(function(file, ...) svglite::svglite(file, ...))
    NULL
  }
                                             
  # ---- FAST ontology/obsolete cache (preserves logic, avoids per-gene select()) ----
  .build_go_cache <- function(go_vec, exclude_set = NULL, verbose = FALSE) {
    go_vec <- as.character(go_vec)

    raw_list <- lapply(go_vec, .extract_go_ids, exclude_set = exclude_set)
    all_ids <- unique(unlist(raw_list, use.names = FALSE))
    all_ids <- all_ids[!is.na(all_ids) & nzchar(all_ids)]

    if (!length(all_ids)) {
      return(list(raw_list = raw_list, id2ont = character(0)))
    }

    # Drop obsolete ONCE
    obs <- AnnotationDbi::mget(all_ids, GO.db::GOOBSOLETE, ifnotfound = NA)
    is_obs <- vapply(
      obs,
      function(x) is.character(x) && length(x) >= 1L && nzchar(x[1L]),
      logical(1)
    )
    ids2 <- all_ids[!is_obs]
    if (!length(ids2)) {
      return(list(raw_list = raw_list, id2ont = character(0)))
    }

    # Lookup ontology ONCE
    tab <- try(
      if (isTRUE(verbose)) {
        AnnotationDbi::select(
          GO.db::GO.db,
          keys = ids2,
          keytype = "GOID",
          columns = "ONTOLOGY"
        )
      } else {
        suppressMessages(
          AnnotationDbi::select(
            GO.db::GO.db,
            keys = ids2,
            keytype = "GOID",
            columns = "ONTOLOGY"
          )
        )
      },
      silent = TRUE
    )

    if (inherits(tab, "try-error") || is.null(tab) || !nrow(tab)) {
      return(list(raw_list = raw_list, id2ont = character(0)))
    }

    tab <- tab[!is.na(tab$GOID) & !is.na(tab$ONTOLOGY), c("GOID", "ONTOLOGY"), drop = FALSE]
    # If select() returns duplicates, keep the first ONTOLOGY per GOID deterministically
    tab <- tab[order(tab$GOID), , drop = FALSE]
    tab <- tab[!duplicated(tab$GOID), , drop = FALSE]
    id2ont <- stats::setNames(as.character(tab$ONTOLOGY), as.character(tab$GOID))

    list(raw_list = raw_list, id2ont = id2ont)
  }

  # Build gene2GO mapping, ontology-aware and GO.db-aware (cached)
  .make_gene2go_cached <- function(ids, go_vec, ont, exclude_set = NULL, verbose = FALSE) {
    ids <- as.character(ids)
    go_vec <- as.character(go_vec)
    ont <- as.character(ont)[1L]

    cache <- .build_go_cache(go_vec, exclude_set = exclude_set, verbose = verbose)
    raw_list <- cache$raw_list
    id2ont <- cache$id2ont

    raw_n <- sum(lengths(raw_list))

    filt_list <- lapply(raw_list, function(x) {
      if (!length(x) || !length(id2ont)) return(character(0))
      o <- id2ont[x]
      x[!is.na(o) & o == ont]
    })
    filt_n <- sum(lengths(filt_list))

    keep <- lengths(filt_list) > 0L
    mapping <- stats::setNames(filt_list[keep], ids[keep])

    if (isTRUE(verbose)) {
      message(
        sprintf(
          "  [go_enrichment] %s: GO tokens=%d, kept_after_filter=%d, genes_with_go=%d/%d",
          ont, raw_n, filt_n, sum(keep), length(ids)
        )
      )
    }

    list(
      mapping = mapping,
      genes = ids[keep],
      raw_token_count = raw_n,
      kept_token_count = filt_n
    )
  }                                                            
                                                            
  base_exclude_set <- .build_exclude_set(exclude_gos)

  .run_one <- function(d, side, ont, comp, comp_dir, runTest_args = NULL) {
    ont <- as.character(ont)
    if (length(ont) != 1L) ont <- ont[[1L]]

    id_col <- if (side == "query") "query_id" else "subject_id"
    go_col <- if (side == "query") "q_go" else "s_go"
    if (!(go_col %in% names(d))) return(NULL)
    d[[go_col]] <- as.character(d[[go_col]])

    d <- .apply_filter(d)
    if (!nrow(d)) return(NULL)

    if (drop_rows_without_go) {
      d <- d[!is.na(d[[go_col]]) & d[[go_col]] != "", , drop = FALSE]
    }
    if (!nrow(d)) return(NULL)

    pos_ids <- d[[id_col]][d$dNdS > pos_threshold]

    exclude_set_run <- base_exclude_set
    if (isTRUE(exclude_descendants) &&
        !is.null(exclude_set_run) &&
        length(exclude_set_run)) {
      seeds <- names(exclude_set_run)
      seeds_exp <- .expand_descendants(
        seeds,
        ont,
        depth = exclude_descendants_depth,
        limit = exclude_descendants_limit
      )
      exclude_set_run <- .build_exclude_set(seeds_exp)
    }

    g2 <- .make_gene2go_cached(
      d[[id_col]],
      d[[go_col]],
      ont = ont,
      exclude_set = exclude_set_run,
      verbose = verbose
    )
    gene2GO <- g2$mapping
    if (!length(gene2GO)) return(NULL)

    universe <- unique(names(gene2GO))
    pos_ids <- intersect(pos_ids, universe)

    geneList <- factor(as.integer(universe %in% pos_ids))
    names(geneList) <- universe

    tgd <- methods::new(
      "topGOdata",
      ontology = ont,
      allGenes = geneList,
      annot = topGO::annFUN.gene2GO,
      gene2GO = gene2GO,
      nodeSize = as.integer(node_min)
    )

    if (is.null(runTest_args)) runTest_args <- list()
    protected <- c("object", "algorithm", "statistic")
    bad <- intersect(names(runTest_args), protected)
    if (length(bad)) {
      warning(
        "[dndsR::go_enrichment] Ignoring arguments in ... that conflict with core topGO::runTest params: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
      runTest_args[bad] <- NULL
    }

    run_args <- list(
      object = tgd,
      algorithm = algorithm,
      statistic = statistic
    )
    res <- do.call(topGO::runTest, c(run_args, runTest_args))

    raw_scores <- topGO::score(res)
    all_terms <- names(raw_scores)
    if (!length(all_terms)) return(NULL)

    if (!is.null(exclude_set_run) && length(all_terms)) {
      all_terms <- setdiff(all_terms, names(exclude_set_run))
    }
    if (!length(all_terms)) return(NULL)

    term_genes <- topGO::genesInTerm(tgd, all_terms)
    annotated <- vapply(term_genes, length, integer(1))
    sig_counts <- vapply(
      term_genes,
      function(gs) {
        sum(
          as.integer(as.character(geneList[gs])) == 1L,
          na.rm = TRUE
        )
      },
      integer(1)
    )

    keep_idx <- annotated <= node_max
    all_terms <- all_terms[keep_idx]
    annotated <- annotated[keep_idx]
    sig_counts <- sig_counts[keep_idx]
    raw_p <- unname(raw_scores[all_terms])

    cols_want <- c("TERM")
    if (isTRUE(include_definition)) cols_want <- c("TERM", "DEFINITION")

    term_map <- try(
      if (isTRUE(verbose)) {
        AnnotationDbi::select(
          GO.db::GO.db,
          keys = all_terms,
          keytype = "GOID",
          columns = cols_want
        )
      } else {
        suppressMessages(
          AnnotationDbi::select(
            GO.db::GO.db,
            keys = all_terms,
            keytype = "GOID",
            columns = cols_want
          )
        )
      },
      silent = TRUE
    )
    
    if (inherits(term_map, "try-error") || is.null(term_map) || !nrow(term_map)) {
      term_name <- rep(NA_character_, length(all_terms))
      term_def <- rep(NA_character_, length(all_terms))
    } else {
      term_map <- term_map[match(all_terms, term_map$GOID), , drop = FALSE]
      term_name <- term_map$TERM
      term_def <- if ("DEFINITION" %in% names(term_map)) {
        term_map$DEFINITION
      } else {
        NA_character_
      }
    }

    valid_idx <- which(is.finite(raw_p) & raw_p < 1 & sig_counts > 0)
    p_adj <- rep(1, length(raw_p))
    if (length(valid_idx)) {
      p_adj[valid_idx] <- if (p_adjust == "BH") {
        stats::p.adjust(raw_p[valid_idx], method = "BH")
      } else {
        raw_p[valid_idx]
      }
    }

    n_pos <- sum(as.integer(as.character(geneList)) == 1L)
    if (n_pos == 0L) return(NULL)
    n_all <- length(geneList)
    enrichment <- (sig_counts / n_pos) / (annotated / n_all)

    out <- data.frame(
      GO_ID = all_terms,
      term_name = term_name,
      term_def = if (isTRUE(include_definition)) term_def else NULL,
      annotated = annotated,
      significant = sig_counts,
      raw_p = raw_p,
      p_adj = p_adj,
      enrichment = enrichment,
      label = ifelse(
        is.na(term_name) | !nzchar(term_name),
        all_terms,
        paste0(all_terms, " ", term_name)
      ),
      side = side,
      ontology = ont,
      comparison = comp,
      stringsAsFactors = FALSE
    )

    out_file <- file.path(
      comp_dir,
      sprintf(
        "%s_%s_GO_%s_topGO_%s_%s.tsv",
        comp,
        if (side == "query") "q" else "s",
        ont,
        algorithm,
        statistic
      )
    )

    utils::write.table(
      out[order(out$p_adj, -out$enrichment, out$GO_ID), ],
      file = out_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    if (isTRUE(make_plots)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning(
        "Skipping plots: package 'ggplot2' is required for plotting.",
        call. = FALSE
      )
    } else {
                plt <- out[order(out$p_adj, -out$enrichment, out$GO_ID), , drop = FALSE]
        plt <- utils::head(plt, top_n)

        # Explicit ordering for the y-axis: top = most significant (smallest p_adj)
        plt$is_significant <- is.finite(plt$p_adj) & !is.na(plt$p_adj) & (plt$p_adj <= alpha)

        # Build a stable ordered factor for y
        plt$y_lab <- plt$label
        plt <- plt[order(plt$p_adj, -plt$enrichment, plt$GO_ID), , drop = FALSE]
        plt$y_lab <- factor(plt$y_lab, levels = rev(plt$y_lab))

        # Boundary line between significant and non-significant (only if both exist in top_n)
        sig_idx <- which(plt$is_significant)
        cut_y <- NULL
        if (length(sig_idx) > 0L && length(sig_idx) < nrow(plt)) {
          cut_y <- max(sig_idx) + 0.5
        }

        base_family <- .pick_sans_family()

        gg <- ggplot2::ggplot(
          plt,
          ggplot2::aes(
            x = enrichment,
            y = y_lab,
            size = significant,
            color = p_adj
          )
        ) +
          ggplot2::geom_point() +
          .padj_scale(alpha) +
          ggplot2::labs(
            x = "Enrichment (pos/bg)",
            y = sprintf("GO %s (%s)", ont, side),
            size = "# pos"
          ) +
          ggplot2::theme_minimal(base_size = 12, base_family = base_family)

        if (!is.null(cut_y)) {
          gg <- gg +
            ggplot2::geom_hline(yintercept = cut_y, linetype = "dashed", linewidth = 0.6) +
            ggplot2::annotate(
              "text",
              x = -Inf,
              y = cut_y,
              label = sprintf("adj p \u2264 %.3g", alpha),
              hjust = -0.05,
              vjust = -0.4,
              size = 3.5
            ) +
            ggplot2::coord_cartesian(clip = "off")
        }

        out_svg <- sub("\\.tsv$", "_topN.svg", out_file)
        dev_fun <- .svg_device()

        if (!is.null(dev_fun)) {
          ggplot2::ggsave(out_svg, gg, device = dev_fun, width = 11, height = 9)
        } else {
          ggplot2::ggsave(out_svg, gg, width = 11, height = 9)
        }
      }
    }
    out_file
  }

  .run_comp_side <- function(comp, comp_dir, side) {
    in_file <- file.path(comp_dir, paste0(comp, "_dnds_annot.tsv"))
    if (!file.exists(in_file)) {
      warning("Missing annotated file: ", in_file, call. = FALSE)
      return(character(0))
    }
    d <- utils::read.table(
      in_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )

    paths <- character(0)
    for (ont in ontologies) {
      p <- .run_one(d, side, ont, comp, comp_dir, runTest_args = topgo_dots)
      if (!is.null(p)) paths <- c(paths, p)
    }
    if (!length(paths)) message("No GO enrichments produced for ", comp, " side=", side)
    paths
  }
                                                            
  # ---- batch vs single ----
  if (!is.null(comparison_file)) {
    df <- .read_comparisons(comparison_file)

    jobs <- expand.grid(
      i = seq_len(nrow(df)),
      side = sides,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    run_one_job <- function(k) {
      i <- jobs$i[k]
      side <- jobs$side[k]
      comp <- df$comparison_name[i]
      comp_dir <- file.path(output_dir, comp)
      dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
      # Clean, non-jumbled progress line in the main console:
      message(sprintf("[go_enrichment] %s %s (logging to %s/%s_%s_topGO.log)",
                      comp, side, comp_dir, comp, side))
      log_file <- file.path(comp_dir, sprintf("%s_%s_topGO.log", comp, side))
      
      .dndsr_with_log(
        log_file = log_file,
        tag = "go_enrichment",
        header = c(
          sprintf("[go_enrichment] comp=%s side=%s", comp, side),
          sprintf("[go_enrichment] pid=%d threads=%d", Sys.getpid(), threads)
        ),
        expr = {
          .run_comp_side(comp, comp_dir, side)
        }
      )
    }

    if (threads > 1L && .Platform$OS.type != "windows") {
      outs_list <- parallel::mclapply(
        seq_len(nrow(jobs)),
        run_one_job,
        mc.cores = threads,
        mc.preschedule = FALSE
      )
    } else {
      if (threads > 1L && .Platform$OS.type == "windows") {
        warning("[go_enrichment] threads>1 requested on Windows; running sequentially.", call. = FALSE)
      }
      outs_list <- lapply(seq_len(nrow(jobs)), run_one_job)
    }

    outs <- unlist(outs_list, use.names = FALSE)
    message("All topGO enrichments complete.")
    return(invisible(outs))
  }
                                                            
  # ---- single mode ----
  if (is.null(dnds_annot_file)) {
    stop(
      "Provide either comparison_file (batch) OR dnds_annot_file (single).",
      call. = FALSE
    )
  }

  comp_dir <- dirname(dnds_annot_file)
  comp_name <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file))

  d <- utils::read.table(
    dnds_annot_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  run_one_side <- function(sd) {
    message(sprintf("[go_enrichment] %s %s (logging to %s/%s_%s_topGO.log)",
                    comp_name, sd, comp_dir, comp_name, sd))
    
    log_file <- file.path(comp_dir, sprintf("%s_%s_topGO.log", comp_name, sd))
    
    .dndsr_with_log(
      log_file = log_file,
      tag = "go_enrichment",
      header = c(
        sprintf("[go_enrichment] comp=%s side=%s", comp_name, sd)
      ),
      expr = {
        paths <- character(0)
        for (ont in ontologies) {
          p <- .run_one(d, sd, ont, comp_name, comp_dir, runTest_args = topgo_dots)
          if (!is.null(p)) paths <- c(paths, p)
        }
        paths
      }
    )

  }

  if (threads > 1L && length(sides) > 1L && .Platform$OS.type != "windows") {
    outs_list <- parallel::mclapply(
      sides,
      run_one_side,
      mc.cores = min(threads, length(sides)),
      mc.preschedule = FALSE
    )
  } else {
    if (threads > 1L && .Platform$OS.type == "windows") {
      warning("[go_enrichment] threads>1 requested on Windows; running sequentially.", call. = FALSE)
    }
    outs_list <- lapply(sides, run_one_side)
  }
                                                                                                                    
  outs <- unlist(outs_list, use.names = FALSE)
  message("topGO enrichment complete for: ", dnds_annot_file)
  invisible(outs)                                                            
}
