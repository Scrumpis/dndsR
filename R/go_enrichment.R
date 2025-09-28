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
#'   (e.g., "q_seqname == s_seqname").
#' @param drop_rows_without_go If TRUE (default), universe is genes with ≥1 GO for that side.
#' @param node_min Minimum term size (topGO nodeSize) to keep (default 10).
#' @param node_max Maximum term size (filtered after test; default Inf).
#' @param p_adjust "BH"|"none" for multiple-testing correction across tested nodes.
#'   (Note: topGO already de-correlates via its algorithm; BH is still common.)
#' @param make_plots If TRUE, writes a top-N bubble plot per result.
#' @param top_n Number of rows to plot.
#' @param exclude_gos Optional character vector of GO IDs to exclude globally
#'   from both positives and background before enrichment (e.g., over-broad terms).
#' @param include_definition If TRUE (default), append TERM definition from GO.db.
#'
#' @return Invisibly: vector of output TSV paths (batch) or a data.frame list (single if make_plots=FALSE).
#' @export
go_enrichment <- function(dnds_annot_file = NULL,
                          comparison_file = NULL,
                          output_dir = getwd(),
                          sides = c("query","subject"),
                          ontologies = c("BP","MF","CC"),
                          algorithm = c("weight01","elim","classic","weight"),
                          statistic = c("fisher","ks"),
                          pos_threshold = 1,
                          max_dnds = 10,
                          filter_expr = NULL,
                          drop_rows_without_go = TRUE,
                          node_min = 10,
                          node_max = Inf,
                          p_adjust = c("BH","none"),
                          make_plots = TRUE,
                          top_n = 20,
                          exclude_gos = NULL,
                          include_definition = TRUE) {

  if (!requireNamespace("topGO", quietly = TRUE) || !requireNamespace("GO.db", quietly = TRUE)) {
    stop("Please install Bioconductor packages 'topGO' and 'GO.db'.")
  }
  # Ensure GO.db is ATTACHED (topGO sometimes expects term envs on search path)
  if (!("package:GO.db" %in% search())) {
    suppressPackageStartupMessages(require("GO.db"))
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Please install Bioconductor package 'AnnotationDbi'.")
  }

  sides      <- match.arg(sides,      c("query","subject"), several.ok = TRUE)
  ontologies <- match.arg(ontologies, c("BP","MF","CC"),    several.ok = TRUE)
  algorithm  <- match.arg(algorithm)
  statistic  <- match.arg(statistic)
  p_adjust   <- match.arg(p_adjust)

  # --- helpers ---
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
    df2 <- .read_ws(x, header_try = FALSE); stopifnot(ncol(df2) >= 5)
    names(df2)[1:5] <- req; df2[, req, drop = FALSE]
  }
  .apply_filter <- function(d) {
    keep <- !is.na(d$dNdS) & d$dNdS < max_dnds
    if (!is.null(filter_expr) && nzchar(filter_expr)) {
      ok <- try(eval(parse(text = filter_expr), envir = d, enclos = parent.frame()), silent = TRUE)
      if (!inherits(ok, "try-error")) keep <- keep & isTRUE(as.vector(ok))
    }
    d[keep, , drop = FALSE]
  }

  # Build an exclusion set: named logical vector (GOID=TRUE)
  .build_exclude_set <- function(excl) {
    if (is.null(excl) || !length(excl)) return(NULL)
    x <- unique(as.character(excl))
    stats::setNames(rep(TRUE, length(x)), x)
  }
  exclude_set <- .build_exclude_set(exclude_gos)

  # Split a vector of "GO:xxxxxxx;GO:yyyyyyy" strings into unique GO IDs,
  # dropping excluded GO IDs if provided.
  .split_go_vec <- function(x, exclude_set = NULL) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) return(list())
    toks <- strsplit(x, ";", fixed = TRUE)
    lapply(toks, function(v) {
      ids <- unique(grep("^GO:\\d{7}$", v, value = TRUE))
      if (!is.null(exclude_set) && length(ids)) {
        # OLD (buggy): ids <- ids[!isTRUE(exclude_set[ids])]
        ex <- as.logical(exclude_set[ids])
        ex[is.na(ex)] <- FALSE
        ids <- ids[!ex]
      }
      ids
    })
  }

  # Build gene2GO mapping after optional exclusion
  .make_gene2go <- function(ids, go_col, exclude_set = NULL) {
    spl <- .split_go_vec(go_col, exclude_set)
    keep <- lengths(spl) > 0L
    if (!any(keep)) return(list(mapping = list(), genes = character(0)))
    genes   <- ids[keep]
    mapping <- setNames(spl[keep], genes)
    list(mapping = mapping, genes = genes)
  }

  # Some setups of topGO look for envs like GOBPTerm; ensure stubs exist so runTest() is happy
  .ensure_GOdb_term_envs <- function() {
    pkg_env <- as.environment("package:GO.db")
    mk_env <- function(onto) {
      nm <- paste0("GO", onto, "Term")  # e.g., GOBPTerm
      if (exists(nm, inherits = TRUE)) return(invisible())
      # Provide a minimal env that declares IDs for this ontology
      all_keys <- AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
      ont_map  <- AnnotationDbi::select(GO.db::GO.db, keys = all_keys,
                                        columns = "ONTOLOGY", keytype = "GOID")
      ids <- unique(ont_map$GOID[ont_map$ONTOLOGY == onto])
      e <- new.env(hash = TRUE, parent = emptyenv())
      for (id in ids) assign(id, TRUE, envir = e)
      assign(nm, e, envir = .GlobalEnv)
    }
    mk_env("BP"); mk_env("MF"); mk_env("CC")
  }
  .ensure_GOdb_term_envs()

  .run_one <- function(d, side, ont, comp, comp_dir) {
    prefix <- if (side == "query") "q_" else "s_"
    id_col <- if (side == "query") "query_id" else "subject_id"
    go_col <- paste0(prefix, "go")
    if (!(go_col %in% names(d)) || !is.character(d[[go_col]])) return(NULL)

    d <- .apply_filter(d)
    if (!nrow(d)) return(NULL)

    pos_ids <- d[[id_col]][d$dNdS > pos_threshold]

    # Keep rows that have at least *some* GO text; precise empties are filtered in mapping too
    if (drop_rows_without_go) d <- d[!is.na(d[[go_col]]) & d[[go_col]] != "", , drop = FALSE]
    if (!nrow(d)) return(NULL)

    # Build mapping with exclusions applied
    g2 <- .make_gene2go(d[[id_col]], d[[go_col]], exclude_set = exclude_set)
    gene2GO <- g2$mapping
    if (!length(gene2GO)) return(NULL)

    universe <- unique(names(gene2GO))
    geneList <- factor(as.integer(universe %in% pos_ids))
    names(geneList) <- universe

    tgd <- methods::new("topGOdata",
                        ontology = ont,
                        allGenes = geneList,
                        annot    = topGO::annFUN.gene2GO,
                        gene2GO  = gene2GO,
                        nodeSize = as.integer(node_min))

    res <- topGO::runTest(tgd, algorithm = algorithm, statistic = statistic)

    # All tested terms & raw p-values
    raw_scores <- topGO::score(res)
    all_terms  <- names(raw_scores)
    if (!is.null(exclude_set) && length(all_terms)) {
      all_terms <- setdiff(all_terms, names(exclude_set))
    }
    if (!length(all_terms)) return(NULL)

    # Compute Annotated & Significant without GenTable()
    term_genes <- topGO::genesInTerm(tgd, all_terms)
    annotated  <- vapply(term_genes, length, integer(1))
    sig_counts <- vapply(term_genes, function(gs) sum(as.integer(as.character(geneList[gs])) == 1L, na.rm = TRUE),
                         integer(1))

    # Optional upper cap
    keep_idx <- annotated <= node_max
    all_terms <- all_terms[keep_idx]
    annotated <- annotated[keep_idx]
    sig_counts <- sig_counts[keep_idx]
    raw_p <- unname(raw_scores[all_terms])

    # Term name (+ definition if requested)
    cols_want <- c("TERM")
    if (isTRUE(include_definition)) cols_want <- c("TERM","DEFINITION")
    term_map <- try(
      AnnotationDbi::select(GO.db::GO.db, keys = all_terms, keytype = "GOID", columns = cols_want),
      silent = TRUE
    )
    if (inherits(term_map, "try-error")) {
      term_name <- rep(NA_character_, length(all_terms))
      term_def  <- rep(NA_character_, length(all_terms))
    } else {
      # Ensure 1:1 in GOID order
      term_map <- term_map[match(all_terms, term_map$GOID), , drop = FALSE]
      term_name <- term_map$TERM
      term_def  <- if ("DEFINITION" %in% names(term_map)) term_map$DEFINITION else NA_character_
    }

    # BH (optional) and enrichment
    p_adj <- if (p_adjust == "BH") stats::p.adjust(raw_p, method = "BH") else raw_p
    n_pos <- sum(as.integer(as.character(geneList)) == 1L)
    n_all <- length(geneList)
    enrichment <- (sig_counts / n_pos) / (annotated / n_all)

    out <- data.frame(
      GO_ID       = all_terms,
      term_name   = term_name,
      term_def    = if (isTRUE(include_definition)) term_def else NULL,
      annotated   = annotated,
      significant = sig_counts,
      raw_p       = raw_p,
      p_adj       = p_adj,
      enrichment  = enrichment,
      label       = ifelse(is.na(term_name) | !nzchar(term_name),
                           all_terms, paste0(all_terms, " — ", term_name)),
      side        = side,
      ontology    = ont,
      comparison  = comp,
      stringsAsFactors = FALSE
    )

    out_file <- file.path(comp_dir, sprintf("%s_%s_GO_%s_topGO_%s_%s.tsv",
                                            comp, if (side=="query") "q" else "s",
                                            ont, algorithm, statistic))
    utils::write.table(out[order(out$p_adj, -out$enrichment, out$GO_ID), ],
                       file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

    if (make_plots && requireNamespace("ggplot2", quietly = TRUE)) {
      plt <- out[order(out$p_adj, -out$enrichment), ]
      plt <- utils::head(plt, top_n)
      gg <- ggplot2::ggplot(plt, ggplot2::aes(x = enrichment,
                                              y = stats::reorder(label, -p_adj),
                                              size = significant,
                                              color = p_adj)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient(low = "red", high = "blue") +
        ggplot2::labs(x = "Enrichment (pos/bg)", y = sprintf("GO %s (%s)", ont, side),
                      size = "# pos", color = "adj p") +
        ggplot2::theme_minimal(base_size = 12)
      ggplot2::ggsave(sub("\\.tsv$", "_topN.svg", out_file), gg, width = 11, height = 9)
    }
    out_file
  }

  .run_comp <- function(comp, comp_dir) {
    in_file <- file.path(comp_dir, paste0(comp, "_dnds_annot.tsv"))
    if (!file.exists(in_file)) { warning("Missing annotated file: ", in_file); return(character(0)) }
    d <- utils::read.table(in_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    paths <- character(0)
    for (sd in sides) for (ont in ontologies) {
      p <- .run_one(d, sd, ont, comp, comp_dir)
      if (!is.null(p)) paths <- c(paths, p)
    }
    if (!length(paths)) message("No GO enrichments produced for ", comp)
    paths
  }

  # ---- batch vs single ----
  if (!is.null(comparison_file)) {
    df <- .read_comparisons(comparison_file)
    outs <- character(0)
    for (i in seq_len(nrow(df))) {
      comp <- df$comparison_name[i]
      comp_dir <- file.path(output_dir, comp)
      dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
      outs <- c(outs, .run_comp(comp, comp_dir))
    }
    message("All topGO enrichments complete.")
    return(invisible(outs))
  }

  if (is.null(dnds_annot_file)) stop("Provide either comparison_file (batch) OR dnds_annot_file (single).")
  comp_dir  <- dirname(dnds_annot_file)
  comp_name <- sub("_dnds_annot\\.tsv$", "", basename(dnds_annot_file))
  d <- utils::read.table(dnds_annot_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  outs <- character(0)
  for (sd in sides) for (ont in ontologies) {
    p <- .run_one(d, sd, ont, comp_name, comp_dir)
    if (!is.null(p)) outs <- c(outs, p)
  }
  message("topGO enrichment complete for: ", dnds_annot_file)
  invisible(outs)
}
