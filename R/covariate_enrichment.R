#' Term enrichment adjusted for gene length and dS
#'
#' Tests whether GO and InterPro terms hold more genes with `dNdS > pos_threshold` than their
#' genes' lengths and dS predict. Short genes and genes with small dS pass that threshold far more
#' often, so Fisher's test, as used by `go_enrichment()`, `ipr_enrichment()`, and
#' `term_enrichment()`, can call a term significant because of the kind of gene it holds. See the
#' calibration evidence in `docs/validation/enrichment-confounding.md`.
#'
#' For each universe of annotated genes, a logistic regression gives each gene's probability of
#' being positive from natural splines of log length and log dS, 3 degrees of freedom each, and no
#' probability is allowed below a tenth of the universe's positive rate. Without that floor the test is
#' badly anticonservative when few genes are positive. A term's p-value is the upper tail of the
#' Poisson-binomial distribution of its positive genes, the sum of independent Bernoulli variables
#' with those probabilities. BH then runs once across all tested terms of a term type, in
#' all requested ontologies. A term is tested if at least `min_total` genes carry it, which does not
#' depend on the outcome.
#'
#' GO terms are propagated to all their ancestors, as `go_enrichment(algorithm = "classic")` does.
#' This command does not use topGO's conditional algorithms.
#'
#' This adjusts for two known confounders, and does not make a term evidence of selection. Terms remain
#' leads to test by other means. It treats genes as independent, and treats fitted probabilities as known.
#'
#' @param dnds_annot_file Path to one `<comparison>_dnds_annot.tsv`.
#' @param comparison_file Batch mode: a comparison file, as for the other commands.
#' @param output_dir Root folder with per-comparison folders (batch mode).
#' @param terms Term types to test, among `"go"` and `"ipr"`. A type with no column is skipped.
#' @param sides Sides to test, among `"query"` and `"subject"`.
#' @param ontologies GO branches, among `"BP"`, `"MF"`, and `"CC"`.
#' @param pos_threshold `dNdS > pos_threshold` defines a positive gene (default 1).
#' @param max_dnds Drops rows with `dNdS >= max_dnds` or undefined `dNdS` (default 10).
#' @param filter_expr Optional logical expression evaluated in the table.
#' @param min_total A term is tested if at least this many genes in the universe carry it (default 10).
#' @param min_pos A term is written if at least this many of its genes are positive (default 2). It is
#'   applied after correction, and changes which rows are written, not their p-values.
#' @param threads Maximum parallel workers across comparisons and sides (forked; default 1).
#' @param verbose If `TRUE`, reports what was tested.
#'
#' @return Invisibly, the paths of the TSV files written. Columns: `term`, `term_type`, `ontology`,
#'   `term_name`, `universe`, `pos_total`, `total_count`, `pos_count`, `expected` (positive genes
#'   predicted from length and dS, after the floor), `obs_over_expected`, `fisher_p` (without adjustment, for
#'   comparison), `p_value`, `p_adj` (BH), `n_tested`, and `model`.
#' @export
covariate_enrichment <- function(dnds_annot_file = NULL,
                                 comparison_file = NULL,
                                 output_dir = getwd(),
                                 terms = c("go", "ipr"),
                                 sides = c("query", "subject"),
                                 ontologies = c("BP", "MF", "CC"),
                                 pos_threshold = 1,
                                 max_dnds = 10,
                                 filter_expr = NULL,
                                 min_total = 10,
                                 min_pos = 2,
                                 threads = 1,
                                 verbose = TRUE) {
  terms <- match.arg(terms, choices = c("go", "ipr"), several.ok = TRUE)
  sides <- match.arg(sides, choices = c("query", "subject"), several.ok = TRUE)
  ontologies <- match.arg(ontologies, choices = c("BP", "MF", "CC"), several.ok = TRUE)
  threads <- max(1L, suppressWarnings(as.integer(threads)), na.rm = TRUE)
  if (is.null(dnds_annot_file) == is.null(comparison_file)) {
    stop("Provide either comparison_file (batch) or dnds_annot_file (single).", call. = FALSE)
  }
  files <- if (!is.null(dnds_annot_file)) {
    dnds_annot_file
  } else {
    comparisons <- utils::read.table(comparison_file, header = FALSE, sep = "", quote = "\"", comment.char = "",
                                     stringsAsFactors = FALSE, strip.white = TRUE, blank.lines.skip = TRUE)
    names_found <- setdiff(comparisons[[1]], "comparison_name")
    file.path(output_dir, names_found, paste0(names_found, "_dnds_annot.tsv"))
  }
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) {
    stop("Annotated table not found: ", paste(missing_files, collapse = ", "), ". Run append_annotations first.", call. = FALSE)
  }
  jobs <- expand.grid(file = files, side = sides, stringsAsFactors = FALSE)
  run <- function(i) {
    .covariate_enrichment_side(jobs$file[i], jobs$side[i], terms, ontologies, pos_threshold, max_dnds, filter_expr,
                               min_total, min_pos, verbose)
  }
  out <- if (threads > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = threads, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(nrow(jobs)), run)
  }
  failed <- vapply(out, inherits, logical(1), "try-error")
  if (any(failed)) stop(paste(unlist(out[failed]), collapse = "\n"), call. = FALSE)
  invisible(unlist(out, use.names = FALSE))
}

#' Upper tail of a Poisson-binomial distribution
#'
#' Exact, by dynamic programming, and absorbing at `k`, so the cost is length(p) * k.
#'
#' @param p Success probability of each independent trial.
#' @param k Number of successes.
#' @return `P(S >= k)`.
#' @keywords internal
.poisson_binomial_upper <- function(p, k) {
  if (k <= 0L) return(1)
  if (k > length(p)) return(0)
  f <- c(1, numeric(k))
  for (q in p) {
    top <- f[k + 1L] + f[k] * q
    f[-1L] <- f[-1L] * (1 - q) + f[-(k + 1L)] * q
    f[1L] <- f[1L] * (1 - q)
    f[k + 1L] <- top
  }
  min(1, f[k + 1L])
}

#' Test one universe of genes for term enrichment, adjusted for length and dS
#'
#' @param positive Logical, whether each gene is positive.
#' @param len,ds Gene length and dS, positive numbers.
#' @param gene_terms List with each gene's terms.
#' @param min_total Smallest term to test.
#' @return Data frame with a row per tested term, not corrected.
#' @keywords internal
.covariate_test_universe <- function(positive, len, ds, gene_terms, min_total) {
  n_genes <- length(positive); n_pos <- sum(positive)
  if (n_pos == 0L || n_pos == n_genes) {
    probability <- rep(n_pos / n_genes, n_genes)
    model <- "none: all genes have the same outcome"
  } else {
    fit <- suppressWarnings(stats::glm(positive ~ splines::ns(log(len), 3) + splines::ns(log(ds), 3), family = stats::binomial()))
    probability <- pmax(stats::fitted(fit), 0.1 * n_pos / n_genes)
    model <- "logistic: ns(log length, 3) + ns(log dS, 3); floor 0.1 * positive rate"
    if (!isTRUE(fit$converged)) model <- paste(model, "(did not converge)")
  }
  members <- split(rep(seq_along(gene_terms), lengths(gene_terms)), unlist(gene_terms, use.names = FALSE))
  members <- members[lengths(members) >= min_total]
  if (!length(members)) return(NULL)
  rows <- lapply(members, function(i) {
    k <- sum(positive[i])
    c(total_count = length(i), pos_count = k, expected = sum(probability[i]),
      fisher_p = stats::phyper(k - 1L, n_pos, n_genes - n_pos, length(i), lower.tail = FALSE),
      p_value = .poisson_binomial_upper(probability[i], k))
  })
  m <- do.call(rbind, rows)
  data.frame(term = names(members), universe = n_genes, pos_total = n_pos, m, model = model, stringsAsFactors = FALSE, row.names = NULL)
}

#' @noRd
.covariate_enrichment_side <- function(path, side, terms, ontologies, pos_threshold, max_dnds, filter_expr, min_total, min_pos, verbose) {
  d <- utils::read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
  prefix <- if (side == "query") "q" else "s"
  len_col <- paste0(prefix, "_len")
  needed <- c("dS", "dNdS", len_col)
  if (!all(needed %in% names(d))) {
    stop("covariate_enrichment needs columns ", paste(setdiff(needed, names(d)), collapse = ", "), " in ", path, call. = FALSE)
  }
  d <- .filter_dnds(d, filter_expr = filter_expr, max_dnds = max_dnds)
  ds <- suppressWarnings(as.numeric(d$dS)); len <- suppressWarnings(as.numeric(d[[len_col]]))
  usable <- is.finite(ds) & ds > 0 & is.finite(len) & len > 0
  d <- d[usable, , drop = FALSE]; ds <- ds[usable]; len <- len[usable]
  positive <- as.numeric(d$dNdS) > pos_threshold
  base <- sub("_dnds_annot\\.tsv$", "", basename(path))
  split_terms <- function(x) {
    x <- as.character(x)
    lapply(strsplit(ifelse(is.na(x), "", x), "[;,|]"), function(t) unique(t[nzchar(t) & t != "NA"]))
  }
  written <- character(0)
  test_universe <- function(gene_terms) {
    has <- lengths(gene_terms) > 0L
    if (!any(has)) return(NULL)
    .covariate_test_universe(positive[has], len[has], ds[has], gene_terms[has], min_total)
  }
  for (type in terms) {
    col <- paste0(prefix, "_", type)
    if (!col %in% names(d)) next
    gene_terms <- split_terms(d[[col]])
    if (!any(lengths(gene_terms) > 0L)) next
    res <- if (type == "go") {
      if (!requireNamespace("GO.db", quietly = TRUE)) {
        stop("Package 'GO.db' is required for GO terms. Install it from Bioconductor.", call. = FALSE)
      }
      ids <- unique(unlist(gene_terms, use.names = FALSE))
      maps <- list(BP = GO.db::GOBPANCESTOR, MF = GO.db::GOMFANCESTOR, CC = GO.db::GOCCANCESTOR)
      parts <- lapply(ontologies, function(o) {
        known <- intersect(ids, AnnotationDbi::keys(maps[[o]]))
        ancestors <- AnnotationDbi::mget(known, maps[[o]])
        closure <- stats::setNames(lapply(seq_along(known), function(i) setdiff(c(known[i], ancestors[[i]]), "all")), known)
        # A term with no ancestors in this branch is a branch root, and belongs here as well.
        roots <- c(BP = "GO:0008150", MF = "GO:0003674", CC = "GO:0005575")[[o]]
        if (roots %in% ids) closure[[roots]] <- roots
        r <- test_universe(lapply(gene_terms, function(t) unique(unlist(closure[t], use.names = FALSE))))
        if (!is.null(r)) r$ontology <- o
        r
      })
      do.call(rbind, parts)
    } else {
      r <- test_universe(gene_terms)
      if (!is.null(r)) r$ontology <- NA_character_
      r
    }
    if (is.null(res) || !nrow(res)) next
    res$p_adj <- stats::p.adjust(res$p_value, "BH")
    res$n_tested <- nrow(res)
    res$term_type <- toupper(type)
    res$obs_over_expected <- ifelse(res$expected > 0, res$pos_count / res$expected, NA_real_)
    res$term_name <- .covariate_term_names(res$term, type)
    n_tested <- nrow(res)
    if (max(res$pos_total) < 30L) {
      message(sprintf("[covariate_enrichment] %s %s %s: at most %d positive genes in a universe. Power is low, so few or no terms is not evidence of absence.",
                      base, side, toupper(type), max(res$pos_total)))
    }
    res <- res[res$pos_count >= min_pos, , drop = FALSE]
    res <- res[order(res$p_value, res$term), c("term", "term_type", "ontology", "term_name", "universe", "pos_total", "total_count", "pos_count",
                                               "expected", "obs_over_expected", "fisher_p", "p_value", "p_adj", "n_tested", "model"), drop = FALSE]
    out <- file.path(dirname(path), sprintf("%s_%s_%s_covariate_enrichment.tsv", base, prefix, toupper(type)))
    utils::write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    if (isTRUE(verbose)) {
      message(sprintf("[covariate_enrichment] %s %s %s: %d terms significant at p_adj <= 0.05 of %d tested",
                      base, side, toupper(type), sum(res$p_adj <= 0.05), n_tested))
    }
    written <- c(written, out)
  }
  written
}

#' @noRd
.covariate_term_names <- function(ids, type) {
  if (type == "go") {
    found <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(ids), columns = "TERM", keytype = "GOID"))
    return(found$TERM[match(ids, found$GOID)])
  }
  entries <- system.file("extdata", "interpro_entry.list.tsv", package = "dndsR")
  if (!nzchar(entries)) return(rep(NA_character_, length(ids)))
  e <- utils::read.delim(entries, quote = "", stringsAsFactors = FALSE)
  e$ENTRY_NAME[match(ids, e$ENTRY_AC)]
}
