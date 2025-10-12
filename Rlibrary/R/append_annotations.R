#' Append annotations from query & subject GFFs to dN/dS results (ordered outputs + simple customs)
#'
#' @param dnds_file Path to a dN/dS TSV (single mode).
#' @param query_gff,subject_gff Paths to GFF3 files (single mode).
#' @param output_file Optional output path (single mode).
#' @param comparison_file Path to whitespace-delimited file (tabs/spaces; header or not)
#'   with columns: comparison_name, query_fasta, query_gff, subject_fasta, subject_gff.
#'   If provided, batch mode is used and outputs are written beside each dN/dS file.
#' @param output_dir Root directory containing per-comparison folders (batch mode).
#' @param custom Optional comma-separated patterns like "REX{5},TOM{8}" meaning
#'   PREFIX followed by exactly N digits; produces q_prefix/s_prefix columns.
#' @param threads Integer threads to use for per-ID annotation. If NULL/NA, uses
#'   options(dndsR.threads) or DNDSR_THREADS env; defaults to 1.
#' @param overwrite Logical; in batch mode, skip a comparison if output exists unless TRUE. Default FALSE.
#' @return In single mode: a data.frame. In batch mode: (invisibly) vector of output paths.
#' @export
append_annotations <- function(dnds_file = NULL,
                               query_gff = NULL,
                               subject_gff = NULL,
                               output_file = NULL,
                               comparison_file = NULL,
                               output_dir = getwd(),
                               custom = NULL,
                               threads = NULL,
                               overwrite = FALSE) {

  ## ---------- lightweight verbose helper ----------
  vmsg <- function(...) {
    if (isTRUE(getOption("dndsR.verbose")) || isTRUE(getOption("dndsR.cli.verbose"))) {
      message(...)
    }
  }

  ## ---------- thread resolver & parallel lapply ----------
  .resolve_threads <- function(th) {
    if (is.null(th) || is.na(th)) {
      opt <- getOption("dndsR.threads", NULL)
      if (is.null(opt)) {
        env <- Sys.getenv("DNDSR_THREADS", "")
        if (nzchar(env) && grepl("^[0-9]+$", env)) as.integer(env) else 1L
      } else as.integer(opt)
    } else {
      as.integer(th)
    }
  }

  plapply <- function(X, FUN, threads = 1L, ...) {
    threads <- max(1L, as.integer(threads))
    if (.Platform$OS.type != "windows" && threads > 1L) {
      return(parallel::mclapply(X, FUN, mc.cores = threads, mc.preschedule = TRUE, ...))
    } else {
      return(lapply(X, FUN, ...))
    }
  }

  ## ---------- readers ----------
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
    if (ncol(df2) < 5) stop("comparison_file must have 5 cols or a header with: ", paste(req, collapse = ", "))
    names(df2)[1:5] <- req
    df2[, req, drop = FALSE]
  }

  # ---- Fast GFF reader with pure-R fallbacks (no shell) ----
  .read_gff <- function(path) {
    start_time <- proc.time()[["elapsed"]]
    ok9 <- function(x) is.data.frame(x) && ncol(x) >= 9L

    use_fread <- requireNamespace("data.table", quietly = TRUE)

    if (use_fread) {
      # We'll try three modes with comment skipping (pure R):
      # 1) tab, 2) auto, 3) single-space. All keep only columns 1:9.
      modes <- list(
        list(sep = "\t", tag = "fread(tab,comment)"),
        list(sep = "auto", tag = "fread(auto,comment)"),
        list(sep = " ", tag = "fread(space,comment)")
      )
      for (mode in modes) {
        g <- try(data.table::fread(
          file = path,
          sep = mode$sep,
          header = FALSE,
          quote = "",
          data.table = FALSE,
          showProgress = FALSE,
          fill = TRUE,
          comment = "#",          # ignore inline comments; full-line '#' become empty -> skipped
          skipEmptyLines = TRUE,
          select = 1:9,           # read only first 9 columns
          col.names = paste0("V", 1:9),
          colClasses = c("character","character","character",
                         "integer","integer","character",
                         "character","character","character"),
          na.strings = c("", "NA", ".")
        ), silent = TRUE)
        if (!inherits(g, "try-error") && ok9(g)) {
          names(g)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attributes")
          et <- proc.time()[["elapsed"]] - start_time
          vmsg(sprintf("GFF read via %s in %.3fs; rows=%s: %s",
                       mode$tag, et, format(nrow(g), big.mark=","), basename(path)))
          return(g)
        }
      }
    }

    # Base fallback(s): first tab (strict), then any whitespace (last resort).
    for (sep in c("\t", "")) {
      g <- try(utils::read.table(path, sep = sep, header = FALSE, quote = "",
                                 stringsAsFactors = FALSE, comment.char = "#", fill = TRUE,
                                 colClasses = c("character","character","character",
                                                "integer","integer","character",
                                                "character","character","character")),
               silent = TRUE)
      if (!inherits(g, "try-error") && ok9(g)) {
        names(g)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attributes")
        et <- proc.time()[["elapsed"]] - start_time
        vmsg(sprintf("GFF read via base (sep='%s') in %.3fs; rows=%s: %s",
                     if (sep=="") "any whitespace" else "\\t",
                     et, format(nrow(g), big.mark=","), basename(path)))
        if (sep == "") {
          vmsg("Warning: parsed with sep=''; if attribute values contain unescaped spaces, they may be split.")
        }
        return(g)
      }
    }

    stop("GFF appears malformed or unreadable: ", path)
  }

  ## ---------- term patterns & normalizers ----------
  .default_term_patterns <- function() {
    list(
      GO      = "GO:\\d{7}",
      IPR     = "(?:InterPro:)?IPR\\d{6}",
      KEGG    = "(?:KEGG:)?K\\d{5}|map\\d{5}",
      PANTHER = "PTHR\\d{5,6}(?:\\:SF\\d+)?",
      PFAM    = "PF\\d{5}",
      TIGRFAM = "TIGR\\d{5}",
      COG     = "COG\\d{4}"
    )
  }
  .default_normalizers <- function() {
    list(IPR = function(x) sub("^InterPro:", "", x, ignore.case = TRUE))
  }

  # Parse comma-separated customs like "REX{5},TOM{8}" -> named list: REX="REX\\d{5}", TOM="TOM\\d{8}"
  .parse_custom <- function(custom) {
    if (is.null(custom) || !nzchar(custom)) return(NULL)
    toks <- strsplit(custom, ",", fixed = TRUE)[[1]]
    toks <- trimws(toks)
    out <- list()
    for (tk in toks) {
      m <- regexec("^([A-Za-z][A-Za-z0-9_]*)\\{(\\d+)\\}$", tk)
      mm <- regmatches(tk, m)[[1]]
      if (length(mm) == 3) {
        key <- toupper(mm[2])  # label (we'll tolower for column names later)
        n   <- as.integer(mm[3])
        out[[key]] <- paste0(mm[2], "\\d{", n, "}")  # pattern: PREFIX followed by exactly n digits
      } else {
        warning("Ignoring custom token (expects PREFIX{N}): ", tk)
      }
    }
    if (!length(out)) NULL else out
  }

  ## Extract terms from ONE attributes string (key-agnostic)
  .extract_terms_from_string <- function(attr_string, patterns, normalizers, case_insensitive = TRUE) {
    out <- structure(vector("list", length(patterns)), names = names(patterns))
    if (isTRUE(is.na(attr_string)) || !nzchar(attr_string)) return(out)
    s <- attr_string
    for (lab in names(patterns)) {
      rx <- patterns[[lab]]
      hits <- unlist(regmatches(s, gregexpr(rx, s, perl = TRUE,
                                            ignore.case = case_insensitive)))
      if (length(hits)) {
        if (!is.null(normalizers[[lab]])) hits <- normalizers[[lab]](hits)
        out[[lab]] <- unique(hits[nzchar(hits)])
      }
    }
    out
  }

  ## Build graph + pre-extract OWN terms and OWN attribute strings that contain any term
  .precompute_maps_and_terms <- function(gff_df, patterns, normalizers) {
    children  <- new.env(hash = TRUE, parent = emptyenv())  # parent -> child IDs
    own_terms <- new.env(hash = TRUE, parent = emptyenv())  # id -> list(term vectors)
    own_attrs <- new.env(hash = TRUE, parent = emptyenv())  # id -> character vector of filtered attribute strings

    any_rx <- paste(sprintf("(?:%s)", unname(patterns)), collapse = "|")

    .union_term_lists <- function(a, b) {
      labs <- unique(c(names(a), names(b)))
      out <- structure(vector("list", length(labs)), names = labs)
      for (lab in labs) {
        va <- if (!is.null(a[[lab]])) a[[lab]] else character(0)
        vb <- if (!is.null(b[[lab]])) b[[lab]] else character(0)
        out[[lab]] <- unique(c(va, vb))
      }
      out
    }

    for (i in seq_len(nrow(gff_df))) {
      attr <- gff_df$attributes[i]
      pieces <- strsplit(attr, ";", fixed = TRUE)[[1]]
      kv <- strsplit(pieces, "=", fixed = TRUE)
      key <- vapply(kv, `[`, "", 1L); val <- vapply(kv, function(x) if (length(x) > 1) x[2] else "", "")
      names(val) <- key
      id  <- unname(val["ID"])
      par <- unname(val["Parent"])

      if (!is.na(id) && nzchar(id)) {
        # children edges (skip self-parent to avoid trivial cycles)
        if (!is.na(par) && nzchar(par)) {
          for (p in strsplit(par, ",", fixed = TRUE)[[1]]) {
            if (identical(p, id)) next
            cur <- if (exists(p, envir = children, inherits = FALSE)) get(p, envir = children) else character(0)
            assign(p, unique(c(cur, id)), envir = children)
          }
        }
        # own terms (union across duplicate IDs if any)
        new_terms <- .extract_terms_from_string(attr, patterns, normalizers)
        if (exists(id, envir = own_terms, inherits = FALSE)) {
          old <- get(id, envir = own_terms)
          assign(id, .union_term_lists(old, new_terms), envir = own_terms)
        } else {
          assign(id, new_terms, envir = own_terms)
        }
        # own attributes: keep ONLY rows that contain any recognized term
        if (grepl(any_rx, attr, perl = TRUE, ignore.case = TRUE)) {
          if (exists(id, envir = own_attrs, inherits = FALSE)) {
            assign(id, unique(c(get(id, envir = own_attrs), attr)), envir = own_attrs)
          } else {
            assign(id, attr, envir = own_attrs)
          }
        } else if (!exists(id, envir = own_attrs, inherits = FALSE)) {
          assign(id, character(0), envir = own_attrs)
        }
      }
    }
    list(children = children, own_terms = own_terms, own_attrs = own_attrs)
  }

  .union_lists <- function(lst) {
    labs <- unique(unlist(lapply(lst, names)))
    out <- structure(vector("list", length(labs)), names = labs)
    for (lab in labs) {
      vals <- unlist(lapply(lst, function(x) if (!is.null(x[[lab]])) x[[lab]] else character(0)), use.names = FALSE)
      out[[lab]] <- unique(vals[nzchar(vals)])
    }
    out
  }

  ## Memoized subtree unions WITH CYCLE GUARDS
  .subtree_terms <- local({
    cache <- new.env(hash = TRUE, parent = emptyenv())
    INPROG <- new.env(hash = TRUE, parent = emptyenv())

    function(id, children_env, own_terms_env) {
      if (is.null(id) || !nzchar(id)) return(structure(list(), names = character(0)))
      if (exists(id, envir = cache, inherits = FALSE)) return(get(id, envir = cache))
      if (exists(id, envir = INPROG, inherits = FALSE)) return(structure(list(), names = character(0)))

      assign(id, TRUE, envir = INPROG)
      on.exit({ if (exists(id, envir = INPROG, inherits = FALSE)) rm(list = id, envir = INPROG) }, add = TRUE)

      own  <- if (exists(id, envir = own_terms_env, inherits = FALSE)) get(id, envir = own_terms_env) else structure(list(), names = character(0))
      kids <- if (exists(id, envir = children_env,  inherits = FALSE)) get(id, envir = children_env)  else character(0)

      res <- if (length(kids)) .union_lists(c(list(own), lapply(kids, .subtree_terms, children_env = children_env, own_terms_env = own_terms_env))) else own
      assign(id, res, envir = cache); res
    }
  })

  .subtree_attrs <- local({
    cache <- new.env(hash = TRUE, parent = emptyenv())
    INPROG <- new.env(hash = TRUE, parent = emptyenv())

    function(id, children_env, own_attrs_env) {
      if (is.null(id) || !nzchar(id)) return(character(0))
      if (exists(id, envir = cache, inherits = FALSE)) return(get(id, envir = cache))
      if (exists(id, envir = INPROG, inherits = FALSE)) return(character(0))

      assign(id, TRUE, envir = INPROG)
      on.exit({ if (exists(id, envir = INPROG, inherits = FALSE)) rm(list = id, envir = INPROG) }, add = TRUE)

      own  <- if (exists(id, envir = own_attrs_env, inherits = FALSE)) get(id, envir = own_attrs_env) else character(0)
      kids <- if (exists(id, envir = children_env,  inherits = FALSE)) get(id, envir = children_env)  else character(0)

      res <- if (length(kids)) unique(c(own, unlist(lapply(kids, .subtree_attrs, children_env = children_env, own_attrs_env = own_attrs_env), use.names = FALSE))
) else own
      assign(id, res, envir = cache); res
    }
  })

  ## Annotate a set of IDs -> data.frame(ID, ATTR, <terms...>)
  .annotate_id_set <- function(ids, children_env, own_terms_env, own_attrs_env, labs, collapse_sep = ";", threads = 1L) {
    to_row <- function(id) {
      terms <- .subtree_terms(id, children_env, own_terms_env)
      attrs <- .subtree_attrs(id, children_env, own_attrs_env)
      row <- c(
        ID   = id,
        ATTR = if (length(attrs)) paste(unique(attrs), collapse = collapse_sep) else NA_character_
      )
      for (L in labs) {
        v <- terms[[L]]
        row[[L]] <- if (!is.null(v) && length(v)) paste(v, collapse = collapse_sep) else NA_character_
      }
      row
    }
    rows <- plapply(ids, to_row, threads = threads)
    df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    rownames(df) <- NULL
    df
  }

  ## ---------- single comparison ----------
  .annotate_single <- function(dnds_path, q_gff, s_gff, out_path = NULL, threads = 1L) {
    if (!file.exists(dnds_path)) stop("dN/dS file not found: ", dnds_path)
    if (!file.exists(q_gff))     stop("query_gff not found: ", q_gff)
    if (!file.exists(s_gff))     stop("subject_gff not found: ", s_gff)

    vmsg(sprintf("Reading dN/dS: %s", dnds_path))
    dnds <- utils::read.table(dnds_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")

    if (!all(c("query_id","subject_id") %in% names(dnds))) {
      stop("dN/dS file must have columns 'query_id' and 'subject_id'.")
    }

    patterns    <- .default_term_patterns()
    customs     <- .parse_custom(custom)
    if (!is.null(customs)) patterns[names(customs)] <- unname(customs)
    normalizers <- .default_normalizers()

    std_order <- c("IPR","GO","KEGG","PANTHER","PFAM","TIGRFAM","COG")
    labs <- c(std_order[std_order %in% names(patterns)],
              setdiff(names(patterns), std_order))

    vmsg("Precomputing query-side maps…")
    q_maps <- .precompute_maps_and_terms(.read_gff(q_gff), patterns, normalizers)
    vmsg("Precomputing subject-side maps…")
    s_maps <- .precompute_maps_and_terms(.read_gff(s_gff), patterns, normalizers)

    q_ids <- unique(dnds$query_id)
    s_ids <- unique(dnds$subject_id)
    vmsg(sprintf("Annotating IDs (q=%d, s=%d) with threads=%d…", length(q_ids), length(s_ids), threads))

    # Split threads across sides to avoid oversubscription
    tq <- max(1L, threads %/% 2L)
    ts <- max(1L, threads - tq)

    if (.Platform$OS.type != "windows" && threads > 1L) {
      j1 <- parallel::mcparallel(.annotate_id_set(q_ids, q_maps$children, q_maps$own_terms, q_maps$own_attrs, labs, threads = tq), silent = TRUE)
      j2 <- parallel::mcparallel(.annotate_id_set(s_ids, s_maps$children, s_maps$own_terms, s_maps$own_attrs, labs, threads = ts), silent = TRUE)
      res <- parallel::mccollect(list(j1, j2), wait = TRUE)
      q_ann <- res[[1]]; s_ann <- res[[2]]
    } else {
      q_ann <- .annotate_id_set(q_ids, q_maps$children, q_maps$own_terms, q_maps$own_attrs, labs, threads = tq)
      s_ann <- .annotate_id_set(s_ids, s_maps$children, s_maps$own_terms, s_maps$own_attrs, labs, threads = ts)
    }

    present_labs <- labs[vapply(labs, function(L) {
      qa <- q_ann[[L]]; sa <- s_ann[[L]]
      any(!is.na(qa) & qa != "") || any(!is.na(sa) & sa != "")
    }, logical(1))]

    d <- dnds
    mapv <- function(df, col) { v <- df[[col]]; names(v) <- df$ID; v }

    d$q_attributes <- mapv(q_ann, "ATTR")[match(d$query_id,   q_ann$ID)]
    d$s_attributes <- mapv(s_ann, "ATTR")[match(d$subject_id, s_ann$ID)]

    for (L in present_labs) {
      lname <- tolower(L)
      d[[paste0("q_", lname)]] <- mapv(q_ann, L)[match(d$query_id,   q_ann$ID)]
      d[[paste0("s_", lname)]] <- mapv(s_ann, L)[match(d$subject_id, s_ann$ID)]
    }

    if (!is.null(out_path)) {
      vmsg(sprintf("Writing: %s", out_path))
      utils::write.table(d, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
      return(out_path)
    }
    d
  }

  ## ---------- batch vs single ----------
  threads <- .resolve_threads(threads)
  if (threads < 1L || is.na(threads)) threads <- 1L

  if (!is.null(comparison_file)) {
    df <- .read_comparisons(comparison_file)
    vmsg(sprintf("Batch mode: %d comparison(s). Using threads=%d per comparison.", nrow(df), threads))
    out_paths <- character(0)
    for (i in seq_len(nrow(df))) {
      comp      <- df$comparison_name[i]
      comp_dir  <- file.path(output_dir, comp)
      dnds_path <- file.path(comp_dir, paste0(comp, "_dnds.tsv"))
      out_path  <- file.path(comp_dir, paste0(comp, "_dnds_annot.tsv"))

      if (!isTRUE(overwrite) && file.exists(out_path) && file.info(out_path)$size > 0) {
        vmsg(sprintf("[%s] exists; skipping (use overwrite=TRUE to regenerate).", basename(out_path)))
        out_paths <- c(out_paths, out_path)
        next
      }

      vmsg(sprintf("[%s] reading GFFs & appending annotations…", comp))
      out <- .annotate_single(dnds_path, df$query_gff[i], df$subject_gff[i], out_path, threads = threads)
      vmsg(sprintf("[%s] wrote %s", comp, out))
      out_paths <- c(out_paths, out)
    }
    message("All annotation appends complete.")
    return(invisible(out_paths))
  }

  if (is.null(dnds_file) || is.null(query_gff) || is.null(subject_gff)) {
    stop("Provide either comparison_file (batch) OR dnds_file + query_gff + subject_gff (single).")
  }
  .annotate_single(dnds_file, query_gff, subject_gff, output_file, threads = threads)
}
