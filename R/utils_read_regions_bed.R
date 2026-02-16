#' Read and normalize a regions BED-like file
#'
#' Returns a data.frame with columns: seqname, start, end, region_name
#'
#' Robust to headerless BED files: tries header=TRUE, then falls back to header=FALSE.
#'
#' Coordinate conventions:
#'   - regions_coord = "bed0": interpret as BED 0-based, half-open [`start`, `end`),
#'     convert to 1-based, closed [`start + 1`, `end`] for overlap with GFF-like coordinates.
#'   - regions_coord = "gff1": interpret as 1-based, closed [`start`, `end`].
#'
#' Column selectors:
#'   - region_*_col may be NULL (default), a column *name*, or a 1-based *integer index*.
#'     For headerless input, indices (1,2,3,4) are often easiest.
#'
#' @keywords internal
.read_regions_bed <- function(regions_bed,
                              region_seq_col   = NULL,
                              region_start_col = NULL,
                              region_end_col   = NULL,
                              region_name_col  = NULL,
                              regions_coord    = c("bed0", "gff1"),
                              stop_on_invalid  = TRUE,
                              ...) {

  regions_coord <- match.arg(regions_coord)

  if (missing(regions_bed) || is.null(regions_bed) || !nzchar(regions_bed)) {
    stop("regions_bed is required.")
  }
  if (!file.exists(regions_bed)) stop("regions_bed not found: ", regions_bed)

  reg_raw <- try(.read_ws(regions_bed, header_try = TRUE), silent = TRUE)
  if (inherits(reg_raw, "try-error") || ncol(reg_raw) < 3) {
    reg_raw <- .read_ws(regions_bed, header_try = FALSE)
  }
  if (ncol(reg_raw) < 3) stop("regions_bed must have at least 3 columns: seqname, start, end.")

  reg_raw <- .clean_colnames(reg_raw)

  # Resolve a selector that can be NULL, integer index, or column name
  .resolve_col <- function(sel, default_idx) {
    if (is.null(sel)) return(default_idx)

    # numeric index (1-based)
    if (is.numeric(sel) && length(sel) == 1L && is.finite(sel)) {
      idx <- as.integer(sel)
      if (idx >= 1L && idx <= ncol(reg_raw)) return(idx)
      stop("Invalid regions_bed column index: ", sel, " (ncol=", ncol(reg_raw), ").")
    }

    # column name
    if (is.character(sel) && length(sel) == 1L) {
      if (sel %in% names(reg_raw)) return(match(sel, names(reg_raw)))
      stop("Invalid regions_bed column name: '", sel, "'. Available: ",
           paste(names(reg_raw), collapse = ", "))
    }

    stop("Invalid regions_bed column selector (must be NULL, integer index, or column name).")
  }

  # If headerless, read.table will name columns V1, V2, ...; that’s fine.
  seq_idx   <- .resolve_col(region_seq_col,   1L)
  start_idx <- .resolve_col(region_start_col, 2L)
  end_idx   <- .resolve_col(region_end_col,   3L)

  # region_name is optional: default to col4 if present, else "region"
  name_idx <- NULL
  if (!is.null(region_name_col)) {
    name_idx <- .resolve_col(region_name_col, 4L)
  } else if (ncol(reg_raw) >= 4L) {
    name_idx <- 4L
  }

  seqname <- as.character(reg_raw[[seq_idx]])
  start0  <- suppressWarnings(as.numeric(reg_raw[[start_idx]]))
  end0    <- suppressWarnings(as.numeric(reg_raw[[end_idx]]))
  rname   <- if (!is.null(name_idx)) as.character(reg_raw[[name_idx]]) else rep("region", length(seqname))

  bad <- is.na(seqname) | !nzchar(seqname) | is.na(start0) | is.na(end0)
  if (any(bad)) {
    msg <- paste0(
      "regions_bed has invalid rows (missing/NA seqname or non-numeric start/end). ",
      "First few bad row indices: ",
      paste(utils::head(which(bad), 10), collapse = ", ")
    )
    if (stop_on_invalid) stop(msg) else warning(msg)
  }

  # Convert coordinate convention to 1-based closed intervals for downstream overlap
  if (regions_coord == "bed0") {
    start <- start0 + 1
    end   <- end0
  } else {
    start <- start0
    end   <- end0
  }

  bad2 <- !is.finite(start) | !is.finite(end) | start > end
  if (any(bad2)) {
    msg <- paste0(
      "regions_bed has invalid intervals (start > end or non-finite). ",
      "First few bad row indices: ",
      paste(utils::head(which(bad2), 10), collapse = ", ")
    )
    if (stop_on_invalid) stop(msg) else warning(msg)
  }

  # region_name: fill blanks
  if (any(is.na(rname) | !nzchar(rname))) {
    rname[is.na(rname) | !nzchar(rname)] <- "region"
  }

  data.frame(
    seqname     = seqname,
    start       = start,
    end         = end,
    region_name = rname,
    stringsAsFactors = FALSE
  )
}
