# PAML codeml as an engine for branch models: discovery, inputs, runs, and parsing.
# Analysis-level code is in phylo_branch_omega.R.

#' PAML release dndsR was validated with
#' @keywords internal
.paml_validated_version <- "4.10.10"

#' Find codeml
#' @keywords internal
.codeml_path <- function(codeml = NULL) {
  if (!is.null(codeml)) {
    if (!file.exists(codeml)) stop("codeml not found at ", codeml, call. = FALSE)
    return(normalizePath(codeml))
  }
  candidates <- c(getOption("dndsR.codeml"), Sys.getenv("DNDSR_CODEML", ""), Sys.which("codeml"))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (p in candidates) if (file.exists(p)) return(normalizePath(p))
  stop("PAML's codeml was not found. Install it (conda install -c conda-forge -c bioconda paml=", .paml_validated_version,
       "), put 'codeml' on PATH, or give its path with --codeml or DNDSR_CODEML.", call. = FALSE)
}

#' Standard input for codeml: closed, because codeml pauses for Enter on some warnings
#' @keywords internal
.paml_stdin <- function() if (.Platform$OS.type == "windows") "" else "/dev/null"

#' codeml's version, such as "4.10.10"
#' @keywords internal
.codeml_version <- function(bin) {
  d <- tempfile("codeml_version")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  writeLines(c("seqfile = absent.phy", "treefile = absent.nwk", "outfile = absent.mlc", "noisy = 0"), file.path(d, "v.ctl"))
  old <- setwd(d)
  on.exit(setwd(old), add = TRUE)
  out <- tryCatch(suppressWarnings(system2(bin, "v.ctl", stdout = TRUE, stderr = TRUE, stdin = .paml_stdin(), timeout = 30)),
                  error = function(e) character())
  v <- regmatches(out, regexpr("paml version [0-9]+\\.[0-9]+(\\.[0-9]+)?", out))
  if (!length(v)) stop("Could not read the PAML version from codeml at '", bin, "'.", call. = FALSE)
  sub("paml version ", "", v[1])
}

#' Write an alignment in PAML's sequential PHYLIP format, named by taxon
#' @keywords internal
.paml_write_phylip <- function(aln, path) {
  writeLines(c(sprintf(" %d %d", nrow(aln), nchar(aln$sequence[1])), paste0(aln$taxon, "  ", aln$sequence)), path)
}

#' Newick for codeml: tip names only, no internal names, optional `#1` marks
#'
#' @param tree Unrooted `ape` phylo (a trifurcation at the root, as codeml expects).
#' @param marks Logical per node: mark the branch above it `#1`.
#' @return A tree file's lines: the PAML header (`ns 1`) and the Newick string.
#' @keywords internal
.paml_tree_lines <- function(tree, marks = NULL) {
  n_tip <- ape::Ntip(tree)
  if (is.null(marks)) marks <- logical(n_tip + tree$Nnode)
  root <- n_tip + 1L
  kids <- split(tree$edge[, 2], factor(tree$edge[, 1], levels = seq_len(n_tip + tree$Nnode)))
  write_node <- function(v) {
    tag <- if (marks[v]) " #1" else ""
    if (v <= n_tip) return(paste0(tree$tip.label[v], tag))
    inner <- paste(vapply(kids[[v]], write_node, character(1)), collapse = ",")
    if (v == root) paste0("(", inner, ");") else paste0("(", inner, ")", tag)
  }
  c(sprintf(" %d 1", n_tip), write_node(root))
}

#' codeml control file for one model and one start
#'
#' Branch models only (`NSsites = 0`). `model`: 0 one omega, 1 free ratios, 2 two or more
#' ratios with `#1`-marked branches. `fix_omega = 1` with `model = 2` fixes the omega of
#' the marked branches at `omega`.
#' @keywords internal
.paml_ctl <- function(model, omega, fix_omega, kappa, codon_freq, seqfile = "alignment.phy", treefile = "tree.nwk") {
  c(paste("seqfile =", seqfile), paste("treefile =", treefile), "outfile = codeml.mlc",
    "noisy = 0", "verbose = 1", "runmode = 0", "seqtype = 1",
    sprintf("CodonFreq = %d", codon_freq), "clock = 0", "aaDist = 0",
    sprintf("model = %d", model), "NSsites = 0", "icode = 0",
    "fix_kappa = 0", sprintf("kappa = %s", format(kappa)),
    sprintf("fix_omega = %d", fix_omega), sprintf("omega = %s", format(omega)),
    "fix_alpha = 1", "alpha = 0", "Malpha = 0", "ncatG = 1",
    "getSE = 0", "RateAncestor = 0", "Small_Diff = .5e-6", "cleandata = 0", "method = 0")
}

#' codeml's auxiliary output files, which dndsR does not read
#' @keywords internal
.paml_aux_files <- c("rst", "rst1", "rub", "lnf", "2NG.t", "2NG.dN", "2NG.dS", "4fold.nuc", "rates", "codeml.log")

#' Run codeml in its own folder
#' @return A list with `status`, `timed_out`, and `log` (path).
#' @keywords internal
.codeml_run <- function(bin, dir, timeout) {
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  status <- suppressWarnings(system2(bin, "codeml.ctl", stdout = "codeml.log", stderr = "codeml.log",
                                     stdin = .paml_stdin(), timeout = timeout))
  list(status = status, timed_out = identical(as.integer(status), 124L), log = file.path(dir, "codeml.log"))
}

#' Parse codeml's main output for a branch model
#'
#' @param path `codeml.mlc`.
#' @param tips Tip names in the order of the sequence file, which codeml numbers 1, 2, ...
#' @return A list: `lnl`, `np`, `kappa`, `omega` (one value, or background and marked for
#'   model 2), `tree_length`, `warnings`, and `branches` (one row per branch: `tips` (the
#'   tip names below it), `t`, `N`, `S`, `omega`, `dN`, `dS`, `NdN`, `SdS`).
#' @keywords internal
.codeml_parse <- function(path, tips) {
  if (!file.exists(path)) stop("codeml wrote no output")
  x <- readLines(path, warn = FALSE)
  num <- function(pattern) {
    hit <- grep(pattern, x, value = TRUE)
    if (!length(hit)) return(NA_real_)
    as.numeric(regmatches(hit[length(hit)], regexpr("-?[0-9]+\\.?[0-9]*([eE][-+]?[0-9]+)?$", trimws(hit[length(hit)]))))
  }
  lnl_line <- grep("^lnL\\(ntime:", x, value = TRUE)
  if (!length(lnl_line)) stop("codeml output has no likelihood")
  lnl_line <- lnl_line[length(lnl_line)]
  lnl <- as.numeric(sub("^.*\\):\\s*(-?[0-9.]+).*$", "\\1", lnl_line))
  np <- as.integer(sub("^.*np:\\s*([0-9]+)\\).*$", "\\1", lnl_line))
  omega <- num("^omega \\(dN/dS\\) =")
  two <- grep("^w \\(dN/dS\\) for branches:", x, value = TRUE)
  if (length(two)) omega <- as.numeric(strsplit(trimws(sub("^w \\(dN/dS\\) for branches:", "", two[length(two)])), "\\s+")[[1]])
  head <- grep("^\\s*branch\\s+t\\s+N\\s+S\\s+dN/dS\\s+dN\\s+dS\\s+N\\*dN\\s+S\\*dS", x)
  if (!length(head)) stop("codeml output has no branch table")
  rows <- character()
  for (i in seq(head[length(head)] + 1L, length(x))) {
    if (!nzchar(trimws(x[i]))) { if (length(rows)) break else next }
    rows <- c(rows, trimws(x[i]))
  }
  f <- do.call(rbind, strsplit(rows, "\\s+"))
  ends <- do.call(rbind, strsplit(f[, 1], "..", fixed = TRUE))
  parent <- as.integer(ends[, 1]); child <- as.integer(ends[, 2])
  n <- length(tips)
  below <- function(v) if (v <= n) tips[v] else unlist(lapply(child[parent == v], below))
  vals <- apply(f[, -1, drop = FALSE], 2, as.numeric)
  if (is.null(dim(vals))) vals <- matrix(vals, nrow = 1)
  branches <- data.frame(parent = parent, child = child, t = vals[, 1], N = vals[, 2], S = vals[, 3], omega = vals[, 4],
                         dN = vals[, 5], dS = vals[, 6], NdN = vals[, 7], SdS = vals[, 8], stringsAsFactors = FALSE)
  branches$tips <- I(lapply(child, below))
  warn <- grep("check convergence|[Ee]rror|not converge", x, value = TRUE)
  list(lnl = lnl, np = np, kappa = num("^kappa \\(ts/tv\\) ="), omega = omega, tree_length = num("^tree length ="),
       warnings = trimws(warn), branches = branches)
}
