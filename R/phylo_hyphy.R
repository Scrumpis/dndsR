# HyPhy engine layer. Method-agnostic: finding the binary, running one analysis in its own
# folder with an explicit seed and thread count, and reading its JSON. Method-specific
# arguments and parsing live in the registry `.hyphy_methods`.

# HyPhy release this code was validated against. Others run, with a notice.
.hyphy_validated_version <- "2.5.103"

#' Find the HyPhy executable
#'
#' Order: the `hyphy` argument, `options(dndsR.hyphy)`, the `DNDSR_HYPHY` environment
#' variable, then `hyphy` on `PATH`.
#' @keywords internal
.hyphy_path <- function(hyphy = NULL) {
  if (!is.null(hyphy)) {
    if (!file.exists(hyphy)) stop("HyPhy not found at ", hyphy, call. = FALSE)
    return(normalizePath(hyphy))
  }
  candidates <- c(hyphy, getOption("dndsR.hyphy"), Sys.getenv("DNDSR_HYPHY", ""), Sys.which("hyphy"))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (p in candidates) if (file.exists(p)) return(normalizePath(p))
  stop("HyPhy was not found. Install it (conda install -c conda-forge -c bioconda hyphy=", .hyphy_validated_version,
       "), put 'hyphy' on PATH, or give its path with --hyphy or DNDSR_HYPHY. Published dndsR images up to v0.4.0 do not include it.",
       call. = FALSE)
}

#' HyPhy version, such as "2.5.103"
#' @keywords internal
.hyphy_version <- function(bin) {
  out <- tryCatch(suppressWarnings(system2(bin, "--version", stdout = TRUE, stderr = TRUE)), error = function(e) character())
  v <- regmatches(out, regexpr("[0-9]+\\.[0-9]+\\.[0-9]+", out))
  if (!length(v)) stop("Could not read the HyPhy version from '", bin, " --version'.", call. = FALSE)
  v[1]
}

#' Run one HyPhy analysis in its own folder
#'
#' HyPhy writes `errors.log` into its working folder, so each run gets one. The random seed
#' and thread count are always set: in testing on HyPhy 2.5.103, a fixed seed and thread
#' count gave bit-identical results, and anything else did not.
#'
#' @param bin HyPhy executable.
#' @param method HyPhy analysis name, such as "busted".
#' @param args Character vector of keyword arguments, already quoted where needed.
#' @param workdir Folder for the run.
#' @param cpu Threads for HyPhy.
#' @param seed Random seed.
#' @return A list with `status` (exit code), `log` (path), `error` (text, or `NA`), and
#'   `command` (as run).
#' @keywords internal
.hyphy_run <- function(bin, method, args, workdir, cpu, seed) {
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(workdir, "hyphy.log")
  err_file <- file.path(workdir, "errors.log")
  if (file.exists(err_file)) unlink(err_file)
  full <- c(sprintf("CPU=%d", as.integer(cpu)), sprintf("ENV=RANDOM_SEED=%d", as.integer(seed)), method, args)
  old <- setwd(workdir)
  on.exit(setwd(old), add = TRUE)
  status <- suppressWarnings(system2(bin, full, stdout = "hyphy.log", stderr = "hyphy.log"))
  error <- NA_character_
  if (file.exists("errors.log")) {
    e <- trimws(readLines("errors.log", warn = FALSE))
    e <- sub("^Error:", "", e[nzchar(e)])
    error <- substr(paste(utils::head(e, 2), collapse = " "), 1, 300)
  }
  list(status = status, log = log, error = error, command = paste(c(shQuote(bin), full), collapse = " "))
}

#' Read HyPhy JSON and check the keys a method needs
#' @keywords internal
.hyphy_read_json <- function(path, required) {
  if (!file.exists(path)) stop("HyPhy wrote no JSON")
  j <- tryCatch(jsonlite::read_json(path, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(j)) stop("HyPhy JSON could not be parsed")
  missing <- setdiff(required, names(j))
  if (length(missing)) stop("HyPhy JSON lacks: ", paste(missing, collapse = ", "))
  j
}

#' Rate classes from a HyPhy "Rate Distributions" entry
#' @keywords internal
.hyphy_rate_classes <- function(x, prefix) {
  if (is.null(x)) return(list())
  x <- x[order(as.integer(names(x)))]
  out <- list()
  for (k in seq_along(x)) {
    out[[sprintf("%s_omega_%d", prefix, k)]] <- as.numeric(x[[k]]$omega)
    out[[sprintf("%s_prop_%d", prefix, k)]] <- as.numeric(x[[k]]$proportion)
  }
  out
}

#' Registry of supported HyPhy methods
#'
#' Each entry gives the HyPhy analysis name, its settings with defaults equal to HyPhy's,
#' a function building keyword arguments, the JSON keys required, a parser returning a
#' one-row list of results, and the unit of testing ("gene" or "branch").
#' @keywords internal
.hyphy_methods <- list(
  busted = list(
    hyphy_name = "busted",
    unit = "gene",
    defaults = list(srv = "Yes", multiple_hits = "None", rates = 3L, syn_rates = 3L, error_sink = "No",
                    grid_size = 250L, starting_points = 1L, kill_zero_lengths = "Yes", genetic_code = "Universal"),
    choices = list(srv = c("Yes", "No", "Branch-site"), multiple_hits = c("None", "Double", "Double+Triple"),
                   error_sink = c("No", "Yes"), kill_zero_lengths = c("Yes", "No"), genetic_code = "Universal"),
    args = function(s, alignment, tree, output) {
      a <- c("--alignment", shQuote(alignment), "--tree", shQuote(tree), "--branches", "FG",
             "--code", s$genetic_code, "--srv", s$srv, "--rates", s$rates)
      if (s$srv != "No") a <- c(a, "--syn-rates", s$syn_rates)
      c(a, "--multiple-hits", shQuote(s$multiple_hits), "--error-sink", s$error_sink,
        "--grid-size", s$grid_size, "--starting-points", s$starting_points,
        "--kill-zero-lengths", s$kill_zero_lengths, "--output", shQuote(output))
    },
    required = c("test results", "fits", "tested", "branch attributes", "input"),
    parse = function(j) {
      tr <- j[["test results"]]
      unc <- j$fits[["Unconstrained model"]]
      con <- j$fits[["Constrained model"]]
      if (is.null(tr$LRT) || is.null(tr[["p-value"]]) || is.null(unc)) stop("HyPhy JSON lacks the BUSTED test")
      ba <- j[["branch attributes"]][["0"]]
      lengths <- vapply(ba, function(b) if (is.null(b$unconstrained)) NA_real_ else as.numeric(b$unconstrained), numeric(1))
      c(list(
        lrt = as.numeric(tr$LRT),
        p_value = as.numeric(tr[["p-value"]]),
        constrained_fitted = !is.null(con),
        lnl_unconstrained = as.numeric(unc[["Log Likelihood"]]),
        lnl_constrained = if (is.null(con)) NA_real_ else as.numeric(con[["Log Likelihood"]]),
        tree_length = sum(lengths, na.rm = TRUE)
      ), .hyphy_rate_classes(unc[["Rate Distributions"]]$Test, "fg"),
      .hyphy_rate_classes(unc[["Rate Distributions"]]$Background, "bg"))
    }
  ),
  absrel = list(
    hyphy_name = "absrel",
    unit = "branch",
    defaults = list(srv = "No", multiple_hits = "None", genetic_code = "Universal"),
    choices = list(srv = c("Yes", "No"), multiple_hits = c("None", "Double", "Double+Triple"), genetic_code = "Universal"),
    args = function(s, alignment, tree, output) {
      c("--alignment", shQuote(alignment), "--tree", shQuote(tree), "--branches", "FG", "--code", s$genetic_code,
        "--srv", s$srv, "--multiple-hits", shQuote(s$multiple_hits), "--output", shQuote(output))
    },
    required = c("test results", "fits", "tested", "branch attributes", "input"),
    parse = function(j) {
      br <- .absrel_branches(j)
      t <- br[br$tested_by_hyphy & !is.na(br$p_holm), , drop = FALSE]
      if (!nrow(t)) stop("HyPhy JSON has no aBSREL test")
      list(
        n_branches_tested = nrow(t),
        n_branches_holm_005 = sum(t$p_holm <= 0.05),
        # HyPhy's p_holm omits Holm's running maximum, so its smallest value is Hochberg's global
        # test that "some tested branch evolved under episodic positive selection" (valid when the
        # branches' p-values are independent or positively dependent).
        p_value = min(t$p_holm),
        branch_min_p = t$branch[which.min(t$p_holm)],
        lnl_baseline = .json_num(j$fits[["Baseline MG94xREV"]][["Log Likelihood"]]),
        lnl_adaptive = .json_num(j$fits[["Full adaptive model"]][["Log Likelihood"]])
      )
    },
    parse_branches = function(j) .absrel_branches(j),
    map = function(hyp, present, gene_tree, mode, ctx) {
      if (isTRUE(ctx$all_branches)) .phylo_map_all_tested(hyp, present, gene_tree, mode) else .phylo_map_foreground(hyp, present, gene_tree, mode)
    },
    adjust_branches = function(branches, results) .absrel_adjust(branches, results),
    testing = list(
      within_gene = "p_holm: HyPhy's correction across the branches tested in each orthogroup (k-th smallest p times m - k + 1, without Holm's running maximum); p_value (per orthogroup) is its smallest value, and p_adj applies BH to p_value across orthogroups",
      per_lineage = "p_bh_lineage: BH over p_raw, across the tested orthogroups in which that species-tree branch was tested (one family per species-tree branch)"
    )
  ),
  relax = list(
    hyphy_name = "relax",
    unit = "gene",
    defaults = list(models = "All", rates = 3L, srv = "No", multiple_hits = "None", grid_size = 250L,
                    starting_points = 1L, genetic_code = "Universal"),
    choices = list(models = c("All", "Minimal"), srv = c("Yes", "No"), multiple_hits = c("None", "Double", "Double+Triple"),
                   genetic_code = "Universal"),
    args = function(s, alignment, tree, output) {
      c("--alignment", shQuote(alignment), "--tree", shQuote(tree), "--test", "TEST", "--reference", "REFERENCE",
        "--code", s$genetic_code, "--models", s$models, "--rates", s$rates, "--srv", s$srv,
        "--multiple-hits", shQuote(s$multiple_hits), "--grid-size", s$grid_size, "--starting-points", s$starting_points,
        "--output", shQuote(output))
    },
    required = c("test results", "fits", "tested", "input"),
    parse = function(j) {
      tr <- j[["test results"]]
      k <- .json_num(tr[["relaxation or intensification parameter"]])
      if (is.na(k) || is.null(tr[["p-value"]])) stop("HyPhy JSON lacks the RELAX test")
      alt <- j$fits[["RELAX alternative"]]
      sep <- j$fits[["MG94xREV with separate rates for branch sets"]][["Rate Distributions"]]
      sep_omega <- function(set) {
        x <- sep[[sprintf("non-synonymous/synonymous rate ratio for *%s*", set)]]
        if (is.null(x)) NA_real_ else as.numeric(x[[1]][[1]])
      }
      c(list(
        k = k,
        direction = if (k > 1) "intensification" else if (k < 1) "relaxation" else "none",
        # HyPhy bounds k to [0, 50]; a value at a bound is an extreme fit, not a measured size.
        k_at_bound = k <= 1e-6 || k >= 49.99,
        lrt = .json_num(tr$LRT),
        p_value = .json_num(tr[["p-value"]]),
        lnl_alternative = .json_num(alt[["Log Likelihood"]]),
        lnl_null = .json_num(j$fits[["RELAX null"]][["Log Likelihood"]]),
        mg94_omega_test = sep_omega("Test"),
        mg94_omega_reference = sep_omega("Reference")
      ), .hyphy_rate_classes(alt[["Rate Distributions"]]$Test, "test"),
      .hyphy_rate_classes(alt[["Rate Distributions"]]$Reference, "ref"))
    },
    map = function(hyp, present, gene_tree, mode, ctx) .phylo_map_relax(hyp, ctx$hyp_ref, present, gene_tree, mode),
    testing = list(note = "RELAX tests whether omega classes on test branches are scaled by k relative to reference branches: k < 1 relaxation, k > 1 intensification. It is not a test of positive selection.")
  )
)

#' Number from a JSON value, `NA` when absent
#' @keywords internal
.json_num <- function(x) if (is.null(x) || !length(x)) NA_real_ else suppressWarnings(as.numeric(x[[1]]))

#' Per-branch results of aBSREL
#'
#' Every branch has the baseline MG94xREV fit (one omega per branch, no site variation)
#' and the adaptive model's rate classes; tested branches also have the likelihood ratio,
#' HyPhy's uncorrected p-value, and HyPhy's corrected p-value (Holm's step without the running maximum) across the branches
#' tested in this gene.
#' @keywords internal
.absrel_branches <- function(j) {
  ba <- j[["branch attributes"]][["0"]]
  tested <- j$tested[["0"]]
  rows <- lapply(names(ba), function(b) {
    x <- ba[[b]]
    rd <- x[["Rate Distributions"]]
    om <- if (length(rd)) vapply(rd, function(r) as.numeric(r[[1]]), numeric(1)) else numeric()
    pr <- if (length(rd)) vapply(rd, function(r) as.numeric(r[[2]]), numeric(1)) else numeric()
    data.frame(
      branch = b,
      tested_by_hyphy = identical(tolower(as.character(tested[[b]])), "test"),
      baseline_omega = .json_num(x[["Baseline MG94xREV omega ratio"]]),
      baseline_length = .json_num(x[["Baseline MG94xREV"]]),
      adaptive_length = .json_num(x[["Full adaptive model"]]),
      adaptive_dn = .json_num(x[["Full adaptive model (non-synonymous subs/site)"]]),
      adaptive_ds = .json_num(x[["Full adaptive model (synonymous subs/site)"]]),
      rate_classes = .json_num(x[["Rate classes"]]),
      omega_max = if (length(om)) max(om) else NA_real_,
      omega_max_prop = if (length(om)) pr[which.max(om)] else NA_real_,
      lrt = .json_num(x[["LRT"]]),
      p_raw = .json_num(x[["Uncorrected P-value"]]),
      p_holm = .json_num(x[["Corrected P-value"]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Benjamini-Hochberg per species-tree branch, across orthogroups (aBSREL)
#'
#' One family per species-tree branch (`species_tree_branches`): the tested orthogroups in
#' which that branch was tested. Branches that exist only in a gene tree, or that merge
#' several species-tree branches after pruning, form families of their own label.
#' @keywords internal
.absrel_adjust <- function(branches, results) {
  tested_og <- results$orthogroup[results$status == "tested"]
  ok <- branches$orthogroup %in% tested_og & !is.na(branches$p_raw) & !is.na(branches$species_tree_branches)
  branches$p_bh_lineage <- NA_real_
  branches$lineage_family_size <- NA_integer_
  for (g in unique(branches$species_tree_branches[ok])) {
    i <- which(ok & branches$species_tree_branches == g)
    branches$p_bh_lineage[i] <- stats::p.adjust(branches$p_raw[i], method = "BH")
    branches$lineage_family_size[i] <- length(i)
  }
  branches
}

#' Branches HyPhy reports as tested
#' @keywords internal
.hyphy_tested_branches <- function(j) {
  t <- j$tested[["0"]]
  sort(names(t)[vapply(t, function(x) identical(tolower(x), "test"), logical(1))], method = "radix")
}

#' Resolve method settings against the registry
#' @keywords internal
.hyphy_settings <- function(method, given) {
  m <- .hyphy_methods[[method]]
  s <- m$defaults
  for (nm in names(given)) if (!is.null(given[[nm]])) s[[nm]] <- given[[nm]]
  for (nm in names(m$choices)) {
    v <- as.character(s[[nm]])
    hit <- m$choices[[nm]][tolower(m$choices[[nm]]) == tolower(v)]
    if (length(hit) != 1L) {
      stop(sprintf("Invalid %s '%s'. Use one of: %s.", nm, v, paste(m$choices[[nm]], collapse = ", ")), call. = FALSE)
    }
    s[[nm]] <- hit
  }
  for (nm in c("rates", "syn_rates", "grid_size", "starting_points")) {
    if (!is.null(s[[nm]])) {
      v <- suppressWarnings(as.integer(s[[nm]]))
      if (is.na(v) || v < 1L) stop("Invalid ", nm, ": ", s[[nm]], call. = FALSE)
      s[[nm]] <- v
    }
  }
  s
}
