#' Choose one representative transcript per gene
#'
#' Prefers a transcript with a valid CDS, then the longest CDS, then the first transcript
#' name, so the choice does not depend on input order.
#'
#' @param tx Transcript names.
#' @param gene Gene identifier for each transcript.
#' @param length_nt CDS length in nucleotides for each transcript.
#' @param valid Whether each transcript's CDS is usable: length divisible by 3 and no
#'   in-frame stop.
#' @return Data frame with one row per gene: `gene`, `tx`, `length_nt`, `n_isoforms`,
#'   ordered by gene.
#' @keywords internal
.select_longest_isoform <- function(tx, gene, length_nt, valid = length_nt %% 3L == 0L) {
  stopifnot(length(tx) == length(gene), length(tx) == length(length_nt), !anyNA(gene))
  ord <- order(gene, !valid, -length_nt, tx)
  first <- ord[!duplicated(gene[ord])]
  data.frame(
    gene = gene[first],
    tx = tx[first],
    length_nt = length_nt[first],
    n_isoforms = as.integer(table(gene)[gene[first]]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Check coding sequences for frame and in-frame stop problems
#'
#' A terminal stop codon is allowed. Lengths that are not a multiple of 3 are reported
#' without looking for stops, because their reading frame is not defined.
#'
#' @param cds A `Biostrings::DNAStringSet` of CDS.
#' @param genetic_code NCBI genetic code identifier.
#' @return Data frame with `multiple_of_3` and `internal_stop` for each sequence.
#' @keywords internal
.cds_quality <- function(cds, genetic_code = 1L) {
  lens <- Biostrings::width(cds)
  in_frame <- lens %% 3L == 0L
  internal <- rep(NA, length(cds))
  if (any(in_frame)) {
    aa <- as.character(Biostrings::translate(
      cds[in_frame],
      genetic.code = Biostrings::getGeneticCode(as.character(genetic_code)),
      if.fuzzy.codon = "X"
    ))
    internal[in_frame] <- grepl("*", sub("\\*$", "", aa), fixed = TRUE)
  }
  data.frame(multiple_of_3 = in_frame, internal_stop = internal, row.names = NULL)
}

#' Identify the file contents and software behind a dN/dS result
#'
#' Hashes identify inputs because file dates are unreliable on purged filesystems.
#'
#' @param query_file,subject_file Paths to the sequence inputs.
#' @param args Named list of options given to `orthologr::dNdS()`.
#' @return A list suitable for writing as JSON.
#' @keywords internal
.dnds_provenance <- function(query_file, subject_file, args) {
  version_of <- function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }
  diamond <- tryCatch(
    system2("diamond", "version", stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_,
    warning = function(w) NA_character_
  )
  keep <- setdiff(names(args), c("query_file", "subject_file", "comp_cores", "quiet", "print_citation", "clean_folders"))
  list(
    dndsR = version_of("dndsR"),
    orthologr = version_of("orthologr"),
    biostrings = version_of("Biostrings"),
    r = paste(R.version$major, R.version$minor, sep = "."),
    diamond = if (length(diamond)) diamond[1] else NA_character_,
    query_md5 = unname(tools::md5sum(query_file)),
    subject_md5 = unname(tools::md5sum(subject_file)),
    parameters = args[keep]
  )
}

#' Explain why an existing dN/dS result might not match its current inputs
#'
#' @param prov_path Path to a provenance JSON file.
#' @param query_file,subject_file Current sequence inputs.
#' @return `NULL` if the recorded hashes match, otherwise a one-line reason.
#' @keywords internal
.provenance_mismatch <- function(prov_path, query_file, subject_file) {
  if (!file.exists(prov_path)) {
    return("it has no provenance record, so it was made by an older dndsR")
  }
  if (!file.exists(query_file) || !file.exists(subject_file)) {
    return("its sequence inputs are missing")
  }
  old <- jsonlite::read_json(prov_path, simplifyVector = TRUE)
  same <- identical(old$query_md5, unname(tools::md5sum(query_file))) &&
    identical(old$subject_md5, unname(tools::md5sum(subject_file)))
  if (same) NULL else "its sequence inputs have changed"
}

#' Assign genes of significant terms to all of their ancestors
#'
#' Used by the parent-child method: a parent is tested without the genes already explained
#' by a significant descendant, including through a non-significant intermediate term.
#'
#' @param order_terms Terms ordered from most specific to most general.
#' @param sig_terms Terms found significant.
#' @param term2rows List mapping each term to its gene rows.
#' @param parents_map List mapping each term to its parents.
#' @param nodes All terms that need an entry in the result.
#' @return List mapping each term to the rows claimed by its significant descendants.
#' @keywords internal
.propagate_claims <- function(order_terms, sig_terms, term2rows, parents_map, nodes) {
  claimed <- stats::setNames(rep(list(integer(0)), length(nodes)), nodes)
  for (t in order_terms) {
    passed_up <- claimed[[t]]
    if (t %in% sig_terms) passed_up <- union(passed_up, term2rows[[t]])
    for (p in parents_map[[t]]) claimed[[p]] <- union(claimed[[p]], passed_up)
  }
  claimed
}
