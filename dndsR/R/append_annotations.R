#' Append GFF functional annotations to dN/dS results
#'
#' This function parses a GFF file to extract functional annotations (e.g., IPR, GO, Note),
#' maps them to the corresponding CDS features, and appends the annotations to a dN/dS results table.
#'
#' @param dnds_file Path to the dN/dS results file (TSV).
#' @param gff_file Path to the GFF file with functional annotations.
#' @param output_file Optional output file path to write annotated results (TSV format).
#' @return A data.frame with dN/dS values and appended functional annotations.
#' @export
append_annotations <- function(dnds_file, gff_file, output_file = NULL) {
  # Helper: Parse GFF attributes into a named list
  parse_gff_attributes <- function(attr_string) {
    kv_pairs <- strsplit(attr_string, ";")[[1]]
    kv_list <- sapply(kv_pairs, function(kv) {
      kv_split <- strsplit(kv, "=")[[1]]
      if (length(kv_split) == 2) return(setNames(kv_split[2], kv_split[1]))
      return(NULL)
    }, simplify = FALSE)
    return(do.call(c, kv_list))
  }

  # Helper: Build a mapping of ID -> annotation info and parent relationships
  build_id_maps <- function(gff_df) {
    id_to_attrs <- list()
    parent_to_children <- list()
    for (i in seq_len(nrow(gff_df))) {
      row <- gff_df[i, ]
      attrs <- parse_gff_attributes(row$attributes)
      if (is.null(attrs$ID)) next

      # Store functional annotation if present
      has_func_info <- any(grepl("IPR|GO|Note", names(attrs), ignore.case = TRUE))
      if (has_func_info) {
        id_to_attrs[[attrs$ID]] <- attrs
      }

      # Store parent-child relationships
      if (!is.null(attrs$Parent)) {
        parents <- strsplit(attrs$Parent, ",")[[1]]
        for (p in parents) {
          parent_to_children[[p]] <- unique(c(parent_to_children[[p]], attrs$ID))
        }
      }
    }
    return(list(id_to_attrs = id_to_attrs, parent_to_children = parent_to_children))
  }

  # Helper: Trace all descendants of each functional ID to find related CDS features
  collect_functional_annotations <- function(id_to_attrs, parent_to_children) {
    cds_to_annotations <- list()
    for (id in names(id_to_attrs)) {
      attrs <- id_to_attrs[[id]]

      # Traverse to children to find CDS features
      all_descendants <- function(id) {
        result <- c()
        queue <- c(id)
        while (length(queue) > 0) {
          current <- queue[1]
          queue <- queue[-1]
          result <- c(result, current)
          if (!is.null(parent_to_children[[current]])) {
            queue <- c(queue, parent_to_children[[current]])
          }
        }
        return(unique(result))
      }
      descendants <- all_descendants(id)
      for (cds_id in descendants) {
        # Append functional terms per CDS
        cds_to_annotations[[cds_id]] <- unique(c(cds_to_annotations[[cds_id]],
                                                 paste0(names(attrs), "=", attrs)))
      }
    }
    # Collapse annotation list to one string per CDS
    cds_annot_df <- data.frame(CDS = names(cds_to_annotations),
                               annotation = sapply(cds_to_annotations, function(x) paste(x, collapse = ";")),
                               stringsAsFactors = FALSE)
    return(cds_annot_df)
  }

  # Helper: Read GFF
  read_gff <- function(path) {
    gff_raw <- read.table(path, sep = "\t", header = FALSE, comment.char = "#", quote = "", stringsAsFactors = FALSE, fill = TRUE)
    colnames(gff_raw)[1:9] <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
    return(gff_raw)
  }

  # Load dN/dS and GFF
  dnds_df <- read.table(dnds_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  gff_df <- read_gff(gff_file)

  # Build maps and collect annotations
  maps <- build_id_maps(gff_df)
  cds_annot_df <- collect_functional_annotations(maps$id_to_attrs, maps$parent_to_children)

  # Merge
  dnds_df_annot <- merge(dnds_df, cds_annot_df, by.x = "query_id", by.y = "CDS", all.x = TRUE)

  # Write output if needed
  if (!is.null(output_file)) {
    write.table(dnds_df_annot, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  return(dnds_df_annot)
}
