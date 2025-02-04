### Appending Annotations to dN/dS Results ###

# Read IPR data once and set keys for faster join
# Download current interpro IPR term index
download.file(ipr_url, destfile = "entry.list", mode = "wb")
ipr_data <- fread("entry.list", header = TRUE, stringsAsFactors = FALSE, data.table = TRUE)
setnames(ipr_data, c("ENTRY_AC", "ENTRY_TYPE", "ENTRY_NAME"), c("IPR_ID", "Type", "Description"))
setkey(ipr_data, IPR_ID) # see if this actually does anything







# Perform analysis in loop
for (i in 1:nrow(file_of_file_names)) {
  
  # Set the working directory and read paths
  setwd(base_dir)
  comparison_basename <- file_of_file_names[i, 1]
  setwd(comparison_basename)
  subject_gff_path <- normalizePath(file_of_file_names[i, 4])
  query_gff_path <- normalizePath(file_of_file_names[i, 5])
  
  # System message
  message("Appending annotations to dNdS output for: ", comparison_basename)
  message(Sys.time())
  
  # Load dNdS results CSVs
  dnds_file <- read.delim(paste0(comparison_basename, ".tsv"), header = TRUE, # Load dNdS file using comparison_basename.tsv
                          sep = "\t", stringsAsFactors = FALSE)
  message("Loaded dNdS file")
  
  # Skip the header lines (lines starting with "##") and then read the GFFs into a data frame
  # Load subject GFF
  subject_gff <- read.delim(text = grep("^[^#]", readLines(subject_gff_path), value = TRUE), 
                            header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  message("Loaded subject GFF")
  
  # Load query GFF
  query_gff <- read.delim(text = grep("^[^#]", readLines(query_gff_path), value = TRUE), 
                          header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  message("Loaded query GFF")
  
  # Add GFF column names (order is first (s_seqname), second (s_source), ... last (s_attributes) column)
  colnames(subject_gff) <- c("s_seqname", "s_source", "s_type", "s_gene_start", "s_gene_end",
                             "s_score", "s_strand", "s_phase", "s_attributes")
  colnames(query_gff) <- c("q_seqname", "q_source", "q_type", "q_gene_start", "q_gene_end",
                           "q_score", "q_strand", "q_phase", "q_attributes")
  
  
  # Subset GFFs #
  # Extract Parent and IPR values directly in one vectorized operation
  subject_gff_subset <- subject_gff %>%
    mutate(
      s_ID = sub("^ID=", "", str_extract(s_attributes, "ID=[^;\\n]+")),
      s_parent = sub("^Parent=", "", str_extract(s_attributes, "Parent=[^;\\n]+")),
      s_IPR_ID = str_extract(s_attributes, "IPR[0-9]{6}") # Searches for current standard of IPR with 6 digits following. If this ever changes in structure, this pattern could be changed to preserve this scripts function.
    ) %>%
    select(s_seqname, s_ID, s_parent, s_gene_start, s_gene_end, s_IPR_ID) %>% 
    filter(!is.na(s_IPR_ID))
  
  subject_gff_parents <- subject_gff %>%
    mutate(
      s_ID = sub("^ID=", "", str_extract(s_attributes, "ID=[^;\\n]+")),
      s_parent = sub("^Parent=", "", str_extract(s_attributes, "Parent=[^;\\n]+"))
    ) %>%
    select(s_ID, s_parent)
  
  # Extract Parent and IPR values directly in one vectorized operation
  query_gff_subset <- query_gff %>%
    mutate(
      q_ID = sub("^ID=", "", str_extract(q_attributes, "ID=[^;\\n]+")),
      q_parent = sub("^Parent=", "", str_extract(q_attributes, "Parent=[^;\\n]+")),
      q_IPR_ID = str_extract(q_attributes, "IPR[0-9]{6}") # Searches for current standard of IPR with 6 digits following. If this ever changes in structure, this pattern could be changed to preserve this scripts function.
    ) %>%
    select(q_seqname, q_ID, q_parent, q_gene_start, q_gene_end, q_IPR_ID) %>% 
    filter(!is.na(q_IPR_ID))
  
  query_gff_parents <- query_gff %>%
    mutate(
      q_ID = sub("^ID=", "", str_extract(q_attributes, "ID=[^;\\n]+")), # Grabs annotation ID (between ID= and following semi-colon) for joining of dNdS value. Annotation ID = query_id
      q_parent = sub("^Parent=", "", str_extract(q_attributes, "Parent=[^;\\n]+"))
    ) %>%
    select(q_ID, q_parent)
  
  
  # Merge GFFs with dnds output #
  # Query
  dnds_merged <- dnds_file %>%
    left_join(query_gff_parents, by = c("query_id" = "q_ID")) %>%
    left_join(., query_gff_subset, by = c("query_id" = "q_parent")) %>%
    left_join(., query_gff_subset, by = c("q_parent" = "q_ID"), keep = TRUE) %>%
    distinct(query_id, .keep_all = TRUE)  
  #distinct(query_id, q_IPR_ID, .keep_all = TRUE)  
  
  # Count the NA values in each column
  q_na_count_IPRx <- sum(is.na(dnds_merged$q_IPR_ID.x))
  q_na_count_IPRy <- sum(is.na(dnds_merged$q_IPR_ID.y))
  
  # Remove the column with more NA values
  if (q_na_count_IPRx > q_na_count_IPRy) {
    dnds_merged$q_seqname.x <- NULL
    dnds_merged$q_ID.x <- NULL # think about this later. We are removing if x > y but not if y > x because it will be NA if parent (gene) contains 
    dnds_merged$q_gene_start.x <- NULL
    dnds_merged$q_gene_end.x <- NULL
    dnds_merged$q_IPR_ID.x <- NULL
    dnds_merged$q_parent.y <- NULL # x > y means parent contains IPR. Parent from dnds_merged always means parent of mRNA, or gene. Genes won't have parent values.
  } else if (q_na_count_IPRy > q_na_count_IPRx) {
    dnds_merged$q_seqname.y <- NULL
    dnds_merged$q_ID.y <- NULL
    dnds_merged$q_parent.y <- NULL
    dnds_merged$q_gene_start.y <- NULL
    dnds_merged$q_gene_end.y <- NULL
    dnds_merged$q_IPR_ID.y <- NULL
  } else {
    cat("Both columns have the same number of NA entries. No column removed.")
  }
  
  
  # Subject 
  dnds_merged <- dnds_merged %>%
    left_join(subject_gff_parents, by = c("subject_id" = "s_ID")) %>%
    left_join(., subject_gff_subset, by = c("subject_id" = "s_parent")) %>%
    left_join(., subject_gff_subset, by = c("s_parent" = "s_ID"), keep = TRUE) %>%
    distinct(subject_id, .keep_all = TRUE)
  #distinct(subject_id, s_IPR_ID, .keep_all = TRUE) # NEED TO SEE HOW MULTIPLE IPRS PER GENE HANDLED with AT
  
  # Count the NA values in each column
  s_na_count_IPRx <- sum(is.na(dnds_merged$s_IPR_ID.x))
  s_na_count_IPRy <- sum(is.na(dnds_merged$s_IPR_ID.y))
  
  # Remove the column with more NA values
  if (s_na_count_IPRx > s_na_count_IPRy) {
    dnds_merged$s_seqname.x <- NULL
    dnds_merged$s_ID.x <- NULL
    dnds_merged$s_gene_start.x <- NULL
    dnds_merged$s_gene_end.x <- NULL
    dnds_merged$s_IPR_ID.x <- NULL
    dnds_merged$s_parent.y <- NULL
  } else if (s_na_count_IPRy > s_na_count_IPRx) {
    dnds_merged$s_seqname.y <- NULL
    dnds_merged$s_parent.y <- NULL
    dnds_merged$s_ID.y <- NULL
    dnds_merged$s_gene_start.y <- NULL
    dnds_merged$s_gene_end.y <- NULL
    dnds_merged$s_IPR_ID.y <- NULL
  } else {
    cat("Both columns have the same number of NA entries. No column removed.")
  }
  
  # Remove .x and .y from column names
  colnames(dnds_merged) <- gsub("\\.x$|\\.y$", "", colnames(dnds_merged))
  
  # Append IPR functional descriptions to IPR terms
  dnds_merged <- dnds_merged %>%
    left_join(ipr_data, by = c("q_IPR_ID" = "IPR_ID")) %>%
    left_join(., ipr_data, by = c("s_IPR_ID" = "IPR_ID")) %>%
    rename(q_type = Type.x) %>%
    rename(q_IPR_description = Description.x) %>%
    rename(s_type = Type.y) %>%
    rename(s_IPR_description = Description.y)
  
  # Extract chromosome number and append to dnds_merge input *ONLY WORKS IF SEQNAME ENDS WITH CHR NUMBER*
  dnds_merged$s_chr <- gsub("^[^1-9]*([1-9].*)", "\\1", dnds_merged$s_seqname)
  dnds_merged$q_chr <- gsub("^[^1-9]*([1-9].*)", "\\1", dnds_merged$q_seqname)
  
  # Save the appended dataframe
  output_file <- paste0(comparison_basename, "_merged.tsv")
  write.table(dnds_merged, file = output_file, sep = "\t", row.names = FALSE, quote = TRUE)
  
  # Print status messages
  message("Finished appending annotations to dNdS results for: ", comparison_basename)
  message(Sys.time())
}

message("Finished appending annotations to dNdS results for all")
message(Sys.time())  
