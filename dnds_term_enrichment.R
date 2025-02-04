### IPR Enrichment Plot ###

# Optional Analysis Flag
#perform_analysis <- TRUE  # Set this to FALSE to skip analysis

#if (perform_analysis) {


#if (term_type=IPR) {
#
#}

# Establish the base working directoy and move into it
base_dir <- "/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio"
setwd(base_dir)

# Read the file containing file names (tsv or space-delimited file)
file_of_file_names <- read.table("/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio/data/ToL_Form_test_fofn_fastas.txt",
                                 header = FALSE, stringsAsFactors = FALSE)                                 

# Enrichment Function
term_enrichment <- function(data_type, base_dir, comparison_basename, analysis_type) {
  # Set working directory to the comparison folder
  setwd(base_dir)
  setwd(comparison_basename)
  
  # Load the dNdS results file
  dnds_merged <- read.csv(paste0(comparison_basename,"_merged.tsv"), sep = "\t") # Read as TSV
  message("Loaded dNdS merged file for: ", comparison_basename, " (", data_type, ")")
  
  # Extract the first character of data_type ("s" or "q") to avoid long header names
  data_type_sub <- substr(data_type, 1, 1)
  
  # Extract chromosome number and append to dnds_merged input
  chr_col <- paste0(data_type_sub, "_chr")
  seqname_col <- paste0(data_type_sub, "_seqname")
  dnds_merged[[chr_col]] <- gsub("^[^1-9]*([1-9].*)", "\\1", dnds_merged[[seqname_col]]) #do we trust the gsub thing?
  
  # Filter dNdS > 10, NAs, and interchromosomal comparisons
  dnds_merged_filter <- dnds_merged %>%
    filter(!is.na(dNdS) &  # Remove rows with dNdS NA values
             dNdS < 10 &    # Remove rows with dNdS > 10, as these are likely errors
             !!rlang::parse_expr(analysis_type)) #%>% # Set chr==chr, all_vs_all, etc comparison type
  
  # Subset dNdS > 1 (positive selection) genes
  dnds_merged_pos <- dnds_merged_filter %>%
    filter(dNdS > 1)
  
  # Calculate IPR counts
  ipr_ID <- paste0(data_type_sub, "_IPR_ID")
  ipr_desc <- paste0(data_type_sub, "_IPR_description")
  ipr_pos_count <- table(dnds_merged_pos[[ipr_ID]])
  ipr_total_count <- table(dnds_merged_filter[[ipr_ID]])
  
  # Run Fisher's Exact Test
  fisher_pvals <- lapply(names(ipr_pos_count), function(ipr) {
    contingency_table <- matrix(
      c(sum(ipr_pos_count[ipr]), sum(ipr_pos_count),
        sum(ipr_total_count[ipr]) - sum(ipr_pos_count[ipr]),
        sum(ipr_total_count) - sum(ipr_pos_count)),
      nrow = 2
    )
    fisher_results <- fisher.test(contingency_table, alternative = "greater")
    data.frame(
      IPR = ipr,
      p_value = fisher_results$p.value,
      total_count = ipr_total_count[ipr],
      pos_count = sum(ipr_pos_count[ipr]),
      stringsAsFactors = FALSE
    )
  })
  
  fisher_pvals <- do.call(rbind, fisher_pvals)
  
  # Perform Benjamini-Hochberg procedure
  fisher_pvals$adjusted_p_value <- p.adjust(fisher_pvals$p_value, method = "BH")
  
  # Add IPR annotations
  fisher_pvals <- fisher_pvals %>%
    right_join(
      dnds_merged_filter %>% distinct(across(c(!!sym(ipr_ID), !!sym(ipr_desc)))),
      by = c("IPR" = ipr_ID)
    )
  
  # Write output into CSV
  write.csv(fisher_pvals, file = paste0(comparison_basename, "_", data_type_sub, "_IPR_Enrichment.csv"),
            row.names = FALSE)
  
  # Create Top 20 Plot
  top_20 <- fisher_pvals %>%
    arrange(adjusted_p_value, desc(total_count)) %>%
    head(20)
  
  IPR_enrichment <- ggplot(top_20, aes(x = pos_count, y = reorder(!!sym(ipr_desc), -adjusted_p_value),
                                       size = pos_count / total_count, color = adjusted_p_value)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Positive IPR Count", y = "IPR Functional Description",
         size = "Pos/Total\n Count", color = "Adjusted\n p-value") +
    theme(axis.text = element_text(face = "bold", size = 16),
          legend.text = element_text(face = "bold", size = 16),
          legend.title = element_text(face = "bold", size = 20),
          axis.title.x = element_text(face = "bold", size = 20),
          axis.title.y = element_text(face = "bold", size = 20))
  
  plot_file <- paste0(comparison_basename, "_", data_type_sub, "_IPR_enrichment_plot.svg")
  ggsave(plot_file, plot = IPR_enrichment, width = 12, height = 10)
  
  message("Completed analysis for: ", comparison_basename, " (IPR enrichment completed)")
}



# Loop through both query and subject with the dnds_ideogram function
for (i in 1:nrow(file_of_file_names)) {
  comparison_basename <- file_of_file_names[i, 1]
  ipr_enrichment("query", base_dir, comparison_basename, analysis_type)
  ipr_enrichment("subject", base_dir, comparison_basename, analysis_type)
}
