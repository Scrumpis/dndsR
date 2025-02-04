message("Generating ideograms")


# Establish the base working directoy and move into it
base_dir <- "/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio"
setwd(base_dir)

# Read the file containing file names (tsv or space-delimited file)
file_of_file_names <- read.table("/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio/data/ToL_Form_test_fofn_fastas.txt",
                                 header = FALSE, stringsAsFactors = FALSE)                                 


### dNdS Ideogram ###
# ***NOTE: This will only work if your contigs/scaffolds are named as chrs
# with the number and subgenome at the end of their name with no numbers 1-9 leading the chr number
# (i.e., chr01D, Cf1D, BCVADSFG01D would all be printed as "1D" in their respective column
# but AT1H1D would be printed as 1H1D). Query must have a corresponding chr to subject to plot.

# Might want to consider having this as a separate script because placement of legend is not automated

dnds_ideogram <- function(data_type, comparison_basename, base_dir, file_of_file_names, analysis_type, window_size = 300000) {
  # Set working directory to base and then to comparison name
  setwd(base_dir)
  setwd(comparison_basename)
  
  # Load dNdS results file with merged annotation information (dnds_merged)
  dnds_merged <- read.csv(paste0(comparison_basename, "_merged.tsv"), sep = "\t")
  message("Loaded dNdS merged file for: ", comparison_basename, " (", data_type, ")")
  
  # Extract the first character of data_type ("s" or "q") to avoid long header names
  data_type_sub <- substr(data_type, 1, 1)
  #data_type <- "query"
  
  # Extract chromosome numbers
  #dnds_merged[[paste0(data_type_sub, "_chr")]] <- gsub("^[^1-9]*([1-9].*)", "\\1", dnds_merged[[paste0(data_type_sub, "_seqname")]])
  
  # Create dNdS windows
  dnds_windows <- dnds_merged %>%
    filter(!is.na(dNdS) &  # Remove rows with dNdS NA values
             dNdS < 10 &    # Remove rows with dNdS > 10, as these are likely errors
             !!rlang::parse_expr(analysis_type)) %>% # Set chr==chr, all_vs_all, etc comparison type
    
    mutate(ln_dNdS = log(dNdS + 1),          # Log transform values for color contrast in ideogram            
           ln_neg_dNdS = ifelse(dNdS >= 1, NA, log(dNdS + 1)),
           ln_pos_dNdS = ifelse(dNdS < 1, NA, log(dNdS + 1))) %>%
    
    group_by(!!sym(paste0(data_type_sub, "_chr")), 
             start = (!!sym(paste0(data_type_sub, "_gene_start")) %/% window_size) * window_size,
             end = start + window_size) %>%
    
    summarize(total_raw_dNdS = sum(dNdS, na.rm = TRUE),
              total_ln_dNdS = sum(ln_dNdS, na.rm = TRUE),
              total_ln_neg_dNdS = sum(ln_neg_dNdS, na.rm = TRUE),
              total_ln_pos_dNdS = sum(ln_pos_dNdS, na.rm = TRUE),
              n_genes = n(),
              neg_genes = sum(dNdS < 1, na.rm = TRUE),
              pos_genes = sum(dNdS >= 1, na.rm = TRUE)) %>%
    
    mutate(avg_ln_dNdS = total_ln_dNdS / n_genes,
           avg_ln_neg_dNdS = total_ln_neg_dNdS / neg_genes,
           avg_ln_pos_dNdS = total_ln_pos_dNdS / pos_genes)
  
  # Load and process karyotype data
  fasta_path <- normalizePath(file_of_file_names[i, ifelse(data_type == "query", 3, 2)])
  fai_path <- paste0(fasta_path, ".fai")
  
  if (!file.exists(fai_path)) system(paste("samtools faidx", fasta_path)) # if fai does not exist, create it
  fai <- read.csv(fai_path, sep = "\t", header = FALSE)
  
  karyotype <- fai %>%
    mutate(V1 = gsub("^[^1-9]*([1-9].*)", "\\1", V1), Start = 0) %>% # Don't know if I like this
    select(V1, Start, V2) %>%
    rename(Chr = V1, End = V2) %>%
    #arrange(as.numeric(Chr)) <- didnt work, worried about if chr11 or chr11B are used
    arrange(substr(Chr, 2, 2), substr(Chr, 1, 1)) # Groups polyploids by subgenome, then orders by Chr number sorts by second character, then first, so would it make chr 12 go after chr 1?
  # Does the above arrange() need to be optional if polyploid flag used?
  # Does this do what I think its doing? Like if it was a number at the end, how would it sort?
  
  # Adjust window ends to match karyotype ends
  last_rows <- dnds_windows %>%
    group_by(!!sym(paste0(data_type_sub, "_chr"))) %>%
    slice_tail(n = 1) %>% # Select the last row for each chromosome
    mutate(end = karyotype$End[match(!!sym(paste0(data_type_sub, "_chr")), karyotype$Chr)]) %>%
    ungroup() # Ensure that grouping is removed if necessary - do I need this?
  
  # Replace that row with the end value of the karyotype for accurate last window
  dnds_windows <- dnds_windows %>%
    left_join(last_rows %>% 
                select(!!sym(paste0(data_type_sub, "_chr")), start, end),
              by = c(paste0(data_type_sub, "_chr"), "start"),
              suffix = c(".x", ".y")) %>%
    mutate(end = coalesce(end.y, end.x)) %>% # check if some of these are needed
    select(-end.x, -end.y) %>%
    relocate(end, .after = start)
  
  # Subset data for Rideogram
  neg_dnds <- dnds_windows %>%
    select(!!sym(paste0(data_type_sub, "_chr")), start, end, avg_ln_neg_dNdS) %>%
    filter(!is.nan(avg_ln_neg_dNdS)) %>%
    arrange(substr(!!sym(paste0(data_type_sub, "_chr")), 2, 2), # might need to check this if works with Chr11
            substr(!!sym(paste0(data_type_sub, "_chr")), 1, 1), start) %>% # these sort by subgenome I think
    mutate(avg_ln_neg_dNdS = avg_ln_neg_dNdS * -1) %>% # invert values for visualization (so yellow is high in neg, but low in pos)
    rename(Chr = !!sym(paste0(data_type_sub, "_chr")), Start = start, End = end, Value = avg_ln_neg_dNdS)
  
  pos_dnds <- dnds_windows %>%
    select(!!sym(paste0(data_type_sub, "_chr")), start, end, avg_ln_pos_dNdS) %>%
    filter(!is.nan(avg_ln_pos_dNdS)) %>%
    arrange(substr(!!sym(paste0(data_type_sub, "_chr")), 2, 2), 
            substr(!!sym(paste0(data_type_sub, "_chr")), 1, 1), start) %>%
    rename(Chr = !!sym(paste0(data_type_sub, "_chr")), Start = start, End = end, Value = avg_ln_pos_dNdS)
  
  # Generate Ideogram
  ideogram(karyotype = karyotype, overlaid = neg_dnds, 
           label = pos_dnds, label_type = "heatmap", 
           colorset1 = c("#FFFFCC", "#e34a33"), 
           colorset2 = c("#FFFFCC", "#2c7fb8"), Ly = 3)
  
  # Adjust RIdeogram SVG to be poster-ready (bold font and arial)
  # Read the SVG file
  svg <- read_xml("chromosome.svg")
  
  # Find all text elements and modify their attributes
  text_nodes <- xml_find_all(svg, ".//text")
  
  xml_attr(text_nodes, "font-size") <- "14"
  xml_attr(text_nodes, "font-weight") <- "bold"
  
  # Print the text content of all text nodes to identify the ones you need to change
  text_content <- xml_text(text_nodes)
  
  # Manually identified indices for "Low" and "High"
  low_indices <- which(text_content == "Low")
  high_indices <- which(text_content == "High")
  
  # Modify the labels
  xml_text(text_nodes[low_indices]) <- "Neut"
  
  # Check the number of high labels and modify accordingly
  if (length(high_indices) == 2) {
    xml_text(text_nodes[high_indices[1]]) <- "Neg"
    xml_text(text_nodes[high_indices[2]]) <- "Pos"
  } else {
    stop("Unexpected number of 'High' labels found.")
  }
  
  # Find the parent node of the scale bars to place the new title correctly
  parent_node <- xml_parent(text_nodes[high_indices[1]])
  
  # Create a new text element for the legend title
  new_text_node <- xml_add_child(parent_node, "text", "Selection Pressure")
  
  # Set attributes for the legend header
  xml_set_attr(new_text_node, "x", "600")  # Set the legend header x position
  xml_set_attr(new_text_node, "y", "0")  # Set the legend header y position
  xml_set_attr(new_text_node, "font-size", "14")
  xml_set_attr(new_text_node, "font-weight", "bold")
  xml_set_attr(new_text_node, "font-family", "Arial")
  xml_set_attr(new_text_node, "text-anchor", "middle")
  
  # Write the modified SVG back to a file
  write_xml(svg, "chromosome.svg")
  convertSVG("chromosome.svg", device = "png") # COULD HAVE FLAG TO TURN THIS ON, good for adjusting labels
  
  # Rename files based on comparison name
  svg_name <- paste0(comparison_basename, "_", data_type, "_ideogram")
  file.rename("chromosome.svg", paste0(svg_name, ".svg"))
  file.rename("chromosome.png", paste0(svg_name, ".png"))
  
  message("Finished processing ", data_type, " for: ", comparison_basename)
}



# Loop through both query and subject with the dnds_ideogram function
for (i in 1:nrow(file_of_file_names)) {
  comparison_basename <- file_of_file_names[i, 1]
  dnds_ideogram("query", comparison_basename, base_dir, file_of_file_names, analysis_type)
  dnds_ideogram("subject", comparison_basename, base_dir, file_of_file_names, analysis_type)
}



message("Finished generating genome-wide dNdS ideograms for all")
message(Sys.time())  
