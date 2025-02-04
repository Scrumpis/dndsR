### Extract CDS Regions from Each Fasta Using GFFs ###
# Nice that this is in R to keep it one place but could just do bash script for pipeline too.
# Establish the base working directoy and move into it
base_dir <- "/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio"
setwd(base_dir)

# Read the file containing file names (tsv or space-delimited file)
file_of_file_names <- read.table("/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio/data/ToL_Form_test_fofn_fastas.txt",
                                 header = FALSE, stringsAsFactors = FALSE)

# Loop through each row of the file_of_file_names
for (i in 1:nrow(file_of_file_names)) {
  
  # Set Working Directory to Base Level
  setwd(base_dir)
  
  # Create directory for comparison basename if it doesn't exist and change into it
  comparison_basename <- file_of_file_names[i, 1]
  dir.create(comparison_basename, showWarnings = FALSE)
  setwd(comparison_basename)
  
  # System message
  message("Extracting CDS for: ", comparison_basename)
  
  # Set output variable for each run
  output_file <- paste0(comparison_basename, ".tsv")
  
  # Check if the output file already exists, and skip if it does
  if (file.exists(output_file)) {
    message("Skipping ", comparison_basename, " comparison as ", output_file, " already exists.")
    next  
  }
  
  # Grab FASTA paths from input file
  subject_fasta_path <- normalizePath(file_of_file_names[i, 2]) 
  query_fasta_path <- normalizePath(file_of_file_names[i, 3])
  
  # Extract the full file path from the current row
  subject_gff_path <- normalizePath(file_of_file_names[i, 4])
  query_gff_path <- normalizePath(file_of_file_names[i, 5])
  
  
  
  # Extract CDS with gffread #
  # Create a temporary output files (prob want to change this so files are printed to dir and can be kept if desired with flag)
  # figure out if this is actually necessary
  subject_CDS <- file.path(getwd(), "subject_CDS.fasta")
  query_CDS <- file.path(getwd(), "query_CDS.fasta")
  
  
  # Need something to replace ?s with .s, or leave up to end user.
  # Extract subject CDS with gffread
  system(paste("gffread -x", subject_CDS, "-g", subject_fasta_path, subject_gff_path))
  message("Extracted subject CDS")
  
  # Extract query CDS
  system(paste("gffread -x", query_CDS, "-g", query_fasta_path, query_gff_path))
  message("Extracted query CDS")
  
  # System message
  message("Finished CDS extraction for: ", comparison_basename)
  message(Sys.time())
  
}
