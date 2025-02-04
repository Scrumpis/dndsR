# Establish the base working directoy and move into it
base_dir <- "/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio"
setwd(base_dir)

# Read the file containing file names (tsv or space-delimited file)
file_of_file_names <- read.table("/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio/data/ToL_Form_test_fofn_fastas.txt",
                                 header = FALSE, stringsAsFactors = FALSE)                                 

# OrthologR Settings #
# Default settings listed below
# For more information on options, type help("dNdS") or go to Orthologr github
aligner <- "diamond"
sensitivity_mode <- "fast"
aligner_path <- NULL
seq_type <- "cds"
format <- "fasta"
ortho_detection <- "RBH" # perform DIAMOND best reciprocal hit orthology inference
delete_corrupt_cds <- TRUE # coding sequences that cannot be divided by 3 (triplets) will be removed
store_locally <- FALSE
cdd.path <- NULL
aligner_params <- NULL
eval <- "1E-5"
ortho_path <- NULL
aa_aln_type <- "pairwise" # perform pairwise global alignments of AA seqs 
aa_aln_tool <- "NW" # using Needleman-Wunsch
aa_aln_path <- NULL
aa_aln_params <- NULL
codon_aln_tool <- "pal2nal" # perform codon alignments using Pal2Nal
kaks_calc_path <- NULL
dnds_est.method <- "Comeron" # use Comeron's method for dN/dS inference
comp_cores <- 10 # number of compute cores
quiet <- TRUE
clean_folders <- FALSE
print_citation <- TRUE


# Set options to prevent scientific notation
options(scipen = 999)


# Load libraries
library(orthologr)



# Loop through each row of the file_of_file_names
for (i in 1:nrow(file_of_file_names)) {
  
  # Set Working Directory to Base Level
  setwd(base_dir)
  
  # Create directory for comparison basename if it doesn't exist and change into it
  comparison_basename <- file_of_file_names[i, 1]
  dir.create(comparison_basename, showWarnings = FALSE)
  setwd(comparison_basename)
  
  # System message
  message("Running OrthologR for: ", comparison_basename)
  
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
  
  
  
  ### Orthologr dN/dS Calculation ###
  
  # Set output variable for each run
  #output_file <- paste0(comparison_basename, ".tsv")
  
  # Check if the output file already exists, and skip if it does
  if (file.exists(output_file)) {
    message("Skipping ", comparison_basename, " comparison as ", output_file, " already exists.")
    next  
  }
  
  # Print status message
  message("Running Orthologr for comparison:", comparison_basename)
  message(Sys.time())
  
  # Run the dNdS function for the current comparison
  dnds_basename <- dNdS(query_file      = query_CDS,
                        subject_file    = subject_CDS,
                        aligner = aligner,
                        sensitivity_mode = sensitivity_mode,
                        aligner_path = aligner_path,
                        seq_type = seq_type,
                        format = format,
                        ortho_detection = ortho_detection,
                        delete_corrupt_cds = delete_corrupt_cds,
                        store_locally = store_locally,
                        cdd.path = cdd.path,
                        aligner_params = aligner_params,
                        eval = eval,
                        ortho_path = ortho_path,
                        aa_aln_type = aa_aln_type,
                        aa_aln_tool = aa_aln_tool,
                        aa_aln_path = aa_aln_path,
                        aa_aln_params = aa_aln_params,
                        codon_aln_tool = codon_aln_tool,
                        kaks_calc_path = kaks_calc_path,
                        dnds_est.method = dnds_est.method,
                        comp_cores = comp_cores,
                        quiet = quiet,
                        clean_folders = clean_folders,
                        print_citation = print_citation
  )
  
  # Save the result based on the comparison_basename
  #write.csv(dnds_basename, file = output_file, row.names = FALSE) #WRITE AS TSV #############
  write.table(dnds_basename, file = output_file, row.names = FALSE, sep = "\t", quote = FALSE)
  
  # Print completion message
  message("Completed Orthologr for comparison: ", comparison_basename)
  message(Sys.time())
  
}

message("Orthologr Step Complete")
message(Sys.time())
