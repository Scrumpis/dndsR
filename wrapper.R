# Wrapper-Script in R?

# Establish the base working directoy and move into it
base_dir <- "/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio"
setwd(base_dir)

# Read the file containing file names (tsv or space-delimited file)
file_of_file_names <- read.table("/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio/data/ToL_Form_test_fofn_fastas.txt",
                                 header = FALSE, stringsAsFactors = FALSE)
#file_of_file_names <- read.table("/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio/at_test/data/ToL_Form_At_test_fofn_fastas.txt",
#                                 header = FALSE, stringsAsFactors = FALSE)
#ToL_vs_At_fofn.txt
#file_of_file_names <- read.table("/Users/john7932/GitHub/whatsaweed/dnds_pipeline/Rstudio/at_test/data/ToL_vs_At_fofn.txt",
#header = FALSE, stringsAsFactors = FALSE)

# Assign number of cores for dNdS calc
#comp <- cores <- 10

# IPR Stuff
# Download IPR entry list from the FTP server
ipr_url <- "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list" # change if updated

# Optional Analysis Flags TO SKIP ANALYSES
perform_analysis <- TRUE  # Set this to FALSE to skip analysis
# perform_IPR <- TRUE
# dNdS_only <- FALSE # write kill command after orthologr if marked true


# Analysis Filter Options #
# Comparison type for enrichments and ideogram 
#(should maybe be used as dictionary to set analysis type later)
####
compare_all <- "TRUE" # compare all
same_chrs <- "q_chr == s_chr" # make sure chromosome numbers, including subgenome, are identical
#### MAKE SURE THE ABOVE HANDLES THINGS LIKE chr01 vs Chr01 vs chr1 vs CHR1 vs etc ####

#same_chrs_diff_sub <- # make sure chromosome numbers are identical, but subgenomes can differ

# Change this if you wish to adjust how the ideograms are constructed (will it use dnds calculations made only between same chr number or use all comparisons?)
analysis_type <- same_chrs
#analysis_type <- compare_all

# Add s_type == q_type filter? Makes it so IPR groupings are same, both are domain, both are superfam, etc  
# Chr numbers match only, between chrs (so subgenomes can differ for diploid against poly)
# Numbers and letters (subgenome) match, will make only same chrs from same subgenomes match for poly comparisons

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


####### End setting variables section ###### 

# Set options to prevent scientific notation
options(scipen = 999)
