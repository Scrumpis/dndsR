# dndsR

***Pre-release. Still performing tests before publication. Use at your own risk.***  
  
<!-- badges: start -->
[![R-CMD-check](https://github.com/Scrumpis/dndsR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Scrumpis/dndsR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
  
An R package for scalable calculation and analysis of the ratio of non-synonymous (dN) to synonymous (dS) substitutions (dN/dS) between orthologous genes in pairs of genomes. For single or multiple pairwise comparisons of genomes, dndsR...
- Extracts CDS or proteins using input genome FASTAs and gene annotation GFF3s
- Calculates dN/dS
- Visualizes selection pressures across the genomes through an ideogram
- Calculates enrichment of gene functional annotation terms (IPR, GO, etc.) under positive/diversifying selection
- Analyzes broad dN/dS patterns across genomes or regions of genomes (under development)
  
dndsR is primarily built for containerized command line usage but is also a loadable R library for more advanced users.

## Setup
CLI usage with Docker or Singularity is recommended, however, if Orthofinder 2.5.4, Diamond 2.1.14, and R dependencies are present on your system, you can use the library directly in RStudio or similar. See Dockerfile for dependencies and dndsR-test-vignette.Rmd (under development) for usage. The Docker image includes a patched ```orthologr``` version which calls the updated ```pwalign::``` commands instead of the deprecated ```pairwiseAlignment::```, which is required by ```dndsR``` to calcualte dNdS. 

### Clone the Repo
```
git clone https://github.com/Scrumpis/dndsR
cd dndsR
```
### Pull the Image
Pull image from [DockerHub](https://hub.docker.com/r/scrumpis/dndsr) using either Apptainer/Singularity or Docker.  

Singularity
```
singularity pull dndsr.sif docker://scrumpis/dndsr:latest
```
Docker
```
docker pull scrumpis/dndsr:latest
```
### Install
Singularity (can call Apptainer too):
```
tools/dndsr-install install --engine singularity --sif dndsr.sif
```
Docker (on Mac):
```
tools/dndsr-install install --engine docker --image scrumpis/dndsr:latest --docker-platform linux/amd64
```
#### Make dndsR findable and allow tab completion  
Bash:
```
bash source ~/.bashrc
```
zsh (macOS default):
```
source ~/.zshrc
```
### Update (if needed)
Used to update to current repo version if using older version.  
While in the dndsR cloned repo:
```
git pull
dndsr
```

## Command-Line Interface (CLI)/Terminal Usage
Recommended for large-scale analysis. All functions will produce outputs for both the query and subject of a comparison by default. All commands allow single or batch comparisons. Batch mode takes as input a space or tab separated text file (comparison_file) containing: comparison_basename, "query_fasta", "query_gff3", "subject_fasta", "subject_gff3".
### Sample comparison_file
```
Cform_v_Calbum "/path_to/Cformosanum.fasta" "/path_to/Cformosanum.gff" "/path_to/Calbum.fasta" "/path_to/Calbum.gff"
CalbumB_v_CalbumC "/path_to/CalbumB.fasta" "/path_to/CalbumB.gff" "/path_to/CalbumC.fasta" "/path_to/CalbumC.gff"
```
**Use** ```dndsr --help``` **or** ```dndsr <function> --help``` **for more information on any commands.**

## Notes
- Can use special characters or spaces in path if double-quoted in the comparison_file
- All flags are called with ```--```
- ```-``` or ```_``` are recognized flag separators (```--min_pos``` and ```--min-pos``` both work)
- Tab completion for dndsr and subcommands available
- Global short-hand flags (case-insensitive):
  - ```-c``` = ```--comparison_file``` Comparison file for batch mode.
  - ```-t``` = ```--threads``` Workers to assign to a job.
  - ```-o``` = ```--output_dir``` Output directory.
  - ```-w``` = ```--warnings``` Prints R warnings: off|summary|all (default: off).
- If warnings are not being logged while multi-threading, try single-threading to debug.
- All functionality from called packages like orthologr can be called even if not specified in docs.

### 1. split_comparisons.R (optional)
Separate subgenomes, haplotypes, or other patterns into their own fastas and gffs and generates a new comparison_file corresponding to the splits. Generally recommended for polyploid comparisons so best matches occur between the same subgenome. [SubPhaser](https://github.com/zhangrengang/SubPhaser?tab=readme-ov-file) can be used to phase allopolyploids lacking diploid progenitor genomes.
```
dndsr split_comparisons -c comparison_file.txt
```
### 2. Extract CDS or Proteins
Extracts CDS or proteins into a new fasta using the genome.fasta and genome.gff files for each species of each comparison in comparison_file.
```
dndsr extract_cds -c comparison_file.txt -t 8
```
### 3. Calculate dN/dS
Uses [orthologr](https://github.com/drostlab/orthologr) and has access to all orthologr::dNdS functionality. Long runtime. If on cluster, consider submitting through SLURM, PBS, or similar.  
```
dndsr calculate_dnds -c comparison_file.txt -t 80
```
### 4. Append annotations
Appends GFF annotation attributes, functional terms, seqname, start, and end values for both query and subject to dN/dS calculations based on gene_id.
```
dndsr append_annotations -c comparison_file.txt -t 8
```
### 5. Annotation term enrichment
Enrichment of various gene annotation functional terms under positive selection (dN/dS>1).  
  
***If conducting annotation term enrichment, ensure GFF term versions are the same in each comparison (i.e., both annotated with IPR 83.0) to avoid erroneous results***  
  
#### InterPro (IPR) term enrichment:
Enrichment of IPR terms under positive selection. Fisher's Exact Test, or optionally account for parent-child relationships, with multiple testing correction.
```
dndsr ipr_enrichment -c comparison_file.txt -t 8
```
#### Gene Ontology (GO) term enrichment:
[TopGO](https://bioconductor.org/packages/release/bioc/html/topGO.html) enrichment of GO terms under positive selection.
```
dndsr go_enrichment -c comparison_file.txt -t 8
```
#### General term enrichment:
Fisher's Exact Test and multiple testing correction. Tests non-IPR and non-GO terms like KEGG, PANTHER, etc. Optionally receives a custom pattern of interest to test for enrichment. 
```
dndsr term_enrichment -c comparison_file.txt -t 8
```
### 6. Selection pressure ideogram
Custom [RIdeogram](https://cran.r-project.org/web/packages/RIdeogram/) visualization of dN/dS binned values accross a genome in an ideogram.
```
dndsr dnds_ideogram -c comparison_file.txt
```
### 7. Comparative Analysis (under development)
The below is for making comparisons between dN/dS outputs.
Outputs dN/dS distributions and stats between comparisons. Requires regions.bed of interest for analysis (will update to default to regionless, whole genome analysis with regional option).
regions.bed (seq_name, start, end, feature_name (optional)  
Example:
```
Chr01B 71000000 72000000 SG3
Chr01B 72000000 73000000 SG3
Chr01B 73000000 74000000 SG3
```
#### Run regional_dnds_summary
```
dndsr regional_dnds_summary \
-c comparison_file.txt \
--regions-bed regions.bed \
-O .
```
#### Contrast dndsR calculations (i.e., enrichment of regional selection pressures)
contrast_file.txt (new_comparison_name, query_comparison_name, query_genome, subject_comparison, subject_genome).  
The query_ and subject_genome entries are the side of the comparison (query/subject) they belonged to in the referenced comparison. B was query in CalbumBvC and I want to use its dN/dS calculations and gene_ids, so I selected query, where subject would've used C's values.  
Example:
```
BCvsBD CalbumBvC query CalbumUkBvD query
BCvsCD CalbumBvC subject CalbumUkCvD query
```
Run regional_dnds_contrasts
```
dndsr regional_dnds_contrasts \
-c comparison_file.txt \
--regions-bed regions.bed \
--contrast-file contrast_file.txt \
-O .
```


## Future Improvements
- Overhaul regional analyses and gene contrast functions
- Add forest plot function and other grand analyses
- Add --help-more/-hm flag for CLI to expand dependency args (dndsr calculate_dnds -hm gives ?orthologr::dNdS()). 
- Update Rvignette with test dataset
- Ensure docker image contains all dependencies of all optional functions
- Cleanup documentation
- Extract recycled helpers into their own scripts
- dnds_ideogram - Remove PNG convsersion and update SVG headers for in-browser viewing, parallelization, add multi-pattern removal for chr labels
- go_enrichment - topN.svg -> top$N.svg output file naming
- term_enrichment - fix font

## Contributing
dndsR is under active development. Contributions, bug reports, and feature requests are welcome. 
If you are using dndsR in your research and encounter issues or have ideas
for improvement, please feel free to open an issue, start a discussion, or contact Nick Johnson directly.

