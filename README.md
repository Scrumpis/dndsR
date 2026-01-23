# dndsR
<!-- badges: start -->
[![R-CMD-check](https://github.com/Scrumpis/dndsR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Scrumpis/dndsR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
  
**dndsR** is a containerized R package for scalable calculation, analysis, and publication-ready visualization of the ratio of non-synonymous (dN) to synonymous (dS) substitutions (dN/dS) between orthologous genes in pairwise genome comparisons. The dN/dS ratio is commonly used to determine whether genes are conserved by selection (dN/dS < 1), diverging (dN/dS > 1), or evolving neutrally (dN/dS ~ 1).
  
Using genome FASTA and gene annotation GFF3 input files, **dndsR** provides end-to-end functionality to:
- Extract coding sequences (CDS) or proteins
- Calculate dN/dS
- Summarize selection pressures across genomes through ideograms
- Calculate enrichment of gene functional annotation terms (InterPro (IPR), Gene Ontology (GO), etc.) under positive/diversifying selection
- Analyze broad dN/dS patterns across whole genomes or genomic regions
  
**dndsR** is primarily built for containerized command-line usage but is also a loadable R library for more advanced users.

## Table of Contents
* [Setup](#setup)
* [Usage](#usage)
* [Contributing](#contributing)

## Setup
**Containerized command-line interface (CLI) usage** with Docker or Singularity is recommended (follow instructions starting at Clone Repo), however, if the below are present locally, dndsR can be used directly in Rstudio or similar (see dndsR_usage_vignette.Rmd).
- ```Orthofinder 2.5.4```
- ```Diamond 2.1.14```
- Patched ```orthologr``` version which calls ```pwalign::``` commands, not deprecated ```pairwiseAlignment::``` (required to calcualte dN/dS; included in Docker image)
  
Advanced users may also execute R scripts or interactive R sessions directly within the container.
### Clone Repo
```
git clone https://github.com/Scrumpis/dndsR
cd dndsR
```
### Pull Docker Image
Pull the Docker image from [DockerHub](https://hub.docker.com/r/scrumpis/dndsr) using either Singularity/Apptainer or Docker. Replace ```singularity``` with ```apptainer``` if needed in the below commands.  
  
**Performance note:** The dndsR Docker image is built for ```linux/amd64```. Apple-silicon (M-series) Macs will run the image under emulation, which may be significantly slower for large genomes. For production-scale analyses, Linux-based HPC usage is strongly recommended.  
  
Singularity/Apptainer:
```
singularity pull dndsr.sif docker://scrumpis/dndsr:latest
```
Docker:
```
docker pull scrumpis/dndsr:latest
```
### Install
Singularity/Apptainer:
```
tools/dndsr-install install --engine singularity --sif dndsr.sif
```
Docker (example using Mac):
```
tools/dndsr-install install --engine docker --image scrumpis/dndsr:latest --docker-platform linux/amd64
```
### Make dndsR findable and allow tab completion  
Bash:
```
bash source ~/.bashrc
```
zsh (macOS default):
```
source ~/.zshrc
```
### Update (if needed)
Used to update to current repo version if using an older version.  
While in the dndsR cloned repo:
```
git pull
dndsr
```

## Usage
_**Command-line interface (CLI) usage recommended for large-scale analysis. See dndsR_usage_vignette.Rmd for R usage.**_  
  
The below documents a typical CLI workflow with dndsR. All functions will produce outputs for both the query and subject of a comparison by default. All commands allow single or batch comparisons. Batch mode takes as input a headerless space or tab separated text file (comparison_file) containing: ```comparison_basename, "query_fasta", "query_gff3", "subject_fasta", "subject_gff3"```, with a row for each pairwise comparison being made.
### Sample comparison_file
_**Use relative paths from ```output_dir``` for Docker**_
```
Cform_v_Calbum ./path_to/Cformosanum.fasta ./path_to/Cformosanum.gff ./path_to/Calbum.fasta ./path_to/Calbum.gff
CalbumB_v_CalbumC ./path_to/CalbumB.fasta ./path_to/CalbumB.gff ./path_to/CalbumC.fasta ./path_to/CalbumC.gff
```

### Usage notes
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
- If conducting annotation term enrichment, ensure GFF term versions are the same in each comparison (i.e., both annotated with IPR 83.0) to avoid erroneous results.
- Use ```dndsr --help``` or ```dndsr <function> --help``` for more information on usage.

### 1. Split by subgenome or haplotype (optional)
Separates subgenomes, haplotypes, or other chromosome/scaffold naming patterns into their own FASTA and GFF3 files and generates a new ```comparison_file``` corresponding to the splits. Generally recommended for polyploid comparisons so best matches occur between the same subgenome. [SubPhaser](https://github.com/zhangrengang/SubPhaser?tab=readme-ov-file) can be used to phase allopolyploids lacking diploid progenitor genomes.
```
dndsr split_comparisons -c comparison_file.txt
```
### 2. Extract CDS or proteins
Extracts CDS or proteins into new FASTA files using the input genome.fasta and genome.gff files for each species of each comparison in comparison_file.
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
### 7. Comparative analysis (under development)
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
- extract_cds - if fasta or gff input paths are not being parsed while multi-threading, logs do not indicate that. Logs print as if run was fine. But in single-thread mode it states files not found.
- dnds_ideogram - Remove PNG convsersion and update SVG headers for in-browser viewing, parallelization, add multi-pattern removal for chr labels
- go_enrichment - topN.svg -> top$N.svg output file naming
- term_enrichment - fix font
- All term enrichment - add note to documentation  If conducting annotation term enrichment, ensure GFF term versions are the same in each comparison (i.e., both annotated with IPR 83.0) to avoid erroneous results.
- Rvignette - add note  If conducting annotation term enrichment, ensure GFF term versions are the same in each comparison (i.e., both annotated with IPR 83.0) to avoid erroneous results.
- Remove warnings flag and just auto-expand warnings always.

## Contributing
dndsR is under active development. Contributions, bug reports, and feature requests are welcome. 
If you are using dndsR in your research and encounter issues or have ideas
for improvement, please feel free to open an issue, start a discussion, or contact Nick Johnson directly.

