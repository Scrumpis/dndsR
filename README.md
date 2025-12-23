# dndsR

***Pre-release. Still performing tests before publication. Use at your own risk.***  
  
A scalable dN/dS analysis R package. For single or multiple pairwise comparisons of genomes, dndsR...
- Extracts CDS or proteins using input genome FASTAs and gene annotation GFF3s
- Calculates dN/dS
- Visualizes selection pressures across the genomes through an ideogram
- Calculates enrichment of gene functional annotation terms (IPR, GO, etc.) for biological functions under positive/diversifying selection
- Analyzes specified retgions for overall differences in dN/dS (i.e., do Dkmer regions of B and C have dN/dS>1 significantly more frequently than
  
dndsR is built for both containerized command line usage and as a loadable library in R.

## Setup
Note: We currently recommend Docker and Singularity with CLI usage, however, if Orthofinder 2.5.4 and Diamond 2.1.14 are present on your system and R dependencies are installed, you can use the library directly in R. See Dockerfile for dependencies.
  
### Clone the Repo
```
git clone https://github.com/Scrumpis/dndsR
```
### Pull the Image
Pull image from [DockerHub](https://hub.docker.com/r/scrumpis/dndsr)  

Singularity
```
singularity pull dndsr.sif docker://scrumpis/dndsr:latest
```
Docker
```
docker pull scrumpis/dndsr:latest
```

## CLI Usage
Recommended for large-scale analysis. There is also an [R.md vignette](https://github.com/Scrumpis/dndsR/blob/main/dndsR-test-vignette.Rmd) for users who want to work in Rstudio or similar. All commands allow single or batch comparisons. All functions will produce outputs for both the query and subject of a comparison by default.    

### Sample comparison_file
dndsR was built to run batches of comparisons and can optionally take as input a space or tab separated text file (comparison_file) containing: comparison_basename, "query_fasta", "query_gff3", "subject_fasta", "subject_gff3". 
Example:
```
Cform_v_Calbum "/path_to/Cformosanum.fasta" "/path_to/Cformosanum.gff" "/path_to/Calbum.fasta" "/path_to/Calbum.gff"
CalbumB_v_CalbumC "/path_to/CalbumB.fasta" "/path_to/CalbumB.gff" "/path_to/CalbumC.fasta" "/path_to/CalbumC.gff"
```
### 1. split_comparisons.R (optional)
Separate subgenomes, haplotypes, or other patterns into their own fastas and gffs and generates a new comparison_file corresponding to the splits. This is generally recommended for calculations between polyploid species so best matches occur between the same subgenome. If needed, [SubPhaser](https://github.com/zhangrengang/SubPhaser?tab=readme-ov-file) can be used to phase allopolyploids lacking diploid progenitor genomes.
```
singularity exec dndsr.sif ./dndsR-launcher run split_comparisons \
-C comparison_file.txt -v -m subgenome
```
### 2. Extract CDS or Proteins
Extracts CDS or proteins into a new fasta using the genome.fasta and genome.gff files for each species of each comparison in comparison_file.
```
singularity exec dndsr.sif ./dndsR-launcher run extract_cds \
-C comparison_file.txt
```
### 3. Calculate dN/dS
Long runtime. If on cluster, consider submitting through SLURM, PBS, or similar.  
```
singularity exec dndsr.sif ./dndsR-launcher run calculate_dnds \
-C comparison_file.txt -t 80
```
### 4. Append annotations
Appends GFF annotation attributes, functional terms, seqname, start, and end values for both query and subject to dN/dS calculations based on gene_id.
```
singularity exec dndsr.sif ./dndsR-launcher run append_annotations \
-C comparison_file.txt -O . -v -t 8
```
### 5. Annotation term enrichment
Enrichment of IPR terms under positive selection. Comparable to topGO in function. Handles parent-child relationships of IPR terms.
#### ipr_enrichment
```
singularity exec dndsr.sif ./dndsR-launcher run ipr_enrichment \
-C comparison_file.txt -t 8 -v -O .
```
#### go_enrichment
TopGO enrichment of GO terms under positive selection.
```
singularity exec dndsr.sif ./dndsR-launcher run go_enrichment \
-C comparison_file.txt -t 8 -v -O .
```
#### term_enrichment
General term enrichment. Fisher's Exact Test and multiple testing correction. Looks for non-IPR or GO terms like KEGG, PANTHER, etc. Optionally takes as input a custom pattern of interest to test for enrichment. 
```
singularity exec dndsr.sif ./dndsR-launcher run term_enrichment \
-C comparison_file.txt -t 8 -v -O .
```
### 6. Selection pressure ideogram
Visualizes dN/dS binned values accross a genome in an ideogram.
```
singularity exec dndsr.sif ./dndsR-launcher run dnds_ideogram \
-C comparison_file.txt -t 8 -v -O .
```

### 7. Comparative Analysis
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
singularity exec dndsr.sif ./dndsR-launcher run regional_dnds_summary \
-C comparison_file.txt \
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
singularity exec dndsr.sif ./dndsR-launcher run regional_dnds_contrasts \
-C comparison_file.txt \
--regions-bed regions.bed \
--contrast-file contrast_file.txt \
-O .
```


## Notes
Can use special characters or spaces in path if quoted in the comparison_file


## Future Improvements
- ggplot::aes_string deprecated
- Update Rvignette with test dataset
- Cleanup documentation
- DBUS warning is container related, can be ignored, and will be updated
- Add forest plot function and other grand analyses
- Ensure docker image contains all dependencies of all optional functions
- Verify these are unusused and rebuild image: Biogenerics, IRanges, S4Vectors, optparse, readr
- Change regional analysis name. Can focus regions, but can do whole genome summaries too.
  - Change default to regions_bed = NULL
- Check if below is needed anymore
### Install dndsR library
#### Singularity
```
singularity exec dndsr.sif ./dndsR-launcher install
```
