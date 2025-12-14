# dndsR

***Functionality seems stable. Still performing tests before publication. Use at your own risk.***  
  
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
chmod 774 -R dndsR
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
Recommended for large-scale analysis. There is also an [R.md vignette](https://github.com/Scrumpis/dndsR/blob/main/dndsR-test-vignette.Rmd) for users who want to work in Rstudio or similar. All commands allow single or batch comparisons.  

### Sample comparison_file

### 1. split_comparisons.R
Separate subgenomes, haplotypes, or other patterns into their own fastas and gffs to prevent erroneous dN/dS analysis. If subgenomes are left unphased, the best matches will be a mix of homeolog and ortholog comparisons. If needed, [SubPhaser](https://github.com/zhangrengang/SubPhaser?tab=readme-ov-file) can be used to phase allopolyploids lacking diploid progenitor genomes.
```
singularity exec dndsr.sif ./dndsR-launcher run split_comparisons \
-C data/CheFo_vs_CheAl_full_fofn.txt -v -m subgenome
```
### 2. Extract CDS or Proteins
```
singularity exec dndsr.sif ./dndsR-launcher run extract_cds \
-C data/CheFo_vs_CheAl_full_fofn_split.txt
```
### 3. Calculate dN/dS
Long runtime. If on cluster, consider submitting through SLURM, PBS, or similar.  
```
singularity exec dndsr.sif ./dndsR-launcher run calculate_dnds \
-C data/CheFo_vs_CheAl_full_fofn_split.txt -t 80
```
### 4. Append annotations
```
singularity exec dndsr.sif ./dndsR-launcher run append_annotations \
-C data/CheFo_vs_CheAl_full_fofn_split.txt -O . -v -t 8
```
### 5. Annotation term enrichment
#### ipr_enrichment
```
singularity exec dndsr.sif ./dndsR-launcher run ipr_enrichment \
-C data/CheFo_vs_CheAl_full_fofn_split.txt -t 8 -v -O .
```
#### go_enrichment
```
singularity exec dndsr.sif ./dndsR-launcher run go_enrichment \
-C data/CheFo_vs_CheAl_full_fofn_split.txt -t 8 -v -O .
```
#### term_enrichment
```
singularity exec dndsr.sif ./dndsR-launcher run term_enrichment \
-C data/CheFo_vs_CheAl_full_fofn_split.txt -t 8 -v -O .
```
### 6. Selection pressure ideogram
```
singularity exec dndsr.sif ./dndsR-launcher run dnds_ideogram \
-C data/CheFo_vs_CheAl_full_fofn_split.txt -t 8 -v -O .
```

### 7. Regional Analysis
Regional summary between comparisons 
```
singularity exec dndsr.sif ./dndsR-launcher run regional_dnds_summary \
  -C data/CheFo_vs_CheAl_full_fofn_split_mod.txt \
  --regions-bed data/Calbum_dkmer_regions.txt \
  --contrast-file data/Calbum_contrast_file.txt \
  -O .
```


dndsR has a wrapper script () to run batches of comparisons which takes as input a space or tab separated text file containing: comparison_basename, query_fasta, query_gff3, subject_fasta, subject_gff3. 
Example:
```
CformvsCalbumIWGC "/Users/john7932/GitHub/dndsR/tests/full/Chenopodium_formosanum_chrs.fasta" "/Users/john7932/GitHub/dndsR/tests/full/CheFo_v3.gff" "/Users/john7932/GitHub/dndsR/tests/full/Chenopodium_album.genome_v2_chrs.fasta" "/Users/john7932/GitHub/dndsR/tests/full/CheAl_v01.0.note_fixed.gff"
```

## Notes
Can use special characters or spaces in path if quoted in the fofn

## Future Improvements
- Rebuild container with biomartr?
- Append annotation - seems to currently single thread each linearlly, if can't get more threads per job, then use available threads to run each in parallel. Ideally, precompute query and subject in parallel, wait until both are done, then move until next step. This might be a little complex though so at least give all comparisons a thread to run in parallel.
- Update Rvignette with test dataset
- Cleanup documentation
- Check if below is needed anymore
### Install dndsR library
#### Singularity
```
singularity exec dndsr.sif ./dndsR-launcher install
```
