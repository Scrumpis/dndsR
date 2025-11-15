# dndsR

***Functionality seems stable. Still performing tests before publication. Use at your own risk***  
  
A scalable dN/dS analysis R package. For single or multiple pairwise comparisons of genomes, dndsR...
- Extracts CDS or proteins using input genome FASTAs and gene annotation GFF3s
- Calculates dN/dS
- Visualizes selection pressures across the genomes through an ideogram
- Calculates enrichment of gene functional annotation terms (IPR, GO, etc.) for biological functions under positive/diversifying selection

Improvements:
- Rebuild container with biomartr
- Add fonts to inst

## Setup
Note: We recommend Docker and Singularity usage for either CLI or RStudio usage, however, if you have a working Conda env or other way to provide the dependencies, Docker and Singularity are not required. You can install directly into RStudio as shown below (note line once we make it).
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

### Install dndsR library
Singularity
```
singularity exec dndsr.sif ./dndsR-launcher install
```

## Data preparation
Separate subgenomes or phased haplotypes into their own genome, protein, or CDS fastas prior to calculating dN/dS. [Orthologr](https://github.com/drostlab/orthologr), which calculates dN/dS, outputs the best match for every query/subject CDS/protein comparison. If subgenomes or haplotypes are left unphased, the best matches will be a mix of homeolog and ortholog comparisons, leading to spurious results. We have had success using [SubPhaser](https://github.com/zhangrengang/SubPhaser?tab=readme-ov-file) to phase subgenomes of allopolyploids with unavalaiable diploid progenitor genomes.
For subgenome-phased genome FASTAs, a BASH command can be used to extract each into their own FASTAs. This also works to separate haplotypes in phased genomes.


## Usage
The below documents the command line interface (CLI) workflow, which is recommended for large-scale analysis. There is also an [R.md vignette](https://github.com/Scrumpis/dndsR/blob/main/dndsR-test-vignette.Rmd) for users who would to work in Rstudio or similar. All commands allow single or batch comparisons.  

### 1. split_comparisons.R
Separate subgenomes, haplotypes, or other patterns into their own fastas and gffs to prevent erroneous dN/dS analysis. If needed, [SubPhaser](https://github.com/zhangrengang/SubPhaser?tab=readme-ov-file) can be used to phase allopolyploids.
#### Singularity
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
Recommended to run on cluster, slurm script available.  
Need shims first? Test again without.
```
bash scripts/make_shims.sh shims
export PATH="/mnt/gs21/scratch/john7932/dndsR/shims:$PATH"
```
Singularity
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


dndsR has a wrapper script () to run batches of comparisons which takes as input a space or tab separated text file containing: comparison_basename, query_fasta, query_gff3, subject_fasta, subject_gff3. 
Example:
```
CformvsCalbumIWGC "/Users/john7932/GitHub/dndsR/tests/full/Chenopodium_formosanum_chrs.fasta" "/Users/john7932/GitHub/dndsR/tests/full/CheFo_v3.gff" "/Users/john7932/GitHub/dndsR/tests/full/Chenopodium_album.genome_v2_chrs.fasta" "/Users/john7932/GitHub/dndsR/tests/full/CheAl_v01.0.note_fixed.gff"
```
The wrapper will accomplish the below tasks, which can also be run piecewise as their own callable functions from the dndsR package.
1. Extract CDS sequence from

## Notes
Can use special characters or spaces in path if quoted in the fofn
