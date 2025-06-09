# dndsR
A scalable dN/dS analysis R package. dndsR was designed with non-computational scientists in mind, aimed to be easy to run with a single command.

## Installation
Docker? Conda? Both? Need orthofinder and diamond2 in conda env

## Data preparation
Separate subgenomes or phased haplotypes into their own genome, protein, or CDS fastas prior to calculating dN/dS. [Orthologr](https://github.com/drostlab/orthologr), which calculates dN/dS, outputs the best match for every query/subject CDS/protein comparison. If subgenomes or haplotypes are left unphased, the best matches will be a mix of homeolog and ortholog comparisons, leading to spurious results. We have had success using [SubPhaser](https://github.com/zhangrengang/SubPhaser?tab=readme-ov-file) to phase subgenomes of allopolyploids with unavalaiable diploid progenitor genomes.

provide a BASH script () and R helper function () to assist in subgenome or haplotype FASTA seperation.

## Usage
Start-up
To run, open Rstudio from a conda env containing OrthoFinder

dndsR has a wrapper script () to run batches of comparisons which takes as input a space or tab separated text file containing: comparison_basename, query_fasta, query_gff3, subject_fasta, subject_gff3. 
Example:
```
CformvsCalbumIWGC "/Users/john7932/GitHub/dndsR/tests/full/Chenopodium_formosanum_chrs.fasta" "/Users/john7932/GitHub/dndsR/tests/full/CheFo_v3.gff" "/Users/john7932/GitHub/dndsR/tests/full/Chenopodium_album.genome_v2_chrs.fasta" "/Users/john7932/GitHub/dndsR/tests/full/CheAl_v01.0.note_fixed.gff"
```
The wrapper will accomplish the below tasks, which can also be run piecewise as their own callable functions from the dndsR package.
1. Extract CDS sequence from

## Notes
Can use special characters or spaces in path if quoted in the fofn
