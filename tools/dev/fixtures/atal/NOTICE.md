# Third-party data: not covered by dndsR's license

The FASTA and GFF3 files in this folder are excerpts of public NCBI RefSeq records, cut to small regions around chosen genes and otherwise unchanged:

- *Arabidopsis thaliana*, assembly GCF_000001735.4 (TAIR10.1), annotation from TAIR/Araport as distributed by NCBI RefSeq.
- *Arabidopsis lyrata*, assembly GCF_000004255.2 (v.1.0), from the U.S. Department of Energy Joint Genome Institute and NCBI RefSeq.

They are not dndsR's work, and the MIT license in `LICENSE.md` does not apply to them. dndsR claims no rights in them and applies no license to them. NCBI states that it places no restrictions on the use or distribution of the data in its molecular databases, and that the original submitters may claim patent, copyright, or other intellectual property rights in some data, which NCBI cannot assess or waive ("Policies and Disclaimers", section "Molecular Data Usage"). They are redistributed here, in this small extract, only to test dndsR, on those terms and with this attribution. To use the data for anything else, obtain it from NCBI.

`genes.tsv` (the gene list), `expected/` (outputs of dndsR), and the scripts in the parent folder are dndsR's. How the excerpt was made, with checksums of the source files: `../README.md`.
