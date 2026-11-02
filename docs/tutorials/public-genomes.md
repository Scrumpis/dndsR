# A full analysis of public genomes

From downloads to a term list, with *Arabidopsis thaliana* and *A. lyrata* RefSeq assemblies, and GO and InterPro terms from UniProt. This assumes you have completed [getting started](getting-started.md), so that `dndsr` works. It also needs `curl` and `python3`.

It downloads about 600 MB, and writes about 1.5 GB. Most of the time is `calculate_dnds`. Measured times and memory are at the end. On a cluster, run it through your scheduler.

## 1. Get the genomes and terms

Run this from the folder holding your `dndsR` clone.

```
mkdir -p public_demo/inputs
base=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF
curl -sSf -o public_demo/inputs/At.fasta.gz $base/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
curl -sSf -o public_demo/inputs/At.gff.gz $base/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz
curl -sSf -o public_demo/inputs/Al.fasta.gz $base/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_genomic.fna.gz
curl -sSf -o public_demo/inputs/Al.gff.gz $base/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_genomic.gff.gz
gunzip -f public_demo/inputs/*.gz
curl -sSf -o public_demo/inputs/At_uniprot.tsv "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000006548&fields=accession,xref_geneid,xref_interpro,go_id&format=tsv"
```

The UniProt table changes with each release, so your term counts can differ a little. Keep a hash of what you downloaded.

## 2. Add terms to the query genome

dndsR reads `GO:` and `IPR` terms from any GFF3 attribute. For your own genomes, use any annotation tool that writes them, and the same database version on both genomes. These RefSeq genes carry NCBI GeneIDs, which UniProt lists, so this script can join them.

```
cd public_demo
awk -F'\t' '$3=="region" && /genome=(chloroplast|mitochondrion)/{drop[$1]=1} !($1 in drop)' inputs/At.gff > inputs/At_nuclear.gff
awk -F'\t' '$3=="region" && /genome=(chloroplast|mitochondrion)/{drop[$1]=1} !($1 in drop)' inputs/Al.gff > inputs/Al_nuclear.gff
python3 ../dndsR/tools/annotate-from-uniprot.py inputs/At_uniprot.tsv inputs/At_nuclear.gff inputs/At_annot.gff
python3 ../dndsR/tools/check-first-exon-phase.py inputs/At_annot.gff
python3 ../dndsR/tools/check-first-exon-phase.py inputs/Al_nuclear.gff
echo "AtAl inputs/At.fasta inputs/At_annot.gff inputs/Al.fasta inputs/Al_nuclear.gff" > comparison_file.txt
```

The first two lines remove organelles, which use other genetic codes, and whose trans-spliced genes can stop `extract_cds`. The last script counts partial genes, which dndsR does not trim to frame. Here there should be few.

## 3. Calculate dN/dS

Use as many threads as you have, with about 1 GB of memory each.

```
dndsr extract_cds -c comparison_file.txt -t 8
dndsr calculate_dnds -c comparison_file.txt -t 32
dndsr append_annotations -c comparison_file.txt -t 8
```

Read the report at the end of `calculate_dnds`: how many pairs are undefined, below dS 0.02, above dS 0.8, and above 3 times the median dS.

## 4. Test terms

```
dndsr covariate_enrichment -c comparison_file.txt -t 2 --sides query
dndsr go_enrichment -c comparison_file.txt -t 2 --sides query
dndsr ipr_enrichment -c comparison_file.txt -t 2 --sides query
```

`AtAl/AtAl_q_IPR_covariate_enrichment.tsv` and `AtAl/AtAl_q_GO_covariate_enrichment.tsv` are the lists to report from, because they account for gene length and dS. Compare them with the other two to see what that removes. See [interpretation](../interpretation.md) and [the evidence](../validation/enrichment-confounding.md).

## 5. Summarize

```
dndsr dnds_summary -c comparison_file.txt
dndsr dnds_ideogram -c comparison_file.txt
```

## What to expect

RESULTS_PLACEHOLDER
