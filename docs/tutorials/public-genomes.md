# A full analysis of public genomes

From downloads to a term list, with *Arabidopsis thaliana* and *A. lyrata* RefSeq assemblies, and GO and InterPro terms from UniProt. This assumes you have completed [getting started](getting-started.md), so that `dndsr` works. It also needs `curl` and `python3`.

It takes about 2.1 GB of disk once unpacked. Most of the time is `calculate_dnds`. Measured times and memory are at the end. On a cluster, run it through your scheduler.

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

Measured on one Linux node with 32 CPUs, in a clean home folder, with `tools/dev/test-public-tutorial.sb`:

| Step | Time |
|---|---:|
| Downloads, unpacking, and input preparation | 40 s |
| `extract_cds` | 50 s |
| `calculate_dnds`, 32 threads | 13 min |
| `append_annotations` | 7 min |
| `covariate_enrichment` | 12 s |
| `go_enrichment` | 40 s |
| `ipr_enrichment` | 2 min |
| `dnds_summary` | 5 s |
| `dnds_ideogram` | 5 min |
| **Whole job** | **30 min, peak memory 26 GB, 2.1 GB on disk** |

`check-first-exon-phase.py` reports 0 of 48,147 *A. thaliana* transcripts and 6 of 39,076 *A. lyrata* transcripts with a first CDS segment of nonzero phase.

`calculate_dnds` reports 22,145 pairs, median dS 0.141, with dN/dS undefined in 55, dS below 0.02 in 24, and dS above 0.8 in 242. Of these pairs, 206 have dN/dS above 1. Significant terms at `p_adj` 0.05, query side:

| Command | GO terms | InterPro terms |
|---|---:|---:|
| `covariate_enrichment` | 1 of 3,324 | 19 of 1,477 |
| `go_enrichment` | 4 of 3,324 | not applicable |
| `ipr_enrichment` | not applicable | 19 of 6,749 |

`ipr_enrichment` and `covariate_enrichment` each list 19 InterPro terms, of which 14 are shared. `covariate_enrichment` tests fewer because, by default, it requires a term to have at least 10 genes, where `ipr_enrichment` requires 2.

**Your numbers will differ somewhat.** UniProt changes between releases, and this download has no version. This run's `At_uniprot.tsv` began `84f39c5f3baefa7e` by SHA-256. With another public source of the same species' terms, and all but one pair's dN/dS the same, [the evidence](../validation/enrichment-confounding.md#how-stable-a-term-list-is) gives 0 GO and 16 InterPro terms. The 16 are all among these 19, and the 3 extra InterPro terms and the GO term here all have `p_adj` above 0.03. Terms close to the threshold move with the annotation, so do not read much into them.

These are leads about labels. They are not evidence that any gene or function is under selection: see [interpretation](../interpretation.md). For a publication, use the [checklist](../publication-protocol.md).
