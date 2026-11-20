# Real-genome regression fixture: *A. thaliana* and *A. lyrata*

A small, permanent piece of two public genomes, used by `tools/dev/regression.sh` to run the whole pipeline (`extract_cds`, `calculate_dnds`, `append_annotations`) in minutes instead of hours. It is regression evidence: it shows that a change leaves results as they were. It is not correctness evidence, which comes from the independent checks in `tests/` and `docs/validation/`. How it fits the other tiers: [docs/dev/TESTING.md](../../../docs/dev/TESTING.md).

## Contents of `atal/`

| File | What |
|---|---|
| `genes.tsv` | The 542 *A. thaliana* and 549 *A. lyrata* genes chosen, and the stratum that chose each. This list defines the fixture |
| `Athaliana.fasta.gz`, `Athaliana.gff3.gz` | 531 contigs, 2.51 Mb of real sequence, 707 genes and 18,030 features |
| `Alyrata.fasta.gz`, `Alyrata.gff3.gz` | 544 contigs, 2.59 Mb, 592 genes and 11,895 features |
| `expected/outputs.sha256` | SHA-256 of every pipeline output except logs and provenance: the baseline |
| `expected/AtAl_dnds.tsv` | The baseline dN/dS table, so a failing run can say which pairs and values changed |

## Source and licence

NCBI RefSeq assemblies, the same ones used in `docs/tutorials/public-genomes.md` and `docs/validation/public-genome.md`:

| Species | Assembly | Files, with NCBI MD5 verified on download | SHA-256 of the downloaded `.gz` |
|---|---|---|---|
| *A. thaliana* | GCF_000001735.4 (TAIR10.1) | `GCF_000001735.4_TAIR10.1_genomic.fna.gz`, `..._genomic.gff.gz` | `a182f3c71662973d...`, `d6a3aa3c00df97cf...` |
| *A. lyrata* | GCF_000004255.2 (v.1.0) | `GCF_000004255.2_v.1.0_genomic.fna.gz`, `..._genomic.gff.gz` | `2ff53c7b052e69b3...`, `1aab0db0d78c0b90...` |

Full hashes are printed by the build command below. NCBI states that it places no restrictions on the use or distribution of the data in its databases, while noting that original submitters may claim rights in some data. The excerpt is redistributed here on that basis, with this attribution. That is an assumption the project maintainer should confirm before a release that ships this directory; it is not in the R package (`tools/dev` is in `.Rbuildignore`) or the image (`.dockerignore` excludes FASTA).

## How genes were chosen

`select-atal-genes.py` read the public full-genome run of v0.4.0 against `9a5c758` (byte-identical to v0.4.0) and sampled with a fixed seed, at most 25 pairs per stratum unless noted, each chosen gene with its partner:

- dS below 0.02, 0.02-0.1 (40), 0.1-0.3 (60), 0.3-0.8 (40), above 0.8, and above 3 times the median (often not orthologs);
- dN/dS above 1, dN/dS undefined, dN = 0;
- reliability flags: fewer than 100 codons, fewer than 3 synonymous differences;
- minus strand, single exon, 10 or more exons, 4 or more isoforms, CDS under 300 nt or over 4,500 nt;
- three gene symbols each shared by two paired *A. thaliana* genes (see below);
- per species, 30 genes with a valid CDS and no reciprocal partner, and up to 12 each that `extract_cds` excludes for an internal stop or a length not a multiple of 3.

Organelle sequences are left out, as the public tutorial removes them before dndsR. `build-atal-fixture.py` then cuts each chosen gene with 1,000 bases on each side, merges windows that overlap, and keeps every feature whose top-level gene lies wholly inside a window, so neighbours, isoforms, UTRs, pseudogenes, and non-coding genes stay as annotated. Both scripts use the Python standard library and no dndsR code.

**Repeated gene symbols are deliberate.** Bioconductor names a gene by its `Name` when every `Name` in the file is unique, and by its `ID` otherwise. The whole *A. thaliana* GFF repeats 168 symbols, so dndsR names genes by locus tag (`AT1G01010`). A subset without a repeat is named by symbol instead (`NAC001`), and `append_annotations` then misses those genes (201 of 502 in a test build): an open issue in `docs/dev/DECISIONS.md`. The shared symbols keep the fixture on the whole-genome path.

## Evidence that it is faithful

At `9a5c758`, v0.4.0 image, 8 threads (results outside Git):

- All 495 full-genome pairs whose two genes were chosen are found in the fixture, and the 503 fixture pairs that are also full-genome pairs have identical dN, dS, and dN/dS.
- 5 more pairs form only in the fixture: with fewer competitors, reciprocal best hits pair genes whose full-genome best hit was left out. This is expected of the method, and is a limit of any subset.
- v0.4.0 (`957cfe1`) and `9a5c758` give byte-identical outputs on the fixture, apart from `AtAl_dnds_reliability.tsv`, which v0.4.0 does not write.
- A rebuild from the same downloads is byte-identical.
- The pipeline takes about 2 minutes with 8 threads.

## Rebuild

Only needed to change the fixture, with a recorded reason. Download the four files above into a folder outside the clone, check them against NCBI's `md5checksums.txt`, then:

```
python3 -I tools/dev/fixtures/build-atal-fixture.py tools/dev/fixtures/atal/genes.tsv \
  <dir>/GCF_000001735.4_TAIR10.1_genomic.fna.gz <dir>/GCF_000001735.4_TAIR10.1_genomic.gff.gz \
  <dir>/GCF_000004255.2_v.1.0_genomic.fna.gz <dir>/GCF_000004255.2_v.1.0_genomic.gff.gz tools/dev/fixtures/atal
sha256sum <dir>/*.gz
tools/dev/regression.sh --update     # new baseline, in the pinned image
```

Re-selecting genes (`select-atal-genes.py <full-genome run dir> <At gff.gz> <Al gff.gz> genes.tsv`) needs a full-genome run, and changes the fixture.
