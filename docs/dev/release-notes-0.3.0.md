# v0.3.0 release notes: draft for approval

Paste into the GitHub release when v0.3.0 is approved. No emails, pinned issues, or Docker Hub notices are planned.

---

**dndsR 0.3.0 corrects serious errors in earlier pre-releases. Results from v0.1.0 and v0.2.0 should be recomputed.**

v0.1.0 and v0.2.0 were pre-releases that had not been scientifically validated. Validating dndsR for this release found the following defects in them. Each one ran without an error, and its output looked plausible.

| Defect | Occurs when | Consequence |
|---|---|---|
| Exons joined in the wrong order | Multi-exon genes on the minus strand, with the default `extract_cds` | Out-of-frame sequence. dN, dS, and dN/dS are wrong, and can give false dN/dS above 1 |
| Isoforms joined into a chimera | Genes with more than one annotated isoform, with the default `extract_cds` | Sequences that do not exist, with the same effects |
| Annotations from the wrong genome | A gene ID also appears in a genome read earlier in the same run, at one thread | Wrong terms and coordinates |
| Annotations missing | RefSeq-style IDs such as `gene-AT1G01010` | Empty annotations and coordinates |
| First region lost | A regions file with no header | One region left out of `dnds_summary` and `dnds_contrast` |
| Enrichment errors | InterPro parent-child method; terms with no known type; a filter that could not be evaluated | Wrong tests, or a run that was silently unfiltered |
| Multiple-testing correction after a filter | `term_enrichment` and `ipr_enrichment`, at their defaults | With no true association, 59–74% and 95–100% of runs reported a significant term. Now 1–3%, and at most 5.3% for any `ipr_enrichment` table |
| Annotations cut at `#` | `append_annotations` and `dnds_ideogram`, when a GFF attribute contains `#` | Pairs silently lost terms: 73 on *A. thaliana* TAIR10.1 |

In *A. thaliana* TAIR10.1, all 5,655 multi-exon minus-strand genes and all 10,695 multi-isoform genes had incorrect sequences, as judged against gffread. Comparing *A. thaliana* with *A. lyrata*, public genomes, the old code reported 1,503 pairs with dN/dS above 1, and the corrected code reports 206. Evidence: `docs/validation/public-genome.md`.

Enrichment, ideograms, summaries, and contrasts built on these values are affected too, even where some genes in a table were correct. Single-exon genes and plus-strand genes with one isoform had correct sequences.

**Rerun affected analyses:** `docs/migrating-to-0.3.md`. **Evidence:** `docs/validation/`.

**Also new:** `dndsr example` and `dndsr verify_example` to check an installation; QC tables and provenance records; tutorials for the CLI and R; an image containing dndsR itself; a portable Agent Skill. Gene-level `dnds_state_contrast` is disabled, because it did not give meaningful results. Failures now stop with an error rather than looking like success.

**Known limitations:** `docs/validation/known-issues.md`. Enrichment is corrected for multiple testing, not for gene length, dS, or annotation completeness, and `ipr_enrichment` corrects each entry type separately by default. A dN/dS above 1 is a candidate, not proof of positive selection: see `docs/interpretation.md`.

**Image:** `scrumpis/dndsr:v0.3.0`, digest `<from the publish workflow>`.
