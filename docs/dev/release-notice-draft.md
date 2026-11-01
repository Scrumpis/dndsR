# Correctness notices: for the maintainer's approval

Nothing here has been published, and no one has been contacted. These use only public evidence: synthetic data and *A. thaliana*.

## Proposed locations

Historical releases and tags are not changed.

| Where | What | Why |
|---|---|---|
| **GitHub issue, pinned**, titled as below | Full notice | Visible on the repository without changing a release, and linkable |
| **Pull request to `main`: a README banner** | Short notice, linking to the issue | The README on the default branch is what visitors see. Needs your authorization to go to `main` |
| **Docker Hub repository description** for `scrumpis/dndsr` | Short notice | People using the images may never open GitHub |
| Direct message to known users | Email text | Only you know who they are |

Adding a line to the v0.1.0 and v0.2.0 release pages would also reach people, but it changes those pages, so it is left out unless you want it. The tags themselves would be unchanged.

## Full notice

**Title: Correctness notice: dN/dS results from dndsR v0.1.0 and v0.2.0 can be wrong, please recompute**

dndsR v0.1.0 and v0.2.0, and commits on `main` and `dev` before the fixes in the upcoming v0.3.0, have defects that can produce incorrect results. They run without errors, and their output looks plausible.

**Defects**

1. **Coding sequences assembled incorrectly.** By default, `extract_cds` joined multi-exon genes on the minus strand in the wrong exon order, and joined all isoforms of a gene into one chimeric sequence. In *A. thaliana* TAIR10.1, all 5,655 such minus-strand genes and all 10,695 multi-isoform genes were wrong, as judged against gffread.
2. **Annotations from the wrong genome.** At one thread, a gene ID shared with a genome read earlier in the same run took that genome's terms and coordinates.
3. **Missing annotations.** With RefSeq-style IDs such as `gene-AT1G01010`, `append_annotations` and `dnds_ideogram` found no annotations, without an error.
4. **Lost region.** A regions file with no header lost its first region in `dnds_summary` and `dnds_contrast`.
5. **Enrichment errors.** These affected the InterPro parent-child method and the handling of terms with no known type. A filter that failed could also run as if unfiltered.

**What is affected**

- **dN and dS, and the dN/dS ratio** for any pair involving a multi-exon minus-strand gene or a multi-isoform gene. Out-of-frame sequence can produce false dN/dS above 1. In a designed test, a gene with a true dN/dS of about 0.11 was reported as 1.6. Other effects include high dS, undefined ratios, and missed orthologs.
- **Everything computed from those values:** enrichment, ideograms, summaries, and contrasts, even when other genes in the same table were correct.
- Single-exon genes, and plus-strand genes with one isoform, had correct sequences.

**What to do**

Do not rely on results from these versions, and please recompute them. A corrected version, v0.3.0, is in validation and is not yet released. The fixes are on the `dev` branch. Rerunning instructions: `docs/migrating-to-0.3.md`. Evidence: `docs/validation/`.

If these results appear in a manuscript, talk, or dataset, please reassess those conclusions once recomputed. Questions are welcome here.

## Short notice

> **Correctness notice:** results from dndsR v0.1.0 and v0.2.0 can be wrong, including false signals of positive selection, because of errors in coding-sequence assembly. Please recompute them with a corrected version: [link to issue].

## Email to known users

> Subject: Please recompute dndsR results: errors found in v0.1.0 and v0.2.0
>
> While validating dndsR, I found errors in how v0.1.0 and v0.2.0 assembled coding sequences. They affect multi-exon genes on the minus strand and genes with several isoforms, so they apply to nearly every genome, and they can produce false signals of positive selection. Details: [link to issue].
>
> If you used dndsR, please do not rely on those results. A corrected v0.3.0 is being validated, and I will let you know when it is released. If any of your results are in a manuscript, talk, or dataset, I'm happy to help work out what changes. I'm sorry for the trouble.

## v0.3.0 release notes, when released

> **dndsR 0.3.0 corrects serious errors. Please rerun analyses from earlier versions.** See NEWS.md for the full list, `docs/migrating-to-0.3.md` for how to rerun, and `docs/validation/` for evidence, including comparison with gffread, PAML, and KaKs_Calculator.
>
> New: `dndsr example` and `dndsr verify_example` to check an installation, QC tables and provenance records, tutorials, and a portable Agent Skill.
>
> Gene-level `dnds_state_contrast` is disabled, because it did not give meaningful results. Known limitations: `docs/validation/known-issues.md`.
