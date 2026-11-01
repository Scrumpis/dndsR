# Draft notices: not published

For the maintainer to review and post. Nothing here has been published, and no one has been contacted.

## For the v0.1.0 and v0.2.0 release pages

> **Warning: do not use this release for new analyses. Rerun earlier results with 0.3.0 or later.**
>
> This version assembles coding sequences incorrectly for multi-exon genes on the minus strand and for genes with several isoforms, which can give false signals of positive selection. On tested genomes, every such gene was affected, about half of orthologous pairs changed, and most genes reported with dN/dS above 1 were artifacts. Enrichment results changed almost entirely.
>
> Details and how to rerun: `docs/migrating-to-0.3.md` and `docs/validation/real-genomes.md`.

## v0.3.0 release notes

> **dndsR 0.3.0 corrects serious errors. Please rerun analyses made with earlier versions.**
>
> The main fix: `extract_cds` could assemble exons in the wrong order or combine isoforms. See NEWS.md for the full list, `docs/migrating-to-0.3.md` to rerun, and `docs/validation/` for the evidence, which includes comparisons with gffread, PAML, and KaKs_Calculator on real genomes.
>
> New in this release: `dndsr example` and `dndsr verify_example` to check an installation; QC tables and provenance records; tutorials; and a portable Agent Skill.
>
> Known limitations: `docs/validation/known-issues.md`.

## For known users or collaborators

> Subject: Please rerun dndsR analyses: errors found in versions before 0.3.0
>
> While validating dndsR, I found errors in how earlier versions assembled coding sequences. They affect most genomes: multi-exon genes on the minus strand, and genes with several isoforms. In the genomes tested, most genes reported with dN/dS above 1 were artifacts, and enrichment results changed almost entirely.
>
> If you used dndsR, please do not rely on those results, and rerun with version 0.3.0: [link]. If any of them are in a manuscript, talk, or dataset, I'm happy to help work out what changes. I'm sorry for the trouble.

## For the maintainer's own work

Existing analyses under the Chenopodium project were produced with the affected code: its *C. formosanum* table has the same dN/dS columns as a fresh run of the old code. Figures, tables, and enrichment summaries built on them should be regenerated before they are used or submitted, including anything already shared.
