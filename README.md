# dndsR
<!-- badges: start -->
[![R-CMD-check](https://github.com/Scrumpis/dndsR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Scrumpis/dndsR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> **Correctness notice:** v0.1.0 and v0.2.0 were unvalidated pre-releases with defects that can give wrong dN/dS, including false signals of positive selection. Please recompute any results from them. [What was affected and how to rerun](docs/migrating-to-0.3.md).

**dndsR** calculates and analyzes dN/dS between orthologous genes in pairwise genome comparisons, starting from genome FASTA and GFF3 files. It is an R package you can call from R, and a command-line tool that runs in a container. It is meant to be easy to use, for researchers and for AI agents alike.

| Step | Command | R function |
|---|---|---|
| Split subgenomes or haplotypes (optional) | `dndsr split_comparisons` | `split_comparisons()` |
| Extract coding sequences | `dndsr extract_cds` | `extract_cds()` |
| Calculate dN/dS | `dndsr calculate_dnds` | `calculate_dnds()` |
| Add annotations | `dndsr append_annotations` | `append_annotations()` |
| InterPro, GO, or other term enrichment | `dndsr ipr_enrichment`, `go_enrichment`, `term_enrichment` | `ipr_enrichment()`, `go_enrichment()`, `term_enrichment()` |
| Enrichment adjusted for gene length and dS | `dndsr covariate_enrichment` | `covariate_enrichment()` |
| Genome ideogram | `dndsr dnds_ideogram` | `dnds_ideogram()` |
| Summaries, optionally by region | `dndsr dnds_summary` | `dnds_summary()` |
| Contrast comparisons | `dndsr dnds_contrast`, `dnds_state_contrast` | `dnds_contrast()`, `dnds_state_contrast()` |

## Quick start

On Linux with Singularity or Apptainer (replace `singularity` with `apptainer` if needed). This takes a few minutes and a 1.3 GB download. The commands install into your home directory.

```
git clone --branch v0.3.0 https://github.com/Scrumpis/dndsR
cd dndsR
singularity pull dndsr.sif docker://scrumpis/dndsr@sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831
tools/dndsr-install install --engine singularity --sif dndsr.sif
export PATH="$HOME/.local/bin:$PATH"

dndsr example
cd dndsr_example
dndsr extract_cds -c comparison_file.txt
dndsr calculate_dnds -c comparison_file.txt -t 2
dndsr verify_example
```

The last command checks your results against independently established values and should end with `All 14 example checks passed.`

## Documentation

- [Getting started](docs/tutorials/getting-started.md): install, run the example, and read the results
- [Using dndsR from R](docs/tutorials/using-r.md)
- [Analyzing your own genomes](docs/tutorials/your-own-data.md)
- [Output files](docs/output-files.md)
- [Interpreting results](docs/interpretation.md): read this before drawing conclusions
- [Installation options](INSTALL.md): containers, conda, and R
- [How dndsR computes its results](docs/methods.md)
- [Validation evidence](docs/validation/README.md) and [known issues](docs/validation/known-issues.md)
- [Changes between versions](NEWS.md)

Use `dndsr --help` or `dndsr <command> --help` for every option.

## For AI agents

[skills/dndsr/SKILL.md](skills/dndsr/SKILL.md) is a portable [Agent Skill](https://agentskills.io) with tested commands, interpretation guidance, and limits. Claude Code and Codex find it automatically in this repository. Elsewhere, copy or link `skills/dndsr` to your agent's skills folder.

## Contributing

Bug reports and feature requests are welcome on the [issue tracker](https://github.com/Scrumpis/dndsR/issues). Developers should start with [CLAUDE.md](CLAUDE.md) and the [development runbook](docs/dev/RUNBOOK.md), whether or not they use AI tools.
