---
name: dndsr
description: Installs and runs dndsR, an R package and CLI for pairwise dN/dS analysis of orthologous genes from genome FASTA and GFF3 files, with annotation, term enrichment, ideograms, and regional or contrast summaries. Use when a user wants to calculate, summarize, or interpret dN/dS ratios with dndsR, or set it up with Singularity, Apptainer, Docker, or R.
license: MIT
compatibility: Needs bash and git, plus Singularity, Apptainer, or Docker, or an R installation with dndsR's dependencies. Only Linux amd64 with Singularity has been tested. Cloning and pulling the image need network access unless both are already present.
metadata:
  software-version: "0.3.0"
  skill-status: "core workflow tested on a designed example; not validated on real genomes"
---

# dndsR

dndsR estimates and analyzes dN/dS between orthologous genes in pairwise genome comparisons. It is both an R package and a command-line tool.

What has been tested is recorded in [references/verification.md](references/verification.md). Tell the user when you rely on something that has not been tested. Open problems are in [references/gaps.md](references/gaps.md).

## 1. Check the environment

```
command -v apptainer singularity docker Rscript
command -v dndsr && dndsr --help
```

If `dndsr --help` lists commands, go to step 3. Without containers, use the tested conda route in [INSTALL.md](../../INSTALL.md). If R has dndsR (`Rscript -e 'packageVersion("dndsR")'`), every command below is also a function: see [Using dndsR from R](../../docs/tutorials/using-r.md).

## 2. Install

Tell the user first: this downloads 1.3 GB, and writes to `~/.config/dndsr`, `~/.dndsr/rlib`, and `~/.local/bin`. The image holds dependencies; dndsR is installed from the clone, which must stay in place. Replace `singularity` with `apptainer` if that is what is available (untested).

```
git clone --branch dev https://github.com/Scrumpis/dndsR
cd dndsR
singularity pull dndsr.sif docker://scrumpis/dndsr@sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831
tools/dndsr-install install --engine singularity --sif dndsr.sif
export PATH="$HOME/.local/bin:$PATH"
```

Pin digests rather than tags, because tags can move. To see published tags: `curl -s "https://hub.docker.com/v2/repositories/scrumpis/dndsr/tags?page_size=20"`. The Docker route is in [INSTALL.md](../../INSTALL.md) and is untested.

## 3. Verify with the example

Always do this after installing. It takes about a minute.

```
dndsr example
cd dndsr_example
dndsr extract_cds -c comparison_file.txt
dndsr calculate_dnds -c comparison_file.txt -t 2
dndsr verify_example
```

Success ends with `All 14 example checks passed.` Any `[FAIL]` line, or a nonzero exit code, means the installation does not behave as expected: report it rather than continuing. This confirms the installation works. It does not validate results on other data.

## 4. Run a real analysis

Input is a comparison file, one line per comparison: `name query.fasta query.gff3 subject.fasta subject.gff3`. Paths are relative to the current folder, which is visible in the container. For data elsewhere, reinstall with `--bind <folder>` and use absolute paths.

```
dndsr extract_cds -c comparison_file.txt -t 8
dndsr calculate_dnds -c comparison_file.txt -t 40
dndsr append_annotations -c comparison_file.txt -t 8
dndsr term_enrichment -c comparison_file.txt -t 8
dndsr ipr_enrichment -c comparison_file.txt -t 8
dndsr go_enrichment -c comparison_file.txt -t 8
dndsr dnds_ideogram -c comparison_file.txt
dndsr dnds_summary -c comparison_file.txt --regions-bed regions.bed
```

`calculate_dnds` can take many hours on whole genomes, so on a cluster use the scheduler. Commands and options: [references/commands.md](references/commands.md). Output columns: [docs/output-files.md](../../docs/output-files.md).

Check `<comparison>/<genome>_CDS_qc.tsv` for exclusions, and the number of pairs in `<comparison>_dnds.tsv`, before interpreting results.

## 5. Interpret honestly

Read [docs/interpretation.md](../../docs/interpretation.md) before explaining results. The essentials:

- **dN/dS above 1 is a candidate, not proof of positive selection.** dndsR applies a threshold to a point estimate. In simulations of 600-codon neutral genes, up to 35% exceeded 1 with the default method.
- Look at dN and dS as well as the ratio. Small dS makes ratios unstable, and large dS (above roughly 1 to 2) suggests saturation. dS = 0 gives NA.
- The default Comeron method tends to give lower ratios than true values in simulation. Evidence: [docs/validation/estimators.md](../../docs/validation/estimators.md).
- Reciprocal best hits can mispair paralogs, especially in polyploids. Annotation errors can look like selection.
- Enrichment tests treat genes as independent.
- `dnds_state_contrast --level gene` is experimental. Do not base conclusions on it.
- Results from dndsR before 0.3.0 should not be reused. In tested genomes, most of their dN/dS calls above 1 were artifacts. See [docs/migrating-to-0.3.md](../../docs/migrating-to-0.3.md).

Report the version, commit, image digest, command, and the number of excluded genes. Parameters and input hashes are in `<comparison>_dnds_provenance.json`.

## 6. Troubleshooting

- `command not found: dndsr`: add `~/.local/bin` to `PATH`, then run `hash -r`.
- `Unknown command` or `Unknown option`: dndsR stops on typos. Run `dndsr --help` or `dndsr <command> --help`.
- `Sequence files not found`: run `extract_cds` first, with the same `-o`.
- `jobs failed`: the message gives each error, and logs are in the comparison folder. Retry with `-t 1`, and check memory if a worker died.
- Enrichment writes no table: a small dataset may lose terms to `--max-prop` (default 0.2), `--min-pos`, or `--min-total`. Read the log before relaxing them.
- `dndsr doctor` and `tools/dndsr-install doctor` check the installation.

## Citation

No DOI exists yet. Use [CITATION.cff](../../CITATION.cff) with the exact commit, and also cite orthologr, DIAMOND, and Comeron (1995).
