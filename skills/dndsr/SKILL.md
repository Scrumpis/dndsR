---
name: dndsr
description: Installs and runs dndsR, an R package and CLI for pairwise dN/dS analysis of orthologous genes from genome FASTA and GFF3 files, with annotation, term enrichment, ideograms, and regional or contrast summaries. Use when a user wants to calculate, summarize, or interpret dN/dS ratios with dndsR, or set it up with Singularity, Apptainer, or Docker.
license: MIT
compatibility: Needs bash, git, and Singularity, Apptainer, or Docker. The Linux amd64 image is the only environment tested. Pulling the image and cloning the repository need network access unless both are already present.
metadata:
  software-version: "0.3.0.000"
  skill-status: "initial; installation and CLI discovery verified, analyses not yet verified"
---

# dndsR

dndsR estimates and analyzes dN/dS between orthologous genes in pairwise genome comparisons. Its main output, the dN/dS ratio, compares nonsynonymous and synonymous substitution rates.

Everything marked **unverified** below comes from the repository's documentation or code and has not been run to completion. Tell the user when you rely on it. Evidence for verified steps is in [references/verification.md](references/verification.md). Known gaps are in [references/gaps.md](references/gaps.md).

## 1. Check the environment

```
command -v apptainer singularity docker
command -v dndsr && dndsr --help
```

If `dndsr --help` lists commands, skip to step 3. If no container engine exists, see "Native installation" in section 5.

## 2. Install

The runtime image holds dependencies only. dndsR itself is installed from a git clone into an R library in the home directory, so the clone must stay in place.

```
git clone https://github.com/Scrumpis/dndsR
cd dndsR
```

The published image is on Docker Hub as `scrumpis/dndsr`, tag `v0.2.0`, digest `sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831`. List current tags with:

```
curl -s "https://hub.docker.com/v2/repositories/scrumpis/dndsr/tags?page_size=20"
```

Pin the digest for reproducibility, because tags can move. Apptainer works the same way when you replace `singularity` with `apptainer`, but only Singularity has been tested.

```
singularity pull dndsr.sif docker://scrumpis/dndsr@sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831
tools/dndsr-install install --engine singularity --sif dndsr.sif
```

The installer writes to `~/.config/dndsr`, `~/.dndsr/rlib`, and `~/.local/bin`, adds shell completions, and may edit `~/.zshrc`. Tell the user before running it. If `~/.local/bin` is not on `PATH`, it prints the line to add.

Docker is **unverified**:

```
docker pull scrumpis/dndsr@sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831
tools/dndsr-install install --engine docker --image scrumpis/dndsr@sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831 --docker-platform linux/amd64
```

Confirm the installation. Expect 11 commands:

```
tools/dndsr-install doctor
dndsr --help
```

## 3. Run commands

Use `dndsr <command> --help` for the arguments of any command. The commands, with their README-documented purpose, are in [references/commands.md](references/commands.md).

The documented order is below. **All of these are unverified.** Most accept a whitespace-delimited comparison file with `comparison_name query_fasta query_gff subject_fasta subject_gff` on each line. See the references for the exact format.

```
dndsr split_comparisons -c comparison_file.txt
dndsr extract_cds -c comparison_file.txt -t 8
dndsr calculate_dnds -c comparison_file.txt -t 80
dndsr append_annotations -c comparison_file.txt -t 8
dndsr ipr_enrichment -c comparison_file.txt -t 8
dndsr go_enrichment -c comparison_file.txt -t 8
dndsr term_enrichment -c comparison_file.txt -t 8
dndsr dnds_ideogram -c comparison_file.txt
```

`calculate_dnds` can be slow. On a cluster, submit it through the scheduler. Paths in the comparison file must be visible inside the container: with Docker the README asks for paths relative to the output directory, and with Singularity you may need `DNDSR_BINDS` in `~/.config/dndsr/config.env`.

## 4. Interpret results carefully

No correctness audit has been completed, and no expected-result example exists yet. Do not present dndsR results as validated.

From the code, `calculate_dnds` wraps `orthologr::dNdS()` with DIAMOND best hits and the Comeron method by default.

General dN/dS background, not specific to dndsR:

- A ratio below 1 is consistent with purifying selection, near 1 with neutral evolution, and above 1 with positive selection. These readings are not proof: a whole-gene ratio averages over sites and lineages.
- Very small or saturated dS makes ratios unstable. Check dN and dS, not only the ratio.
- Misassigned orthologs, poor alignments, and annotation or frame errors can create false signals. Best-hit pairing is a heuristic, especially in polyploids.
- Enrichment comparisons need matching annotation versions between genomes, as the README warns.

Report the dndsR commit, image digest, command, and parameters with any result.

## 5. Native installation and troubleshooting

Native installation is **unverified**. It needs the dependencies in `DESCRIPTION`, DIAMOND, OrthoFinder, and the patched orthologr. `INSTALL.md` is partial and describes itself as temporary. Prefer the container.

- `command not found: dndsr`: add `~/.local/bin` to `PATH`, then run `hash -r`.
- `Unknown command` or `Unknown option`: dndsR stops rather than ignoring a typo. Run `dndsr --help` or `dndsr <command> --help`.
- `jobs failed`: a parallel job failed. The message includes each error; logs are in the comparison folder. Retry with `--threads 1` if the cause is unclear, and check memory if a worker died.
- The image tag is older than this code. Report that mismatch if a command fails on a missing dependency.

## Citation

No DOI or formal citation exists yet. Cite the repository URL and exact commit: https://github.com/Scrumpis/dndsR
