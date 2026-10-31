# dndsR

dndsR is a proper R package, callable from R and intended to be fit for CRAN one day, and also a containerized command-line tool. It serves people and AI agents equally, with ease of use as a primary goal. It does pairwise dN/dS calculation, annotation, enrichment, and visualization. Roadmap: audit `dev`, merge to `main`, release v0.3.0, then v1.0.0 for publication. Current state: [docs/dev/STATUS.md](docs/dev/STATUS.md). How to inspect and resume work: [docs/dev/RUNBOOK.md](docs/dev/RUNBOOK.md).

## Repository facts

Verified by reading the code on `dev`:

- `R/<command>.R` holds exported analysis functions. `R/cli_<command>.R` holds thin wrappers whose `target` attribute names the function.
- `R/cli_main.R` finds commands by looking for `cli_*` functions in the namespace. `R/cli_parse.R` parses arguments. `inst/scripts/dndsr` calls `dndsR::cli_main()`.
- `calculate_dnds()` wraps `orthologr::dNdS()`. Its defaults are `aligner = "diamond"`, `sensitivity_mode = "fast"`, and `dnds_method = "Comeron"`.
- `tools/Dockerfile` builds a runtime image with dependencies only. dndsR is not in the image: `tools/dndsr-install` installs it from the clone into a host R library and creates a `dndsr` shim that runs it in the container.
- The image needs a patched orthologr (`inst/patches/orthologr-pwalign.patch`).
- There is no `tests/` directory, and `.github/workflows/R-CMD-check.yaml` runs only for `main`/`master`.

Everything else, including whether the science is correct, is unaudited. Say which claims you verified and which are assumptions.

## Two kinds of agent

- Development agents (`.claude/agents/`) audit and change this repository. They never ship to users.
- The end-user skill (`skills/dndsr/`) teaches external AI agents to run dndsR. It is a first-class deliverable and is released with the software.

## Non-negotiable rules

Details are in `.claude/rules/`.

- Commit only to `dev`. **Agents do not push**: finish with the exact command for the maintainer, normally `git push origin dev`.
- Never touch `main`, tags, releases, GitHub settings, or published images without explicit authorization for that specific action. Never rewrite shared history.
- Never silently change an algorithm, numerical definition, default, output meaning, or public interface. Such changes need approval, a regression test, and documentation.
- A software change is incomplete until the skill and tutorials have been reviewed and updated.
- A feature is complete only when a human following the docs and a suitably equipped agent can both use it.
- Prefer the simplest solution. New dependencies and new infrastructure need a reason.

## Development

- Set up a clone once: `tools/dev/setup.sh`.
- Keep dependencies declared in `DESCRIPTION`, and in `tools/Dockerfile` when the image needs them. The HPC default R lacks many dependencies, so run package code in the pinned container (see the runbook).
- Regenerate `man/` and `NAMESPACE` with roxygen2 rather than editing them by hand. roxygen2 is not in the current image.
- Before committing, run what applies: `tools/dev/test-safeguards.sh`, `Rscript --vanilla tools/dev/validate-skill.R`, and package tests or `R CMD check`.

## Evidence

- Report what was actually run, with its command and result. Never describe an unrun check as passing, and state failures plainly.
- Scientific tests need independently derived expected values, not values copied from dndsR's own output.
- Record decisions and open questions in `docs/dev/DECISIONS.md`, and end each substantive session by updating `docs/dev/STATUS.md` and `docs/dev/SESSION_LOG.md`.

## Files

- Never commit transcripts, credentials, private paths, containers, temporary output, or bulky results. `tools/dev/check-public-safe.sh` checks staged files.
- Put temporary files outside the repository. Inspect a target before deleting or overwriting it.
- Do not delete code, docs, tests, or results without a documented reason and explicit approval.

## Stop and ask a human before

Changing scientific behavior or a public interface; adding a major dependency; deleting tracked content; any action on `main`, tags, releases, image registries, Zenodo, or GitHub settings; resolving a scientific question where the evidence is ambiguous.
