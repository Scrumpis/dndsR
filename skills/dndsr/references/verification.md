# Verification record

Only checks actually run are listed. Rerun and update these entries when the relevant software changes.

## 2026-10-31: installation and CLI discovery

Environment: Linux amd64 HPC node, singularity-ce 4.1.2, no Docker. Code: `dev` at `c8b534f` with bootstrap infrastructure added. The installer ran with a temporary `HOME`.

| Check | Result |
|---|---|
| Docker Hub lists tags `v0.2.0` and `latest` with one digest | `sha256:726657e929b0a026dc5b81fc3576c56e6cee98d30ec4e6f94992fa2a11c55831` |
| Pull by digest | Succeeded, 1.3 GB SIF |
| Image contents | R 4.5.2, DIAMOND 2.1.14, OrthoFinder 2.5.4, orthologr 0.4.2, pwalign 1.6.0, topGO 2.62.0, IHW 1.38.0, qvalue 2.42.0, RIdeogram 0.2.2 |
| dndsR in image | Absent, as expected from the Dockerfile |
| `tools/dndsr-install install --engine singularity --sif ...` | Exit 0 |
| `tools/dndsr-install doctor` | Reports the engine, SIF, library, and shim |
| `dndsr --help`, `dndsr __commands` | 11 commands listed |
| `dndsr calculate_dnds --help` | Prints R documentation |
| `dndsr dnds_contrasts --help` | Prints `Unknown command` and usage |

## 2026-10-31: skill discovery

| Check | Result |
|---|---|
| Claude Code 2.1.316, new session in the clone | Reports the `dndsr` skill as available |
| Codex CLI 0.165.0, `codex exec` in the clone | Finds `.agents/skills/dndsr/SKILL.md` |

These show discovery only. No agent has yet completed a workflow using the skill.

## Not verified

Any analysis run to completion; all outputs; Docker; Apptainer; native installation; macOS; offline use; correctness of any numerical result.
