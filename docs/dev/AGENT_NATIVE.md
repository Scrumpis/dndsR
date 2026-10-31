# Agent-native development policy

This is a permanent project policy.

## What dndsR is

- **A real R package.** Every analysis is a documented function that can be called from R. Package changes should keep it fit for CRAN one day: clean `R CMD check`, optional dependencies guarded, nothing written outside the paths a user gives.
- **A command-line tool with containers.** The CLI is the recommended route for large analyses, and it must not be the only one.
- **Accessible to people and AI agents.** Ease of use comes first for both.

The R and CLI interfaces have the same capabilities. Tutorials cover both.

dndsR is also a portable Agent Skill that gives compatible agents the knowledge to run it correctly. A person must never need an AI agent, and an agent must never need undocumented project knowledge. A feature is complete only when both can use it.

The per-change checklist is in `.claude/rules/skill-maintenance.md`. No model is trained or fine-tuned; expertise lives in version-controlled instructions, references, scripts, and tests.

## Skill layout

| Path | Role |
|---|---|
| `skills/dndsr/` | Canonical source, following the Agent Skills specification |
| `.claude/skills/dndsr` | Symlink for Claude Code project discovery |
| `.agents/skills/dndsr` | Symlink for Codex and other `.agents/skills` readers |
| `docs/tutorials/` | Human tutorials using the same commands |
| `.claude/agents/` | Internal development agents, not part of the skill |

Claude Code and Codex document symlinked skill directories, so a symlink avoids a second copy. Git on Windows may check symlinks out as plain files; `skills/dndsr/` remains usable directly there. Users outside a clone can copy or link `skills/dndsr` to `~/.claude/skills/dndsr` or `~/.agents/skills/dndsr`.

A later phase can add a Claude Code plugin manifest (`.claude-plugin/`) that packages `skills/`. Phase 0 leaves it out because no release exists to link it to and the marketplace format has not been checked.

## Release linkage

Each release should tie together one git tag, one `software-version` in the skill, one image digest, and eventually one Zenodo DOI. Commands use digests because tags can move. No image, release, or DOI is created until the maintainer authorizes it.

## Validation strategy

| Layer | What | When |
|---|---|---|
| Static | `tools/dev/validate-skill.R`: frontmatter, symlinks, links, commands exist in the CLI, no `:latest`, one digest, version matches `DESCRIPTION` | Pre-push; run by hand |
| Change gate | `commit-msg` hook: software commits must change `skills/` or carry a `Skill-Impact:` trailer | Every commit |
| Evidence | `references/verification.md` records what ran; unrun steps are marked unverified | Each change |
| Workflow | Run skill commands in a clean `HOME` with the pinned image | Before release |
| End to end | Fresh-agent test below | Before release |

The static validator proves consistency only. It does not show that a workflow works or that a result is correct.

## End-to-end agent test (planned)

Requires the validated example dataset that does not exist yet.

1. Start from a new container or temporary `HOME` with only a public clone, a container engine, and an agent CLI. Record the OS, engine, agent version, and network access.
2. Give the agent a plain request, with no project hints: "Use the dndsr skill in this repository to analyze this example pair and explain the result."
3. Check from its transcript that it loaded the skill, chose the version and digest matching the checkout, and installed or accessed the executable.
4. Check that it ran the documented example.
5. Compare the output with expected values established independently of dndsR, with a stated tolerance, using a script rather than the agent's own claim.
6. Check that its explanation covers the meaning of the result and its limitations, using a written rubric.
7. Run it on Claude Code and Codex, and store the results and transcripts locally. Record verified and untested paths in `references/verification.md`.
