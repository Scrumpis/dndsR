# Channels between agents and collaborators

A **channel** is a two-way, persistent thread between the agents of two people who own different parts of this repository. It is not a handoff: a handoff carries findings one way and is deleted when decided ([handoffs](../handoffs/README.md), [rules](../../../.claude/rules/handoffs.md)). A channel stays, so both sides can answer, quote, and track what changed.

Why it exists: when someone outside the repository's ownership boundary contributes (for example a collaborator's agent working on `R/mk.R`), findings and answers otherwise pass through the maintainer, in order, by hand. A channel puts them in one file both sides can read at the start of a session and write to at the end.

Rules: [.claude/rules/channels.md](../../../.claude/rules/channels.md).

## Open channels

| File | Between | Covers | Started |
|---|---|---|---|
| [mcdonald-kreitman.md](mcdonald-kreitman.md) | The dndsR integration agent and the MK maintainer's agents | MK methods (`R/mk.R`, `tests/testthat/test-mk.R`, `docs/mcdonald-kreitman.md`, MK sections of tutorials, skill, methods); the MK columns of the comparison file are shared | 2026-11-14 (US Pacific) |
