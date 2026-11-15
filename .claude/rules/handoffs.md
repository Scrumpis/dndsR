# Handoffs between agents and contributors

Some areas have an owner other than you. The McDonald-Kreitman code (`R/mk.R`, `tests/testthat/test-mk.R`, `docs/mcdonald-kreitman.md`) belongs to its author, the MK maintainer. A handoff is how a finding about someone else's area reaches its owner without a second agent editing that code.

`docs/dev/handoffs/` holds **open handoffs only**: one Markdown file each, starting with this header.

```
Status: open
To: <owner role, for example "MK maintainer">
From: <who wrote it>, at <commit>
Covers: <files or functions>
Decision needed from: <owner, project maintainer, or both>
```

## At the start of a session

- List `docs/dev/handoffs/`, ignoring `README.md`. Read every handoff addressed to you, and every one whose `Covers` overlaps the files you are about to change.
- A handoff gives evidence and options, not decisions. Do not implement its recommendations until the owner or project maintainer has decided, and do not re-review it unless asked.
- If a handoff is not addressed to you, do not change the files it covers.

## Writing one

- Put findings about another owner's area in a handoff, not in their code. Give the evidence and how to reproduce it, why it matters, the options and their tradeoffs, and a recommendation labelled as one.
- Use roles, not names. Keep it public-safe, like everything committed.

## Closing one

Only the project maintainer closes a handoff, or the owner with the maintainer's decision. Once it is decided, in one commit:

1. Record the decision in `docs/dev/DECISIONS.md`, citing the handoff's file name and the commit where it was last present.
2. Delete the handoff file. Git history keeps it, and deleting it is the documented reason this rule gives. Nothing decided stays in the folder, so no agent reviews it again.

A handoff that is only partly decided stays open, with the decided points marked `Decided: <summary>` at the top of their section.
