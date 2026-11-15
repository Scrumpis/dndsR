# Agent instructions

Project instructions for all agents are in [CLAUDE.md](CLAUDE.md) and [.claude/rules/](.claude/rules/). Read them before changing anything.

Before you start, sync `dev` with GitHub as [git-policy](.claude/rules/git-policy.md) describes (use the remote that works on this machine, and merge rather than rebase commits that a frozen study or document cites), so work and messages that landed since your last session are in front of you. Then run `tools/dev/setup.sh`, which reinstalls the git hooks from the synced sources. Ownership of each area is in [CLAUDE.md](CLAUDE.md). Then read [docs/dev/channels/](docs/dev/channels/) and [docs/dev/handoffs/](docs/dev/handoffs/): they hold open findings and messages addressed to the owners of particular areas. Read any addressed to you or covering files you will change ([channels rules](.claude/rules/channels.md), [handoff rules](.claude/rules/handoffs.md)).

The end-user skill is in [skills/dndsr/](skills/dndsr/SKILL.md).
