# Human-native and agent-native development

Every workflow must be easy to discover, execute, understand, reproduce, and automate, for a person without AI and for an agent without undocumented knowledge. Full policy: `docs/dev/AGENT_NATIVE.md`.

For every relevant software change:

1. Implement it.
2. Test its behavior.
3. Record the successful commands, inputs, outputs, and relevant failures.
4. Update `skills/dndsr/` and `docs/tutorials/`, using the same commands in both.
5. Add or update automated validation where practical.
6. Run the documented workflow and confirm it works.
7. Commit software, tests, skill, and docs together on `dev`.

Rules for the skill:

- `skills/dndsr/` is the only source. `.claude/skills/dndsr` and `.agents/skills/dndsr` are symlinks to it.
- Document only what has been verified, and record it in `skills/dndsr/references/verification.md`. Mark anything else as unverified, or leave it out.
- Use immutable image digests in commands. Do not assume Docker, internet access, specific paths, or elevated privileges.
- Keep `metadata.software-version` in `SKILL.md` equal to `Version` in `DESCRIPTION`; review the skill before changing it.
