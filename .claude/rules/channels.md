# Channel rules

A channel is one Markdown file in `docs/dev/channels/`, between two named sides. It is how agents on both sides communicate without either editing the other's code.

## Header

Every channel file starts with a table naming the two sides, the files it covers, and who decides what. After that, messages only.

## Messages

Append, never overwrite. One message is:

```
## <n>. <side> to <side> · <YYYY-MM-DD> · commit <sha or "unpushed">
Status: open | answered | decided: <one line> | withdrawn
Decision needed from: <side> | maintainer | none
```

Then the message body: findings with evidence, what changed, what is asked. Keep it public-safe, like everything committed: no unpublished data, no credentials, no machine-specific paths.

- Reply in the same file, below the message, with an incremented number.
- The side that asked sets `Status:` in its own message; the side that answered sets `answered` on the question it addresses, or `decided: ...` once the maintainer has ruled.
- Never delete another side's message. Trim only your own superseded text, and say what you removed.
- Code changes belong in commits, not in the channel: cite the commit and describe the effect in one or two lines.
- A decision that changes behaviour is recorded in [DECISIONS.md](../DECISIONS.md), citing the channel and the message number. The channel keeps the discussion; DECISIONS keeps the ruling.

## At the start of a session

List `docs/dev/channels/`. Read every channel whose `Covers` overlaps the files you are about to change, and every message addressed to your side that is still `open`.

## Closing

The channel is closed by the project maintainer when the collaboration ends: move the decided points to `DECISIONS.md`, delete the channel file, and drop its row from the README. Git history keeps the text.
