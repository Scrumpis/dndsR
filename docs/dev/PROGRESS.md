# Progress log

Durable checkpoints for resuming autonomous work. Newest first. Evidence is in `docs/validation/`. The plan is at the end.

## Environment

- Runtime files live outside the clone, in a `dndsr-runtime` folder beside it: image, `SHA256SUMS`, R library, reference tools (PAML 4.10.10 and KaKs_Calculator 2.0.1 in a micromamba env). Recreate them as described in `RUNBOOK.md` if purged.
- Scratch is purged after 45 days, and the maintainer may recursively touch files. **File dates mean nothing here: identify everything by content hash.**
- The maintainer pushes. Agents report the push command.

## Checkpoints

### Exon order, isoforms, and QC
Fixed the CDS assembly bug (`docs/validation/cds-extraction.md`), added the designed example, independent NG86 reference, first tests, input hashes, and hash-based audit archiving.
