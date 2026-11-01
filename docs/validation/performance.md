# Thread scaling

How much each step gains from more threads, on the public *A. thaliana* and *A. lyrata* pair: 22,146 pairs, one comparison, with InterPro and GO terms on the query only. Code `59faf30`, image `dndsr_v0.2.0_sha256-726657e9`. Each row is one run in a fresh process under Slurm, with `dndsr <step> -c comparison.txt -t <threads>`, so times are not averaged. Memory is Slurm's `MaxRSS` for that step.

| Step | Time by threads | Peak memory | What limits it |
|---|---|---|---|
| `extract_cds` | 71 s at 1, 43 s at 8 | Not recorded | One worker per genome, so at most 2 per comparison |
| `calculate_dnds` | 36 min at 8, 18 min at 16, 12 min at 32 | 15 GB at 16, 29 GB at 32 | Scales well |
| `append_annotations` | 15.4 min at 1, 17.2 min at 8 | 17 GB at 8 | **No gain.** 87% of the time is GFF parsing before the threaded part |
| `ipr_enrichment` | 5.3 min at 1, 4.6 min at 8 | 0.5 GB | One worker per side and comparison. 92% of the time is one per-row term filter |
| `go_enrichment` | 37 s at 1, 41 s at 8 | 1.2 GB | One worker per side and comparison |
| `term_enrichment` | 3 s at 1 and at 8 | Small | Fast already |
| `dnds_summary` | 2 s; no threads option | Small | Fast already |
| `dnds_ideogram` | 5 min, of which about 4 min is PNG conversion; no threads option | 1 GB without PNG | `RIdeogram::convertSVG` |

**Thread count did not change results.** `AtAl_dnds.tsv` is byte-identical at 8, 16, and 32 cores, and the annotated table and all enrichment tables are identical at 1 and 8 threads. `tests/testthat/test-threads.R` checks this on the example.

**Memory.** `calculate_dnds` needs about 1 GB per core when run from the command line. Calling it from an R session after other steps is riskier, because every worker is forked from that session: at 32 cores, in a session that had already run the whole workflow, it was killed at the 96 GB limit, and the error was only `error reading from connection`. Start a fresh session, or use fewer cores, for large runs from R.

## Where the gains are

Ordered by measured time saved on one comparison at 32 CPUs, where the times are `append_annotations` 17 min, `calculate_dnds` 12 min, `dnds_ideogram` 5 min, and `ipr_enrichment` 5 min.

1. **`append_annotations`:** parse GFF lines in parallel, and each genome separately.
2. **`ipr_enrichment`:** the term filter repeats a named-vector lookup for each gene. Filtering once, without that lookup, may help more than threads. This has not been tried.
3. **`dnds_ideogram`:** run sides side by side, and make PNG conversion optional in the CLI docs.
4. **`extract_cds` and the enrichment commands:** also split work within one genome, by ontology or term type.

None of this should change a number, and each change must pass the identical-output check. These are not done: `R/` is frozen while `dev` is the source of the v0.3.0 pull request.
