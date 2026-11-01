# Annotation leak on real data

## Setup

`Calbum_BvC` and `CalbumUk_BvC` are subgenome comparisons from two assemblies whose annotations share 80,767 gene IDs. Both were run in one `append_annotations` batch, using the same existing dN/dS tables and read-only inputs, with four settings: old code (`c8b534f`) and current code (`fedd59e`), each at 1 and 8 threads.

## Rows differing from current code at one thread

| Comparison | Version, threads | Rows | Query start | Query InterPro | Query GO | Subject GO |
|---|---|---:|---:|---:|---:|---:|
| `Calbum_BvC` (first in batch) | Old, 1 | 15,992 | 0 | 0 | 0 | 0 |
| | Old, 8 | 15,992 | 0 | 0 | 0 | 0 |
| | Current, 8 | 15,992 | 0 | 0 | 0 | 0 |
| `CalbumUk_BvC` (second in batch) | **Old, 1** | 15,777 | **7,915** | **7,471** | **6,299** | **6,598** |
| | Old, 8 | 15,777 | 0 | 0 | 0 | 0 |
| | Current, 8 | 15,777 | 0 | 0 | 0 | 0 |

## Conclusions

- **Confirmed in old code, at one thread:** about half of the second comparison's rows took coordinates and terms from the first assembly. Results depended on thread count.
- **Current code gives the same result at 1 and 8 threads.**
- **Old code at 8 threads did not show the leak,** because each side runs in its own process.
- **GO and InterPro extraction are otherwise unchanged:** apart from the leak, current and old code agree on every row, which also confirms the `%2C` regression is fixed.
- **The maintainer's existing annotation tables were not affected by this leak:** coordinates and InterPro terms in `Calbum_BvC` (10,578 rows) and `CalbumUk_BvC` (10,763 rows) match fresh results on every row. They are still affected by the CDS bug, and they also hold fewer rows than the dN/dS tables used here, whose origin was not investigated.

Aggregate table: `annotation_leak.tsv`.
