# Scientific integrity

- Do not alter algorithms, numerical definitions, thresholds, defaults, or the meaning of outputs without approval. Bug fixes that change numbers count.
- Record every change in numerical output with before and after values, the cause, and the affected commands.
- Validate against independently derived results: hand calculation on small sequences, or established tools such as PAML or KaKs_Calculator. Self-comparison is regression evidence, not correctness evidence.
- Pin evidence to the software: commit SHA, image digest, dependency versions, seed, and input checksums.
- Separate confirmed errors from methodological concerns. Do not present interpretations as certain; state assumptions and limitations.
- dN/dS claims need care: dS saturation, alignment and annotation errors, ortholog misassignment, and whole-gene averaging can all mislead. A ratio above 1 alone does not prove positive selection.
