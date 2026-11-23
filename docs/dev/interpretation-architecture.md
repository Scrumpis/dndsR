# Evidence-grounded interpretation: minimal architecture, validation plan, and release scope (proposal)

Status: **proposal for the maintainer's review (2026-11-23); nothing implemented.** Deliverables 3-5 of the maintainer's message of 2026-11-23. The aim is more accurate, transparent and appropriately limited conclusions, not more conservative defaults.

## Principle

Three layers, each testable on its own, and the first two complete without any language model:

```
A. Execution          B. Diagnostics                      C. Interpretation (optional)
engines, versioned -> machine-readable evidence record -> report from templates (no LLM)
configurations        (every item: trigger, value,        and, optionally, an agent that
results tables        basis, status)                      explains and drafts text, citing
                                                          only items in A and B and references
```

## A. Statistical execution (exists)

What dndsR already does: pinned engines (HyPhy 2.5.103, PAML 4.10.10, KaKs_Calculator, MAFFT), settings in `_provenance.json`, deterministic tables, seeds and thread counts recorded. Nothing changes here; backward compatibility is kept. The only addition B needs is a stable identifier per run (a hash of the inputs, settings and versions) so later records can point at exactly one result.

## B. Deterministic diagnostics: the evidence record

One JSON file per run, `<run>_evidence.json`, written beside the results; a schema versioned with dndsR. Each item:

| Field | Meaning |
|---|---|
| `id` | stable identifier, for example `busted.fg_window_flagged` |
| `scope` | run, orthogroup, branch or term |
| `value` | the measured quantity (for example the number of flagged windows) |
| `trigger` | the rule that sets the status, stated exactly (for example `dissimilar_windows > 0 on a foreground taxon`) |
| `status` | `known_failure` (validated to cause wrong inferences), `validated_ok` (inside a tested and passed domain), `unvalidated_concern` (plausible, untested), `not_assessable` (the data needed are absent) |
| `basis` | the evidence: a dndsR validation report and table, or a literature reference with its verification status |
| `consequence` | what a user should not conclude, in one sentence, from the matrix ([interpretation-limits-matrix.md](interpretation-limits-matrix.md)) |

Rules: an item exists only if its trigger is computable from dndsR outputs and its basis is recorded; `unvalidated_concern` items are labelled as such and never counted as QC failures; absence of a `known_failure` item never means a result is valid (the record says so in its header); the record never alters a p-value or a call. Most items already exist as columns or warnings (`qc_flags`, `dissimilar_windows`, `alignment_check`, calibration flags, `n_starts`, `omega_at_upper_bound`, enrichment notices); B gathers them with their basis and status.

**The ordinary report (no LLM).** `dndsr report <run>` renders A and B into a Markdown or HTML page with fixed templates: what the test asks, the result, the evidence items that apply to each significant or notable row, and the matrix's "not justified" sentences for that method. Deterministic and testable byte for byte. This is the version every scientist gets.

## C. Optional agentic interpretation

An agent (through the dndsR skill) reads A, B and the methods' reference notes, and writes: plain-language explanations, uncertainties, qualified biological interpretation, draft Results and Discussion wording, and suggested further checks.

Hard constraints, each checkable mechanically:

1. Every factual sentence carries a citation to an A row, a B item `id`, or a reference key; a checker rejects any sentence whose cited value differs from the source (numbers re-read from the files).
2. The agent may not create, remove or restate the status of a B item, change a statistical outcome, or cite its own judgement as evidence.
3. A warning is never written as proof that a result is false; `unvalidated_concern` is never written as a failure.
4. When A and B lack what a claim needs, the text says "not established by these results" rather than guessing.
5. The generated text is stored separately from A and B, with the model name, version, prompt and skill version, and is never written into result tables.

## Is this better than conventional reporting? (critical assessment)

- **What B alone adds** is most of the reliability gain: explicit triggers with evidence, the same for every user. It is cheap, reproducible and testable. It should come first, whatever happens to C.
- **What C could add** is applying the matrix to a specific result: noticing that a "significant" gene has a flagged foreground window and a single short branch, or that a claim of adaptation in a draft overreaches. Whether it does so better than the deterministic report is an empirical question (below), not an assumption.
- **Risks of C:** fluent text that sounds more certain than the evidence; selective emphasis; fabricated references; variation between runs and models. The constraints above reduce but do not remove them, which is why C stays optional and every claim traceable.

## Validation plan (prioritized)

1. **B's correctness** (required before any release use): unit tests that each trigger fires exactly on constructed inputs; that every item's basis resolves to an existing report or reference; regression tests on the validation studies' outputs (for example, every study B corrupted gene gets `busted.fg_window_flagged`; realigned-indel genes do not, and the record says the flag cannot detect them).
2. **The deterministic report**: snapshot tests, and a user test with scientists not involved in development reading reports for clean, compromised and ambiguous analyses; outcome: wrong conclusions drawn, against the same results without the report.
3. **C, if pursued**: the benchmark in [methods-paper-concept.md](methods-paper-concept.md): known-truth synthetic analyses (clean positives and negatives, corrupted alignments, realigned indels, confounded enrichment, underpowered tests), adversarial cases where a significant result invites an unjustified biological claim, and expert-reviewed realistic examples; primary endpoint the rate of unsupported conclusions, secondary the loss of valid conclusions; comparison with B's report alone; development and evaluation cases separate; several runs per case and several models.

Demonstrating reduced misleading conclusions requires comparing readers (people or agents) given A only, A+B report, and A+B+C, on cases with known truth, with conclusions scored blind against a rubric fixed in advance.

## Release scope

| For the next release (essential) | Optional, future |
|---|---|
| Finish the running validation studies and act on them | The agentic layer C |
| Correct the masking-power wording (forensics) and the BUSTED-E description | The methods paper |
| The limitations matrix reflected in docs and skill ("not justified" sentences per method) | `dndsr report` HTML rendering |
| An evidence record (B) for the phylogenetic commands, built from existing flags, with tests | Evidence records for all commands |
| | Reader studies |

B can be scoped small for the release: the existing flags and warnings with their basis and status, in one JSON file, no new statistics.
