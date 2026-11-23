# Methods-paper concept: does evidence-grounded interpretation improve evolutionary inference? (prospective)

Status: **prospective research concept for the maintainer's review (2026-11-23); not authorized for implementation.** Deliverable 7 of the maintainer's message of 2026-11-23. Written to be falsifiable: the paper is worth writing whether the hypothesis holds or not, but only if the evaluation is rigorous.

## Central hypothesis (and its null)

H1: Giving the interpreter (a person or an agent) a deterministic, evidence-linked diagnostic record reduces unsupported conclusions from evolutionary selection analyses, without removing valid ones; and an agent reading that record adds a further, measurable reduction beyond the record alone.

H0 for each step: no difference in the rate of unsupported conclusions, or a reduction bought only by indiscriminate over-qualification (loss of valid conclusions).

The second step is the one most likely to fail; the design isolates it.

## Positioning (literature reviewed in this session; preprints flagged)

- **General and workflow agents:** Biomni (bioRxiv 2025, self-reported gains), BioMaster, AutoBA, BioAgents (Sci Rep 2025), CellAgent. These automate analysis; their evaluations score completion and answers, not the validity of the scientific conclusion drawn from correct output.
- **Benchmarks:** BixBench (arXiv 2503.00096; open answers graded by LLMs), BioAgent Bench (ICML 2026; robustness to corrupted inputs and decoys exposes weaknesses), FlowBench (bioRxiv preprint 2026; separates planning, fault recovery, interpretation and output fidelity; reports data-grounded interpretation unsolved and fixes that "force a clean exit while leaving the underlying data invalid"), PromptBio-Bench, OmicsBench (preprints).
- **Closest framing:** "Evaluating Agentic Bioinformatics through Function, Evidence, and Validation" (arXiv 2607.27556, preprint): across 109 systems, explicit verification and provenance-linked evidence are rare; it proposes judging agents by inspectable trajectories. It is a framework and review, not an empirical test in one domain with ground truth.

**Defensible novelty** (none of which is the combination of tools by itself): (1) diagnostics whose every trigger is tied to simulation evidence about *when the inference is wrong* (dndsR's validation studies), with an explicit status taxonomy (known failure, validated, unvalidated concern, not assessable); (2) a ground-truth benchmark of *interpretive validity* in molecular evolution, where the statistics can be correct and the conclusion still unsupported; (3) an ablation that separates the value of the record from the value of the agent, across models.

## Arms

1. Conventional: engine outputs and standard documentation.
2. Conventional plus deterministic QC and the template report (layer B, no LLM).
3. Generic LLM given arm 1's outputs only.
4. Evidence-grounded agent: arm 2's information plus the agent (layer C).
5. Ablation: arm 4 without the methodological reference notes.

The same statistical results and engines in every arm; arms 2, 4 and 5 receive identical deterministic information. Readers: several language models (from at least two providers) and, for a subset, human analysts at two experience levels, blinded to the study's hypothesis.

## Benchmark

- **Synthetic, known truth (primary):** from dndsR's simulators: clean null and selected genes; corrupted foreground segments; realigned indels; multinucleotide changes; confounded enrichment; underpowered designs; RELAX on 16 taxa. Truth: whether selection is present, whether the test is valid in that configuration, and what claims are supported.
- **Adversarial:** a valid significant result paired with an unsupported prompt ("does this show adaptation to altitude?"); significant results only on flagged genes; nonsignificant results in underpowered designs.
- **Realistic, expert-reviewed:** public genome analyses (dndsR's public pairs and study E if run), each judged by at least two evolutionary geneticists with an adjudication rule.
- Development cases and evaluation cases are disjoint, and the evaluation set is frozen before any arm is tuned.

## Endpoints (fixed before evaluation)

Primary: rate of unsupported conclusions per case (false declaration of selection where the test is invalid; claims of adaptation or causality; absence-of-evidence written as evidence of absence), scored blind against a rubric. Co-primary guardrail: retention of valid conclusions (a reduction in the primary endpoint achieved by losing valid conclusions does not count). Secondary: recognition of violated assumptions; calibration of stated uncertainty; fidelity to outputs (every number checked mechanically); fabricated diagnostics or references; traceability (share of claims with a resolvable citation); reproducibility across repeated runs.

Analysis: cases as the unit, mixed models with case and reader random effects, repeated agent runs per case; pre-specified contrasts: arm 2 against 1 (value of the record), arm 4 against 2 (value of the agent beyond the record), arm 4 against 3 (value of grounding), arm 4 against 5 (value of references); generalization by method, data source and model.

## Main figures (planned)

1. Architecture: execution, evidence record, optional interpretation, with provenance links.
2. The diagnostic evidence base: for each known-failure trigger, the simulation result that established it.
3. Primary endpoint by arm, with the valid-conclusion guardrail.
4. Ablation: record alone against agent with and without references, by model.
5. Failure analysis: what each arm still gets wrong.

## Minimum convincing evidence

A pre-registered evaluation set of at least a few hundred cases across at least four methods; at least two model families; arm 2 against arm 1 and arm 4 against arm 2 both estimated with intervals; zero tolerance on fabricated numbers enforced mechanically; human-reader subset showing the same direction.

## Why it might not work, and what a negative result would mean

- Deterministic QC (arm 2) may capture almost all of the benefit, because most failures in this domain are detectable by fixed rules. Then the publishable finding is that the evidence record, not the agent, is what matters, which is still useful to the field and argues against agent-first tools.
- Agents may over-qualify, trading valid conclusions for safety; the guardrail endpoint is there to show it.
- Models change quickly; results may not transfer between versions. Reporting by model, and releasing the benchmark, is the defence.
- The benchmark could be built to suit the framework. Separate case authors, expert adjudication, and adversarial cases written by people outside development reduce that.
- Positioning: if a comparable domain-specific ground-truth benchmark appears first, the contribution shrinks to the dndsR-specific evidence base.

**Best positioning, on current knowledge:** a benchmark and evaluation of interpretive validity for evolutionary inference, with dndsR's evidence record as the reference implementation. A framework paper without the empirical comparison would not be competitive.
