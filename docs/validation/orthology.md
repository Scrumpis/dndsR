# Ortholog pairs against OMA

dndsR pairs genes by reciprocal best hits (RBH), a heuristic. How often do they match an established, independent resource, and does it matter?

## Method

`tools/dev/validation/orthology-oma.py` downloads pairwise orthologs from the [OMA browser](https://omabrowser.org) API, for *A. thaliana* against *A. lyrata* and against *Brassica rapa*. OMA's identifiers differ from RefSeq's, so a protein is matched by the MD5 of its sequence, then by accession. A dndsR pair then is:

- **agree:** its subject is among the mapped OMA orthologs of its query.
- **disagree:** its query has mapped OMA orthologs, and its subject is not one.
- **unscored:** its query has none.

Tables: [data/](data/), `orthology_*`.

## Results

| Comparison | dndsR pairs | Agree | Disagree | Unscored | Disagree, of scored | OMA pairs, mapped |
|---|---:|---:|---:|---:|---:|---:|
| *A. thaliana*, *A. lyrata* | 22,146 | 13,840 | 517 | 7,789 | 3.6% | 24,814, 16,366 |
| *A. thaliana*, *B. rapa* | 19,309 | 10,408 | 1,832 | 7,069 | 15.0% | 31,485, 18,303 |

Disagreement rises steeply with dS. These are the shares of scored pairs that disagree, with scored pairs in brackets:

| dS | *A. lyrata* | *B. rapa* |
|---|---:|---:|
| 0.02 to 0.1 | 0.6% (2,179) | none scored |
| 0.1 to 0.3 | 1.2% (11,614) | 17.3% (897) |
| 0.3 to 0.8 | 56.3% (419) | 14.1% (10,834) |
| 0.8 to 1.5 | 91.9% (86) | 24.7% (465) |
| 1.5 and above | 100% (31) | 79.4% (34) |
| Above 3 times the median dS | 84.8% (368) | 64.2% (53) |
| At most 3 times the median dS | 1.4% (13,975) | 14.7% (12,177) |

**Reading this.**

- **A pair's dS, relative to its comparison's median, is a useful flag for a wrong pairing.** High dS at low overall divergence often means paralogs. `calculate_dnds` reports pairs above 3 times the median dS.
- **Whole-genome duplication matters more than divergence alone.** *B. rapa* is triplicated relative to *Arabidopsis*, and 14–17% of its pairs disagree at ordinary dS, where a best hit can be any of several copies. This may not be an RBH error in a biological sense: both resources may hold true homeoologs, and they pick differently. Splitting subgenomes first can help where they are known, and this has not been measured.
- **dN/dS above 1 was not much affected here.** Of 88 scored *A. lyrata* pairs above 1, 3 disagree, near the overall rate. In *B. rapa*, all 3 pairs above 1 are unscored.

## Effect on term lists

Removing the 517 disagreeing *A. lyrata* pairs and retesting, with the same terms:

| Test | Family | Terms, all pairs | Terms, without them | Lost | Gained |
|---|---|---:|---:|---|---|
| `covariate_enrichment`'s method | GO | 0 | 0 | none | none |
| | InterPro | 16 | 13 | IPR005508, IPR034086, IPR052421 | none |
| Fisher's test, all GO ancestors and InterPro | GO | 9 | 10 | none | GO:0098772 |
| | InterPro | 24 | 20 | IPR005508, IPR006501, IPR022364, IPR035513 | none |

So most terms survive this one, but a few depended on pairs OMA does not support. Repeat this sensitivity check for a term list you intend to report, by filtering on dS: see the [publication protocol](../publication-protocol.md).

## Limits

- OMA is an independent resource, not truth. It can be wrong, and it uses its own annotation.
- 35% of *A. lyrata* and 37% of *B. rapa* pairs could not be scored, because the query's OMA orthologs did not map. These pairs are not a random sample.
- Two comparisons, one query genome, one family of plants.
- Pairs were compared at the gene level, for dndsR's chosen isoform.
- OMA's service refuses Python's default client name, so the script sends an identifying one, `dndsR-validation`, and waits a second between pages.
