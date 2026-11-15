# McDonald-Kreitman review: handoff to the MK maintainer

For the author of PR #5 and #6 and their agent. Written by the dndsR integration agent at `5f8e0b8` (`dev`), reviewing `R/mk.R` as merged in `4c6254c`. **Nothing here is approved or decided.** The recommendations are one reviewer's view, for the MK maintainer and the project maintainer to accept, change, or reject. The integration agent will not modify `R/mk.R`, MK tests, or MK methodology docs.

## What the integration agent already changed in MK-owned files

None of these changes counts or statistics. Keep or revise them as you prefer.

| Commit | File | Change |
|---|---|---|
| `b306639` | `R/mk.R` | `setNames(` → `stats::setNames(` (R CMD check note) |
| `b306639` | `DESCRIPTION` | `GenomicRanges` and `IRanges` added to `Imports` (R CMD check warning; both were already hard dependencies of GenomicFeatures) |
| `cfb2a59` | `docs/mcdonald-kreitman.md`, `docs/tutorials/your-own-data.md` | The single-diploid statement was wrong (see 3); singleton removal is now described as a choice, and `AC>=1 && AC<AN` is given for keeping singletons. The documented population command is unchanged |
| `cfb2a59` | `R/calculate_dnds.R` (roxygen only), `man/calculate_dnds.Rd`, `NEWS.md` | The headerless column 6/7 interpretation is documented (see 4) |

CI passes on every platform at `c78abfa`.

## 1. Divergence class depends on which species is `snp_side`

**Observed.** `.mk_classify_codon(a, b)` classifies each step by the degeneracy of the changed position in the codon the step *starts from*. `.mk_pair_walk()` always starts from the SNP side's codon (`ref_seq`). For some codon pairs the position is 0-fold in one codon and 2-fold in the other, so the same fixed difference counts as `Dn` or is excluded, depending on `snp_side`.

**Evidence** (script A below, genetic code only):

- 526 ordered single-difference pairs of sense codons; for 38 (7.2%) the class depends on direction, always 0-fold one way and mid-fold the other, never 4-fold. Examples: AGA (R) → GGA (G), CGA (R) → GGA (G), ATA (I) → CTA (L), mainly first positions of the six-fold Arg and Leu codons and their neighbours.
- No step classed 4-fold changes the amino acid, and none classed 0-fold is synonymous: the classes are internally correct.
- On 1,472 public *A. thaliana*–*A. lyrata* orthologs (script B; codon alignments from public NCBI genomes): reading from *A. thaliana* gives Dn 33,693.7, reading from *A. lyrata* 33,624.6 (0.2%). Ds is identical. But Dn changes in 965 genes (66%), by up to 8 steps in one gene.

**Why it matters.** Per-gene MK tables, alpha, NI, and Fisher p-values change when the polymorphism data come from the other species, for reasons unrelated to biology. The pooled effect is small; the per-gene effect is not negligible for genes with few differences.

**Options.**

| Option | Effect | Cost |
|---|---|---|
| a. Keep as is, document that classes are read from the SNP side | No code change | Per-gene direction dependence remains |
| b. Symmetric divergence: a step counts as 0-fold (4-fold) only if the position is 0-fold (4-fold) in both its start and end codon; otherwise mid-fold | Removes direction dependence for divergence | Dn −4.5% on the sample above (32,163 vs 33,694); Ds unchanged; 843 genes change |
| c. Strict site classes: a position counts as 0-fold only if it is 0-fold in the reference codon and in every codon reached by changing that position; 4-fold likewise | Symmetric, allele-free, so the same rule serves divergence, VCF, and BED polymorphisms | Similar or slightly larger loss than b (not measured) |

**Recommendation (not a decision): c,** because it applies the same class definition to `Pn`/`Ps` and `Dn`/`Ds`, which the MK test needs, works without alleles (BED), and removes the dependence on `snp_side`. If only b is adopted, polymorphisms stay classified from the reference codon alone, so the two halves of the table use slightly different definitions.

**Files:** `R/mk.R` `.mk_syn_alternatives()`, `.mk_classify_codon()`, `.mk_poly_counts()`, `.mk_pair_walk()`.

**Suggested tests:** (i) for every single-difference sense pair, `.mk_classify_codon(a, b)` and `.mk_classify_codon(b, a)` give the same class; (ii) the bundled example gives identical `_mk.tsv` with `snp_side` swapped when the SNP file holds no sites (divergence only); (iii) hand-derived classes for AGA↔GGA and ATA↔CTA.

## 2. Strict 0-fold/4-fold versus allele-aware counting

**Observed.** PR #6 counts only 0-fold and 4-fold positions, for both VCF and BED, and reports 2- and 3-fold positions in `snps_midfold` and `div_steps_midfold`. PR #5 classified VCF polymorphisms and all divergence by the actual amino acid change, and BED polymorphisms by degeneracy, so the two input formats used different classes.

**Evidence** (script B, divergence only, same 1,472 genes): 31.5% of divergence steps are mid-fold and excluded. Allele-aware counting (PR #5's divergence rule) gives Ds 60,884 (+71%) and Dn 40,400 (+20%) with nothing excluded. In script A, 38 of 134 synonymous and 59 of 392 nonsynonymous single steps fall at mid-fold positions.

**Convention versus consequence.**

- *Convention:* both are established. The original MK test counts replacement and synonymous changes by their actual effect; genome-scale and asymptotic MK analyses commonly use 0-fold and 4-fold sites, with 4-fold as the closest-to-neutral class.
- *Statistical consequence:* the test needs `Pn`/`Ps` and `Dn`/`Ds` defined the same way. PR #6 does that for both formats; PR #5 did not for BED. Excluding mid-fold sites costs power, so more genes have `Ps = 0` or `Ds = 0` and an undefined alpha or NI.
- *Biological consequence:* 2-fold synonymous sites (mostly third-position transitions) are under weak selection in some taxa (codon usage), so leaving them out keeps the synonymous class closer to neutral. Allele-aware counting mixes them in.

**Recommendation (not a decision):** keep strict 0-fold/4-fold as the primary method (with option c above). If power matters, consider an opt-in allele-aware mode for VCF only, with its own output columns, so tables from the two modes cannot be confused. BED cannot support allele-aware counts.

## 3. VCF prefiltering

**Observed and already corrected in docs (`cfb2a59`):** the documented filter `-i 'AC>=2 && AC<AN'` was said to keep "only heterozygous sites" for one diploid sample. With `AN = 2` it keeps nothing.

**Evidence:** BCFtools 1.18 on synthetic VCFs (recorded in `skills/dndsr/references/verification.md`):

- Three diploids: `AC>=2 && AC<AN` keeps only the common SNP; `AC>=1 && AC<AN` also keeps the singleton.
- One diploid: `AC>=2 && AC<AN` keeps nothing; `AC>=1 && AC<AN` and `GT="het"` keep the heterozygous site.
- Both runs: `-f PASS -v snps -m2 -M2` dropped the LowQual, multiallelic, and indel records as intended.

**Open question for you:** `AC>=2` removes singletons, a frequency filter, while the Limits section says dndsR applies none. The docs now present it as the user's choice.

**Recommendation (not a decision):** keep dndsR itself free of SNP filtering, and keep presenting singleton removal as an explicit, reported choice in the docs, rather than adding an option to the package.

## 4. Headerless comparison files

**Observed (from code, not run end to end):** `R/calculate_dnds.R` (near line 249) and `R/split_comparisons.R` now read headerless columns 6 and 7 as `snp_file` and `snp_side`. Before PR #5, extra columns were ignored. A headerless file with any extra column now reaches `.mk_validate_snp_args()` only after that comparison's dN/dS is computed (near line 535), and stops with "snp_file is set … but snp_side is missing", which does not mention the column interpretation. No silent wrong result is possible.

**Options:** (a) keep, as now documented in NEWS and `?calculate_dnds`; (b) validate MK columns for every row before any dN/dS work, with a message naming the headerless column rule; (c) read columns 6 and 7 as MK only when column 7 is `query` or `subject`, otherwise ignore them as before (restores the old behavior, at the cost of guessing).

**Recommendation (not a decision): b.** It fails fast and clearly, without guessing. It is a public-interface change, so it needs the project maintainer's approval. `extract_cds` and the other commands still ignore the extra columns.

## Not investigated

Fisher's exact test on rounded path-averaged `Dn`/`Ds` (`dn_int`, `ds_int`): R's `round()` rounds halves to even, so 0.5 becomes 0. This matters only for genes with very few differences, and was not examined.

## Reproduction

Run inside the dndsR image (or any R with Biostrings). Script A needs nothing else. Script B needs codon alignments (FASTA, one file per gene, names like `At|AT1G01010`); the run above used the `prepared/alignments` folder from `phylo_prepare` on the public five-species benchmark (`docs/validation/phylo-calibration-pilot.md`). Any pairwise codon alignments with that naming work.

```
git show 4c6254c:R/mk.R > mk_current.R     # PR #6
git show 6225bae:R/mk.R > mk_pr5.R         # PR #5, for the allele-aware rule
Rscript repro_asymmetry.R mk_current.R
Rscript repro_rules.R mk_current.R mk_pr5.R <alignment dir> At Al 1500
```

Script A, `repro_asymmetry.R`:

```r
source(commandArgs(TRUE)[1])
b <- c("A", "C", "G", "T")
cod <- as.vector(outer(outer(b, b, paste0), b, paste0))
cod <- cod[vapply(cod, .mk_aa, "") != "*"]
cls <- function(x) c("4fold", "0fold", "mid")[which(x > 0)]
pairs <- do.call(rbind, lapply(cod, function(a) {
  o <- cod[vapply(cod, function(c2) sum(strsplit(a, "")[[1]] != strsplit(c2, "")[[1]]) == 1L, logical(1))]
  if (length(o)) cbind(a, o)
}))
fwd <- apply(pairs, 1, function(p) cls(.mk_classify_codon(p[1], p[2])))
rev <- apply(pairs, 1, function(p) cls(.mk_classify_codon(p[2], p[1])))
cat("ordered single-difference sense codon pairs:", nrow(pairs), "\n")
cat("class depends on direction:", sum(fwd != rev), "\n")
print(table(read_from_first = fwd, read_from_second = rev))
syn <- apply(pairs, 1, function(p) .mk_aa(p[1]) == .mk_aa(p[2]))
cat("4-fold but amino acid changes:", sum(fwd == "4fold" & !syn), "; 0-fold but synonymous:", sum(fwd == "0fold" & syn), "\n")
cat("synonymous single steps:", sum(syn), "of which mid-fold:", sum(syn & fwd == "mid"), "\n")
cat("nonsynonymous single steps:", sum(!syn), "of which mid-fold:", sum(!syn & fwd == "mid"), "\n")
```

Expected: 526 pairs; 38 direction-dependent (19 + 19, 0-fold ↔ mid); 0 and 0 misclassified; 134 synonymous (38 mid-fold); 392 nonsynonymous (59 mid-fold).

Script B, `repro_rules.R`:

```r
a <- commandArgs(TRUE); n_max <- if (length(a) >= 6) as.integer(a[6]) else 1500L
e6 <- new.env(); sys.source(a[1], e6); e5 <- new.env(); sys.source(a[2], e5)
sym_codon <- function(x, y) {   # a step counts as 4-/0-fold only if so in its start AND end codon
  xa <- strsplit(x, "")[[1]]; xb <- strsplit(y, "")[[1]]; tot <- c(0, 0, 0); ok <- 0
  for (p in e6$.mk_permutations(which(xa != xb))) {
    cur <- xa; r <- c(0, 0, 0); good <- TRUE
    for (pos in p) {
      d1 <- e6$.mk_syn_alternatives(paste(cur, collapse = ""), pos); cur[pos] <- xb[pos]
      s2 <- paste(cur, collapse = ""); if (e6$.mk_aa(s2) == "*") { good <- FALSE; break }
      d2 <- e6$.mk_syn_alternatives(s2, pos)
      k <- if (d1 == 3 && d2 == 3) 1 else if (d1 == 0 && d2 == 0) 2 else 3; r[k] <- r[k] + 1
    }
    if (good) { tot <- tot + r; ok <- ok + 1 }
  }
  if (ok) tot / ok else c(0, 0, 0)
}
rules <- list(current_fromA = function(x, y) e6$.mk_classify_codon(x, y),
              current_fromB = function(x, y) e6$.mk_classify_codon(y, x),
              symmetric = sym_codon,
              allele_aware_pr5 = function(x, y) c(e5$.mk_classify_codon(x, y), 0))
walk <- function(s, t, f) { out <- c(0, 0, 0)
  for (i in seq_len(nchar(s) %/% 3)) { x <- substr(s, 3 * i - 2, 3 * i); y <- substr(t, 3 * i - 2, 3 * i)
    if (x == y || !grepl("^[ACGT]{3}$", x) || !grepl("^[ACGT]{3}$", y) || e6$.mk_aa(x) == "*" || e6$.mk_aa(y) == "*") next
    out <- out + f(x, y) }
  out }
fs <- list.files(a[3], "\\.fasta$", full.names = TRUE); set.seed(1); fs <- sample(fs, min(n_max, length(fs)))
res <- list()
for (f in fs) { x <- Biostrings::readDNAStringSet(f); nm <- sub("[|].*", "", names(x))
  if (!all(c(a[4], a[5]) %in% nm)) next
  s <- as.character(x[[match(a[4], nm)]]); t <- as.character(x[[match(a[5], nm)]])
  res[[length(res) + 1]] <- sapply(rules, function(r) walk(s, t, r)) }
tot <- Reduce(`+`, res); rownames(tot) <- c("Ds (4-fold or syn)", "Dn (0-fold or nonsyn)", "excluded (2-/3-fold)")
cat("gene pairs:", length(res), "\n"); print(round(tot, 1))
dn <- sapply(res, function(m) m[2, ])
cat("genes whose Dn changes with direction:", sum(abs(dn[1, ] - dn[2, ]) > 1e-9), "; largest change:", max(abs(dn[1, ] - dn[2, ])), "steps\n")
cat("genes whose Dn differs, symmetric vs current_fromA:", sum(abs(dn[3, ] - dn[1, ]) > 1e-9), "\n")
```

Expected on the benchmark sample (1,472 pairs): Ds 35,666.3 under the three degeneracy rules and 60,884.1 allele-aware; Dn 33,693.7 (from *A. thaliana*), 33,624.6 (from *A. lyrata*), 32,163.3 (symmetric), 40,399.9 (allele-aware); 965 genes change with direction, at most 8 steps; 843 differ under the symmetric rule. The walk skips gapped, non-ACGT, and stop codons but has no SNP exclusion, so it measures divergence classes only.
