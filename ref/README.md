# `ref/` — externally sourced reference data

Files here are **not computed by the pipeline**. They hold facts that come from an outside
authority and would otherwise have to be typed into a table by hand (and so could be typed
wrong, or invented). Keeping them as committed, sourced files means every published number
has a traceable origin and can be re-checked.

Rule: nothing in this directory may be filled in from memory. Every row carries the URL it
came from and the date it was retrieved.

## `pgs_catalog_metadata.tsv`

Per-score metadata for the 12 polygenic scores in `MOMI_PANEL`, retrieved from the PGS
Catalog on **2026-07-18**. Consumed by `05_prs/deliv_S1_prs_manifest.R` [B04], which joins
it onto the panel to produce Supp Table S1.

| Column | Meaning |
|---|---|
| `id` | PGS Catalog accession |
| `n_snps_in_score` | variants in the **published** score ("Number of Variants" on the catalog page) |
| `method` | score development method as stated by the catalog |
| `PMID` | PubMed ID of the **score's own** publication — *not* of the underlying discovery GWAS |
| `discovery_N` | total N of "Source of Variant Associations (GWAS)"; `NA` where the catalog reports no new discovery GWAS |
| `training_N` | total N of "Score Development/Training"; `NA` where not separately reported |
| `discovery_ancestry` | ancestry composition string as displayed |
| `publication` | first author, journal, year |
| `catalog_url` | page the row was read from |
| `retrieved` | retrieval date |

### Two distinctions that matter when reading this file

**`n_snps_in_score` is the published count, not what we recovered.** The number of variants
actually present in MOMI genotype data is computed separately by B04 from the `.sscore`
`DENOM` column. The ratio of the two is the score recovery rate, and that ratio is the point
of Table S1 — a published 7.4M-variant score of which we recover 5.4M is a different
instrument from the one the authors validated.

## `methods_citations.tsv`

Methodological rules that govern a **design decision** in the pipeline, each with the verbatim
sentence it rests on, the caveats that qualify it, and a DOI. Added because a design choice
defended by "this is standard practice" is not defensible — it has to be attributable.

Currently two entries, both from the Burgess MR guidelines (Wellcome Open Res 2019;4:186 v3,
DOI 10.12688/wellcomeopenres.15555.3, PMID 32760811):

- **MR adjusts for age only** (not the observational confounder set). The guidelines recommend
  "only including as covariates age, sex, genomic principal components of ancestry, and
  technical covariates", because further adjustment risks conditioning on a mediator or
  inducing collider bias. ⚠️ PCs *would* also be legitimate, but no adequate PCA exists — the
  only one on the cluster covers 15% of the sample (`docs/paths_reference.md`), so
  **population stratification is uncontrolled in our MR arm and that is a stated limitation**.
- **First stage and reduced form must use the same covariates**, or the Wald ratio is
  incoherent.

⚠️ Read the `caveat` column before quoting. The first is phrased "in general, we recommend" —
guidance, not a prohibition — and BMI/education/gravidity are *our* examples of covariates
falling outside the recommended set, not the paper's.

## `external_mr_estimates.tsv`

The external comparator arm of **T4 / F4**, the triangulation figure. Retrieved 2026-07-19
from Morales-Berstein F et al., *BMC Medicine* 2026;24:2 (DOI 10.1186/s12916-025-04548-3) —
a wide-angled two-sample MR of maternal blood pressure on pregnancy and perinatal outcomes.

Units are **OR per 10 mmHg**, matching our Wald ratio, so the numbers drop straight into T4.
Headline rows: PTB 1.12 (1.06–1.17) · LBW 1.33 (1.26–1.41) · SGA 1.16 (1.07–1.27).

### ⚠️ Two things to state in the manuscript, not bury

**1. This is not an independent instrument.** Their exposure GWAS is Keaton et al.
(N = 1,028,980), which is the *same* GWAS our EUR score PGS004603/4604 is derived from. So
our PRS-MR and their two-sample MR share a genetic instrument. Agreement between them
demonstrates that the effect **transports** from European to South-Asian and African
pregnancy cohorts — which is exactly what this paper is for — but it is *not* two
independent lines of genetic evidence, and "both MRs agree" must not be allowed to carry
more weight than it earns. The genuinely independent arm of the triangle is the adjusted
observational estimate.

**2. The mirror-image outcomes are a free validity check.** The same paper reports HBW 0.76,
LGA 0.87 and post-term 0.94 — the opposites of LBW, SGA and PTB, all moving the other way.
That internal consistency supports a growth-restriction/shortened-gestation mechanism rather
than a chance pattern, and the explicit nulls (stillbirth 1.00, miscarriage 1.00) are useful
negative comparators.

**Fillable:** DBP estimates. The paper says DBP results are "broadly similar to SBP, although
estimated with higher imprecision" but per-outcome values live in Additional File 3,
Supplementary Table 6, which has not been retrieved. T4's DBP panel needs them.

---

## `pgs_catalog_metadata.tsv` (continued)

**`discovery_N` and `training_N` are not interchangeable.** Some scores rest on a new
discovery GWAS (Keaton: 1,028,980 European). Others perform no new GWAS at all and instead
re-weight existing scores in a training sample — both South Asian scores (PGS004830,
PGS004758) are PRSmix/PRSmixPlus fitted in 28,752 / 28,743 British South Asians from Genes
& Health, which is why their `discovery_N` is `NA`. Reporting the PRSmix training N in a
column headed "GWAS N" would misrepresent what those scores are. The manuscript must
describe them as mixture scores, not as a South Asian GWAS.
