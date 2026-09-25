# genetic_pcs — within-cohort joint-platform principal components

The regression models of this study (the portability models in `analysis/03`, the
Mendelian-randomization models in `analysis/04` and `analysis/05`) adjust for the first five
genetic principal components computed **within each cohort, jointly across the two genotyping
technologies** (low-pass whole-genome sequencing and the Global Screening Array). This folder
computes them, checks them for technical artefacts, and also runs the separate, descriptive
cross-cohort analysis shown in Supplementary Figure S3, with the pairwise Hudson FST.

It runs in the secure environment that holds the genotypes, with PLINK 2 and R
(`data.table`); `bcftools` is used when present.

## Run it

```sh
cp genetic_pcs/config_pcs.example.sh genetic_pcs/config_pcs.sh   # edit the paths
Rscript run_all.R 01                                           # the cleaned-sample identifier lists
bash genetic_pcs/run_genetic_pcs.sh all                        # five cohorts, then the cross-cohort step, then the checks
```

`bash genetic_pcs/run_genetic_pcs.sh cohort <cohort>` runs one cohort; `cross` and `assemble`
run the later stages once all five cohorts have finished.

## What it does

| script | step |
|---|---|
| `select_records.R` | One genotype record per woman of the cleaned sample: her low-pass sequencing record if she has one that passed quality control, otherwise her array record. A woman genotyped on both technologies contributes one record only; her array record is used for nothing but the within-woman agreement check. The unrelated reference: all women of the cohort minus relatives, removing, while any pair with kinship ≥ 0.0884 (second degree or closer) remains, the woman with most remaining relatives (ties broken by a seeded random order). Each other woman gets a relatedness weight κ = min(1, Σ 2 × kinship to reference women). |
| `cohort_genotypes.sh` | Variant quality control on each technology among the records used: autosomes; biallelic A/C/G/T single-nucleotide variants; the 43 long-range high-linkage-disequilibrium regions in `highld_regions_hg38.txt` (GRCh38) excluded; variant missingness ≤ 5%; minor-allele frequency ≥ 0.05; low-pass imputation R² ≥ 0.8. Then the shared variant set (`match_snps.R`: position matches, strand-ambiguous and allele-mismatched variants removed, REF/ALT swaps aligned), the cross-technology frequency check (`frequency_concordance.R`: variants whose aligned allele frequencies differ by more than 0.10 removed), LD pruning (`--indep-pairwise 200 50 0.2` in the low-pass records), the array genotypes at the pruned variants aligned to the low-pass alleles (dosages from the re-imputed array VCF with `array_dosage.sh` and `array_dosage_align.R`, or hard calls), and export in blocks of variants. |
| `joint_pca.R`, `pca_lib.R` | The relationship matrix G = ZZ′/m accumulated over the blocks, with Z standardised by the reference allele frequencies of both technologies pooled; the axes are the eigenvectors of the reference block; reference women take their in-sample scores and every other record is projected, P = G<sub>xr</sub> U λ<sup>−1/2</sup>, and corrected for projection shrinkage by d = 1 − ḡ(1 − κ)/λ. Diagnostics: variance of each component explained by the record's platform, standardised differences between platforms, the within-woman cross-technology agreement, duplicates, missingness, Tracy–Widom statistics and cross-validation folds. |
| `cross_cohort_genotypes.sh`, `cross_cohort_pca.R` | The variants passing the within-cohort filters in all five cohorts with identical alleles, pruned in each cohort in turn and thinned at random to at most `JPCA_CROSS_MAXSNP`; the descriptive cross-cohort PCA on the union of the five unrelated references (written as 2-D bins for the figure); pairwise Hudson FST as a ratio of averages over variants, with sample sizes counted in chromosomes and a delete-one-block jackknife over 5-Mb windows. |
| `assemble_joint_pcs.R` | Writes the joint principal-component file and applies the quality gates below. |
| `helpers/` | Identifier normalisation, the Tracy–Widom statistic and the genotyping-technology variable, shared with `analysis/`. |

## Quality gates

The joint principal components may be used only if `aggregate/joint_pca_decision.tsv` reports
PASS. It reports STOP when, in any cohort: fewer than 20,000 variants enter the within-cohort
PCA (10,000 for the cross-cohort set); the cross-technology allele-frequency correlation is below
0.98, or more than 1% of shared variants differ in frequency by more than 0.10 (checked on the
shared variants and again on the genotypes that enter the PCA); in dosage mode, the array dosages
agree with the array hard calls in fewer than 98% of genotypes; in a cohort with at least 30 women
genotyped on both technologies, the median same-woman cross-technology relationship is below 0.80;
the record's platform explains more than 50% of the variance of any of the first five components
(cohorts with at least five array-record women), the standardised mean difference between the
platforms exceeds 1.0 (both groups at least 30 women), or leaving out the women genotyped on both
technologies moves a component's platform difference by more than 0.50 SD (at least 30 such women); a woman enters
a fit twice, or two different women share a genotype; any woman of the cleaned sample lacks the
five components; or the genotyping-technology variable differs from the one the analysis uses.
Smaller departures are reported as warnings.

## Outputs

Everything is written below `JPCA_OUT`:

* `participant_level/joint_pcs.rds` (and `.tsv`) — one row per woman: `IID`, `cohort`,
  `technology`, the `record` used, whether she is `in_reference`, her relatedness weight `kappa`,
  and `PC1`…`PC10`. Participant-level: it stays in the secure environment. `config.R` points
  `pc_file` at it.
* `aggregate/` — the decision table and the diagnostics above, per cohort and combined, and the
  cross-cohort results (`aggregate/cross/`). The 2-D bins of the cross-cohort figure are close to
  individual-level at the extremes of the distribution and are not published; the published
  aggregate values of Supplementary Figure S3 are the counts per cohort and technology.
