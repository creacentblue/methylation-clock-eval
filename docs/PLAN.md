[English](PLAN.md) · [简体中文](PLAN.zh-CN.md)

# Analysis plan

Goal: learn how to evaluate an aging clock without fooling yourself, using public data and computation only. This is a learning and
reproducibility project — not a discovery, and not medical advice.

## 1. Steps (see `scripts/`)

1. **Download** (`01_download.sh`): GEO processed matrices, no IDAT; resumable, `gzip -t` check, sha256 recorded.
2. **Clock CpG lists** (`02_clock_features.py`, `02b_cell_reference.py`): features and metadata for each clock from pyaging; a 450K cell-type reference from purified leukocytes (Reinius 2012).
3. **Parse** (`03_parse.py`): chunked streaming, keep only needed CpGs; QC over all cg probes (per-sample NA fraction, mean beta); parse age/sex.
4. **Run clocks** (`04_run_clocks.py`): pyaging `predict_age`; missing-CpG fractions are computed independently (pyaging substitutes a reference value or 0 and only warns); stray NA betas are mean-imputed per CpG.
5. **Evaluate** (`05_evaluate.py`): MAE, median AE, RMSE, R², Pearson/Spearman, mean signed error; by age bin and sex; plots. DunedinPACE and phenotypic-age clocks are not scored as age predictors.
6. **Optional own model** (`06_train_elasticnet.py`): split first → imputer/scaler fitted on train only → ElasticNetCV → freeze → external test on GSE87571.
7. **Sample audit** (`07_sample_identity.py`) and **batch/dedup metrics** (`09_subset_metrics.py`) for GSE87571.
8. **Cell composition** (`08a_estimate_cells.py`, `08_cell_composition.py`): Houseman projection; acceleration vs cell proportions; cell-adjusted comparison and rankings.

## 2. Training-overlap check (evidence gathered 2026-11-11 from papers / pyaging metadata)

| Clock | GSE40279 | GSE87571 | Evidence |
|---|---|---|---|
| Hannum 2013 | **yes** | no_evidence | model derived (N=482) and validated (N=174) inside GSE40279 |
| Horvath 2013 | **yes** | no_evidence | Genome Biol 14:R115: training data set 3 is the Hannum et al. whole-blood data (exact samples not verified); GSE87571 was deposited in 2016 |
| Lin 2016 | **yes** | no_evidence | Aging 2016: the 99-CpG model was trained on the 656 Hannum blood profiles; GSE87571 is not among its 12 validation sets |
| ZhangEN 2019 | likely_yes | no_evidence | Genome Med 2019, Table 1: GSE40279 is one of 14 pooled cohorts (final model trained on 13,566 pooled samples); GSE87571 is not in Table 1 |
| DNAmPhenoAge | no_evidence | no_evidence | pyaging notes: coefficients fitted on InCHIANTI (other validation cohorts in the paper not checked one by one) |
| DunedinPACE | no_evidence | no_evidence | trained in the Dunedin cohort (EPIC) |
| SkinAndBlood, IntrinClock | unknown | unknown | not verified |

`no_evidence` ≠ independent; it only means no evidence was found. Next: check each paper's supplementary dataset list.

## 3. Known pitfalls

1. **Training leakage**: testing on a cohort a clock was trained on is only a reproduction check (Hannum, Horvath, Lin on GSE40279 here). For own models split first, then select features and scale.
2. **Between-cohort batch / preprocessing differences**: the two cohorts were normalised differently and clock coefficients assume a particular preprocessing, so global offsets appear. No BMIQ, probe-type or batch correction was done.
3. **Imputation of missing CpGs**: pyaging fills absent CpGs with a reference value/0, silently changing predictions. Here ≤ 0.26% of any clock's CpGs are missing, so the effect should be small; PC clocks and EPIC-designed clocks (DunedinPACE) need their own coverage check on 450K. Stray NA betas are filled with the cohort mean per CpG.
4. **Blood cell composition confounding**: part of age-related CpG change is white-cell proportion change. A post-hoc analysis is in the README: Hannum, DNAmPhenoAge and DunedinPACE acceleration correlate consistently with estimated CD4T/granulocyte proportions (|r| about 0.2–0.33); cells explain ~13% of Hannum acceleration variance and ~1% for IntrinClock. But the cell estimates are biased (CD8T r = 0.66 vs the depositor's estimates; PBMC granulocyte estimate 0.12 where ~0 is expected) and there are no subtypes (naive/memory T cells), so "epigenetic age reflects cell-intrinsic ageing" still cannot be concluded.
5. **Clock acceleration ≠ lifespan evidence**: a predicted-minus-chronological age gap needs a longitudinal cohort with outcomes before one can talk about risk; there are no outcomes here. A "rejuvenated" readout is not proof of being biologically younger (it can be cell composition or calibration).
6. **Regression to the mean and age range**: most clocks read low in the very old and high in the young; GSE87571 includes ages 14–30 while Hannum/Horvath mainly target adults, so do not pool age ranges when judging bias (`results/residuals_by_age_bin.csv`).
7. **Sample count and independence (audited)**: GSE87571's 732 GEO records = Johansson 2013's 421 + 311 added samples (Ahsan 2017, same NSPHS cohort); 3 records have no age/genotype/smoking, leaving 729 individuals. No repeat individuals were found (104 genotype variants + correlation across ~22k CpGs), but the population is fairly closed and relatives were not removed. The two `Sample_Group` batches differ in clock bias by up to ~2–2.5 years (with different age distributions, no batch correction) — a warning for cross-batch comparison.
8. **"Recalibrated" vs "raw" error are different questions**: raw-RMSE and recalibrated-RMSE rankings are almost uncorrelated (Kendall τ 0.07 / 0.20). Say whether you compare "out-of-the-box error" or "information content after removing offset/slope".

## 4. Suggested next steps (best value first)

1. **Verify training lists**: resolve "unknown" for SkinAndBlood and IntrinClock, and make Horvath's Additional file 1 GEO ids machine-readable so the `training_overlap` column is checkable.
2. **Make cell estimates more trustworthy**: renormalise reference and both cohorts jointly (e.g. noob/BMIQ or joint quantile normalisation); add subtypes (naive/memory CD4T/CD8T); sensitivity analysis around the weak CD8T estimate; recompute the acceleration–cell table and see which associations persist.
3. **Add a truly independent third cohort**: a public whole-blood 450K/EPIC cohort absent from every clock's training list (e.g. a recent GEO deposit); repeat `metrics.csv` and the cell analysis to see whether rankings and cell sensitivity are stable.
4. **Handle batch and relatedness in GSE87571**: use `Sample_Group` and slide/chip ids (`Slide`, `chip:ass` in the depositor table) as covariates; estimate kinship from genotypes and drop close relatives.
5. **Harden the own model**: nested CV (scaling inside CV), different panels (clock-union vs all probes vs stable probes), age-stratified splits, confidence intervals; optionally compare adding cell proportions as covariates; aim to understand generalisation, not to chase MAE.
