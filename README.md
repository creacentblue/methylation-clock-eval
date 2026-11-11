# methylation-clock-eval

[English](README.md) · [简体中文](README.zh-CN.md)

Cross-cohort evaluation of published DNA-methylation aging clocks on two public whole-blood Illumina 450K cohorts
(GEO GSE40279 and GSE87571), done with a "don't fool yourself" attitude: every clock/cohort pair carries a
training-overlap flag, missing-CpG fraction, bias breakdowns, a sample-count audit of GSE87571, a blood-cell-composition
analysis, and an optional frozen ElasticNet external test. Everything is scripted and rerunnable.

> **Not medical advice. Not a clinical tool.** This is an educational/research reproducibility project that uses public,
> de-identified GEO data. Nothing here may be used to estimate anyone's biological age for personal or clinical decisions,
> and clock age acceleration is **not** evidence about lifespan, health, or the effect of any intervention.

## Motivation

Epigenetic clocks are easy to run and easy to misread. This repository is a worked example of the unglamorous parts:
checking whether a "validation" cohort was actually in a clock's training data, how many CpGs were silently missing,
how much of the error is a calibration offset versus noise, whether the numbers hold across two cohorts and two
batches inside one cohort, and how much blood cell composition explains.

## What it does / does not claim

**Does:** reproduce published clocks through `pyaging` on public data; report error metrics, bias by age bin and sex; flag
training overlap with sourced evidence; audit the GSE87571 sample count; estimate cell composition and quantify its link to
age acceleration; train and externally test a simple ElasticNet.

**Does not:** propose a new clock; claim any clock is "best" (rankings change with the cohort and with calibration);
make clinical, mortality or intervention claims; prove independence of any evaluation cohort (`no_evidence` only means no evidence
of overlap was found).

## Quickstart

Requirements: Python ≥ 3.11 (developed on 3.13), ~4.3 GB free disk for the downloaded GEO files and intermediates,
a machine like the development one (15 GB RAM; peak memory was not measured). The pinned environment in
`requirements.txt` includes PyTorch (pyaging backend); a CPU-only PyTorch build is enough.

```bash
git clone <this repository> && cd methylation-clock-eval     # or unpack the archive
./run_all.sh                    # creates .venv, downloads data (~4 GB), runs every step
# afterwards:  SKIP_DOWNLOAD=1 ./run_all.sh     (reuse data/raw)     or     make all
pytest                          # lightweight tests; need no downloaded data
```

Measured on the development machine (8 cores, 15 GB RAM, warm disk cache): all analysis steps after the download took
about 3 minutes wall-clock (`SKIP_DOWNLOAD=1 ./run_all.sh`). Download time depends on your connection (≈ 4 GB from NCBI GEO).
Results were regenerated twice and the key tables were byte-identical. Pinned versions: `requirements.txt` (direct) and
`requirements.lock.txt` (full `pip freeze`).

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 1 | `scripts/01_download.sh` | Download GEO processed matrices (resumable, `gzip -t` check, sha256 in `results/raw_files.sha256`) |
| 2 | `scripts/02_clock_features.py` | CpG lists and metadata (incl. pyaging citation text) for the evaluated clocks |
| 2b | `scripts/02b_cell_reference.py` | 450K cell-type reference from Reinius 2012 purified leukocytes (GSE35069), minfi-style probe selection |
| 3 | `scripts/03_parse.py` | Stream-parse the matrices; keep only needed CpGs; parse age/sex; QC over all probes |
| 4 | `scripts/04_run_clocks.py` | Run clocks with pyaging; compute missing-CpG fractions; mean-impute stray NA betas |
| 5 | `scripts/05_evaluate.py` | Metrics, age-bin and sex breakdowns, plots |
| 6 | `scripts/06_train_elasticnet.py` | Optional: ElasticNetCV on GSE40279, frozen, tested on GSE87571 |
| 7 | `scripts/07_sample_identity.py` | GSE87571 732-vs-421 audit: repeat individuals, `Sample_Group`, dedup flag |
| 8a | `scripts/08a_estimate_cells.py` | Houseman constrained projection → 6 cell proportions |
| 8 | `scripts/08_cell_composition.py` | Acceleration vs cell proportions; cell-adjusted comparison; rankings |
| 9 | `scripts/09_subset_metrics.py` | GSE87571 metrics by dedup / batch (`Sample_Group`) / adults |

Shared code: `scripts/config.py` (clock list, training-overlap table with evidence notes), `scripts/metrics_utils.py`,
`scripts/cell_deconv.py`. Tests: `tests/`. Plan and pitfalls: `docs/PLAN.md`. References: `docs/REFERENCES.md`.

## Data sources and attribution

Nothing from GEO is redistributed here; scripts download it at run time. Please cite the source papers if you use the data.

| Data | Source |
|---|---|
| GEO GSE40279 — 656 whole-blood 450K profiles, ages 19–101 | Hannum et al. 2013, *Mol Cell* 49:359–367, doi:10.1016/j.molcel.2012.10.016 |
| GEO GSE87571 — 732 records / 729 individuals, ages 14–94 (NSPHS, Northern Sweden) | Johansson et al. 2013, *PLoS ONE* 8:e67378, doi:10.1371/journal.pone.0067378 (421-sample batch); Ahsan et al. 2017, *PLoS Genet* 13:e1007005, doi:10.1371/journal.pgen.1007005 (additional samples; additional sample characteristics) |
| GEO GSE35069 — purified blood cell types, 6 male donors (cell reference) | Reinius et al. 2012, *PLoS ONE* 7:e41361, doi:10.1371/journal.pone.0041361 |
| Clock implementations and weights | pyaging — de Lima Camillo 2024, *Bioinformatics* 40:btae200, doi:10.1093/bioinformatics/btae200 (v0.6.2) |

Each clock's original paper (as given in pyaging's metadata), plus the Houseman 2012 and minfi 2014 method papers, are listed in
[`docs/REFERENCES.md`](docs/REFERENCES.md). GEO data and the Reinius reference keep their own terms; the MIT license covers this repository's code and text only.

## Key results

All numbers below come from files in `results/` produced by the scripts (not typed by hand).

### 1. Clock performance (`results/metrics.csv`)

MAE/RMSE in years; bias = predicted − chronological (positive = predicted too old); "CpGs missing" = clock CpGs absent from the
cohort matrix. *Training overlap*: `yes` = documented in the clock's paper (then the result is only a reproduction check),
`likely_yes` / `no_evidence` / `unknown` as defined in `scripts/config.py` — `no_evidence` is **not** proof of independence.

| Cohort | Clock | n | MAE | RMSE | R² | r | bias | CpGs missing | Training overlap |
|---|---|---|---|---|---|---|---|---|---|
| GSE40279 | horvath2013 | 656 | 4.77 | 6.33 | 0.815 | 0.918 | -2.34 | 0/353 | yes |
| GSE40279 | hannum | 656 | 5.51 | 6.78 | 0.788 | 0.946 | +4.75 | 0/71 | yes |
| GSE40279 | skinandblood | 656 | 4.36 | 5.79 | 0.846 | 0.940 | -2.45 | 1/391 | unknown |
| GSE40279 | zhangen | 656 | 4.00 | 5.38 | 0.867 | 0.964 | -0.04 | 1/514 | likely_yes |
| GSE40279 | lin | 656 | 4.12 | 5.34 | 0.868 | 0.932 | -0.00 | 0/99 | yes |
| GSE40279 | intrinclock | 656 | 3.69 | 4.93 | 0.888 | 0.944 | -0.56 | 0/380 | unknown |
| GSE40279 | dnamphenoage* | 656 | 9.04 | 11.05 | 0.436 | 0.852 | -7.89 | 0/513 | no_evidence |
| GSE40279 | dunedinpace† | 656 | – | – | – | 0.135 | – | 16/20000 | no_evidence |
| GSE87571 | horvath2013 | 729 | 3.90 | 5.05 | 0.942 | 0.971 | +0.33 | 0/353 | no_evidence |
| GSE87571 | hannum | 729 | 5.13 | 6.14 | 0.914 | 0.975 | +4.01 | 0/71 | no_evidence |
| GSE87571 | skinandblood | 729 | 2.90 | 3.63 | 0.970 | 0.986 | +0.61 | 0/391 | unknown |
| GSE87571 | zhangen | 729 | 6.29 | 8.24 | 0.845 | 0.991 | -3.23 | 0/514 | no_evidence |
| GSE87571 | lin | 729 | 5.67 | 7.21 | 0.881 | 0.960 | +3.12 | 0/99 | no_evidence |
| GSE87571 | intrinclock | 729 | 3.21 | 4.17 | 0.960 | 0.982 | -1.32 | 0/380 | unknown |
| GSE87571 | dnamphenoage* | 729 | 9.34 | 10.78 | 0.735 | 0.968 | -8.90 | 0/513 | no_evidence |
| GSE87571 | dunedinpace† | 729 | – | – | – | 0.593 | – | 0/20000 | no_evidence |

\* DNAmPhenoAge predicts a phenotypic age, so age-scale error and bias are not meaningful; read correlation only.
† DunedinPACE predicts the *pace* of aging, not an age in years (EPIC-trained): only its Pearson r with age is shown.

Observations: (i) Hannum overshoots by about +4–5 years even on its own training cohort (MAE 5.5; 3.8 after removing the mean offset).
pyaging's Hannum has a zero intercept, i.e. the published coefficients as released; I did not check the original calibration, so this is not
an exact reproduction of the paper's error. (ii) Biases depend on age (`results/residuals_by_age_bin.csv`): e.g. Horvath/SkinAndBlood read
young at 80+, ZhangEN reads very young under 30 (about −12 to −14 years). (iii) Rankings differ between cohorts; ZhangEN has the highest
correlation on GSE87571 (0.991) but a large raw MAE, i.e. a calibration problem rather than a ranking problem. (iv) Missing CpGs are at most 0.26% per clock.

### 2. Optional frozen ElasticNet (`results/own_model_metrics.csv`)

Trained only on GSE40279 (80/20 split made before any other step; imputer and scaler fitted on the training part only), then frozen and applied to GSE87571.
Panel = 21,739 CpGs of the evaluated clocks present in both cohorts (fixed in advance); 581 non-zero coefficients.

| Evaluation set | n | MAE | RMSE | R² | r | bias |
|---|---|---|---|---|---|---|
| GSE40279 train (in-sample) | 524 | 1.13 | 1.45 | 0.990 | 0.996 | -0.00 |
| GSE40279 held-out 20% | 132 | 3.39 | 5.02 | 0.884 | 0.942 | +0.51 |
| GSE87571 external (all) | 729 | 2.85 | 3.65 | 0.970 | 0.986 | +0.37 |
| GSE87571 external (age>=19) | 645 | 2.82 | 3.63 | 0.963 | 0.982 | +0.75 |

The held-out score is optimistic because the panel itself was chosen by published studies, some of which used GSE40279; GSE87571 has a wider and
younger age range, so its errors are not comparable with GSE40279's.

### 3. GSE87571: 732 samples vs N = 421 (`scripts/07_sample_identity.py`)

Resolved. The depositor's supplementary sample characteristics contain `Sample_Group`: 421 samples in group 2012 (median age 44, matching Johansson 2013's N = 421, 14–94)
and 311 in group 2013 → 732. (The meaning of `Sample_Group` is my inference from these numbers; GEO does not define it.) The extra samples come from the same NSPHS
cohort (Karesuando 2006: 683 records; Soppero 2009: 49), analysed in Ahsan et al. 2017, which reports 729 samples after QC. Three GEO records lack age, genotype and
smoking information (GSM2334015, GSM2334342, GSM2334619); they are not identifiable individuals and are excluded (729 left). I could not confirm whether they are
the paper's control samples.

Repeat-individual check: 104 genotype variants compared across all 267,546 sample pairs gave a smallest mean dosage difference of 0.109 (median 0.472; a repeat
should be near 0); centered methylation correlation across about 22,000 CpGs peaked at 0.736. With thresholds of genotype difference < 0.05 or correlation > 0.90,
**0 duplicate pairs** were found, so deduplicated metrics equal the 729-sample metrics. The closest genotype pairs differ in sex or age and look like relatives in a fairly closed
population; relatives were not removed. The two batches differ systematically (e.g. Hannum bias +4.96 in group 2012 vs +2.73 in group 2013) while also differing in age
distribution (median 44 vs 53.5) and being uncorrected for batch (`results/metrics_gse87571_subsets.csv`).

### 4. Blood cell composition and age acceleration (`scripts/08*.py`)

Cell proportions: Houseman constrained projection (w ≥ 0, Σw ≤ 1) on a 450K reference built from Reinius 2012 (550 probes chosen minfi-style; 10 absent in GSE40279).
Agreement with the depositor's own estimates on GSE87571: Pearson r 0.99 (Gran), 0.96 (NK), 0.93 (B), 0.92 (Mono), 0.92 (CD4T), **0.66 (CD8T)**. On Reinius PBMC samples
(not in the reference) the granulocyte estimate is 0.12 where it should be near 0, so treat estimates as relative covariates, not absolute proportions. pyaging's built-in EPIC
deconvolution was not used because ~57% of its 600 reference CpGs do not exist on 450K (diagnostic comparison in `results/cell_estimates_agreement_gse87571.csv`).

Definitions: **AA** = residual of `predicted ~ age`; **IEAA** = residual of AA on 5 cell proportions (6th dropped, proportions sum to ~1). After cell adjustment, MAE on chronological age is
no longer the right target, so the comparison uses recalibrated 10-fold-CV RMSE (`age ~ predicted`, and `age ~ predicted + cells`) and SD(AA) vs SD(IEAA); recalibration uses each cohort's own
labels, so it is a diagnostic, not a deployable clock.

(a) Pearson correlation of age acceleration with estimated cell proportion (GSE40279 / GSE87571):

| Clock | CD8T (40279 / 87571) | CD4T (40279 / 87571) | NK (40279 / 87571) | Bcell (40279 / 87571) | Mono (40279 / 87571) | Gran (40279 / 87571) |
|---|---|---|---|---|---|---|
| horvath2013 | +0.11 / +0.18 | -0.18 / -0.13 | +0.07 / +0.24 | -0.08 / -0.08 | +0.09 / -0.16 | +0.05 / -0.09 |
| hannum | +0.04 / -0.04 | -0.33 / -0.29 | +0.05 / +0.23 | -0.15 / -0.16 | +0.08 / -0.01 | +0.21 / +0.12 |
| skinandblood | +0.01 / +0.11 | -0.07 / -0.01 | +0.07 / +0.18 | -0.09 / -0.12 | +0.06 / -0.13 | +0.02 / -0.09 |
| zhangen | -0.07 / +0.05 | -0.08 / -0.00 | -0.02 / +0.14 | -0.10 / -0.09 | +0.07 / -0.07 | +0.10 / -0.07 |
| lin | +0.10 / +0.13 | -0.18 / -0.05 | -0.02 / +0.15 | -0.09 / -0.06 | +0.07 / -0.13 | +0.10 / -0.07 |
| intrinclock | -0.07 / +0.02 | +0.05 / +0.01 | -0.05 / +0.01 | -0.00 / -0.01 | +0.00 / -0.11 | +0.02 / +0.01 |
| dnamphenoage | -0.11 / -0.01 | -0.27 / -0.26 | -0.08 / +0.02 | -0.14 / -0.18 | +0.15 / -0.06 | +0.27 / +0.24 |
| dunedinpace | -0.02 / -0.01 | -0.24 / -0.30 | -0.06 / -0.06 | -0.04 / -0.06 | +0.11 / +0.08 | +0.19 / +0.25 |

Hannum, DNAmPhenoAge and DunedinPACE accelerations correlate consistently with CD4T (negative) and granulocytes (positive) in both cohorts; IntrinClock has the smallest correlations (|r| ≤ 0.11).
Horvath/SkinAndBlood/Lin correlate with NK in GSE87571 only (+0.15 to +0.24), not in GSE40279 — not a stable finding (the cause was not investigated).

(b) Effect of removing cell composition:

| Cohort | Clock | SD(AA) | SD(IEAA) | R² cells→AA | RMSE raw | RMSE recal. (CV) | RMSE recal.+cells (CV) |
|---|---|---|---|---|---|---|---|
| GSE40279 | horvath2013 | 5.06 | 4.94 | 0.043 | 6.33 | 5.86 | 5.86 |
| GSE40279 | hannum | 4.24 | 3.97 | 0.125 | 6.78 | 4.78 | 4.63 |
| GSE40279 | skinandblood | 4.19 | 4.16 | 0.018 | 5.79 | 5.03 | 5.06 |
| GSE40279 | zhangen | 4.74 | 4.69 | 0.022 | 5.38 | 3.92 | 3.94 |
| GSE40279 | lin | 4.98 | 4.88 | 0.038 | 5.34 | 5.36 | 5.33 |
| GSE40279 | intrinclock | 4.33 | 4.32 | 0.007 | 4.93 | 4.86 | 4.83 |
| GSE40279 | dnamphenoage | 6.96 | 6.56 | 0.111 | 11.05 | 7.74 | 7.66 |
| GSE40279 | dunedinpace | 0.13 | 0.12 | 0.073 | – | – | – |
| GSE87571 | horvath2013 | 5.02 | 4.74 | 0.109 | 5.05 | 4.99 | 4.86 |
| GSE87571 | hannum | 4.50 | 4.18 | 0.140 | 6.14 | 4.66 | 4.47 |
| GSE87571 | skinandblood | 3.58 | 3.46 | 0.064 | 3.63 | 3.54 | 3.48 |
| GSE87571 | zhangen | 3.78 | 3.72 | 0.031 | 8.24 | 2.86 | 2.85 |
| GSE87571 | lin | 6.40 | 6.24 | 0.049 | 7.21 | 5.85 | 5.80 |
| GSE87571 | intrinclock | 3.93 | 3.90 | 0.011 | 4.17 | 3.96 | 3.95 |
| GSE87571 | dnamphenoage | 5.84 | 5.57 | 0.089 | 10.78 | 5.25 | 5.13 |
| GSE87571 | dunedinpace | 0.11 | 0.10 | 0.110 | – | – | – |

Cells explain the most acceleration variance for Hannum (12.5% / 14.0%) and the least for IntrinClock (0.7% / 1.1%). **Rankings change far more with calibration than with cells:** raw-RMSE vs recalibrated-RMSE
rankings have Kendall τ of only 0.07 (GSE40279) and 0.20 (GSE87571) — e.g. ZhangEN goes from last to first on GSE87571 — while adding cells on top of recalibration leaves the ranking unchanged (τ = 1.00 in both cohorts).
The SD(AA)→SD(IEAA) ranking has one swap (Hannum/SkinAndBlood, GSE40279). Using the depositor's cell estimates for GSE87571 gives the same conclusions. Full tables: `results/accel_*.csv`, `results/clock_rankings_before_after_cells.csv`.

## Limitations

- Processed GEO betas were used as published; no BMIQ/noob renormalisation, no batch correction, and the two cohorts were preprocessed differently.
- Training overlap evidence comes from papers/pyaging metadata and is incomplete (`unknown` for SkinAndBlood and IntrinClock); `no_evidence` is not independence.
- Cell estimates are biased (see CD8T and PBMC checks); no cell subtypes. Associations are not causal.
- Relatives in GSE87571 were not removed; samples are not fully independent. Batch (`Sample_Group`) and age distribution are confounded.
- Not run: PC clocks, GrimAge (needs sex/age inputs), ZhangBLUP (~320k CpGs).
- Age acceleration is not lifespan evidence, and a "younger" clock readout does not mean someone is biologically younger.

## Reproducibility

- Pinned environment: `requirements.txt` (+ `requirements.lock.txt`); results regenerated twice with byte-identical key tables (`metrics.csv`, `own_model_metrics.csv`, `accel_adjustment_summary.csv`, `metrics_gse87571_subsets.csv`).
- Randomness is seeded (`random_state=0`).
- Disk: ~4.1 GB raw downloads + ~0.1 GB processed under `data/` (git-ignored). The development virtualenv with CUDA-enabled PyTorch wheels was ~5.9 GB; a CPU-only PyTorch install is much smaller.
- pyaging downloads clock weights from Hugging Face at first use (network needed; set `HF_TOKEN` only if you hit rate limits).
- Per-sample outputs (`results/predictions.csv`, `results/cell_props_*.csv`) are regenerated locally and git-ignored; aggregate tables and plots are tracked.

## Repository layout

```
scripts/   pipeline steps and shared helpers        tests/   lightweight tests (no data needed)
docs/      PLAN.md (pitfalls, next steps), REFERENCES.md        results/  aggregate tables and plots (per-sample tables git-ignored)
run_all.sh / Makefile   one-command rerun            requirements*.txt   pinned environments
```

## License

MIT for the code and documentation (see `LICENSE`). Data and third-party tools keep their own terms and are not redistributed.
Citing: see `CITATION.cff` (author name is a placeholder).
