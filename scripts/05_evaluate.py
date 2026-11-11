"""Metrics (results/metrics.csv), age-bin and sex breakdowns, plots."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from metrics_utils import core
from config import CLOCKS, COHORTS, AGE_BINS, DEFAULT_OVERLAP, TRAINING_OVERLAP, RES, PLOTS

pred = pd.read_csv(RES / "predictions.csv")
miss = pd.read_csv(RES / "missing_cpgs.csv").set_index(["cohort", "clock"])


rows, bins_rows, sex_rows = [], [], []
for cohort in COHORTS:
    d = pred[pred.cohort == cohort]
    d = d[d.age.notna()]
    for clock, kind in CLOCKS.items():
        dd = d[d[clock].notna()]
        m = core(dd.age.to_numpy(), dd[clock].to_numpy())
        flag, note = TRAINING_OVERLAP.get((clock, cohort), DEFAULT_OVERLAP)
        if kind == "pace":   # not a years-scale age: report only correlation with age
            m = dict(n=m["n"], pearson_r=m["pearson_r"], spearman_rho=m["spearman_rho"], slope_pred_on_age=m["slope_pred_on_age"])
        rows.append(dict(cohort=cohort, clock=clock, clock_kind=kind,
                         age_scale_metrics_meaningful=(kind == "age"), **m,
                         frac_cpgs_missing=miss.loc[(cohort, clock), "frac_cpgs_missing"],
                         n_cpgs_missing=miss.loc[(cohort, clock), "n_missing_cpgs"],
                         n_cpgs_total=miss.loc[(cohort, clock), "n_features"],
                         training_overlap=flag, training_overlap_note=note))
        if kind == "pace":
            continue
        dd = dd.assign(err=dd[clock] - dd.age, ae=(dd[clock] - dd.age).abs(),
                       age_bin=pd.cut(dd.age, AGE_BINS, right=False))
        for b, g in dd.groupby("age_bin", observed=True):
            bins_rows.append(dict(cohort=cohort, clock=clock, age_bin=str(b), n=len(g), MAE=g.ae.mean(), bias=g.err.mean(),
                                  RMSE=np.sqrt((g.err ** 2).mean())))
        for s, g in dd.groupby("sex"):
            sex_rows.append(dict(cohort=cohort, clock=clock, sex=s, n=len(g), MAE=g.ae.mean(), bias=g.err.mean(),
                                 RMSE=np.sqrt((g.err ** 2).mean()), pearson_r=stats.pearsonr(g.age, g[clock])[0]))
M = pd.DataFrame(rows)
M.to_csv(RES / "metrics.csv", index=False, float_format="%.4f")
pd.DataFrame(bins_rows).to_csv(RES / "residuals_by_age_bin.csv", index=False, float_format="%.3f")
pd.DataFrame(sex_rows).to_csv(RES / "residuals_by_sex.csv", index=False, float_format="%.3f")

# ---- plots
age_clocks = [c for c, k in CLOCKS.items() if k in ("age", "phenotypic_age")]
fig, ax = plt.subplots(len(COHORTS), len(age_clocks), figsize=(3.1 * len(age_clocks), 3.2 * len(COHORTS)), squeeze=False)
fig2, ax2 = plt.subplots(len(COHORTS), len(age_clocks), figsize=(3.1 * len(age_clocks), 3.2 * len(COHORTS)), squeeze=False)
for i, cohort in enumerate(COHORTS):
    d = pred[(pred.cohort == cohort) & pred.age.notna()]
    for j, clock in enumerate(age_clocks):
        r = M[(M.cohort == cohort) & (M.clock == clock)].iloc[0]
        a = ax[i, j]; a.scatter(d.age, d[clock], s=5, alpha=.5)
        lim = [0, 110]; a.plot(lim, lim, "k--", lw=.8); a.set_xlim(lim); a.set_ylim(-15, 115)
        a.set_title(f"{clock}\n{cohort}  MAE={r.MAE:.1f} r={r.pearson_r:.2f}\noverlap: {r.training_overlap}", fontsize=8)
        a.set_xlabel("chronological age"); a.set_ylabel("predicted")
        b = ax2[i, j]; b.scatter(d.age, d[clock] - d.age, s=5, alpha=.5); b.axhline(0, color="k", lw=.8)
        b.set_title(f"{clock} / {cohort}\nbias={r.bias_mean_signed_error:.1f}", fontsize=8)
        b.set_xlabel("chronological age"); b.set_ylabel("predicted - chronological")
fig.tight_layout(); fig.savefig(PLOTS / "pred_vs_age.png", dpi=130)
fig2.tight_layout(); fig2.savefig(PLOTS / "residual_vs_age.png", dpi=130)
fig3, ax3 = plt.subplots(1, len(COHORTS), figsize=(4.5 * len(COHORTS), 3.8), squeeze=False)
for i, cohort in enumerate(COHORTS):
    d = pred[(pred.cohort == cohort) & pred.age.notna()]
    r = M[(M.cohort == cohort) & (M.clock == "dunedinpace")].iloc[0]
    ax3[0, i].scatter(d.age, d.dunedinpace, s=5, alpha=.5)
    ax3[0, i].set_title(f"DunedinPACE / {cohort}  r={r.pearson_r:.2f}", fontsize=9)
    ax3[0, i].set_xlabel("chronological age"); ax3[0, i].set_ylabel("DunedinPACE (pace, not years)")
fig3.tight_layout(); fig3.savefig(PLOTS / "dunedinpace_vs_age.png", dpi=130)
cols = ["cohort", "clock", "n", "MAE", "median_AE", "RMSE", "R2", "pearson_r", "bias_mean_signed_error", "frac_cpgs_missing", "training_overlap"]
print(M[cols].to_string(float_format=lambda x: f"{x:.3f}"))
