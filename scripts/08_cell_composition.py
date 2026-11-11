"""Blood cell composition vs clock age acceleration.

Cell proportions (main): Houseman-style constrained projection onto a 450K reference built from Reinius 2012 purified leukocytes (GSE35069; 02b/08a),
6 cell types CD8T, CD4T, NK, Bcell, Mono, Gran. Sensitivity (GSE87571 only): the depositor's own estimates (GEO additional sample characteristics).
Diagnostic only: pyaging's EPIC-based DeconvoluteBlood* (57% of its 600 reference CpGs are absent on 450K) -- agreement table only, not used for conclusions.

Definitions (per cohort x clock, chronological-age clocks only for the ranking):
  AA   = age acceleration = residual of  pred ~ a + b*age  (OLS, within cohort; removes the age-dependent bias/compression)
  IEAA = intrinsic-style acceleration = residual of  AA ~ cell proportions (5 cell types; the 6th is dropped because proportions sum to ~1)
  (a) correlation of AA with each cell proportion
  (b) variance of AA explained by cells (R2), SD(AA) vs SD(IEAA)
  Because MAE on chronological age is not the right target after adjusting for cells, we compare rankings on:
     RMSE_raw        : plain RMSE of pred vs age (as in metrics.csv)
     RMSE_cal        : 10-fold CV RMSE of  age ~ a + b*pred                (linear recalibration; removes offset/slope only)
     RMSE_cal_cells  : 10-fold CV RMSE of  age ~ a + b*pred + cells        (does adding cell composition improve age recovery?)
     SD_AA / SD_IEAA : spread of acceleration before/after cell adjustment (lower = less non-age variation)
  Recalibration uses the evaluation labels (cross-validated within the cohort) so it is a diagnostic, not a deployable clock."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold, cross_val_predict
from config import CLOCKS, COHORTS, CELLS, PROC, RES, PLOTS

pred = pd.read_csv(RES / "predictions.csv")
cells = pd.read_csv(RES / "cell_props_houseman450k.csv")
py_cells = pd.read_csv(RES / "cell_props_pyaging.csv")
ex = pd.read_csv(PROC / "GSE87571_extra.csv", index_col=0)
DEPOSITOR = {c: c for c in CELLS}   # depositor columns CD8T CD4T NK Bcell Mono Gran

# ---------- cell estimates: QC / agreement
qc_rows = []
for cohort in COHORTS:
    c = cells[cells.cohort == cohort]
    for cell in CELLS:
        qc_rows.append(dict(cohort=cohort, source="houseman450k_reinius", cell=cell, mean=c[cell].mean(), sd=c[cell].std(), min=c[cell].min(), max=c[cell].max()))
    qc_rows.append(dict(cohort=cohort, source="houseman450k_reinius", cell="SUM_of_6", mean=c[CELLS].sum(axis=1).mean(), sd=c[CELLS].sum(axis=1).std(),
                        min=c[CELLS].sum(axis=1).min(), max=c[CELLS].sum(axis=1).max()))
pd.DataFrame(qc_rows).to_csv(RES / "cell_props_summary.csv", index=False, float_format="%.4f")
c87 = cells[cells.cohort == "GSE87571"].set_index("sample")
p87 = py_cells[py_cells.cohort == "GSE87571"].set_index("sample")
agree = []
for cell, dep in DEPOSITOR.items():
    b = ex.loc[c87.index, dep]
    for name, tab in (("houseman450k_reinius (this project)", c87), ("pyaging_EPIC_deconvolution (diagnostic)", p87)):
        a = tab.loc[c87.index, cell]
        agree.append(dict(cell=cell, estimate=name, pearson_r=stats.pearsonr(a, b)[0], spearman_rho=stats.spearmanr(a, b)[0],
                          mean_estimate=a.mean(), mean_depositor=b.mean(), mean_diff=(a - b).mean()))
pd.DataFrame(agree).to_csv(RES / "cell_estimates_agreement_gse87571.csv", index=False, float_format="%.4f")
print(pd.DataFrame(agree).round(3).to_string())

# ---------- build analysis frames: (cohort, source) -> DataFrame[sample, age, clocks..., cell cols]
frames = {}
for cohort in COHORTS:
    p = pred[pred.cohort == cohort].set_index("sample")
    c = cells[cells.cohort == cohort].set_index("sample")[CELLS]
    d = p.join(c); d = d[d.age.notna()]
    if cohort == "GSE87571":
        d = d.join(ex[["keep_dedup"]]); d = d[d.keep_dedup]
    frames[(cohort, "houseman450k_reinius")] = d
    if cohort == "GSE87571":
        dd = d.copy()
        for cell, dep in DEPOSITOR.items():
            dd[cell] = ex.loc[dd.index, dep].astype(float)
        frames[(cohort, "depositor_houseman")] = dd

ADJ = [c for c in CELLS if c != "Gran"]            # 5 cells in regressions (proportions are ~compositional)
kf = KFold(10, shuffle=True, random_state=0)


def resid(y, X):
    return y - LinearRegression().fit(X, y).predict(X)


def cv_rmse(X, y):
    return float(np.sqrt(np.mean((y - cross_val_predict(LinearRegression(), X, y, cv=kf)) ** 2)))


corr_rows, adj_rows = [], []
for (cohort, src), d in frames.items():
    age = d.age.to_numpy()
    for cell in CELLS:   # correlation of the cell with chronological age
        corr_rows.append(dict(cohort=cohort, cell_source=src, clock="(chronological age)", cell=cell, pearson_r=stats.pearsonr(age, d[cell])[0],
                              spearman_rho=stats.spearmanr(age, d[cell])[0], n=len(d)))
    for clock, kind in CLOCKS.items():
        pr = d[clock].to_numpy()
        AA = resid(pr, age.reshape(-1, 1))
        for cell in CELLS:
            r, pv = stats.pearsonr(AA, d[cell])
            corr_rows.append(dict(cohort=cohort, cell_source=src, clock=clock, cell=cell, pearson_r=r, pearson_p=pv,
                                  spearman_rho=stats.spearmanr(AA, d[cell])[0], n=len(d)))
        Xc = d[ADJ].to_numpy()
        IEAA = resid(AA, Xc)
        r2 = 1 - np.var(IEAA) / np.var(AA)
        row = dict(cohort=cohort, cell_source=src, clock=clock, clock_kind=kind, n=len(d), SD_AA=AA.std(), SD_IEAA=IEAA.std(), R2_AA_on_cells=r2)
        if kind != "pace":
            row.update(RMSE_raw=float(np.sqrt(np.mean((pr - age) ** 2))), RMSE_cal=cv_rmse(pr.reshape(-1, 1), age),
                       RMSE_cal_cells=cv_rmse(np.column_stack([pr, Xc]), age))
        adj_rows.append(row)
C = pd.DataFrame(corr_rows); C.to_csv(RES / "accel_cell_correlations.csv", index=False, float_format="%.4f")
A = pd.DataFrame(adj_rows)

# rankings among the 6 chronological-age clocks
age_clocks = [c for c, k in CLOCKS.items() if k == "age"]
rank_rows = []
for (cohort, src), g in A.groupby(["cohort", "cell_source"]):
    g = g[g.clock.isin(age_clocks)].set_index("clock")
    rk = pd.DataFrame({m: g[m].rank(method="min").astype(int) for m in ["RMSE_raw", "RMSE_cal", "RMSE_cal_cells", "SD_AA", "SD_IEAA"]})
    rk.insert(0, "cell_source", src); rk.insert(0, "cohort", cohort)
    rank_rows.append(rk.reset_index())
    for m1, m2 in [("RMSE_raw", "RMSE_cal"), ("RMSE_cal", "RMSE_cal_cells"), ("SD_AA", "SD_IEAA"), ("RMSE_raw", "SD_IEAA")]:
        A.loc[(A.cohort == cohort) & (A.cell_source == src), f"kendall_tau[{m1}|{m2}]"] = stats.kendalltau(g[m1], g[m2])[0]
R = pd.concat(rank_rows)
A.to_csv(RES / "accel_adjustment_summary.csv", index=False, float_format="%.4f")
R.to_csv(RES / "clock_rankings_before_after_cells.csv", index=False)

# ---------- plots
for src_label, srcs in (("pyaging", [("GSE40279", "houseman450k_reinius"), ("GSE87571", "houseman450k_reinius")]),):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, key in zip(axs, srcs):
        sub = C[(C.cohort == key[0]) & (C.cell_source == key[1]) & (C.clock != "(chronological age)")]
        M = sub.pivot(index="clock", columns="cell", values="pearson_r").loc[list(CLOCKS), CELLS]
        im = ax.imshow(M.to_numpy(), cmap="RdBu_r", vmin=-0.6, vmax=0.6)
        ax.set_xticks(range(len(CELLS))); ax.set_xticklabels(CELLS); ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=8)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, f"{M.iloc[i, j]:.2f}", ha="center", va="center", fontsize=7)
        ax.set_title(f"{key[0]}: corr(age acceleration, est. cell proportion)", fontsize=9)
    fig.colorbar(im, ax=axs, shrink=0.8); fig.savefig(PLOTS / "accel_cell_corr_heatmap.png", dpi=130, bbox_inches="tight")
fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for ax, cohort in zip(axs, COHORTS):
    g = A[(A.cohort == cohort) & (A.cell_source == "houseman450k_reinius") & A.clock.isin(age_clocks)].set_index("clock")
    x = np.arange(len(g)); w = 0.27
    ax.bar(x - w, g.RMSE_raw, w, label="RMSE raw"); ax.bar(x, g.RMSE_cal, w, label="RMSE recalibrated (CV)"); ax.bar(x + w, g.RMSE_cal_cells, w, label="recalibrated + cells (CV)")
    ax.set_xticks(x); ax.set_xticklabels(g.index, rotation=30, fontsize=8); ax.set_title(cohort); ax.set_ylabel("RMSE (years)")
axs[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(PLOTS / "rmse_raw_vs_recalibrated_vs_cells.png", dpi=130)

pd.set_option("display.width", 220)
print(A[["cohort", "cell_source", "clock", "SD_AA", "SD_IEAA", "R2_AA_on_cells", "RMSE_raw", "RMSE_cal", "RMSE_cal_cells"]].round(3).to_string())
print(R.to_string())
