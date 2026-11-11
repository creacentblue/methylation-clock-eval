"""Estimate blood cell proportions for both cohorts with the Reinius-based 450K reference (02b) + Houseman constrained projection.
Outputs results/cell_props_houseman450k.csv, results/cell_reference_coverage.csv, results/cell_method_sanity.csv"""
import numpy as np, pandas as pd
from scipy import stats
from cell_deconv import project
from config import COHORTS, PROC, RES

ref = pd.read_csv(RES / "cell_reference_450k_reinius.csv", index_col=0)
order = list(ref.columns)              # CD8T CD4T NK Bcell Mono Gran
out, cov = [], []
for cohort in COHORTS:
    beta = pd.read_pickle(PROC / f"{cohort}_beta_needed.pkl")
    probes = [p for p in ref.index if p in beta.columns]
    cov.append(dict(cohort=cohort, ref_probes=len(ref), present_in_cohort=len(probes), frac_missing=1 - len(probes) / len(ref)))
    Y = beta[probes]
    Y = Y.fillna(Y.mean()).to_numpy(float)
    P = project(Y, ref.loc[probes].to_numpy(float))
    out.append(pd.DataFrame(P, columns=order, index=beta.index).assign(cohort=cohort).reset_index(names="sample"))
pd.concat(out).to_csv(RES / "cell_props_houseman450k.csv", index=False, float_format="%.5f")
pd.DataFrame(cov).to_csv(RES / "cell_reference_coverage.csv", index=False)
print(pd.DataFrame(cov))

# sanity: Reinius whole blood / PBMC (not in reference); PBMC should have ~0 granulocytes, WB mostly Gran
wb = pd.read_csv(PROC / "reinius_wb_pbmc_selected.csv", index_col=0)
P = project(wb.T.loc[:, ref.index].to_numpy(float), ref.to_numpy(float))
san = pd.DataFrame(P, columns=order, index=wb.columns)
san["sum"] = san.sum(axis=1)
san["group"] = [i.split("_")[0] for i in san.index]
sanity = san.groupby("group").mean().round(3)
sanity.to_csv(RES / "cell_method_sanity.csv")
print(sanity)
