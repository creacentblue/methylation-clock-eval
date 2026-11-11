"""Run pyaging clocks on each cohort's beta matrix. Writes results/predictions.csv and results/missing_cpgs.csv.
Missing-feature fraction is computed here from the clock's CpG list vs. the cohort matrix (pyaging itself silently/with a warning
substitutes its reference value or 0 for absent CpGs). Remaining NaN betas within present CpGs are mean-imputed (per CpG, over the cohort)."""
import io, json, contextlib, warnings
import numpy as np, pandas as pd
import pyaging as pya
from config import CLOCKS, COHORTS, DECONV, PROC, RES

warnings.filterwarnings("ignore")
feats = json.load(open(PROC / "clock_features.json"))
needed = sorted(set(sum(feats.values(), [])))
dfeats = json.load(open(PROC / 'deconv_features.json'))
preds, miss, cells = [], [], []
for cohort in COHORTS:
    beta = pd.read_pickle(PROC / f"{cohort}_beta_needed.pkl")
    meta = pd.read_csv(PROC / f"{cohort}_meta.csv", index_col=0)
    for clock, f in feats.items():
        present = [c for c in f if c in beta.columns]
        sub = beta[present]
        miss.append(dict(cohort=cohort, clock=clock, n_features=len(f), n_missing_cpgs=len(f) - len(present),
                         frac_cpgs_missing=(len(f) - len(present)) / len(f),
                         frac_na_values_in_present_cpgs=float(sub.isna().to_numpy().mean()) if len(present) else np.nan))
    for clock, f in dfeats.items():
        present = [c for c in f if c in beta.columns]
        miss.append(dict(cohort=cohort, clock=clock, n_features=len(f), n_missing_cpgs=len(f) - len(present),
                         frac_cpgs_missing=(len(f) - len(present)) / len(f),
                         frac_na_values_in_present_cpgs=float(beta[present].isna().to_numpy().mean())))
    df = beta.copy()
    df.index.name = None
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        adata = pya.pp.df_to_adata(df, imputer_strategy="mean", verbose=True)
    log = buf.getvalue()
    out = pd.DataFrame(index=beta.index)
    for clock in CLOCKS:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            pya.pred.predict_age(adata, clock, verbose=True)
        log += f"\n===== {clock} =====\n" + buf.getvalue()
        out[clock] = adata.obs[clock].astype(float).values
        print(cohort, clock, "done", flush=True)
    cp = pd.DataFrame(index=beta.index)
    for dname, cell in DECONV.items():
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            pya.pred.predict_age(adata, dname, verbose=True)
        log += f"\n===== {dname} =====\n" + buf.getvalue()
        cp[cell] = adata.obs[dname].astype(float).values
    cells.append(cp.assign(cohort=cohort).reset_index(names="sample"))
    open(RES / f"pyaging_log_{cohort}.txt", "w").write(log)
    out = out.join(meta[["age", "sex"]]).assign(cohort=cohort)
    preds.append(out.reset_index(names="sample"))
pd.concat(preds).to_csv(RES / "predictions.csv", index=False)
pd.concat(cells).to_csv(RES / "cell_props_pyaging.csv", index=False)
pd.DataFrame(miss).to_csv(RES / "missing_cpgs.csv", index=False)
print(pd.DataFrame(miss).to_string())
