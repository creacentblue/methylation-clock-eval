"""Stretch: ElasticNetCV trained on GSE40279 only, frozen, tested on GSE87571 without refitting.
Order: (1) split GSE40279 80/20 (seeded) -> (2) imputer+scaler fit on TRAIN only -> (3) ElasticNetCV on train ->
(4) evaluate internal test + external GSE87571. Feature panel = CpGs of the evaluated clocks present in both cohorts (fixed a priori,
not chosen using any test/external label). Caveat: that panel was itself chosen by published studies (some used GSE40279)."""
import json, numpy as np, pandas as pd, joblib
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from scipy import stats
from config import PROC, RES

SEED = 0
bA = pd.read_pickle(PROC / "GSE40279_beta_needed.pkl"); mA = pd.read_csv(PROC / "GSE40279_meta.csv", index_col=0)
bB = pd.read_pickle(PROC / "GSE87571_beta_needed.pkl"); mB = pd.read_csv(PROC / "GSE87571_meta.csv", index_col=0)
clock_cpgs = set(sum(json.load(open(PROC / 'clock_features.json')).values(), []))   # deconvolution CpGs deliberately excluded
panel = [c for c in bA.columns if c in bB.columns and c in clock_cpgs]
okA = mA.age.notna(); okB = mB.age.notna()
XA, yA = bA.loc[okA, panel], mA.loc[okA, "age"]
XB, yB = bB.loc[okB, panel], mB.loc[okB, "age"]
XA_tr, XA_te, yA_tr, yA_te = train_test_split(XA, yA, test_size=0.2, random_state=SEED)
pipe = make_pipeline(SimpleImputer(strategy="mean"), StandardScaler(),
                     ElasticNetCV(l1_ratio=[0.5, 0.9, 1.0], cv=5, alphas=40, max_iter=5000, n_jobs=-1, random_state=SEED))
pipe.fit(XA_tr.to_numpy(np.float64), yA_tr.to_numpy())
en = pipe[-1]
def met(name, X, y):
    p = pipe.predict(X.to_numpy(np.float64)); r = p - y.to_numpy()
    return dict(eval_set=name, n=len(y), MAE=np.mean(np.abs(r)), median_AE=np.median(np.abs(r)), RMSE=np.sqrt(np.mean(r**2)),
                R2=1 - np.sum(r**2) / np.sum((y - y.mean())**2), pearson_r=stats.pearsonr(y, p)[0], bias_mean_signed_error=np.mean(r))
rows = [met("GSE40279 train (in-sample)", XA_tr, yA_tr), met("GSE40279 held-out 20%", XA_te, yA_te),
        met("GSE87571 external (all)", XB, yB), met("GSE87571 external (age>=19)", XB[yB >= 19], yB[yB >= 19])]
out = pd.DataFrame(rows).assign(n_panel_cpgs=len(panel), n_nonzero_coef=int((en.coef_ != 0).sum()),
                                alpha=en.alpha_, l1_ratio=en.l1_ratio_)
out.to_csv(RES / "own_model_metrics.csv", index=False, float_format="%.4f")
joblib.dump(pipe, PROC / "elasticnet_gse40279.joblib")
print(out.to_string(float_format=lambda x: f"{x:.3f}"))
