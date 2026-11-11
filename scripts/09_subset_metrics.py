"""GSE87571 metrics on (a) deduplicated individual records, (b) Sample_Group 2012 (= Johansson 2013's N=421 batch), (c) Sample_Group 2013
(added samples), (d) adults >=19. Also GSE40279 reference row. Output: results/metrics_gse87571_subsets.csv"""
import numpy as np, pandas as pd
from metrics_utils import core
from config import CLOCKS, PROC, RES

pred = pd.read_csv(RES / "predictions.csv")
ex = pd.read_csv(PROC / "GSE87571_extra.csv", index_col=0)
d = pred[pred.cohort == "GSE87571"].set_index("sample").join(ex[["Sample_Group", "keep_dedup", "year_of_collection"]])
d = d[d.age.notna()]
subsets = {
    "all_with_age (n as in metrics.csv)": d,
    "dedup_individuals": d[d.keep_dedup],
    "Sample_Group_2012 (Johansson 2013 N=421 batch)": d[d.Sample_Group == 2012],
    "Sample_Group_2013 (added samples)": d[d.Sample_Group == 2013],
    "adults_age>=19": d[d.age >= 19],
}
rows = []
for name, dd in subsets.items():
    for clock, kind in CLOCKS.items():
        if kind == "pace":
            from scipy import stats
            rows.append(dict(subset=name, clock=clock, n=len(dd), pearson_r=stats.pearsonr(dd.age, dd[clock])[0]))
        else:
            rows.append(dict(subset=name, clock=clock, **core(dd.age.to_numpy(), dd[clock].to_numpy())))
out = pd.DataFrame(rows)
out.to_csv(RES / "metrics_gse87571_subsets.csv", index=False, float_format="%.4f")
print(out[["subset", "clock", "n", "MAE", "RMSE", "pearson_r", "bias_mean_signed_error"]].to_string(float_format=lambda x: f"{x:.3f}"))
