"""Build a 450K blood cell-type reference (Houseman-style) from Reinius et al. 2012 purified leukocytes (GEO GSE35069, 6 male donors).
Probe selection mimics minfi::pickCompProbes (as used by estimateCellCounts / FlowSorted.Blood.450k):
  one-way F-test across the 6 cell types (p < 1e-8); for each cell type the 50 most hypermethylated and 50 most hypomethylated probes
  (cell type mean minus mean of the other cell types) among the significant ones; union of selected probes = reference.
Cell types: CD8T, CD4T, NK(CD56+), Bcell(CD19+), Mono(CD14+), Gran(sorted granulocytes). Input are the GEO series-matrix betas as deposited
(NOT re-normalized; minfi would quantile-normalize reference and target jointly) -> limitation documented in README."""
import gzip, json
import numpy as np, pandas as pd
from scipy import stats
from config import RAW, PROC, RES

TYPES = {"CD4+": "CD4T", "CD8+": "CD8T", "CD56+": "NK", "CD19+": "Bcell", "CD14+": "Mono", "Gran": "Gran"}
path = RAW / "GSE35069_series_matrix.txt.gz"
titles = None; n_hdr = 0
with gzip.open(path, "rt") as f:
    for i, line in enumerate(f):
        if line.startswith("!Sample_title"):
            titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_geo_accession"):
            gsms = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!series_matrix_table_begin"):
            n_hdr = i + 1; break
lab = {g: TYPES[t.rsplit("_", 1)[0]] for g, t in zip(gsms, titles) if t.rsplit("_", 1)[0] in TYPES}
df = pd.read_csv(path, sep="\t", skiprows=n_hdr, index_col=0, quotechar='"', low_memory=False)
df = df[df.index.astype(str).str.startswith("cg")]
df = df[list(lab)].dropna()
types = np.array([lab[c] for c in df.columns]); order = ["CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"]
X = df.to_numpy(float)
F, p = stats.f_oneway(*[X[:, types == t] for t in order], axis=1)
means = pd.DataFrame({t: X[:, types == t].mean(1) for t in order}, index=df.index)
sel = set()
sig = p < 1e-8
for t in order:
    diff = means[t] - means[[o for o in order if o != t]].mean(axis=1)
    d = diff[sig]
    sel |= set(d.nlargest(50).index) | set(d.nsmallest(50).index)
ref = means.loc[sorted(sel)]
ref.to_csv(RES / "cell_reference_450k_reinius.csv", float_format="%.5f")
json.dump(list(ref.index), open(PROC / "cell_ref_probes.json", "w"))
print("probes in reference matrix:", df.shape[0], "samples:", df.shape[1], "significant (F p<1e-8):", int(sig.sum()), "selected probes:", len(ref))

# sanity-check material: Reinius whole blood (WB) and PBMC samples (not used to build the reference), selected probes only
other = {g: t for g, t in zip(gsms, titles) if t.startswith("WB_") or t.startswith("PBMC_")}
full = pd.read_csv(path, sep="\t", skiprows=n_hdr, index_col=0, quotechar='"', low_memory=False, usecols=["ID_REF"] + list(other))
full = full.loc[full.index.isin(ref.index)]
full.columns = [other[c] for c in full.columns]
full.to_csv(PROC / "reinius_wb_pbmc_selected.csv", float_format="%.5f")
