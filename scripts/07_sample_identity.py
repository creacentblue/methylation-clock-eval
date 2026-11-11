"""Why does GSE87571 have 732 GEO samples while Johansson 2013 reports N=421? Check:
  (1) depositor's additional sample characteristics (Sample_Group, year of collection, genotype dosages, smoking, cell estimates),
  (2) repeat individuals: genotype-dosage identity (104 variants) and near-duplicate methylation profiles (centered correlation),
  (3) records that are not individuals (no age / genotype / smoking).
Outputs: data/processed/GSE87571_extra.csv, results/gse87571_identity_summary.csv, results/gse87571_duplicate_candidates.csv"""
import numpy as np, pandas as pd
from config import RAW, PROC, RES

GENO_DIST_DUP = 0.05      # mean |dosage difference| below which two samples are called the same person
METH_CORR_DUP = 0.90      # centered-beta correlation above which two samples are called near-duplicate profiles

x = pd.read_excel(RAW / "GSE87571_additional_sample_chararcteristics.xlsx", header=2).rename(columns={"geo accession": "gsm"}).set_index("gsm")
x.columns = [c.replace("characteristics: ", "") for c in x.columns]
x = x[x.index.astype(str).str.startswith("GSM")]
cols = list(x.columns)
geno_cols = cols[2:cols.index("Adrenomedullin")]            # layout verified: year_of_collection, smoke, <104 genotype columns>, <protein biomarkers ...>
G = x[geno_cols].apply(pd.to_numeric, errors="coerce").to_numpy()
meta = pd.read_csv(PROC / "GSE87571_meta.csv", index_col=0).join(x[["year_of_collection", "smoke", "Sample_Group", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"]])
meta["has_genotype"] = ~np.isnan(G).all(axis=1)
meta["is_individual_record"] = meta.age.notna() & meta.has_genotype

n = len(G)
D = np.full((n, n), np.nan)
for a in range(n):
    with np.errstate(all="ignore"):
        D[a] = np.nanmean(np.abs(G - G[a]), axis=1)
np.fill_diagonal(D, np.nan)
beta = pd.read_pickle(PROC / "GSE87571_beta_needed.pkl").loc[meta.index]
B = beta.to_numpy(float); B = np.where(np.isnan(B), np.nanmean(B, axis=0), B)
Bc = B - B.mean(0); Bc /= np.linalg.norm(Bc, axis=1, keepdims=True)
C = Bc @ Bc.T; np.fill_diagonal(C, np.nan)

iu = np.triu_indices(n, 1)
gd, mc = D[iu], C[iu]
cand = np.where((gd < GENO_DIST_DUP) | (mc > METH_CORR_DUP))[0]
rows = [dict(gsm_a=meta.index[iu[0][t]], gsm_b=meta.index[iu[1][t]], geno_dist=gd[t], meth_corr=mc[t]) for t in cand]
pd.DataFrame(rows, columns=["gsm_a", "gsm_b", "geno_dist", "meth_corr"]).to_csv(RES / "gse87571_duplicate_candidates.csv", index=False)
# closest pairs regardless of threshold (for the record)
order = np.argsort(np.where(np.isnan(gd), 9, gd))[:10]
closest = pd.DataFrame([dict(gsm_a=meta.index[iu[0][t]], gsm_b=meta.index[iu[1][t]], geno_dist=gd[t], meth_corr=mc[t],
                             age_a=meta.age.iloc[iu[0][t]], age_b=meta.age.iloc[iu[1][t]],
                             sex_a=meta.sex.iloc[iu[0][t]], sex_b=meta.sex.iloc[iu[1][t]]) for t in order])
closest.to_csv(RES / "gse87571_closest_genotype_pairs.csv", index=False, float_format="%.3f")

# dedup flag: drop non-individual records, and the later member of any duplicate pair (none expected)
drop = set(meta.index[~meta.is_individual_record])
for r in rows:
    drop.add(r["gsm_b"])
meta["keep_dedup"] = ~meta.index.isin(drop)
meta.to_csv(PROC / "GSE87571_extra.csv")

summ = {
    "geo_samples": n,
    "samples_with_age_and_genotype (individual records)": int(meta.is_individual_record.sum()),
    "records_without_age_genotype_smoking": int((~meta.is_individual_record).sum()),
    "Sample_Group=2012": int((meta.Sample_Group == 2012).sum()),
    "Sample_Group=2013": int((meta.Sample_Group == 2013).sum()),
    "year_of_collection=2006": int((meta.year_of_collection == 2006).sum()),
    "year_of_collection=2009": int((meta.year_of_collection == 2009).sum()),
    "median_age_Sample_Group_2012": float(meta.loc[meta.Sample_Group == 2012, "age"].median()),
    "n_variants_used_for_identity": len(geno_cols),
    "pairs_compared": int(len(gd)),
    "geno_dist_min_over_pairs": float(np.nanmin(gd)),
    "geno_dist_median": float(np.nanmedian(gd)),
    "meth_centered_corr_max": float(np.nanmax(mc)),
    "meth_centered_corr_median": float(np.nanmedian(mc)),
    "duplicate_pairs_found (geno<%.2f or corr>%.2f)" % (GENO_DIST_DUP, METH_CORR_DUP): len(rows),
    "samples_kept_after_dedup": int(meta.keep_dedup.sum()),
}
pd.Series(summ, name="value").to_csv(RES / "gse87571_identity_summary.csv", header=True)
print(pd.Series(summ).to_string())
print(closest.head(5).to_string())
