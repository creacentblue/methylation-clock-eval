"""Stream-parse GEO processed matrices -> (a) beta matrix restricted to CpGs needed by the clocks,
(b) sample metadata (age, sex, ...), (c) QC summary over ALL probes. Never loads the full matrix into memory."""
import gzip, json, re
import numpy as np, pandas as pd
from config import RAW, PROC, RES

NEEDED = set()          # CpGs to keep; filled by load_needed() (kept out of import time so the module is testable without data)
CLOCK_CPGS = set()


def load_needed():
    """Load clock / deconvolution / cell-reference CpG lists produced by 02_clock_features.py and 02b_cell_reference.py."""
    feats = json.load(open(PROC / "clock_features.json"))
    dfeats = json.load(open(PROC / "deconv_features.json"))
    ref = set(json.load(open(PROC / "cell_ref_probes.json")))      # Reinius-based 450K cell reference probes (02b)
    CLOCK_CPGS.update(sum(feats.values(), []))
    NEEDED.update(CLOCK_CPGS | set(sum(dfeats.values(), [])) | ref)


CHUNK = 20000


def series_meta(path):
    """Parse '!Sample_*' lines of a GEO series matrix -> DataFrame indexed by GSM."""
    cols, chars = {}, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            k, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip('"') for v in vals]
            if k == "!Sample_characteristics_ch1":
                chars.append(vals)
            else:
                cols.setdefault(k[len("!Sample_"):], vals)
    df = pd.DataFrame({k: v for k, v in cols.items() if k in ("title", "geo_accession", "source_name_ch1")})
    for vals in chars:
        key = vals[0].split(":")[0].strip()
        df[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
    return df.set_index("geo_accession")


def process(chunks, sample_ids):
    keep, nprobe = [], 0
    na = np.zeros(len(sample_ids)); s = np.zeros(len(sample_ids)); cnt = np.zeros(len(sample_ids))
    for ch in chunks:
        ch = ch.apply(pd.to_numeric, errors="coerce")
        keep.append(ch.loc[ch.index.isin(NEEDED)].astype(np.float32))   # includes non-cg probes (e.g. 'ch.2.30415474F') some clocks use
        ch = ch[ch.index.astype(str).str.startswith("cg")]               # QC statistics over cg probes only
        X = ch.to_numpy(dtype=np.float64)
        nprobe += X.shape[0]
        na += np.isnan(X).sum(0); s += np.nansum(X, 0); cnt += (~np.isnan(X)).sum(0)
    beta = pd.concat(keep).T          # samples x CpGs
    beta = beta.loc[sample_ids]
    qc = pd.DataFrame({"n_probes_in_matrix": nprobe, "frac_na_all_probes": na / nprobe,
                       "mean_beta_all_probes": s / np.maximum(cnt, 1)}, index=sample_ids)
    return beta, qc


def parse_gse40279():
    meta = series_meta(RAW / "GSE40279_series_matrix.txt.gz")
    meta = meta.rename(columns={"age (y)": "age", "gender": "sex"})
    meta["age"] = pd.to_numeric(meta["age"], errors="coerce")
    meta["sex"] = meta["sex"].map({"F": "F", "M": "M"})
    n_hdr = 0
    with gzip.open(RAW / "GSE40279_series_matrix.txt.gz", "rt") as f:
        for i, line in enumerate(f):
            if line.startswith("!series_matrix_table_begin"):
                n_hdr = i + 1; break
    def chunks():
        for ch in pd.read_csv(RAW / "GSE40279_series_matrix.txt.gz", sep="\t", index_col=0, skiprows=n_hdr,
                              chunksize=CHUNK, quotechar='"', low_memory=False):
            ch = ch[~ch.index.astype(str).str.startswith("!")]
            yield ch
    beta, qc = process(chunks(), list(meta.index))
    return meta, beta, qc


def parse_gse87571():
    meta = series_meta(RAW / "GSE87571_series_matrix.txt.gz")
    meta["name"] = meta["title"].str.split(" ").str[0]          # 'X1 genomic DNA ...' -> 'X1' (column name in suppl matrices)
    meta = meta.rename(columns={"gender": "sex"})
    meta["age"] = pd.to_numeric(meta["age"], errors="coerce")
    meta["sex"] = meta["sex"].map({"Female": "F", "Male": "M"})
    name2gsm = dict(zip(meta["name"], meta.index))
    betas, qcs = [], []
    for fn in ("GSE87571_matrix1of2.txt.gz", "GSE87571_matrix2of2.txt.gz"):
        p = RAW / fn
        hdr = pd.read_csv(p, sep="\t", nrows=0, index_col=0).columns
        # layout verified: for every sample 'Xn' = beta, 'Xn.1' = detection p-value
        bcols = [c for c in hdr if c in name2gsm]
        ids = [name2gsm[c] for c in bcols]
        def chunks():
            for ch in pd.read_csv(p, sep="\t", index_col=0, usecols=["ID_REF"] + bcols, chunksize=CHUNK, low_memory=False):
                ch.columns = [name2gsm[c] for c in ch.columns]
                yield ch
        b, q = process(chunks(), ids)
        betas.append(b); qcs.append(q)
    beta = pd.concat(betas).loc[meta.index]; qc = pd.concat(qcs).loc[meta.index]
    return meta, beta, qc


if __name__ == "__main__":
    load_needed()
    summ = []
    for name, fn in (("GSE40279", parse_gse40279), ("GSE87571", parse_gse87571)):
        print("parsing", name, flush=True)
        meta, beta, qc = fn()
        meta = meta.join(qc)
        meta["needed_cpg_na_frac"] = beta.isna().mean(axis=1)
        meta.to_csv(PROC / f"{name}_meta.csv")
        beta.to_pickle(PROC / f"{name}_beta_needed.pkl")
        present = [c for c in NEEDED if c in beta.columns]
        summ.append(dict(cohort=name, n_samples=len(meta), n_age_missing=int(meta.age.isna().sum()),
                         age_min=meta.age.min(), age_median=meta.age.median(), age_max=meta.age.max(),
                         n_female=int((meta.sex == "F").sum()), n_male=int((meta.sex == "M").sum()),
                         n_probes_in_matrix=int(qc.n_probes_in_matrix.iloc[0]),
                         union_needed_cpgs=len(NEEDED), union_needed_cpgs_present=len(present),
                         clock_cpgs=len(CLOCK_CPGS), clock_cpgs_present=len([c for c in CLOCK_CPGS if c in beta.columns]),
                         median_sample_na_frac_all_probes=float(qc.frac_na_all_probes.median()),
                         max_sample_na_frac_all_probes=float(qc.frac_na_all_probes.max()),
                         mean_beta_all_probes_median=float(qc.mean_beta_all_probes.median())))
        print(summ[-1], flush=True)
    pd.DataFrame(summ).to_csv(RES / "qc_summary.csv", index=False)
