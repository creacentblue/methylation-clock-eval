"""Collect metadata + CpG feature lists for the evaluated pyaging clocks (weights come from the pyaging HF repo)."""
import json, warnings
import pandas as pd
import pyaging as pya
from pyaging.predict._pred_utils import load_clock
from config import CLOCKS, DECONV, PROC, RES

warnings.filterwarnings("ignore")
rows, feats = [], {}
for name, kind in CLOCKS.items():
    m = load_clock(name, verbose=False)
    md = m.metadata
    feats[name] = list(m.features)
    rows.append(dict(clock=name, kind=kind, n_features=len(m.features), year=md.get("year"),
                     tissue="|".join(md.get("tissue", [])), platform="|".join(md.get("platform", [])),
                     predicts="|".join(md.get("predicts", [])), unit="|".join(md.get("unit", [])),
                     doi=md.get("doi"), citation=md.get("citation"), pyaging_notes=md.get("notes")))
pd.DataFrame(rows).to_csv(RES / "clock_info.csv", index=False)
json.dump(feats, open(PROC / "clock_features.json", "w"))
union = sorted(set(sum(feats.values(), [])))
print("clocks:", {k: len(v) for k, v in feats.items()}, "union CpGs:", len(union))

# deconvolution references (kept separate so the ElasticNet panel and clock tables are unaffected)
dfeats = {n: list(load_clock(n, verbose=False).features) for n in DECONV}
json.dump(dfeats, open(PROC / "deconv_features.json", "w"))
print("deconv:", {k: len(v) for k, v in dfeats.items()}, "union:", len(set(sum(dfeats.values(), []))))
