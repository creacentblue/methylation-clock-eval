"""Shared configuration for the methylation-clock-eval project."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
PLOTS = RES / "plots"
for _p in (RAW, PROC, RES, PLOTS):
    _p.mkdir(parents=True, exist_ok=True)

COHORTS = ["GSE40279", "GSE87571"]

# Clocks evaluated (pyaging names, lowercase). kind = 'age' -> predicts years (MAE etc. are meaningful)
# 'pace'/'other' -> NOT a years-scale chronological-age predictor: only correlation with age is reported.
CLOCKS = {
    "horvath2013": "age",
    "hannum": "age",
    "skinandblood": "age",
    "zhangen": "age",
    "lin": "age",
    "intrinclock": "age",
    "dnamphenoage": "phenotypic_age",   # trained on phenotypic age, not chronological age -> MAE not meaningful
    "dunedinpace": "pace",              # pace of aging (EPIC-trained), not an age in years
}

# Training-overlap flags: does the evaluation cohort appear in the clock's training data?
#   yes        : documented in the clock's paper (source noted)
#   likely_yes : cohort is in the paper's pooled cohort list, exact inclusion in the final model not verified
#   no_evidence: not found in the paper's/pyaging's documentation of training data (NOT a proof of independence)
#   unknown    : not checked / could not verify
# Evidence was gathered on 2026-11-11 from the papers/pyaging metadata; see docs/PLAN.md "训练重叠核查".
TRAINING_OVERLAP = {
    ("hannum", "GSE40279"): ("yes", "Hannum 2013: model derived (N1=482) and validated (N2=174) within GSE40279 itself -> reproduction check only"),
    ("hannum", "GSE87571"): ("no_evidence", "GSE87571 registered at GEO in Oct 2016, after Hannum 2013; not listed in the paper"),
    ("horvath2013", "GSE40279"): ("yes", "Horvath 2013 (Genome Biol 14:R115): training data set 3 = whole-blood data of Hannum et al.; exact sample subset not verified"),
    ("horvath2013", "GSE87571"): ("no_evidence", "GSE87571 deposited 2016, after Horvath 2013; Horvath's Additional file 1 table not machine-checked"),
    ("zhangen", "GSE40279"): ("likely_yes", "Zhang 2019 Genome Med Table 1: GSE40279 is one of the 14 pooled training cohorts; final model trained on 13,566 pooled samples"),
    ("zhangen", "GSE87571"): ("no_evidence", "GSE87571 is not in Zhang 2019 Table 1 (14 cohorts)"),
    ("lin", "GSE40279"): ("yes", "Lin 2016 Aging: 99-CpG model trained on the 656 blood profiles of Hannum et al. (GSE40279) -> reproduction check only; bias ~0 here is consistent with this"),
    ("lin", "GSE87571"): ("no_evidence", "GSE87571 is not among the 12 validation datasets listed in Lin 2016"),
    ("dnamphenoage", "GSE40279"): ("no_evidence", "pyaging notes: coefficients fitted on InCHIANTI 450K blood"),
    ("dnamphenoage", "GSE87571"): ("no_evidence", "pyaging notes: coefficients fitted on InCHIANTI 450K blood"),
    ("dunedinpace", "GSE40279"): ("no_evidence", "trained on Dunedin cohort (pyaging notes); EPIC-based"),
    ("dunedinpace", "GSE87571"): ("no_evidence", "trained on Dunedin cohort (pyaging notes); EPIC-based"),
}
DEFAULT_OVERLAP = ("unknown", "training data not verified for this clock/cohort")

AGE_BINS = [0, 30, 40, 50, 60, 70, 80, 120]

# Reference-based blood cell deconvolution (pyaging implementation of Salas 2018 FlowSorted.Blood.EPIC-based constrained projection, 6 cell types).
# Reference was built for EPIC but its CpGs are mostly shared with 450K; coverage is reported per cohort.
DECONV = {
    "deconvolutebloodepicbcell": "Bcell",
    "deconvolutebloodepiccd4tcell": "CD4T",
    "deconvolutebloodepiccd8tcell": "CD8T",
    "deconvolutebloodepicmonocyte": "Mono",
    "deconvolutebloodepicneutrophil": "Gran",
    "deconvolutebloodepicnkcell": "NK",
}
CELLS = ["CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"]   # main analysis: Reinius-450K-reference Houseman projection (Gran = granulocytes)
