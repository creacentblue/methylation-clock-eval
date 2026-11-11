"""Schema checks on the tracked result tables (skipped if a table has not been generated)."""
import pandas as pd
import pytest

import config

METRIC_COLS = {"cohort", "clock", "clock_kind", "n", "MAE", "median_AE", "RMSE", "R2", "pearson_r", "bias_mean_signed_error",
               "frac_cpgs_missing", "training_overlap", "training_overlap_note"}


def _read(name):
    p = config.RES / name
    if not p.exists():
        pytest.skip(f"{name} not generated")
    return pd.read_csv(p)


def test_metrics_schema():
    m = _read("metrics.csv")
    assert METRIC_COLS <= set(m.columns)
    assert set(m.cohort) == set(config.COHORTS) and set(m.clock) == set(config.CLOCKS)
    assert len(m) == len(config.COHORTS) * len(config.CLOCKS)
    assert m.frac_cpgs_missing.between(0, 1).all() and m.pearson_r.between(-1, 1).all()
    assert set(m.training_overlap) <= {"yes", "likely_yes", "no_evidence", "unknown"}
    age = m[m.clock_kind == "age"]
    assert age[["MAE", "RMSE", "R2"]].notna().all().all() and (age.MAE <= age.RMSE).all()


def test_own_model_metrics_schema():
    m = _read("own_model_metrics.csv")
    assert {"eval_set", "n", "MAE", "RMSE", "R2"} <= set(m.columns)
    assert m.eval_set.str.contains("GSE87571 external").any()


def test_rankings_are_permutations():
    r = _read("clock_rankings_before_after_cells.csv")
    for _, g in r.groupby(["cohort", "cell_source"]):
        for col in ["RMSE_raw", "RMSE_cal", "RMSE_cal_cells", "SD_AA", "SD_IEAA"]:
            assert g[col].min() == 1 and g[col].max() <= len(g)      # competition ranking (ties allowed)
