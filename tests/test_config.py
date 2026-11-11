import config

ALLOWED_FLAGS = {"yes", "likely_yes", "no_evidence", "unknown"}


def test_paths_are_inside_repo():
    for p in (config.RAW, config.PROC, config.RES, config.PLOTS):
        assert config.ROOT in p.parents or p == config.ROOT


def test_clock_kinds_and_cells():
    assert set(config.CLOCKS.values()) <= {"age", "phenotypic_age", "pace"}
    assert len(config.CELLS) == 6 and len(set(config.CELLS)) == 6
    assert set(config.DECONV.values()) == set(config.CELLS)


def test_training_overlap_table_is_consistent():
    for (clock, cohort), (flag, note) in config.TRAINING_OVERLAP.items():
        assert clock in config.CLOCKS and cohort in config.COHORTS
        assert flag in ALLOWED_FLAGS and note
    assert config.DEFAULT_OVERLAP[0] in ALLOWED_FLAGS


def test_age_bins_increasing():
    assert config.AGE_BINS == sorted(config.AGE_BINS)
