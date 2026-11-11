import numpy as np
import pytest

from metrics_utils import core


def test_perfect_prediction():
    age = np.array([20.0, 35.0, 50.0, 65.0, 80.0])
    m = core(age, age.copy())
    assert m["MAE"] == 0 and m["RMSE"] == 0 and m["bias_mean_signed_error"] == 0
    assert m["R2"] == pytest.approx(1.0) and m["pearson_r"] == pytest.approx(1.0)


def test_constant_offset_is_bias_not_noise():
    age = np.array([20.0, 35.0, 50.0, 65.0, 80.0])
    m = core(age, age + 3.0)
    assert m["MAE"] == pytest.approx(3.0) and m["bias_mean_signed_error"] == pytest.approx(3.0)
    assert m["MAE_after_removing_bias"] == pytest.approx(0.0, abs=1e-12)
    assert m["pearson_r"] == pytest.approx(1.0)


def test_known_small_example():
    age = np.array([10.0, 20.0, 30.0, 40.0])
    pred = np.array([12.0, 18.0, 33.0, 37.0])      # errors +2, -2, +3, -3
    m = core(age, pred)
    assert m["n"] == 4
    assert m["MAE"] == pytest.approx(2.5)
    assert m["median_AE"] == pytest.approx(2.5)
    assert m["RMSE"] == pytest.approx(np.sqrt((4 + 4 + 9 + 9) / 4))
    assert m["bias_mean_signed_error"] == pytest.approx(0.0)
