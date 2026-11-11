import numpy as np

from cell_deconv import project


def _mixture(seed=0, n_probes=300, k=6, n=20, noise=0.01):
    rng = np.random.default_rng(seed)
    R = rng.uniform(0.05, 0.95, size=(n_probes, k))
    W = rng.dirichlet(np.ones(k) * 2, size=n)
    Y = W @ R.T + rng.normal(0, noise, size=(n, n_probes))
    return R, W, Y


def test_projection_recovers_known_mixture():
    R, W, Y = _mixture()
    P = project(Y, R)
    assert np.abs(P - W).max() < 0.03


def test_projection_respects_constraints():
    R, _, Y = _mixture(seed=1, noise=0.2)           # noisy input must still give a valid composition
    P = project(Y, R)
    assert (P >= -1e-9).all() and (P.sum(axis=1) <= 1 + 1e-6).all()
