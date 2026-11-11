"""Houseman 2012 constrained projection: for each sample y (beta at reference probes) solve  min_w ||y - R w||^2  s.t.  w >= 0, sum(w) <= 1
(same constraint set as minfi::projectCellType, which uses quadprog)."""
import numpy as np
from scipy.optimize import minimize


def project(Y, R):
    """Y: samples x probes, R: probes x cells (same probe order). Returns samples x cells proportions."""
    RtR = R.T @ R
    k = R.shape[1]
    out = np.zeros((Y.shape[0], k))
    cons = [{"type": "ineq", "fun": lambda w: 1.0 - w.sum(), "jac": lambda w: -np.ones(k)}]
    w0 = np.full(k, 1.0 / k)
    for i, y in enumerate(Y):
        Rty = R.T @ y
        res = minimize(lambda w: 0.5 * w @ RtR @ w - Rty @ w, w0, jac=lambda w: RtR @ w - Rty,
                       bounds=[(0, None)] * k, constraints=cons, method="SLSQP", options={"maxiter": 200, "ftol": 1e-12})
        out[i] = res.x
    return out
