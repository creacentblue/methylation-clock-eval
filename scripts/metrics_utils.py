"""Shared metric helpers."""
import numpy as np
from scipy import stats


def core(age, p):
    """Standard error metrics of predicted vs chronological age (inputs: 1-D numpy arrays)."""
    r = p - age
    ss_res, ss_tot = np.sum(r ** 2), np.sum((age - age.mean()) ** 2)
    return dict(n=len(age), MAE=np.mean(np.abs(r)), median_AE=np.median(np.abs(r)), RMSE=np.sqrt(np.mean(r ** 2)),
                R2=1 - ss_res / ss_tot, pearson_r=stats.pearsonr(age, p)[0], spearman_rho=stats.spearmanr(age, p)[0],
                bias_mean_signed_error=np.mean(r),
                MAE_after_removing_bias=np.mean(np.abs(r - r.mean())),   # diagnostic only: uses eval labels to remove the mean offset
                slope_pred_on_age=stats.linregress(age, p).slope)
