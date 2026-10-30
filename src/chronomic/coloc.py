"""Approximate Bayes factor colocalisation (Giambartolomei et al. 2014; coloc.abf).

Single-causal-variant model with default priors p1 = p2 = 1e-4 and p12 = 5e-6
(Wallace 2020 recommendation).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import logsumexp


def sdy_est(varbeta: np.ndarray, maf: np.ndarray, n: np.ndarray | float) -> float:
    """coloc::sdY.est - estimate the trait SD from varbeta, MAF and N."""
    oneover = 1.0 / varbeta
    nvx = 2 * np.asarray(n) * maf * (1 - maf)
    cf = (oneover @ nvx) / (oneover @ oneover)  # regression through the origin nvx ~ oneover
    return float(np.sqrt(cf))


def labf_quant(beta: np.ndarray, se: np.ndarray, sdy: float, prior_sd: float = 0.15) -> np.ndarray:
    v = se ** 2
    w2 = (prior_sd * sdy) ** 2
    r = w2 / (w2 + v)
    z = beta / se
    return 0.5 * (np.log(1 - r) + r * z ** 2)


def _logdiff(a: float, b: float) -> float:
    mx = max(a, b)
    return mx + np.log(np.exp(a - mx) - np.exp(b - mx))


def coloc_abf(l1: np.ndarray, l2: np.ndarray, p1: float = 1e-4, p2: float = 1e-4, p12: float = 5e-6) -> dict:
    s1, s2, s12 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
    lh = np.array([0.0,
                   np.log(p1) + s1,
                   np.log(p2) + s2,
                   np.log(p1) + np.log(p2) + _logdiff(s1 + s2, s12),
                   np.log(p12) + s12])
    pp = np.exp(lh - logsumexp(lh))
    return {f"PP.H{i}": float(pp[i]) for i in range(5)} | {"nsnps": int(len(l1))}


def snp_pp_h4(l1: np.ndarray, l2: np.ndarray) -> np.ndarray:
    s = l1 + l2
    return np.exp(s - logsumexp(s))
