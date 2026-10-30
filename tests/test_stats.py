"""Unit tests for the statistical helpers (run with: python -m pytest tests/)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from chronomic.common import bh_fdr  # noqa: E402
from chronomic.coloc import coloc_abf, labf_quant  # noqa: E402
from chronomic.enrich import poisson_binomial_sf  # noqa: E402


def test_bh_matches_reference():
    p = np.array([0.01, 0.04, 0.03, 0.005])
    # R: p.adjust(c(0.01,0.04,0.03,0.005), "BH") = 0.02 0.04 0.04 0.02
    assert np.allclose(bh_fdr(p), [0.02, 0.04, 0.04, 0.02])


def test_poisson_binomial_equals_binomial():
    from scipy.stats import binom
    probs = np.full(50, 0.1)
    for k in (0, 3, 8, 20):
        assert abs(poisson_binomial_sf(probs, k) - binom.sf(k - 1, 50, 0.1)) < 1e-12


def test_coloc_shared_vs_distinct_signal():
    rng = np.random.default_rng(1)
    n = 200
    se = np.full(n, 0.01)
    b1 = rng.normal(0, 0.01, n)
    b2 = rng.normal(0, 0.01, n)
    b1[100] = b2[100] = 0.1                      # shared causal variant
    shared = coloc_abf(labf_quant(b1, se, 1.0), labf_quant(b2, se, 1.0))
    assert shared["PP.H4"] > 0.9
    b2[100] = 0.0
    b2[20] = 0.1                                 # distinct causal variants
    distinct = coloc_abf(labf_quant(b1, se, 1.0), labf_quant(b2, se, 1.0))
    assert distinct["PP.H3"] > 0.9
