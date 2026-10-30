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


def test_relative_fraction_is_conditional_not_ratio():
    """Phenotype set not nested in the reference set: 5 loci hit the phenotype set, 10 the reference set and
    2 hit both. The relative statistic is the fraction of reference loci that also hit the phenotype set
    (2/10), not the ratio of locus counts (5/10)."""
    from chronomic.enrich import relative_fraction_test
    fo, fn, p = relative_fraction_test(2, 10, np.array([1, 2, 3, 4]), np.array([10, 10, 10, 10]))
    assert np.isclose(fo, 0.2)
    assert np.isclose(fn, 0.25)
    assert np.isclose(p, (1 + 3) / (1 + 4))   # null fractions >= 0.2: 0.2, 0.3, 0.4
    fo, fn, p = relative_fraction_test(0, 0, np.array([1, 2]), np.array([3, 4]))   # empty reference
    assert fo != fo and fn != fn and p != p
