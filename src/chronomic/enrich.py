"""Matched-null (SNPsnap-style) locus-level enrichment testing."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

GENE_BINS = [0, 1, 2, 3, 4, 6, 9, 14, 10**9]


def _qbins(x_ref: np.ndarray, x: np.ndarray, q: int) -> np.ndarray:
    edges = np.unique(np.quantile(x_ref, np.linspace(0, 1, q + 1)))
    edges[0], edges[-1] = -np.inf, np.inf
    return np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)


def _gbins(n: np.ndarray) -> np.ndarray:
    return np.searchsorted(np.array(GENE_BINS), n, side="right") - 1


def match_candidates(obs: pd.DataFrame, pool: pd.DataFrame, ngene_col: str, q: int = 10,
                     min_candidates: int = 50) -> tuple[list[np.ndarray], list[str]]:
    """For every observed locus, the pool indices with matching MAF/LD/gene-count/distance bins.

    Gene count is always matched; the other features are relaxed in turn
    (distance, then LD buddies, then MAF) if fewer than `min_candidates` remain.
    """
    feats = {
        "maf": (_qbins(pool.maf.values, pool.maf.values, q), _qbins(pool.maf.values, obs.maf.values, q)),
        "buddies": (_qbins(pool.n_buddies.values, pool.n_buddies.values, q),
                    _qbins(pool.n_buddies.values, obs.n_buddies.values, q)),
        "dist": (_qbins(np.log1p(pool.nearest_dist.values), np.log1p(pool.nearest_dist.values), q),
                 _qbins(np.log1p(pool.nearest_dist.values), np.log1p(obs.nearest_dist.values), q)),
        "genes": (_gbins(pool[ngene_col].values), _gbins(obs[ngene_col].values)),
    }
    always = ("genes",)
    schemes = [always + ("maf", "buddies", "dist"), always + ("maf", "buddies"), always + ("maf",), always]
    keys = {}
    for sch in schemes:
        pk = np.zeros(len(pool), dtype=np.int64)
        for f in sch:
            pk = pk * 100 + feats[f][0]
        order = np.argsort(pk, kind="stable")
        keys[sch] = (pk, order, pk[order])
    cands, level = [], []
    for i in range(len(obs)):
        for sch in schemes:
            ok = 0
            for f in sch:
                ok = ok * 100 + feats[f][1][i]
            pk, order, sk = keys[sch]
            lo, hi = np.searchsorted(sk, ok, "left"), np.searchsorted(sk, ok, "right")
            if hi - lo >= min_candidates or sch == schemes[-1]:
                cands.append(order[lo:hi])
                level.append("+".join(sch))
                break
    return cands, level


def set_hit_matrix(gene_mat: sparse.csr_matrix, set_ind: np.ndarray) -> np.ndarray:
    """(n_loci x n_sets) boolean: locus contains >= 1 gene of the set."""
    return np.asarray((gene_mat.astype(np.int32) @ set_ind.astype(np.int32)) > 0)


def permute(cands: list[np.ndarray], pool_hits: np.ndarray, n_perm: int, rng: np.random.Generator,
            chunk: int = 1000) -> np.ndarray:
    """Null distribution of the number of loci hitting each set: (n_perm x n_sets)."""
    n_sets = pool_hits.shape[1]
    out = np.zeros((n_perm, n_sets), dtype=np.int32)
    for s in range(0, n_perm, chunk):
        e = min(n_perm, s + chunk)
        draw = np.stack([c[rng.integers(0, len(c), size=e - s)] for c in cands], axis=1)  # (perm x loci)
        for j in range(n_sets):
            out[s:e, j] = pool_hits[draw, j].sum(axis=1)
    return out


def relative_fraction_test(obs_both: int, obs_ref: int, null_both: np.ndarray,
                           null_ref: np.ndarray) -> tuple[float, float, float]:
    """Relative (specificity) test. Among loci containing a gene of a reference set, the fraction that also
    contain a gene of the test set, observed vs the same fraction in matched null loci.

    obs_both / null_both count loci containing a gene of both sets; null arrays must come from the same
    null draws as null_ref. Returns (observed fraction, mean null fraction, upper-tail permutation p).
    """
    fo = obs_both / obs_ref if obs_ref else np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        fn = null_both / null_ref
    fn = fn[np.isfinite(fn)]
    if fo != fo or not fn.size:
        return fo, np.nan, np.nan
    return fo, float(fn.mean()), (1 + (fn >= fo).sum()) / (1 + fn.size)


def summarise(obs_hits: np.ndarray, null: np.ndarray, set_ids: list[str], ref_col: int | None = None,
              obs_both_ref: np.ndarray | None = None, null_both_ref: np.ndarray | None = None) -> pd.DataFrame:
    """Enrichment, empirical p-values and (optionally) enrichment relative to a reference set.

    For the relative columns, obs_both_ref[j] / null_both_ref[:, j] count loci containing a gene of set j and
    a gene of the reference set (see relative_fraction_test).
    """
    if ref_col is not None and (obs_both_ref is None or null_both_ref is None):
        raise ValueError("relative enrichment needs joint hit counts with the reference set")
    n_perm = null.shape[0]
    mean = null.mean(axis=0)
    rows = []
    for j, sid in enumerate(set_ids):
        o = obs_hits[j]
        ge = (null[:, j] >= o).sum()
        le = (null[:, j] <= o).sum()
        p_up = (1 + ge) / (1 + n_perm)
        p_lo = (1 + le) / (1 + n_perm)
        r = {"set_id": sid, "obs_loci": int(o), "exp_loci": mean[j],
             "null_q025": np.quantile(null[:, j], 0.025), "null_q975": np.quantile(null[:, j], 0.975),
             "fold": o / mean[j] if mean[j] > 0 else np.nan,
             "p_enrich": p_up, "p_two_sided": min(1.0, 2 * min(p_up, p_lo))}
        if ref_col is not None and j != ref_col:
            fo, fn, p = relative_fraction_test(obs_both_ref[j], obs_hits[ref_col], null_both_ref[:, j], null[:, ref_col])
            r["frac_of_mendelian_obs"] = fo
            r["frac_of_mendelian_exp"] = fn
            r["p_relative"] = p
        rows.append(r)
    return pd.DataFrame(rows)


def hit_probabilities(cands: list[np.ndarray], pool_hits: np.ndarray) -> np.ndarray:
    """(n_loci x n_sets) probability that a matched null locus hits each set."""
    return np.stack([pool_hits[c].mean(axis=0) for c in cands])


def poisson_binomial_sf(probs: np.ndarray, k: int) -> float:
    """P(X >= k) for X = sum of independent Bernoulli(probs) (exact dynamic programme)."""
    dist = np.zeros(len(probs) + 1)
    dist[0] = 1.0
    for q in probs:
        dist[1:] = dist[1:] * (1 - q) + dist[:-1] * q
        dist[0] *= (1 - q)
    return float(dist[k:].sum())
