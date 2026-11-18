"""
bootstrap.py -- a small, dependency-light bootstrap toolkit.

Only needs numpy (and pandas for convenience). Normal quantiles come from the
standard library, so scipy is not required.

Main pieces
-----------
bootstrap()              generic nonparametric bootstrap of any statistic,
                         optionally resampling within strata (e.g. exposure groups)
jackknife()              leave-one-out values (used for the BCa acceleration)
BootResult.ci()          normal, basic, percentile and BCa intervals
ols()                    ordinary least squares with model-based SEs
case_bootstrap_ols()     resample rows (x, y) together
residual_bootstrap_ols() keep X fixed, resample residuals
monte_carlo_error()      repeat a bootstrap for several B to see how much it moves
"""

from __future__ import annotations

from dataclasses import dataclass, field
from statistics import NormalDist
from typing import Callable, Sequence

import numpy as np
import pandas as pd

_N = NormalDist()


def _z(p):
    """Standard normal quantile, vectorised over p."""
    return np.vectorize(_N.inv_cdf)(np.clip(p, 1e-12, 1 - 1e-12))


def _phi(x):
    """Standard normal CDF, vectorised over x."""
    return np.vectorize(_N.cdf)(x)


def _rng(seed):
    return seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)


# --------------------------------------------------------------------------- #
# Result container
# --------------------------------------------------------------------------- #
@dataclass
class BootResult:
    """Holds the observed statistic t0 and the B bootstrap replicates."""

    t0: np.ndarray              # shape (k,)
    reps: np.ndarray            # shape (B, k)
    names: list[str] = field(default_factory=list)
    jack: np.ndarray | None = None   # shape (n, k), only needed for BCa

    @property
    def B(self) -> int:
        return self.reps.shape[0]

    @property
    def se(self) -> np.ndarray:
        return self.reps.std(axis=0, ddof=1)

    @property
    def bias(self) -> np.ndarray:
        return self.reps.mean(axis=0) - self.t0

    # ---- confidence intervals ------------------------------------------- #
    def ci(self, method: str = "percentile", level: float = 0.95) -> np.ndarray:
        """Return an array of shape (k, 2) with lower/upper limits."""
        a = (1 - level) / 2
        m = method.lower()
        if m == "normal":
            z = _z(1 - a)
            centre = self.t0 - self.bias            # bias-corrected normal
            return np.column_stack([centre - z * self.se, centre + z * self.se])
        if m == "percentile":
            return np.quantile(self.reps, [a, 1 - a], axis=0).T
        if m == "basic":
            lo, hi = np.quantile(self.reps, [a, 1 - a], axis=0)
            return np.column_stack([2 * self.t0 - hi, 2 * self.t0 - lo])
        if m == "bca":
            return self._bca(a)
        raise ValueError(f"unknown method {method!r}")

    def _bca(self, a: float) -> np.ndarray:
        if self.jack is None:
            raise ValueError("BCa needs jackknife values: call bootstrap(..., bca=True)")
        out = np.empty((self.t0.size, 2))
        for j in range(self.t0.size):
            r = self.reps[:, j]
            # bias correction z0: how often replicates fall below t0
            # (ties counted as half, which matters for medians of rounded data)
            prop = (np.sum(r < self.t0[j]) + 0.5 * np.sum(r == self.t0[j])) / r.size
            z0 = _z(prop)
            # acceleration from the jackknife
            th = self.jack[:, j]
            d = th.mean() - th
            denom = 6 * (np.sum(d ** 2) ** 1.5)
            acc = np.sum(d ** 3) / denom if denom > 0 else 0.0
            zs = _z(np.array([a, 1 - a]))
            adj = _phi(z0 + (z0 + zs) / (1 - acc * (z0 + zs)))
            out[j] = np.quantile(r, adj)
        return out

    def summary(self, level: float = 0.95,
                methods: Sequence[str] = ("normal", "basic", "percentile", "bca")) -> pd.DataFrame:
        rows = []
        for m in methods:
            if m == "bca" and self.jack is None:
                continue
            lims = self.ci(m, level)
            for j, nm in enumerate(self.names or [f"t{j}" for j in range(self.t0.size)]):
                rows.append({"term": nm, "method": m, "estimate": self.t0[j],
                             "boot_se": self.se[j], "lower": lims[j, 0], "upper": lims[j, 1]})
        return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Generic bootstrap
# --------------------------------------------------------------------------- #
def _as_vec(x) -> np.ndarray:
    return np.atleast_1d(np.asarray(x, dtype=float))


def resample_index(n: int, rng, strata: np.ndarray | None = None) -> np.ndarray:
    """Indices for one bootstrap sample; if strata given, resample within each."""
    if strata is None:
        return rng.integers(0, n, n)
    idx = []
    for s in np.unique(strata):
        members = np.flatnonzero(strata == s)
        idx.append(rng.choice(members, members.size, replace=True))
    return np.concatenate(idx)


def jackknife(data: pd.DataFrame, statistic: Callable[[pd.DataFrame], object]) -> np.ndarray:
    n = len(data)
    return np.vstack([_as_vec(statistic(data.drop(data.index[i]))) for i in range(n)])


def bootstrap(data: pd.DataFrame,
              statistic: Callable[[pd.DataFrame], object],
              B: int = 2000,
              strata: str | None = None,
              bca: bool = False,
              seed=None,
              names: list[str] | None = None) -> BootResult:
    """
    Nonparametric bootstrap.

    statistic : function taking a DataFrame and returning a scalar or 1-D array.
    strata    : column name to resample within (keeps group sizes fixed).
    bca       : also compute jackknife values so BCa intervals are available.
    """
    rng = _rng(seed)
    data = data.reset_index(drop=True)
    n = len(data)
    s = data[strata].to_numpy() if strata else None
    t0 = _as_vec(statistic(data))
    reps = np.empty((B, t0.size))
    for b in range(B):
        reps[b] = _as_vec(statistic(data.iloc[resample_index(n, rng, s)]))
    jack = jackknife(data, statistic) if bca else None
    return BootResult(t0, reps, names or [], jack)


# --------------------------------------------------------------------------- #
# Linear regression
# --------------------------------------------------------------------------- #
@dataclass
class OLSFit:
    coef: np.ndarray
    se: np.ndarray
    fitted: np.ndarray
    resid: np.ndarray
    names: list[str]
    df_resid: int
    hat: np.ndarray | None = None

    def table(self, level: float = 0.95) -> pd.DataFrame:
        # t-quantile approximated by normal: with ~1550 residual df they agree to 3 dp
        z = _z(1 - (1 - level) / 2)
        return pd.DataFrame({"term": self.names, "estimate": self.coef, "se": self.se,
                             "lower": self.coef - z * self.se, "upper": self.coef + z * self.se})


def design_matrix(df: pd.DataFrame, numeric: Sequence[str], categorical: Sequence[str] = ()) -> pd.DataFrame:
    """Intercept + numeric columns + dummy-coded categoricals (first level dropped)."""
    X = df[list(numeric)].astype(float)
    if categorical:
        dummies = pd.get_dummies(df[list(categorical)].astype("category"), drop_first=True, dtype=float)
        X = pd.concat([X, dummies], axis=1)
    X.insert(0, "Intercept", 1.0)
    return X


def ols(X: np.ndarray | pd.DataFrame, y, names: list[str] | None = None, leverage: bool = False) -> OLSFit:
    if isinstance(X, pd.DataFrame):
        names = names or list(X.columns)
        X = X.to_numpy(float)
    y = np.asarray(y, float)
    n, p = X.shape
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    fitted = X @ beta
    resid = y - fitted
    sigma2 = resid @ resid / (n - p)
    se = np.sqrt(np.diag(XtX_inv) * sigma2)
    hat = np.einsum("ij,jk,ik->i", X, XtX_inv, X) if leverage else None
    return OLSFit(beta, se, fitted, resid, names or [f"x{j}" for j in range(p)], n - p, hat)


def robust_se(X: np.ndarray | pd.DataFrame, fit: OLSFit) -> np.ndarray:
    """HC3 sandwich SE -- a useful model-free yardstick for the case bootstrap."""
    X = np.asarray(X, float)
    XtX_inv = np.linalg.inv(X.T @ X)
    h = np.einsum("ij,jk,ik->i", X, XtX_inv, X)
    u2 = (fit.resid / (1 - h)) ** 2
    meat = X.T @ (X * u2[:, None])
    return np.sqrt(np.diag(XtX_inv @ meat @ XtX_inv))


def case_bootstrap_ols(X, y, B: int = 2000, seed=None, strata=None) -> BootResult:
    """Resample (x_i, y_i) pairs. Only assumes rows are i.i.d."""
    rng = _rng(seed)
    names = list(X.columns) if isinstance(X, pd.DataFrame) else None
    X = np.asarray(X, float); y = np.asarray(y, float)
    n = len(y)
    t0 = np.linalg.lstsq(X, y, rcond=None)[0]
    reps = np.empty((B, X.shape[1]))
    for b in range(B):
        i = resample_index(n, rng, strata)
        reps[b] = np.linalg.lstsq(X[i], y[i], rcond=None)[0]
    return BootResult(t0, reps, names or [])


def residual_bootstrap_ols(X, y, B: int = 2000, seed=None, rescale: bool = True) -> BootResult:
    """
    Keep X fixed; build y* = X b + e*, with e* drawn from the residuals.
    Assumes the mean model is right and errors are i.i.d. (constant variance).
    rescale=True inflates residuals by sqrt(n/(n-p)) to undo their shrinkage.
    """
    rng = _rng(seed)
    names = list(X.columns) if isinstance(X, pd.DataFrame) else None
    X = np.asarray(X, float); y = np.asarray(y, float)
    n, p = X.shape
    proj = np.linalg.pinv(X)                 # (X'X)^-1 X', reused every replicate
    t0 = proj @ y
    fitted = X @ t0
    e = y - fitted
    if rescale:
        e = e * np.sqrt(n / (n - p))
    e = e - e.mean()
    reps = np.empty((B, p))
    for b in range(B):
        reps[b] = proj @ (fitted + rng.choice(e, n, replace=True))
    return BootResult(t0, reps, names or [])


# --------------------------------------------------------------------------- #
# Monte Carlo error
# --------------------------------------------------------------------------- #
def monte_carlo_error(run: Callable[[int, int], BootResult],
                      Bs: Sequence[int] = (50, 100, 200, 500, 1000, 2000, 5000),
                      n_repeats: int = 20,
                      term: int = 0,
                      level: float = 0.95,
                      methods: Sequence[str] = ("percentile",),
                      seed: int = 2026) -> pd.DataFrame:
    """
    Re-run a bootstrap n_repeats times at each B (with different seeds) and record
    the SE and CI limits. run(B, seed) must return a BootResult.
    """
    ss = np.random.SeedSequence(seed)
    rows = []
    for B in Bs:
        for r, child in enumerate(ss.spawn(n_repeats)):
            res = run(B, np.random.default_rng(child))
            row = {"B": B, "repeat": r, "se": res.se[term]}
            for m in methods:
                lo, hi = res.ci(m, level)[term]
                row[f"{m}_lower"], row[f"{m}_upper"] = lo, hi
            rows.append(row)
    return pd.DataFrame(rows)
