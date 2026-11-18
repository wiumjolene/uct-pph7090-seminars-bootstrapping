"""
analysis.py -- Option B (Bootstrapping), Seminars in Epi & Biostats.

Runs the full analysis from start to finish on nhefs_2026.csv:

    python analysis.py ../data/nhefs_2026.csv

All resampling uses fixed seeds (SEED below), so results are reproducible.
The report (report/option_b_bootstrap.qmd) calls run() and formats the output.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import python.bootstrap as bs  # noqa: E402

SEED = 2026
B_MEDIAN = 10_000          # Q1
B_REG = 5_000              # Q2
MC_BS = (50, 100, 200, 500, 1000, 2000, 5000)   # Q3
MC_REPEATS = 20

NUMERIC = ["qsmk", "age", "smokeintensity", "smokeyrs", "wt71"]
CATEGORICAL = ["sex", "race", "education", "exercise", "active"]
CI_METHODS = ["normal", "basic", "percentile", "bca"]


def load(path) -> tuple[pd.DataFrame, int]:
    raw = pd.read_csv(path)
    d = raw.dropna(subset=["wt82_71"]).reset_index(drop=True)   # complete cases
    return d, len(raw)


def median_diff(df: pd.DataFrame) -> float:
    y = df["wt82_71"].to_numpy()
    q = df["qsmk"].to_numpy()
    return np.median(y[q == 1]) - np.median(y[q == 0])


def run(path) -> dict:
    ss = np.random.SeedSequence(SEED)
    seeds = ss.spawn(7)
    d, n_raw = load(path)
    y_all, q_all = d["wt82_71"].to_numpy(), d["qsmk"].to_numpy()
    out = {"d": d, "n_raw": n_raw}

    # ---------------- Q1: difference in medians -------------------------- #
    boot_med = bs.bootstrap(d, median_diff, B=B_MEDIAN, bca=True,
                            seed=np.random.default_rng(seeds[0]))
    boot_med_s = bs.bootstrap(d, median_diff, B=B_MEDIAN, strata="qsmk", bca=True,
                              seed=np.random.default_rng(seeds[1]))
    out["q1"] = {"ordinary": boot_med, "stratified": boot_med_s}

    # ---------------- Q2: regression ------------------------------------- #
    X = bs.design_matrix(d, NUMERIC, CATEGORICAL)
    y = d["wt82_71"]
    fit = bs.ols(X, y)
    iq = list(X.columns).index("qsmk")
    out["q2"] = {
        "X": X, "fit": fit, "iq": iq,
        "hc3": bs.robust_se(X, fit),
        "case": bs.case_bootstrap_ols(X, y, B=B_REG, seed=np.random.default_rng(seeds[2])),
        "case_strat": bs.case_bootstrap_ols(X, y, B=B_REG, seed=np.random.default_rng(seeds[3]),
                                            strata=q_all),
        "resid": bs.residual_bootstrap_ols(X, y, B=B_REG, seed=np.random.default_rng(seeds[4])),
    }

    # ---------------- Q3: Monte Carlo error ------------------------------ #
    t0_med = boot_med.t0

    def run_median(B, rng):          # array version of median_diff, for speed
        n = len(y_all)
        reps = np.empty((B, 1))
        for b in range(B):
            i = rng.integers(0, n, n)
            yy, qq = y_all[i], q_all[i]
            reps[b, 0] = np.median(yy[qq == 1]) - np.median(yy[qq == 0])
        return bs.BootResult(t0_med, reps)

    def run_case(B, rng):
        return bs.case_bootstrap_ols(X, y, B=B, seed=rng)

    out["q3"] = {
        "median": bs.monte_carlo_error(run_median, MC_BS, MC_REPEATS, term=0,
                                       seed=int(seeds[5].generate_state(1)[0])),
        "qsmk": bs.monte_carlo_error(run_case, MC_BS, MC_REPEATS, term=iq,
                                     seed=int(seeds[6].generate_state(1)[0])),
    }
    return out


def mc_summary(mc: pd.DataFrame) -> pd.DataFrame:
    g = mc.groupby("B")
    return pd.DataFrame({
        "se_mean": g.se.mean(),
        "se_cv_pct": 100 * g.se.std() / g.se.mean(),
        "lim_sd": np.maximum(g.percentile_lower.std(), g.percentile_upper.std()),
    })


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "data" / "nhefs_2026.csv"
    res = run(path)
    d = res["d"]
    print(f"Complete cases: {len(d)} of {res['n_raw']}")
    print("\nQ1  difference in medians")
    for k, r in res["q1"].items():
        print(f"  [{k}]"); print(r.summary().round(3).to_string(index=False))
    q2 = res["q2"]; iq = q2["iq"]
    print("\nQ2  qsmk coefficient")
    print(f"  estimate {q2['fit'].coef[iq]:.3f}; model SE {q2['fit'].se[iq]:.3f}; HC3 SE {q2['hc3'][iq]:.3f}")
    for k in ("case", "case_strat", "resid"):
        r = q2[k]; lo, hi = r.ci("percentile")[iq]
        print(f"  {k:<10} SE {r.se[iq]:.3f}  95% CI {lo:.3f} to {hi:.3f}")
    print("\nQ3  Monte Carlo error")
    for k, mc in res["q3"].items():
        print(f"  [{k}]"); print(mc_summary(mc).round(4).to_string())
