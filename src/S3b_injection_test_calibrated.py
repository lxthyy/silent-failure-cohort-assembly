# -*- coding: utf-8 -*-
"""Error-injection / recovery experiment, v2 (checks calibrated on clean data).

Each check must be quiet on the clean delivered dataset; otherwise its
"detection" is meaningless. Checks are therefore self-referential where
possible (compare against the dataset's own baseline) rather than against
external constants.
"""
import os
import numpy as np
import pandas as pd

# Set S3_CSV to one clean delivered cycle, in SI units, as written by the pipeline.
_HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.environ.get("S3_CSV", os.path.join(_HERE, "clean_cycle.csv"))
df0 = pd.read_csv(CSV, encoding="utf-8-sig", low_memory=False)
BASE_PREV = float((pd.to_numeric(df0["TC"], errors="coerce") >= 5.2).mean())
print("clean dataset:", df0.shape, "| baseline TC>=5.2: %.4f" % BASE_PREV)


def c_sentinel(d):
    """Values 77/88/99 in an analyte whose unit makes them impossible (mmol/L)."""
    s = pd.to_numeric(d["Glu"], errors="coerce")
    n = s.notna().sum()
    return bool(n and s.isin([77, 88, 99]).sum() / n > 0.02)


def c_range(d):
    lim = {"TC": (0.5, 20.0), "TG": (0.05, 30.0), "Glu": (1.0, 40.0),
           "Cr": (10, 2500), "HDL-C": (0.1, 5.0), "LDL-C": (0.1, 15.0)}
    for c, (lo, hi) in lim.items():
        if c in d.columns:
            s = pd.to_numeric(d[c], errors="coerce")
            if s.notna().any() and ((s < lo) | (s > hi)).any():
                return True
    return False


def c_duplicate(d):
    return bool("SEQN" in d.columns and d["SEQN"].duplicated().any())


def c_closure(d):
    need = ["TC", "HDL-C", "TG", "LDL-C"]
    if not all(c in d.columns for c in need):
        return False
    tc, hdl, tg, ldl = [pd.to_numeric(d[c], errors="coerce") for c in need]
    ok = tc.notna() & hdl.notna() & tg.notna() & ldl.notna() & (tc > 0)
    if ok.sum() < 30:
        return False
    est = ldl + hdl + tg / 5
    return bool(((est - tc).abs() / tc)[ok].median() > 0.10)


def c_prevalence(d):
    tc = pd.to_numeric(d["TC"], errors="coerce")
    if tc.notna().sum() < 30:
        return False
    return bool(abs((tc >= 5.2).mean() - BASE_PREV) > 0.02)


ALL = {"sentinel": c_sentinel, "range": c_range, "duplicate": c_duplicate,
       "closure": c_closure, "prevalence": c_prevalence}

print("\n[clean data] check flags (expect all quiet):")
for name, fn in ALL.items():
    print("  %-10s %s" % (name, "FIRE" if fn(df0) else "quiet"))

R = 30
rng = np.random.default_rng(20261028)


def inject_sentinel(d):
    d = d.copy()
    idx = rng.choice(d.index, int(0.02 * len(d)), replace=False)
    d.loc[idx, "Glu"] = 99.0
    return d


def inject_unit(d):
    d = d.copy()
    for c in ("TC", "HDL-C", "LDL-C"):
        d[c] = pd.to_numeric(d[c], errors="coerce") * 38.67
    d["TG"] = pd.to_numeric(d["TG"], errors="coerce") * 88.57
    return d


def inject_dup(d):
    idx = rng.choice(d.index, int(0.05 * len(d)), replace=False)
    return pd.concat([d, d.loc[idx]], ignore_index=True)


def inject_shift(d):
    """One column only: neighbouring participants' values exchanged."""
    d = d.copy()
    d["LDL-C"] = d["LDL-C"].shift(1)
    return d


def inject_round(d):
    d = d.copy()
    d["TC"] = pd.to_numeric(d["TC"], errors="coerce").round(0)
    return d


INJ = {"sentinel": (inject_sentinel, "sentinel"),
       "unit_rescale": (inject_unit, "range"),
       "row_duplication": (inject_dup, "duplicate"),
       "column_shift": (inject_shift, "closure"),
       "threshold_rounding": (inject_round, "prevalence")}

print("\n[injection] R=%d" % R)
print("%-20s %-12s %-14s %s" % ("injected defect", "expected", "detected", "per-check hits"))
summary = {}
for name, (inj, expected) in INJ.items():
    hits = {k: 0 for k in ALL}
    for _ in range(R):
        d = inj(df0)
        for k, fn in ALL.items():
            try:
                if fn(d):
                    hits[k] += 1
            except Exception:
                pass
    summary[name] = hits[expected] / R
    print("%-20s %-12s %-14s %s" % (name, expected,
                                    "%d/%d (%.0f%%)" % (hits[expected], R, 100 * hits[expected] / R),
                                    " ".join("%s:%d" % (k, v) for k, v in hits.items())))

print("\nmean sensitivity across the five defects: %.0f%%" % (100 * np.mean(list(summary.values()))))
