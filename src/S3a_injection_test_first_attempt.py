# -*- coding: utf-8 -*-
"""Error-injection / recovery experiment.

Take the clean delivered dataset (one survey cycle), inject known defects with
known ground truth, and ask whether the corresponding check fires — and whether
any check fires on the clean data (false positives).

Outputs a compact table to stdout.
"""
import os
import numpy as np
import pandas as pd

# Set S3_CSV to one clean delivered cycle, in SI units, as written by the pipeline.
_HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.environ.get("S3_CSV", os.path.join(_HERE, "clean_cycle.csv"))
df0 = pd.read_csv(CSV, encoding="utf-8-sig", low_memory=False)
print("clean dataset:", df0.shape)

CHECKS = {}


def c_sentinel(d):
    """share of values in {77,88,99} among an integer-labelled analyte."""
    bad = 0
    for c in ("Glu", "Cr", "TC"):
        s = pd.to_numeric(d[c], errors="coerce")
        n = s.notna().sum()
        if n:
            share = s.isin([77, 88, 99]).sum() / n
            if share > 0.02:
                bad += 1
    return bad > 0


def c_range(d):
    """physiologically impossible values."""
    lim = {"TC": (1.0, 15.0), "TG": (0.1, 20.0), "Glu": (2.0, 35.0),
           "Cr": (10, 2000), "HDL-C": (0.2, 4.5)}
    for c, (lo, hi) in lim.items():
        if c in d.columns:
            s = pd.to_numeric(d[c], errors="coerce")
            if s.notna().any() and ((s < lo) | (s > hi)).any():
                return True
    return False


def c_duplicate(d):
    if "SEQN" in d.columns and d["SEQN"].duplicated().any():
        return True
    return False


def c_closure(d):
    """LDL + HDL + TG/5 = TC, median relative deviation."""
    need = ["TC", "HDL-C", "TG", "LDL-C"]
    if not all(c in d.columns for c in need):
        return False
    tc, hdl, tg, ldl = [pd.to_numeric(d[c], errors="coerce") for c in need]
    ok = tc.notna() & hdl.notna() & tg.notna() & ldl.notna() & (tc > 0)
    if ok.sum() < 30:
        return False
    est = ldl + hdl + tg / 5
    dev = ((est - tc).abs() / tc)[ok].median()
    return dev > 0.10


def c_prevalence(d, base=0.3556):
    if "TC" not in d.columns:
        return False
    tc = pd.to_numeric(d["TC"], errors="coerce")
    if tc.notna().sum() < 30:
        return False
    prev = (tc >= 5.2).sum() / tc.notna().sum()
    return abs(prev - base) > 0.02


ALL = {"sentinel": c_sentinel, "range": c_range, "duplicate": c_duplicate,
       "closure": c_closure, "prevalence": c_prevalence}

# ---------------- false positives on clean data ----------------
print("\n[clean data] check flags:")
fp = {}
for name, fn in ALL.items():
    fp[name] = bool(fn(df0))
    print("  %-10s %s" % (name, "FIRE" if fp[name] else "quiet"))

# ---------------- injections ----------------
R = 30
rng = np.random.default_rng(20261028)


def inject_sentinel(d, r):
    d = d.copy()
    idx = d.index[rng.choice(len(d), int(0.02 * len(d)), replace=False)]
    d.loc[idx, "Glu"] = 99.0
    return d


def inject_unit(d, r):
    d = d.copy()
    d["TC"] = pd.to_numeric(d["TC"], errors="coerce") * 38.67
    d["HDL-C"] = pd.to_numeric(d["HDL-C"], errors="coerce") * 38.67
    d["LDL-C"] = pd.to_numeric(d["LDL-C"], errors="coerce") * 38.67
    return d


def inject_dup(d, r):
    d = d.copy()
    idx = rng.choice(d.index, int(0.05 * len(d)), replace=False)
    return pd.concat([d, d.loc[idx]], ignore_index=True)


def inject_shift(d, r):
    d = d.copy()
    cols = [c for c in ("TC", "HDL-C", "LDL-C", "TG") if c in d.columns]
    block = d[cols].copy()
    n = len(block)
    k = int(0.03 * n)
    start = rng.integers(0, n - k)
    block.loc[block.index[start:start + k], cols] = block.loc[
        block.index[start:start + k], cols].shift(1)
    d[cols] = block
    return d


def inject_round(d, r):
    d = d.copy()
    d["TC"] = pd.to_numeric(d["TC"], errors="coerce").round(0)
    return d


INJ = {
    "sentinel": (inject_sentinel, "sentinel"),
    "unit": (inject_unit, "range"),
    "duplicate": (inject_dup, "duplicate"),
    "shift": (inject_shift, "closure"),
    "rounding": (inject_round, "prevalence"),
}

print("\n[injection] R=%d replicates per defect" % R)
rows = []
for name, (inj, expected) in INJ.items():
    hits = {k: 0 for k in ALL}
    for i in range(R):
        d = inj(df0, i)
        for k, fn in ALL.items():
            try:
                if fn(d):
                    hits[k] += 1
            except Exception:
                pass
    target = hits[expected]
    any_fire = sum(1 for k, v in hits.items() if v > 0)
    rows.append((name, expected, target, target / R * 100,
                 " ".join("%s:%d" % (k, v) for k, v in hits.items())))

print("\n%-10s %-10s %-14s %s" % ("defect", "expected", "detected/R", "per-check hits"))
for r in rows:
    print("%-10s %-10s %-14s %s" % (r[0], r[1], "%d/%d (%.0f%%)" % (r[2], R, r[3]), r[4]))

clean_prev = (pd.to_numeric(df0["TC"], errors="coerce") >= 5.2).mean()
print("\nclean TC>=5.2 prevalence: %.4f" % clean_prev)
for c in ("TC", "HDL-C", "LDL-C", "TG", "Glu", "Cr"):
    s = pd.to_numeric(df0[c], errors="coerce")
    print("  %-6s median=%.3f  distinct=%d" % (c, s.median(), s.nunique()))
