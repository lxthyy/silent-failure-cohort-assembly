# -*- coding: utf-8 -*-
"""
Analysis A: lossy unit conversion (rounding) in the NHANES downloader
Analysis B: fasting-subsample weight (WTSAF2YR) vs examination weight (WTMEC2YR)

Reads CDC XPT directly from the local cache; does NOT call the downloader.
Outputs a plain-text report to the console and a CSV block for the paper.
"""
import os
import numpy as np
import pandas as pd

# Set NHANES_CACHE to the folder holding the CDC XPT files (DEMO_E.xpt, ...).
_HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.environ.get("NHANES_CACHE", os.path.join(_HERE, "nhanes_cache"))
CYCLES = {"E": "2007-2008", "F": "2009-2010", "G": "2011-2012"}


def rx(name, cyc):
    return pd.read_sas(os.path.join(CACHE, f"{name}_{cyc}.xpt"), format="xport")


ANALYTES = {"LBXSCR", "LBXSUA", "LBXSBU", "LBXTC", "LBDHDD", "LBXTR", "LBDLDL",
            "LBXGLU", "LBDGLUSI"}


def clean(s, name=None):
    """Mirror the downloader's missing-value rules.

    The v2.05+ tool exempts the clinical analytes from the 77/88/99 codes
    (NORMAL_VALUE_RANGES), so those values must survive here as well;
    otherwise this script would reproduce the very defect under study.
    """
    s = pd.to_numeric(s, errors="coerce")
    s = s.mask((s.abs() > 1e-100) & (s.abs() < 1e-15))
    big = [77777, 99999, 88888, 66666, 7777, 9999, 8888, 6666]
    small = [777, 999, 888, 77, 99, 88]
    s = s.replace(big, np.nan)
    if name not in ANALYTES:
        s = s.replace(small, np.nan)
    return s


# ----------------------------------------------------------------------------
# load
# ----------------------------------------------------------------------------
frames = []
for cyc in CYCLES:
    demo = rx("DEMO", cyc)[["SEQN", "SDMVSTRA", "SDMVPSU", "WTMEC2YR",
                            "RIDAGEYR", "RIAGENDR"]].copy()
    bio = rx("BIOPRO", cyc)[["SEQN", "LBXSCR", "LBXSUA", "LBXSBU"]]
    tch = rx("TCHOL", cyc)[["SEQN", "LBXTC"]]
    hdl = rx("HDL", cyc)[["SEQN", "LBDHDD"]]
    tri = rx("TRIGLY", cyc)[["SEQN", "LBXTR", "LBDLDL"]]
    glu = rx("GLU", cyc)[["SEQN", "LBXGLU", "LBDGLUSI", "WTSAF2YR", "PHAFSTHR"]]
    bmx = rx("BMX", cyc)[["SEQN", "BMXBMI"]]
    d = demo
    for t in (bio, tch, hdl, tri, glu, bmx):
        d = d.merge(t, on="SEQN", how="left")
    d["cycle"] = CYCLES[cyc]
    frames.append(d)
df = pd.concat(frames, ignore_index=True)
for c in ["LBXSCR", "LBXSUA", "LBXSBU", "LBXTC", "LBDHDD", "LBXTR", "LBDLDL",
          "LBXGLU", "LBDGLUSI", "BMXBMI", "RIDAGEYR"]:
    df[c] = clean(df[c], c)
df["WTMEC2YR"] = clean(df["WTMEC2YR"])
df["WTSAF2YR"] = clean(df["WTSAF2YR"])
df["PHAFSTHR"] = clean(df["PHAFSTHR"])
df["PHAFSTMN"] = np.nan

print("=" * 100)
print("ANALYSIS A — lossy unit conversion (rounding) as implemented in UNIT_CONVERSIONS")
print("=" * 100)

# (variable, label, factor, decimals, partner, partner_label, threshold, op)
SPECS = [
    ("LBXSCR", "Creatinine", 88.4, 0, "LBXSUA", "Uric acid", 88.4, ">"),
    ("LBXSUA", "Uric acid", 59.48, 0, "LBXSCR", "Creatinine", 420.0, ">"),
    ("LBXSBU", "BUN", 0.357, 1, "LBXSCR", "Creatinine", np.nan, ">"),
    ("LBXTC", "Total cholesterol", 0.02586, 2, "LBDLDL", "LDL-C", 5.2, ">="),
    ("LBDHDD", "HDL-C", 0.02586, 2, "LBXTR", "Triglyceride", 1.0, "<"),
    ("LBXTR", "Triglyceride", 0.01129, 2, "LBDHDD", "HDL-C", 1.7, ">="),
    ("LBDLDL", "LDL-C", 0.02586, 2, "LBXTC", "Total cholesterol", 3.4, ">="),
    ("LBXGLU", "Fasting glucose", 0.0555, 1, "BMXBMI", "BMI", 7.0, ">="),
]

rows = []
for var, lab, fac, nd, pvar, plab, thr, op in SPECS:
    x = df[var].dropna()
    exact = x * fac
    rnd = exact.round(nd)
    dev = (rnd - exact)
    n = len(x)
    nuniq_exact = exact.nunique()
    nuniq_rnd = rnd.nunique()
    # correlation with partner
    sub = df[[var, pvar]].dropna()
    pe = np.corrcoef(sub[var] * fac, sub[pvar])[0, 1]
    pr = np.corrcoef((sub[var] * fac).round(nd), sub[pvar])[0, 1]
    spe = sub[var].rank().corr(sub[pvar].rank())
    spr = (sub[var] * fac).round(nd).rank().corr(sub[pvar].rank())
    # threshold flips
    if np.isnan(thr):
        flips = np.nan
    else:
        if op == ">":
            f_exact = (exact > thr)
            f_rnd = (rnd > thr)
        elif op == ">=":
            f_exact = (exact >= thr)
            f_rnd = (rnd >= thr)
        else:
            f_exact = (exact < thr)
            f_rnd = (rnd < thr)
        flips = int((f_exact != f_rnd).sum())
    rows.append(dict(var=var, label=lab, factor=fac, decimals=nd, n=n,
                     mean_exact=exact.mean(), mean_round=rnd.mean(),
                     bias_pct=(rnd.mean() - exact.mean()) / exact.mean() * 100,
                     sd_exact=exact.std(), sd_round=rnd.std(),
                     max_abs_dev=np.abs(dev).max(),
                     rmse=np.sqrt((dev ** 2).mean()),
                     levels_lost_pct=(1 - nuniq_rnd / nuniq_exact) * 100,
                     pearson_exact=pe, pearson_round=pr,
                     spearman_exact=spe, spearman_round=spr,
                     partner=plab, threshold=thr, n_flip=flips,
                     flip_pct=flips / n * 100 if n else np.nan))
A = pd.DataFrame(rows)
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 50)
show = A[["label", "factor", "decimals", "n", "mean_exact", "mean_round", "bias_pct",
          "max_abs_dev", "rmse", "levels_lost_pct", "pearson_exact", "pearson_round",
          "spearman_exact", "spearman_round", "partner", "threshold", "n_flip"]].copy()
for c in ["mean_exact", "mean_round", "bias_pct", "max_abs_dev", "rmse",
          "levels_lost_pct", "pearson_exact", "pearson_round", "spearman_exact",
          "spearman_round"]:
    show[c] = show[c].round(6)
print(show.to_string(index=False))

print()
print("=" * 100)
print("ANALYSIS B — fasting subsample weight vs examination weight (glucose, cycles E/F/G)")
print("=" * 100)

def wmean_var(x, w, strata, psu):
    m = x.notna() & w.notna() & (w > 0)
    x, w, strata, psu = x[m].values, w[m].values, strata[m].values, psu[m].values
    W = w.sum()
    est = (w * x).sum() / W
    u = w * (x - est) / W
    tmp = pd.DataFrame({"h": strata, "j": psu, "u": u})
    U = tmp.groupby(["h", "j"])["u"].sum().reset_index()
    var = 0.0
    for h, g in U.groupby("h"):
        nh = len(g)
        if nh < 2:
            continue
        var += nh / (nh - 1) * ((g["u"] - g["u"].mean()) ** 2).sum()
    return est, np.sqrt(var), int(m.sum()), W


def report(tag, x, w, sub=None):
    d = df if sub is None else df[sub]
    est, se, n, W = wmean_var(d[x], d[w], d["SDMVSTRA"], d["SDMVPSU"])
    print(f"  {tag:<52s} mean={est:8.3f}  SE={se:6.3f}  n={n:6d}  sumW={W:12.0f}")
    return est, se, n

n_all = len(df)
n_glu = df["LBXGLU"].notna().sum()
n_saf = ((df["WTSAF2YR"].notna()) & (df["WTSAF2YR"] > 0)).sum()
print(f"  pooled participants (E/F/G)            : {n_all}")
print(f"  with LBXGLU (fasting subsample)        : {n_glu} ({n_glu/n_all*100:.1f}%)")
print(f"  with WTSAF2YR > 0                      : {n_saf} ({n_saf/n_all*100:.1f}%)")
print(f"  PHAFSTHR >= 8 h among LBXGLU non-missing: "
      f"{int(((df['PHAFSTHR']>=8) & df['LBXGLU'].notna()).sum())}")

print("\n  -- all ages --")
g1 = report("LBXGLU, WTMEC2YR on fasting subsample (as in _p3_weighted.R)", "LBXGLU", "WTMEC2YR")
g2 = report("LBXGLU, WTSAF2YR on fasting subsample (design-correct)", "LBXGLU", "WTSAF2YR")
print(f"  difference: {g1[0]-g2[0]:+.3f} mg/dL ({(g1[0]-g2[0])/g2[0]*100:+.2f}%)")

print("\n  -- adults >= 20 y --")
ad = df["RIDAGEYR"] >= 20
g3 = report("LBXGLU, WTMEC2YR, age>=20", "LBXGLU", "WTMEC2YR", ad)
g4 = report("LBXGLU, WTSAF2YR, age>=20", "LBXGLU", "WTSAF2YR", ad)
print(f"  difference: {g3[0]-g4[0]:+.3f} mg/dL ({(g3[0]-g4[0])/g4[0]*100:+.2f}%)")

print("\n  -- 8-hour fasted only --")
f8 = (df["PHAFSTHR"] >= 8) & (df["RIDAGEYR"] >= 20)
g5 = report("LBXGLU, WTSAF2YR, age>=20 & fasted>=8h", "LBXGLU", "WTSAF2YR", f8)
g6 = report("LBXGLU, WTMEC2YR, age>=20 & fasted>=8h", "LBXGLU", "WTMEC2YR", f8)

print("\n  -- per cycle (adults >= 20) --")
for cyc, name in CYCLES.items():
    m = ad & (df["cycle"] == name)
    e1, s1, n1, _ = wmean_var(df.loc[m, "LBXGLU"], df.loc[m, "WTMEC2YR"],
                              df.loc[m, "SDMVSTRA"], df.loc[m, "SDMVPSU"])
    e2, s2, n2, _ = wmean_var(df.loc[m, "LBXGLU"], df.loc[m, "WTSAF2YR"],
                              df.loc[m, "SDMVSTRA"], df.loc[m, "SDMVPSU"])
    print(f"  {name}: WTMEC {e1:7.2f} (SE {s1:.2f})  WTSAF {e2:7.2f} (SE {s2:.2f})  "
          f"diff {e1-e2:+6.2f} mg/dL   n={n1}")

print("\n  -- SI value already supplied by CDC (LBDGLUSI) vs tool re-derivation --")
sub = df[["LBXGLU", "LBDGLUSI"]].dropna()
tool = (sub["LBXGLU"] * 0.0555).round(1)
print(f"  n={len(sub)}  CDC LBDGLUSI vs tool round(LBXGLU*0.0555,1): "
      f"identical {int((tool == sub['LBDGLUSI']).sum())}/{len(sub)} "
      f"({100*(tool == sub['LBDGLUSI']).mean():.2f}%), "
      f"max |diff| {np.abs(tool - sub['LBDGLUSI']).max():.3f} mmol/L")

A.to_csv(os.path.join(_HERE, "A_rounding.csv"), index=False, encoding="utf-8-sig")
print("\n[saved] A_rounding.csv")
