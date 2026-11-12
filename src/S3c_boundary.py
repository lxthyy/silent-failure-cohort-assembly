# -*- coding: utf-8 -*-
"""
S3c: boundary experiment — substitute one same-unit, same-distribution column for another.
BPXSY2 := BPXSY1 (both systolic, mm Hg, same observed range), 30 replicates.
Tests whether any distribution-level check fires. Expected: none.
Writes S3c_boundary.csv.
"""
import os, sys
import pandas as pd, numpy as np
sys.stdout.reconfigure(encoding="utf-8")

CACHE = os.environ.get("NHANES_CACHE", r"C:\Users\lxddz\Desktop\科研稿件总库\99_外部数据\nhanes_主缓存_286XPT")
OUT = os.environ.get("S3C_OUT", r"C:\Users\lxddz\Desktop\科研稿件总库\07_方法学\30_源锚定审计_SAFI\_审计\S3c_boundary.csv")


def rx(n):
    return pd.read_sas(os.path.join(CACHE, f"{n}.xpt"), format="xport")


demo = rx("DEMO_E")[["SEQN"]]
bpx = rx("BPX_E")[["SEQN", "BPXSY1", "BPXSY2"]]
bp = rx("BIOPRO_E")[["SEQN", "LBXSCH", "LBXSCR", "LBXSGL"]]
tg = rx("TRIGLY_E")[["SEQN", "LBXTR", "LBDLDL"]]
hdl = rx("HDL_E")[["SEQN", "LBDHDD"]]

d = demo.merge(bpx, on="SEQN", how="left").merge(bp, on="SEQN", how="left") \
         .merge(tg, on="SEQN", how="left").merge(hdl, on="SEQN", how="left")
for c in ["BPXSY1", "BPXSY2", "LBXSCH", "LBXSCR", "LBXSGL", "LBXTR", "LBDLDL", "LBDHDD"]:
    d[c] = pd.to_numeric(d[c], errors="coerce")

panel = pd.DataFrame({
    "SEQN": d["SEQN"],
    "SBP1": d["BPXSY1"], "SBP2": d["BPXSY2"],
    "TC": d["LBXSCH"] * 0.02586, "TG": d["LBXTR"] * 0.01129,
    "Glu": d["LBXSGL"] * 0.0555, "Cr": d["LBXSCR"] * 88.4,
    "HDL-C": d["LBDHDD"] * 0.02586, "LDL-C": d["LBDLDL"] * 0.02586,
}).dropna(subset=["SBP2", "SBP1"])
print("panel rows", len(panel))

BASE_PREV = float((panel["TC"] >= 5.2).mean())


def checks(x):
    out = {}
    g = x["Glu"]
    out["sentinel"] = bool(g.notna().sum() and g.isin([77, 88, 99]).sum() / g.notna().sum() > 0.02)
    lim = {"TC": (0.5, 20), "TG": (0.05, 30), "Glu": (1, 40), "Cr": (10, 2500),
           "HDL-C": (0.1, 5), "LDL-C": (0.1, 15)}
    out["range"] = any(((x[c].notna()) & ((x[c] < lo) | (x[c] > hi))).any()
                       for c, (lo, hi) in lim.items() if c in x)
    out["duplicate"] = bool(x["SEQN"].duplicated().any())
    need = ["TC", "HDL-C", "TG", "LDL-C"]
    ok = x[need].notna().all(axis=1) & (x["TC"] > 0)
    est = x["LDL-C"] + x["HDL-C"] + x["TG"] / 5
    out["closure"] = bool(ok.sum() >= 30 and (((est - x["TC"]).abs() / x["TC"])[ok].median() > 0.10))
    out["prevalence"] = bool(x["TC"].notna().sum() >= 30 and abs((x["TC"] >= 5.2).mean() - BASE_PREV) > 0.02)
    return out


print("\nclean:", checks(panel))
rng = np.random.default_rng(20261112)
rows = []
for rep in range(30):
    x = panel.copy()
    x["SBP2"] = x["SBP1"]                      # same-unit, same-distribution substitution
    r = checks(x)
    # distribution overlap: cannot distinguish the two columns
    m1, m2 = x["SBP2"].median(), panel["SBP2"].median()
    for k, v in r.items():
        rows.append(dict(rep=rep, check=k, fired=bool(v),
                         median_substituted=m1, median_original=m2,
                         median_shift=abs(m1 - m2),
                         n=int(x["SBP2"].notna().sum())))
res = pd.DataFrame(rows)
res.to_csv(OUT, index=False, encoding="utf-8-sig")
print("\n[boundary] substitution SBP2 := SBP1, 30 replicates")
print(res.groupby("check")["fired"].sum().to_string())
print("fire rate: %d / %d" % (res["fired"].sum(), len(res)))
print("median shift (substituted vs original SBP2):", res["median_shift"].iloc[0])
print("[saved]", OUT)
