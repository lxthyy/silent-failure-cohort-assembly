# -*- coding: utf-8 -*-
"""
Supplementary Table S6: per-marker aggregate table for institutional test bed B.

For each of 100 laboratory markers, compares the mood-disorder group (major depression
or a mood disorder) with all other rows using two-sided Mann-Whitney tests, in both the
2,578-record extract and the 18,836-record deduplicated extract. Writes aggregate rows
only; no patient-level record is emitted.
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, false_discovery_control as fdc

sys.stdout.reconfigure(encoding="utf-8")
OUT = os.environ.get("S7_OUT", os.path.dirname(os.path.abspath(__file__)))
_HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get("S7_BASE", _HERE)
A = os.environ.get("S7_A", os.path.join(BASE, "心理疾病_全量提取.xlsx"))
B = os.environ.get("S7_B", os.path.join(BASE, "deduplicated_extract.csv"))
REG = os.environ.get("S7_REG", os.path.join(BASE, "signal_registry.json"))
PDW = "血小板体积分布宽度"

def scan(df, diag, tag):
    mood = diag.str.contains("心境") | diag.str.contains("抑郁")
    rows = []
    for m in MARKERS:
        if m not in df.columns:
            continue
        v = pd.to_numeric(df[m], errors="coerce")
        g1 = v[mood].dropna(); g0 = v[~mood].dropna()
        if len(g1) < 10 or len(g0) < 10:
            continue
        try:
            p = float(mannwhitneyu(g1, g0, alternative="two-sided").pvalue)
        except Exception as e:
            print("  skip %s: %s" % (m, e)); continue
        s1, s0 = float(g1.std(ddof=1)), float(g0.std(ddof=1))
        pooled = np.sqrt(((len(g1) - 1) * s1 ** 2 + (len(g0) - 1) * s0 ** 2) / (len(g1) + len(g0) - 2))
        d = (float(g1.mean()) - float(g0.mean())) / pooled if pooled > 0 else np.nan
        rows.append(dict(marker=m, version=tag, n_case=len(g1), n_other=len(g0),
                         median_case=round(float(g1.median()), 2),
                         median_other=round(float(g0.median()), 2),
                         smd=round(float(d), 3) if np.isfinite(d) else np.nan,
                         p=p))
    r = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    r["p_bh"] = fdc(r.p.values, method="bh")
    r["signal"] = r.p < 0.05
    print("  %-10s mood_n=%d  testable=%d  P<0.05=%d  BH<0.05=%d" %
          (tag, int(mood.sum()), len(r), int(r.signal.sum()), int((r.p_bh < 0.05).sum())))
    return r

MARKERS = sorted({k.split("|", 1)[1] for k in json.load(open(REG, encoding="utf-8"))})
print("markers:", len(MARKERS))

a = pd.read_excel(A); ra = scan(a, a["诊断"].astype(str), "misaligned")
b = pd.read_csv(B, low_memory=False); rb = scan(b, b["诊断"].astype(str), "deduplicated")

t = ra[["marker", "n_case", "n_other", "median_case", "median_other", "smd", "p", "p_bh", "signal"]] \
    .merge(rb[["marker", "n_case", "n_other", "median_case", "median_other", "smd", "p", "p_bh", "signal"]],
           on="marker", how="outer", suffixes=("_misaligned", "_deduplicated"))
t["direction_misaligned"] = np.sign(t.smd_misaligned)
t["direction_deduplicated"] = np.sign(t.smd_deduplicated)
t["direction_reversed"] = (t.direction_misaligned * t.direction_deduplicated < 0) & \
                          t.signal_misaligned.fillna(False) & t.signal_deduplicated.fillna(False)
t = t.sort_values("p_misaligned").reset_index(drop=True)

p = os.path.join(OUT, "S6_marker_aggregate.csv")
t.to_csv(p, index=False, encoding="utf-8-sig")
print("\n>>>", p, t.shape)
print("flipped(direction reversed & both significant):", int(t.direction_reversed.sum()))
print("misaligned signals:", int(t.signal_misaligned.sum()), "| deduplicated signals:", int(t.signal_deduplicated.sum()))
print("\nPDW:", t[t.marker == PDW][["median_case_misaligned", "median_other_misaligned", "p_misaligned",
                                  "median_case_deduplicated", "median_other_deduplicated", "p_deduplicated"]].to_string(index=False))
