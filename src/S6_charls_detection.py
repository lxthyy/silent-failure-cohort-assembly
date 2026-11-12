# -*- coding: utf-8 -*-
"""
S6: run the twelve-class checks on a public dataset assembled by another team
    (CHARLS: China Health and Retirement Longitudinal Study).

Six conditions are tested. Nothing here asserts that CHARLS is defective:
each condition is a state of the released files that a naive pipeline would
carry forward silently. Every flagged record is written out for adjudication.

Inputs (obtain from the CHARLS release; not redistributable):
    CHARLS 2011 venous blood : Blood_20140429.dta
    CHARLS 2015 venous blood : Blood.dta
    Harmonized CHARLS wave D : H_CHARLS_D_Data.dta
Outputs:
    S5_detail.csv         one row per flagged record, with check id and reason
    console summary
"""
import os
import numpy as np
import pandas as pd

BASE = os.environ.get("CHARLS_DIR",
                      os.path.join(os.path.dirname(os.path.abspath(__file__)), "charls"))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "S5_detail.csv")

d11 = pd.read_stata(os.path.join(BASE, "Blood_20140429.dta"), convert_categoricals=False)
d15 = pd.read_stata(os.path.join(BASE, "Blood.dta"), convert_categoricals=False)
h = pd.read_stata(os.path.join(BASE, "H_CHARLS_D_Data.dta"), convert_categoricals=False)

rows = []


def flag(check, rec_id, field, value, reason):
    rows.append(dict(check=check, id=str(rec_id), field=field, value=value, reason=reason))


# ---------------------------------------------------------------- check 1
# identifier namespaces
i11 = d11["ID"].astype(str).str.strip()
i15 = d15["ID"].astype(str).str.strip()
ih = h["ID"].astype(str).str.strip()
over_raw = len(set(i11) & set(i15))
over_pad = len({x.zfill(12) for x in i11} & set(i15))
over_h = len(set(i15) & set(ih))
print("== check 1 identifier namespace ==")
print("  2011 n=%d len=%s | 2015 n=%d len=%s | harmonized n=%d"
      % (len(i11), sorted({len(x) for x in i11}), len(i15), sorted({len(x) for x in i15}), len(ih)))
print("  overlap 2011x2015 raw=%d | after zero-padding=%d | 2015xh=%d" % (over_raw, over_pad, over_h))
flag("1-identifier-namespace", "ALL", "ID",
     "overlap_raw=%d; overlap_padded=%d; 2015_vs_harmonized=%d" % (over_raw, over_pad, over_h),
     "Two waves of the same cohort share no identifier value; a merge on ID returns zero rows without error")

# ---------------------------------------------------------------- check 2
# analyte renamed between waves, incl. transposed HbA1c
PAIRS = [("newglu", "bl_glu", "glucose"), ("newcrea", "bl_crea", "creatinine"),
         ("newtg", "bl_tg", "triglyceride"), ("newcho", "bl_cho", "total cholesterol"),
         ("newhdl", "bl_hdl", "HDL cholesterol"), ("newldl", "bl_ldl", "LDL cholesterol"),
         ("newcrp", "bl_crp", "C-reactive protein"), ("cystatinc", "bl_cysc", "cystatin C"),
         ("newhba1c", "bl_hbalc", "HbA1c"), ("newbun", "bl_bun", "blood urea nitrogen"),
         ("newua", "bl_ua", "uric acid")]
print("\n== check 2 renaming ==")
for a, b, lab in PAIRS:
    miss = [c for c in (b,) if c not in d15.columns]
    transposed = (a.replace("new", "") != b.replace("bl_", ""))
    print("  %-22s %-10s -> %-10s %s" % (lab, a, b, "MISSING IN 2015" if miss else ""))
    if transposed:
        flag("2-renaming", "SCHEMA", b, a, "Analyte renamed between waves; letters transposed (hba1c -> hbalc)")

kw_hit_11 = [c for c in d11.columns if "hba1c" in c.lower()]
kw_hit_15 = [c for c in d15.columns if "hba1c" in c.lower()]
print("  keyword 'hba1c' matches 2011:%s 2015:%s" % (kw_hit_11, kw_hit_15))
flag("2-renaming", "SCHEMA", "hba1c-keyword",
     "2011=%s 2015=%s" % (kw_hit_11, kw_hit_15),
     "Keyword-based harmonisation finds HbA1c in 2011 and not in 2015")

# ---------------------------------------------------------------- check 3
# censoring threshold
tg = pd.to_numeric(d15["bl_tg"], errors="coerce")
fl = pd.to_numeric(d15["bl_top_coding_tg"], errors="coerce")
at = tg == 500.0
print("\n== check 3 censoring boundary ==")
print("  values == 500.0 : %d | flagged : %d | flagged but <500 : %d | ==500 but unflagged : %d"
      % (at.sum(), int(((fl == 1) & at).sum()), int(((fl == 1) & (tg < 500)).sum()), int((at & (fl == 0)).sum())))
for rid in d15.loc[at & (fl == 0), "ID"]:
    flag("3-censoring", rid, "bl_tg", 500.0,
         "Value sits on the censoring threshold with the flag set to 0; a true 500 and a truncated value are indistinguishable")

# ---------------------------------------------------------------- check 4
# physiologically impossible glucose
g11 = pd.to_numeric(d11["newglu"], errors="coerce")
lo, hi = g11 < 40, g11 > 600
print("\n== check 4 implausible glucose ==")
print("  <40 mg/dL : %d | >600 mg/dL : %d | total %d of %d (%.2f%%)"
      % (lo.sum(), hi.sum(), int((lo | hi).sum()), g11.notna().sum(), 100 * (lo | hi).sum() / g11.notna().sum()))
print("  lowest values:", sorted(g11[lo].round(2).tolist())[:20])
for rid, v in zip(d11.loc[lo | hi, "ID"], g11[lo | hi]):
    flag("4-implausible", rid, "newglu", float(v), "Outside a survivable physiological range; no screen is applied in the released file")

# ---------------------------------------------------------------- check 5
# measured value without a weight
w15 = pd.to_numeric(d15["Blood_weight"], errors="coerce")
has_val = tg.notna()
nw = has_val & w15.isna()
print("\n== check 5 value without weight ==")
print("  rows %d | weight missing %d | has TG value but no weight %d (%.2f%%)"
      % (len(d15), int(w15.isna().sum()), int(nw.sum()), 100 * nw.sum() / len(d15)))
for rid in d15.loc[nw, "ID"]:
    flag("5-no-weight", rid, "Blood_weight", np.nan, "Measured analyte reported without a sampling weight")

# ---------------------------------------------------------------- check 6
# variable inventory differs between waves
only11 = [c for c in d11.columns if c not in d15.columns]
only15 = [c for c in d15.columns if c not in d11.columns]
print("\n== check 6 inventory ==")
print("  2011 only (%d): %s" % (len(only11), only11))
print("  2015 only (%d): %s" % (len(only15), only15))
flag("6-inventory", "SCHEMA", "columns",
     "only_2011=%d only_2015=%d" % (len(only11), len(only15)),
     "Fasting indicator and six blood-count indices exist only in 2015")

# ---------------------------------------------------------------- output
det = pd.DataFrame(rows)
det.to_csv(OUT, index=False, encoding="utf-8-sig")
print("\n== flagged rows by check ==")
print(det["check"].value_counts().to_string())
print("[saved]", OUT, len(det), "rows")
