# -*- coding: utf-8 -*-
"""
Retrospective application of the automatable checks (Section 3.7).

Runs the checks that can be automated on three institutional datasets and reports,
for each, whether it fired, the key numbers and the elapsed time.
  D1 psychiatric   18,836 x 455
  D2 stroke        wide table
  D3 checkup       11,359
Writes a CSV and a text log to S4_OUT.
"""
import os, re, sys, time
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
_HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get("S4_BASE", _HERE)
OUT = os.environ.get("S4_OUT", os.path.join(BASE, "output")); os.makedirs(OUT, exist_ok=True)

DATASETS = {
    "D1_psychiatric": dict(
        path=os.environ.get("S4_D1", os.path.join(BASE, "D1_psychiatric.csv")),
        reader=pd.read_csv,
        id_cands=["住院号", "病人ID"], age_cands=["年龄"], sex_cands=["性别"],
        date_cands=["采样日期"], diag_cands=["诊断"],
        numeric_probe=["血红蛋白"], period_cands=["年份", "采样日期"]),
    "D2_stroke": dict(
        path=os.environ.get("S4_D2", os.path.join(BASE, "D2_stroke.csv")),
        reader=pd.read_csv,
        id_cands=["样本号", "住院号", "病人ID"], age_cands=["年龄"], sex_cands=["性别"],
        date_cands=["采样日期", "日期"], diag_cands=["诊断"],
        numeric_probe=["白细胞计数", "WBC"], period_cands=["年份", "采样日期", "日期"]),
    "D3_checkup": dict(
        path=os.environ.get("S4_D3", os.path.join(BASE, "D3_checkup.xlsx")),
        reader=pd.read_excel,
        id_cands=["住院号", "病人ID", "体检号", "ID"],
        age_cands=["年龄"], sex_cands=["性别"],
        date_cands=["采样日期", "日期", "体检日期"], diag_cands=["诊断"],
        numeric_probe=["TSH", "促甲状腺激素", "ALT"], period_cands=["年份", "采样日期", "日期"]),
}

L = []
def p(s=""):
    L.append(str(s)); print(s)

def pick(df, cands):
    for c in cands:
        if c in df.columns:
            return c
    return None

def norm_age(s):
    return pd.to_numeric(s.astype(str).str.replace(r"[^0-9.]", "", regex=True), errors="coerce")

def norm_sex(s):
    t = s.astype(str).str.strip()
    out = pd.Series(np.nan, index=t.index, dtype=object)
    out[t.isin(["男", "Male", "M", "1"])] = "M"
    out[t.isin(["女", "Female", "F", "2"])] = "F"
    return out

rows = []
def rec(ds, check, name, applicable, fired, detail, sec):
    rows.append(dict(dataset=ds, check=check, name=name, applicable=applicable,
                     fired=fired, detail=detail, seconds=round(sec, 2)))
    p("  [%s] %-12s fired=%-3s (%.2fs) | %s" % (check, name[:12], fired, sec, detail))

p("=" * 100)
p("Retrospective application across three datasets")
p("=" * 100)

loaded = {}
for dname, cfg in DATASETS.items():
    t0 = time.time()
    try:
        df = cfg["reader"](cfg["path"], low_memory=False, on_bad_lines="skip") \
            if cfg["reader"] is pd.read_csv else cfg["reader"](cfg["path"])
    except Exception as e:
        p("\n[%s] read failed: %s" % (dname, e)); continue
    dt = time.time() - t0
    loaded[dname] = (df, cfg)
    p("\n[%s] rows %d, cols %d, read in %.2fs" % (dname, df.shape[0], df.shape[1], dt))

for dname, (df, cfg) in loaded.items():
    p("\n" + "-" * 100)
    p("Dataset %s (n=%d, %d variables)" % (dname, len(df), df.shape[1]))
    p("-" * 100)
    cid = pick(df, cfg["id_cands"]); cage = pick(df, cfg["age_cands"])
    csex = pick(df, cfg["sex_cands"]); cdate = pick(df, cfg["date_cands"])
    cdiag = pick(df, cfg["diag_cands"]); cper = pick(df, cfg["period_cands"])

    # ---- Check 1 ----
    t0 = time.time()
    f = []
    if cid and cage:
        g = df.groupby(cid)[cage].nunique()
        f.append(("age-conflict IDs", int((g > 1).sum())))
    if cid and csex:
        g = df.groupby(cid)[csex].nunique()
        f.append(("sex-conflict IDs", int((g > 1).sum())))
    if cdiag:
        n_ts = int(df[cdiag].astype(str).str.contains(
            r"\d{4}[/\-]\d{1,2}[/\-]\d{1,2}|\d{1,2}:\d{2}:\d{2}", regex=True).sum())
        f.append(("date-like diag values", n_ts))
    if cdate:
        d = pd.to_datetime(df[cdate], errors="coerce")
        f.append(("invalid dates", int(d.isna().sum())))
    if cage:
        a = norm_age(df[cage])
        f.append(("implausible ages", int(((a < 0) | (a > 110)).sum())))
    sec = time.time() - t0
    total = sum(v for _, v in f)
    rec(dname, "1", "Field-row alignment", "Yes",
        "Yes" if total > 0 else "No",
        "; ".join("%s=%d" % (k, v) for k, v in f), sec)

    # ---- Check 3 ----
    if cdiag:
        t0 = time.time()
        s = df[cdiag].astype(str)
        n_uniq = int(s.nunique()); n_blank = int(((s == "nan") | (s.str.strip() == "")).sum())
        n_q = int(s.str.contains(r"待查|？|\?", regex=True).sum())
        n_c = int(s.str.contains(",").sum())
        sec = time.time() - t0
        fired = "Yes" if (n_uniq > 200 or n_blank > 0 or n_q > 0) else "No"
        rec(dname, "3", "Outcome composition", "Yes", fired,
            "distinct=%d, blank=%d, provisional=%d, composite=%d" % (n_uniq, n_blank, n_q, n_c), sec)
    else:
        rec(dname, "3", "Outcome composition", "No (no diagnosis field)", "n/a",
            "dataset has no diagnosis/outcome column", 0.0)

    # ---- Check 4 ---- definition sensitivity (not applicable: no repeat definition here)
    rec(dname, "4", "Definition sensitivity", "Conditional", "n/a",
        "requires a re-specifiable exposure/outcome definition; demonstrated on D1 only", 0.0)

    # ---- Check 5 ----
    cprobe = pick(df, cfg["numeric_probe"])
    if cprobe:
        t0 = time.time()
        v = pd.to_numeric(df[cprobe], errors="coerce").dropna()
        m1 = float(v.mean()); m2 = float(np.sum(v.values) / len(v.values))
        s1 = float(v.std(ddof=1)); s2 = float(np.sqrt(np.sum((v.values - m2) ** 2) / (len(v) - 1)))
        sec = time.time() - t0
        ok = abs(m1 - m2) < 1e-9 and abs(s1 - s2) < 1e-9
        rec(dname, "5", "Independent re-computation", "Yes", "No" if ok else "Yes",
            "%s mean %.6f vs %.6f, SD %.6f vs %.6f, agree=%s" % (cprobe, m1, m2, s1, s2, ok), sec)
    else:
        rec(dname, "5", "Independent re-computation", "No", "n/a", "no probe numeric column", 0.0)

    # ---- Check 6 ----
    if cper:
        t0 = time.time()
        if cper in ("年份",) or df[cper].dtype.kind in "iu":
            strata = df[cper].astype(str)
        else:
            strata = pd.to_datetime(df[cper], errors="coerce").dt.year.astype("Int64").astype(str)
        if cprobe:
            v = pd.to_numeric(df[cprobe], errors="coerce")
            tab = pd.DataFrame({"s": strata, "v": v}).dropna().groupby("s")["v"].median()
            tab = tab[tab.index.isin([x for x in tab.index if x not in ("nan", "<NA>")])]
            if len(tab) >= 2:
                spread = float(tab.max() - tab.min())
                sec = time.time() - t0
                rec(dname, "6", "Period stratification", "Yes",
                    "Yes" if spread > 0.5 else "No",
                    "%s median by period %s, spread=%.3f" % (
                        cprobe, dict(tab.round(2)), spread), sec)
            else:
                rec(dname, "6", "Period stratification", "Yes", "n/a", "fewer than 2 periods", time.time() - t0)
        else:
            rec(dname, "6", "Period stratification", "No", "n/a", "no probe numeric column", 0.0)
    else:
        rec(dname, "6", "Period stratification", "No", "n/a", "no period column", 0.0)

# ---- Check 2, D1 only (a second file exists) ----
p("\n" + "-" * 100)
p("Check 2  cohort identity (needs a second file)")
p("-" * 100)
DS_B = os.environ.get("S4_B", os.path.join(BASE, "D1_replication.xlsx"))
if "D1_psychiatric" in loaded and os.path.exists(DS_B):
    t0 = time.time()
    A = loaded["D1_psychiatric"][0]
    B = pd.read_excel(DS_B)
    def nn(s):
        return s.astype(str).str.replace(r"\s", "", regex=True).str.replace(r"[^\u4e00-\u9fa5A-Za-z]", "", regex=True)
    a = pd.DataFrame({"n": nn(A["姓名"]), "a": norm_age(A["年龄"]), "s": norm_sex(A["性别"])}).dropna(subset=["n", "a"])
    b = pd.DataFrame({"n": nn(B["姓名"]), "a": norm_age(B["年龄"]), "s": norm_sex(B["性别"])}).dropna(subset=["n", "a"])
    ka = set(map(tuple, a[["n", "a", "s"]].dropna().values))
    kb = set(map(tuple, b[["n", "a", "s"]].dropna().values))
    inter = ka & kb
    ovB = 100 * len(inter) / max(len(kb), 1)
    sec = time.time() - t0
    rec("D1_psychiatric", "2", "Cohort identity", "Yes", "Yes" if ovB > 10 else "No",
        "A triplets=%d, B triplets=%d, shared=%d (%.1f%% of B)" % (len(ka), len(kb), len(inter), ovB), sec)
for d in ["D2_stroke", "D3_checkup"]:
    rec(d, "2", "Cohort identity", "No (single file)", "n/a", "no second file for this dataset", 0.0)

# ---- Check 7 ----
for d in loaded:
    rec(d, "7", "Extrapolation audit", "Manual", "n/a", "requires manual claim-to-coefficient mapping", 0.0)

# ===== summary =====
p("\n" + "=" * 100)
p("[Summary: check x dataset]")
p("=" * 100)
T = pd.DataFrame(rows)
piv = T.pivot_table(index=["check", "name"], columns="dataset", values="fired", aggfunc="first")
p(piv.to_string())
p("\n[Elapsed time, seconds]")
tt = T.groupby(["check", "name"])["seconds"].sum().round(2)
p(tt.to_string())
p("\nTotal: %.2f s (three datasets, all checks)" % T["seconds"].sum())

T.to_csv(os.path.join(OUT, "S4_retrospective_application.csv"),
         index=False, encoding="utf-8-sig", quoting=1)
with open(os.path.join(OUT, "S4_retrospective_application.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(L))
p("\n>>> " + os.path.join(OUT, "S4_retrospective_application.csv"))
