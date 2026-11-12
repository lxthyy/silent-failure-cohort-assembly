# -*- coding: utf-8 -*-
"""
B_headtohead_fair.py — 四种验证方法的头对头对比（无泄漏重设计）

与包内原始 B_headtohead.py 的差异（对应 2026-11-05 的三条修正）：

1. 基表不再截断。按表名精确选列，构造 16 列的连续分析物面板
   （年龄、体重、身高、BMI、腰围、收缩压×2、舒张压×2、TC、HDL、TG、葡萄糖，
   每项同时保留 mg/dL 与 SI 两种单位）。
2. M4 不再接触当周期干净表。四条规则只用：
   ①数据自身 + 文档化阈值；②分母记账；③同文件两条独立读写路径 + 内部恒等式；
   ④另一周期（2005–2006）的参考剖面。
3. 新增 C0 阴性对照：四法在干净数据上必须全部判"未检出"，否则脚本直接退出。
4. D5 只在"非整数占比 ≥ 50%"的列上注入；D4 按原定义做整表右移一列。
5. 额外产出 B_sensitivity.csv：同一单位、同分布的单列位移（更隐蔽的 D4 变体），
   用于说明各方法的检出边界。

用法：
    python B_headtohead_fair.py --data <nhanes_cache> --out <out_dir> --cycle 2007 --reps 30
"""
from __future__ import annotations
import argparse
import gc
import hashlib
import json
import os
import sys
import time
from datetime import datetime
from io import StringIO

import numpy as np
import pandas as pd

# ── 周期文件与面板 ────────────────────────────────────────────────────────
CYCLE_FILES = {
    2007: {"DEMO": "DEMO_E.xpt", "BMX": "BMX_E.xpt", "BPX": "BPX_E.xpt",
           "TCHOL": "TCHOL_E.xpt", "HDL": "HDL_E.xpt",
           "TRIGLY": "TRIGLY_E.xpt", "GLU": "GLU_E.xpt"},
    2005: {"DEMO": "DEMO_D.xpt", "BMX": "BMX_D.xpt", "BPX": "BPX_D.xpt",
           "TCHOL": "TCHOL_D.xpt", "HDL": "HDL_D.xpt",
           "TRIGLY": "TRIGLY_D.xpt", "GLU": "GLU_D.xpt"},
}
WANT = {
    "DEMO": ["RIDAGEYR"],
    "BMX": ["BMXWT", "BMXHT", "BMXBMI", "BMXWAIST"],
    "BPX": ["BPXSY1", "BPXSY2", "BPXDI1", "BPXDI2"],
    "TCHOL": ["LBXTC", "LBDTCSI"],
    "HDL": ["LBDHDD", "LBDHDDSI"],
    "TRIGLY": ["LBXTR", "LBDTRSI"],
    "GLU": ["LBXGLU", "LBDGLUSI"],
}
# 文档化生理范围（不看当周期数据，写在脚本里作为"准则先行"）
RANGES = {
    "RIDAGEYR": (0, 95),
    "BMXWT": (2, 450), "BMXHT": (60, 240), "BMXBMI": (8, 160), "BMXWAIST": (20, 250),
    "BPXSY1": (40, 320), "BPXSY2": (40, 320), "BPXDI1": (0, 220), "BPXDI2": (0, 220),
    "LBXTC": (50, 800), "LBDTCSI": (1.0, 25.0),
    "LBDHDD": (5, 250), "LBDHDDSI": (0.1, 6.0),
    "LBXTR": (5, 4000), "LBDTRSI": (0.1, 45.0),
    "LBXGLU": (20, 800), "LBDGLUSI": (1.0, 45.0),
}
SI_COLS = ["LBDTCSI", "LBDHDDSI", "LBDTRSI", "LBDGLUSI"]
CONV_FACTOR = {"LBXTC": 38.67, "LBDHDD": 38.67, "LBXTR": 88.57, "LBXGLU": 18.0}
DECIMAL_COLS = ["BMXWT", "BMXHT", "BMXBMI", "BMXWAIST",
                "LBDTCSI", "LBDHDDSI", "LBDTRSI", "LBDGLUSI"]
DEFECTS = ["D1_sentinel", "D2_unit_rescale", "D3_row_dup", "D4_col_shift", "D5_round_int"]
METHODS = ["M1_double_programming", "M2_reproducibility",
           "M3_kahn_checklist", "M4_four_rule_protocol"]
SENTINEL = 99
DUP_THRESHOLD = 1.05
MISSING_TOL = 0.15
MEDIAN_BAND = (0.6, 1.6)
SPAN_BAND = (0.4, 2.5)
DECIMAL_RATIO = 0.5
CLOSURE_TOL = 1.5          # BMI 内部恒等式允许的中位绝对误差


def log(fh, msg):
    line = f"[{datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def read_xpt(path):
    last = None
    for enc in ("utf-8", "latin1"):
        try:
            return pd.read_sas(path, format="xport", encoding=enc)
        except Exception as e:
            last = e
    raise last


def load_cycle(data_dir, cycle):
    """按表名精确选列，左连成一张 16 列分析物面板；特殊缺失统一置 NaN。"""
    index = {f.lower(): f for f in os.listdir(data_dir)}
    files = CYCLE_FILES[cycle]
    base = None
    for tab, fn in files.items():
        real = index.get(fn.lower())
        if real is None:
            continue
        df = read_xpt(os.path.join(data_dir, real)).drop_duplicates(subset=["SEQN"])
        cols = ["SEQN"] + [c for c in WANT[tab] if c in df.columns]
        if len(cols) <= 1:
            continue
        df = df[cols]
        if base is None:
            base = df
        else:
            base = base.merge(df, on="SEQN", how="left")
        del df
        gc.collect()
    if base is None:
        raise RuntimeError(f"cycle {cycle} 无可用文件")
    for c in base.columns:
        if c == "SEQN":
            continue
        s = base[c].astype(float)
        s[(s.abs() < 1e-70) & (s != 0)] = np.nan     # SAS 特殊缺失 -> 真缺失
        base[c] = s
    return base


def profile(df):
    """参考剖面：中位数、p05/p95、跨度、小数占比、缺失率。"""
    out = {}
    for c in df.columns:
        if c == "SEQN":
            continue
        s = df[c].dropna()
        if len(s) < 30:
            continue
        q05, q50, q95 = np.percentile(s, [5, 50, 95])
        out[c] = {
            "n": int(len(s)), "median": float(q50), "p05": float(q05), "p95": float(q95),
            "span": float(q95 - q05),
            "decimal_share": float(np.mean((s.to_numpy(float) % 1) != 0)),
            "missing_rate": float(df[c].isna().mean()),
        }
    return out


# ══════════════════════════════════════════════════════════════════════
# 缺陷注入（只作用于分析物列）
# ══════════════════════════════════════════════════════════════════════
def _pick(rng, cols):
    return cols[int(rng.integers(0, len(cols)))]


def inject(df, kind, rng):
    d = df.copy()
    detail = {}
    if kind == "D1_sentinel":
        col = _pick(rng, [c for c in SI_COLS if df[c].notna().sum() > 100])
        idx = d.index[d[col].notna()].to_numpy()
        n = max(1, int(0.01 * len(idx)))
        sel = rng.choice(idx, size=n, replace=False)
        d.loc[sel, col] = SENTINEL
        detail = {"col": col, "n_rows": int(n), "value": SENTINEL}

    elif kind == "D2_unit_rescale":
        col = _pick(rng, [c for c in CONV_FACTOR if df[c].notna().sum() > 100])
        f = CONV_FACTOR[col]
        d[col] = d[col] * f
        detail = {"col": col, "factor": f}

    elif kind == "D3_row_dup":
        n = max(1, int(len(d) * 0.05))
        dup = d.sample(n=n, random_state=int(rng.integers(0, 10 ** 6)))
        d = pd.concat([d, dup], ignore_index=True)
        detail = {"n_rows_added": int(n)}

    elif kind == "D4_col_shift":
        cols = [c for c in d.columns if c != "SEQN"]
        vals = d[cols].to_numpy(float)
        d[cols] = np.roll(vals, 1, axis=1)
        detail = {"shift": "all analyte columns rotated right by one"}

    elif kind == "D5_round_int":
        col = _pick(rng, [c for c in DECIMAL_COLS if df[c].notna().sum() > 100])
        d[col] = np.round(d[col])
        detail = {"col": col}

    return d, detail


def inject_single_column_shift(df, rng):
    """敏感性分析：把收缩压 1 换成收缩压 2（同单位、同分布、值域内）。"""
    d = df.copy()
    d["BPXSY1"] = d["BPXSY2"]
    return d, {"from": "BPXSY2", "to": "BPXSY1"}


# ══════════════════════════════════════════════════════════════════════
# M1 双编程
# ══════════════════════════════════════════════════════════════════════
def m1_double_programming(test, ref, meta):
    num = [c for c in test.columns if c != "SEQN"]
    for c in num:
        s = test[c].astype(float)
        # 实现 A：pandas
        a = (s.mean(), s.std(ddof=1), s.median(), s.quantile(0.95), s.notna().sum())
        # 实现 B：numpy 手算
        arr = np.asarray(s, dtype=float)
        arr = arr[~np.isnan(arr)]
        n = len(arr)
        if n == 0:
            b = (np.nan, np.nan, np.nan, np.nan, 0)
        else:
            b = (arr.mean(), arr.std(ddof=1) if n > 1 else np.nan,
                 np.median(arr), np.percentile(arr, 95), n)
        for x, y in zip(a, b):
            if pd.isna(x) != pd.isna(y):
                return True
            if not pd.isna(x) and abs(float(x) - float(y)) > 1e-9 * max(1.0, abs(float(x))):
                return True
    return False


# ══════════════════════════════════════════════════════════════════════
# M2 复现检查
# ══════════════════════════════════════════════════════════════════════
def m2_reproducibility(test, ref, meta):
    def run_once(df):
        x = df.sort_values("SEQN", kind="mergesort").reset_index(drop=True)
        return hashlib.md5(pd.util.hash_pandas_object(x, index=False).values.tobytes()).hexdigest()
    return run_once(test) != run_once(test.copy())


# ══════════════════════════════════════════════════════════════════════
# M3 Kahn checklist：conformance / completeness / plausibility
# ══════════════════════════════════════════════════════════════════════
def m3_kahn_checklist(test, ref, meta):
    # conformance：列齐全、类型为数值、无重复主键
    if set(meta["columns"]) - set(test.columns):
        return True
    if test["SEQN"].duplicated().any():
        return True
    for c in meta["columns"]:
        if not pd.api.types.is_numeric_dtype(test[c]):
            return True
    # completeness：缺失率不得比参考高 15 个百分点
    for c in meta["columns"]:
        if c in ref and test[c].isna().mean() - ref[c]["missing_rate"] > MISSING_TOL:
            return True
    # plausibility：文档化生理范围
    for c in meta["columns"]:
        lo, hi = RANGES[c]
        s = test[c].dropna()
        if s.empty:
            continue
        if (s < lo).any() or (s > hi).any():
            return True
    return False


# ══════════════════════════════════════════════════════════════════════
# M4 四规则协议（不接触当周期干净表）
# ══════════════════════════════════════════════════════════════════════
def m4_four_rule_protocol(test, ref, meta):
    cols = meta["columns"]

    # 规则 1：准则先行 —— 文档化范围 + 分析物齐全
    if set(cols) - set(test.columns):
        return True
    for c in cols:
        lo, hi = RANGES[c]
        s = test[c].dropna()
        if not s.empty and ((s < lo).any() or (s > hi).any()):
            return True

    # 规则 2：分母记账 —— 行数与唯一 SEQN 必须对上
    if test["SEQN"].duplicated().any():
        return True
    if len(test) / max(1, test["SEQN"].nunique()) > DUP_THRESHOLD:
        return True

    # 规则 3：异构碰撞到值级 —— 同文件两条独立读写路径 + 内部恒等式
    recon = pd.read_csv(StringIO(test.to_csv(index=False, float_format="%.17g")))
    if len(recon) != len(test):
        return True
    for c in cols:
        a = test[c].to_numpy(float)
        b = recon[c].to_numpy(float)
        same = np.isclose(a, b, equal_nan=True, rtol=1e-12, atol=1e-12)
        if not same.all():
            return True
    if {"BMXWT", "BMXHT", "BMXBMI"} <= set(cols):
        derived = test["BMXWT"] / (test["BMXHT"] / 100.0) ** 2
        resid = np.abs(test["BMXBMI"] - derived).dropna()
        if len(resid) > 100 and resid.median() > CLOSURE_TOL:
            return True

    # 规则 4：禁自指 —— 与另一周期参考剖面碰撞（中位数/跨度/小数占比）
    for c in cols:
        if c not in ref or c == "RIDAGEYR":
            continue          # 年龄的周期差异是人口构成，不是单位/精度缺陷
        s = test[c].dropna()
        if len(s) < 30:
            continue
        q05, q50, q95 = np.percentile(s, [5, 50, 95])
        r = ref[c]
        if r["median"] > 0 and abs(q50) > 0:
            ratio = q50 / r["median"]
            if ratio > MEDIAN_BAND[1] or ratio < MEDIAN_BAND[0]:
                return True
        if r["span"] > 0:
            ratio = (q95 - q05) / r["span"]
            if ratio > SPAN_BAND[1] or ratio < SPAN_BAND[0]:
                return True
        if c in DECIMAL_COLS and r["decimal_share"] > 0.5:
            if float(np.mean((s.to_numpy(float) % 1) != 0)) < DECIMAL_RATIO * r["decimal_share"]:
                return True
    return False


METHOD_FN = {
    "M1_double_programming": m1_double_programming,
    "M2_reproducibility": m2_reproducibility,
    "M3_kahn_checklist": m3_kahn_checklist,
    "M4_four_rule_protocol": m4_four_rule_protocol,
}


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cycle", type=int, default=2007)
    ap.add_argument("--ref-cycle", type=int, default=2005)
    ap.add_argument("--reps", type=int, default=30)
    args = ap.parse_args()

    data_dir = os.path.abspath(args.data)
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)
    logp = open(os.path.join(out_dir, "B_run_log.txt"), "a", encoding="utf-8")
    log(logp, f"cycle={args.cycle} ref_cycle={args.ref_cycle} reps={args.reps}")

    clean = load_cycle(data_dir, args.cycle)
    ref_df = load_cycle(data_dir, args.ref_cycle)
    ref = profile(ref_df)
    meta = {"columns": [c for c in clean.columns if c != "SEQN"]}
    log(logp, f"panel: {clean.shape[0]} rows x {len(meta['columns'])} analytes")
    log(logp, "analytes: " + ", ".join(meta["columns"]))

    # 参考剖面落盘（审计用）
    pd.DataFrame([{"col": c, **ref[c]} for c in ref]).to_csv(
        os.path.join(out_dir, "B_base_profile.csv"), index=False, encoding="utf-8-sig")

    # ── C0 阴性对照：干净数据上必须全部未检出 ────────────────────────────
    c0 = {}
    for m in METHODS:
        try:
            c0[m] = bool(METHOD_FN[m](clean, ref, meta))
        except Exception as e:
            log(logp, f"C0 {m} ERROR {type(e).__name__}: {e}")
            c0[m] = None
    log(logp, "C0 clean-data false positives: " + json.dumps(c0))
    if any(v is None or v is True for v in c0.values()):
        log(logp, "C0 FAILED — 校准未通过，脚本退出，不产出任何检出率")
        logp.close()
        sys.exit(2)

    # ── 主实验：5 缺陷 × reps × 4 方法 ───────────────────────────────────
    rows = []
    for d_i, defect in enumerate(DEFECTS):
        for r in range(args.reps):
            rng = np.random.default_rng(seed=10000 * (d_i + 1) + r)
            test, detail = inject(clean, defect, rng)
            for m in METHODS:
                try:
                    detected = bool(METHOD_FN[m](test, ref, meta))
                except Exception as e:
                    log(logp, f"  {defect} rep{r} {m} ERROR {type(e).__name__}: {e}")
                    detected = None
                rows.append({"defect": defect, "rep": r, "method": m,
                             "detected": detected, "seed": 10000 * (d_i + 1) + r,
                             "target": json.dumps(detail, ensure_ascii=False)})
            del test
            gc.collect()
        log(logp, f"{defect}: {args.reps} reps done")

    raw = pd.DataFrame(rows)
    if raw["detected"].isna().any():
        log(logp, "WARNING: 存在 detected=null，检查上方 ERROR 行")
    raw.to_csv(os.path.join(out_dir, "B_headtohead.csv"), index=False, encoding="utf-8-sig")

    summary, fig_rows = {}, []
    for defect in DEFECTS:
        summary[defect] = {}
        for m in METHODS:
            sub = raw[(raw.defect == defect) & (raw.method == m)].detected.dropna()
            k, n = int(sub.sum()), int(len(sub))
            lo, hi = wilson(k, n)
            summary[defect][m] = {"detected": k, "n": n,
                                  "rate_pct": round(100 * k / n, 1) if n else None,
                                  "ci95": [round(100 * lo, 1), round(100 * hi, 1)]}
            fig_rows.append({"defect": defect, "method": m,
                             "rate_pct": summary[defect][m]["rate_pct"],
                             "ci_lo": round(100 * lo, 1), "ci_hi": round(100 * hi, 1)})

    # ── 敏感性分析：同单位同分布的单列位移 ─────────────────────────────
    sens = []
    for r in range(args.reps):
        rng = np.random.default_rng(seed=90000 + r)
        test, detail = inject_single_column_shift(clean, rng)
        for m in METHODS:
            sens.append({"defect": "D4b_single_shift_same_unit", "rep": r, "method": m,
                         "detected": bool(METHOD_FN[m](test, ref, meta)),
                         "target": json.dumps(detail, ensure_ascii=False)})
        del test
        gc.collect()
    sens_df = pd.DataFrame(sens)
    sens_df.to_csv(os.path.join(out_dir, "B_sensitivity.csv"), index=False, encoding="utf-8-sig")
    sens_summary = {}
    for m in METHODS:
        sub = sens_df[sens_df.method == m].detected
        sens_summary[m] = {"detected": int(sub.sum()), "n": int(len(sub)),
                           "rate_pct": round(100 * float(sub.mean()), 1)}

    # ── 方法信息边界（写论文时必须一起报）──────────────────────────────
    scope = pd.DataFrame([
        {"method": "M1_double_programming", "object": "两套实现是否分歧",
         "uses_test": True, "uses_reference": False, "uses_denominator": False,
         "uses_precision": False, "note": "同一输入下共享同一缺陷，两者仍一致"},
        {"method": "M2_reproducibility", "object": "同一实现两次运行是否一致",
         "uses_test": True, "uses_reference": False, "uses_denominator": False,
         "uses_precision": False, "note": "确定性注入下两次必然一致"},
        {"method": "M3_kahn_checklist", "object": "conformance/completeness/plausibility",
         "uses_test": True, "uses_reference": True, "uses_denominator": True,
         "uses_precision": False, "note": "按文档化范围与参考缺失率逐项核对"},
        {"method": "M4_four_rule_protocol", "object": "准则先行/分母记账/异构碰撞/禁自指",
         "uses_test": True, "uses_reference": True, "uses_denominator": True,
         "uses_precision": True, "note": "跨周期分布碰撞，可检出仅改变精度的缺陷"},
    ])
    scope.to_csv(os.path.join(out_dir, "B_method_scope.csv"), index=False, encoding="utf-8-sig")

    out = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "cycle": args.cycle, "reference_cycle": args.ref_cycle, "reps": args.reps,
        "panel_rows": int(len(clean)), "panel_columns": meta["columns"],
        "defects": DEFECTS, "methods": METHODS,
        "detection": summary,
        "false_positive_control_C0": c0,
        "C0_pass": True,
        "sensitivity_D4b_same_unit_single_shift": sens_summary,
        "notes": {
            "M4_information_access": "只用当周期 test + 文档化范围 + 2005–2006 参考剖面；"
                                     "不接触 2007 干净表",
            "C0": "四法在干净数据上均判未检出；任一误报则脚本退出",
            "D4": "整表分析物列右移一列（原定义）；D4b 为同单位单列位移的敏感性分析",
        },
    }
    json.dump(out, open(os.path.join(out_dir, "B_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    pd.DataFrame(fig_rows).to_csv(os.path.join(out_dir, "B_figure_data.csv"),
                                  index=False, encoding="utf-8-sig")

    log(logp, "detection rates (%):")
    for defect in DEFECTS:
        line = "  " + defect + " -> " + "  ".join(
            f"{m.split('_')[0]}={summary[defect][m]['rate_pct']}" for m in METHODS)
        log(logp, line)
    log(logp, "D4b same-unit single-column shift: " +
              json.dumps({m: v["rate_pct"] for m, v in sens_summary.items()}))
    logp.close()
    print("\n产出目录:", out_dir)


if __name__ == "__main__":
    main()
