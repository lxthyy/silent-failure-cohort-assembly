# -*- coding: utf-8 -*-
"""
checks_calibrated.py — 12 类静默失效的**校准版**共享检查（任务 A）

与包内原始 common_checks.py 的差异（每一处都在 A_adjudication.csv 里给出理由）：

C1  只在"该哨兵值落在这列自身非哨兵分布的 p1–p99 之内"时计为碰撞。
    这样 77/88/99 在问卷编码（拒答/不知道）里不会被误判，而在年龄 77 岁、
    脉搏 88 次/分、血糖 99 mg/dL 这类真实测量上会被正确计为碰撞。
C2  同时读两遍（pandas / pyreadstat），统计 SAS 特殊缺失被物化成
    5.3976e-79 或 0.0 的单元格数，以及两个读取器给出不同值的单元格数。
C3  保留原判据（行数/唯一 SEQN > 1.05），但把结果标记为"结构性一对多风险"，
    因为它命中的是处方表这类事件级文件，本身不是数据错误。
C4  SEQN 的 dtype / 小数 / 位数在模块间是否漂移（跨模块汇总）。
C7  区分"NHANES 规范的 LBX 惯用单位 + LBD SI 单位配对"与"非标准重名"；
    前者只记为命名风险，不计为缺陷。
C8  同时含全样本权重与子样本权重的模块；另外记录子样本权重的特殊缺失数。
C9  短名是长名前缀且短名近空、长名近满时，再检查两列数值是否同一变量；
    仅在确为同一变量时计为缺陷（排除 LBXTO5 / LBXTO5IN 这类解释码）。
C10 跨周期同名变量的 5%–95% 分位跨度比 > 10 且中位数比 > 5 时计为单位漂移。

作者：Xin Li 校准版   |   生成：2026-11-05
"""
from __future__ import annotations
import re
from typing import Any

import numpy as np
import pandas as pd

# ── 哨兵值 ───────────────────────────────────────────────────────────────
SENTINELS = (77, 88, 99, 777, 888, 999, 7777, 8888, 9999, 77777, 88888, 99999)
SENTINEL_SET = set(SENTINELS)

# ── 权重列 ───────────────────────────────────────────────────────────────
FULL_W = ("WTINT", "WTMEC")
SUB_W = ("WTSAF", "WTSB", "WTSV", "WTSH", "WTSA", "WTSC", "WTPH", "WTDRD")

# ── 已知一人多行的事件级文件（本身不是缺陷，只作结构性风险登记）─────────
ONE_TO_MANY_OK = ("RXQ_RX", "RXQ_RX2", "RXQ_RX3", "DSQTOT", "DSQIDS",
                  "DSQCON", "DSQFOOD", "DSQDRINK", "DSQNUM")

# ── 单位后缀（用于 C7 重名归并）───────────────────────────────────────────
UNIT_SUFFIX = ("SI", "MGDL", "MMOLL", "UMOLL", "NMOLL", "PMOLL", "GDL", "UGL",
               "MGL", "NGL", "PGL", "IUL", "MIUL", "UIUL", "MEQL", "MMOL",
               "UMOL", "KUL", "FL", "PCT", "G", "MG")
# 解释/格式码后缀：不是同一分析物的重复列
CODE_SUFFIX = ("IN", "LC", "IND", "FLAG")


def norm_col(c: Any) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(c).upper())


def numeric_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]


def is_identifier(s: pd.Series, name: str) -> bool:
    """排除 SEQN 这类标识列：它们不是测量值，不能参与哨兵码/单位漂移判定。"""
    if norm_col(name) in ("SEQN", "ID", "IDCODE"):
        return True
    v = s.dropna()
    return bool(len(v) > 30 and v.nunique() == len(v)
                and pd.api.types.is_numeric_dtype(v)
                and not norm_col(name)[:2] in ("LB", "BM", "BP", "DR", "UR"))


# ── 连续测量列白名单（只在这些列上判 C1/C10，问卷编码一律不进入）──────────
MEAS_PREFIX = ("BMX", "BPX", "LBX", "LBD", "LB2", "LBD2", "LBXS", "LBDS",
               "URX", "URD", "DR1T", "DR2T", "WT", "INDF")
MEAS_EXACT = ("RIDAGEYR", "RIDAGEMN", "RIDAGEEX", "DMDHRAGE")
# 明确排除的问卷/编码前缀（防止未来白名单被误扩）
QUESTIONNAIRE_PREFIX = (
    "CDQ", "DIQ", "MCQ", "KIQ", "OCQ", "RXQ", "SMQ", "PAQ", "ALQ", "DPQ",
    "HUQ", "SLQ", "OHX", "OHQ", "DBQ", "IND", "DMD", "RID", "SIA", "FIA",
    "MIA", "AIA", "DMQ", "WHQ", "OSQ", "PUQ", "FSQ", "CSQ", "AGQ", "AUQ",
    "CVQ", "DEQ", "ECQ", "FNQ", "HCQ", "HIQ", "HOQ", "IMQ", "INQ", "PFQ",
    "RHQ", "SPQ", "SUQ", "VTQ", "BPQ", "DED", "FED", "PAD", "SLD", "SMD",
    "KIQ", "WHD", "DID", "SIA",
)


def is_measurement(name: str) -> bool:
    """连续测量列：白名单前缀命中、且不落在问卷编码前缀里。"""
    n = norm_col(name)
    if n in MEAS_EXACT:
        return True
    if any(n.startswith(p) for p in QUESTIONNAIRE_PREFIX):
        return False
    if _analyte_key(name) is None:      # LBXTO5IN 这类解释码
        return False
    return any(n.startswith(p) for p in MEAS_PREFIX)


def _cycle_of(module: str) -> str:
    """从模块名取出周期标签：DEMO_E -> E, 1999_LAB02 -> 1999, ALQ_B -> B."""
    m = re.match(r"^(\d{4})_", module)
    if m:
        return m.group(1)
    m = re.search(r"_([A-Z])\.xpt$", module, re.I)
    if m:
        return m.group(1).upper()
    return "A"


def _non_sentinel_stats(s: pd.Series) -> dict | None:
    """用非哨兵值估计该列的"正常"分布，避免哨兵码把分位数拉高。"""
    v = s.dropna().to_numpy(dtype=float, copy=False)
    if v.size < 30:
        return None
    v = v[~np.isin(v, SENTINELS)]
    if v.size < 30:
        return None
    q01, q05, q50, q95, q99 = np.percentile(v, [1, 5, 50, 95, 99])
    return {"q01": float(q01), "q05": float(q05), "q50": float(q50),
            "q95": float(q95), "q99": float(q99),
            "span": float(q95 - q05), "n_ref": int(v.size)}


# ══════════════════════════════════════════════════════════════════════
# C1  Missing-code collision（限连续测量：哨兵值落在本列正常分布之内）
# ══════════════════════════════════════════════════════════════════════
def check_1_missing_code_collision(df: pd.DataFrame) -> dict:
    evidence, total = [], 0
    for c in numeric_cols(df):
        if not is_measurement(c):
            continue
        s = df[c].dropna()
        if s.nunique() < 5 or is_identifier(s, c):
            continue
        st = _non_sentinel_stats(s)
        if st is None:
            continue
        hits = {}
        for v in SENTINELS:
            n = int((s == v).sum())
            if n > 0 and st["q01"] <= v <= st["q99"]:
                hits[v] = n
        if hits:
            n_hit = int(sum(hits.values()))
            total += n_hit
            evidence.append({
                "col": str(c), "n": n_hit, "codes": hits,
                "share_pct": round(100 * n_hit / len(s), 3),
                "non_sentinel_range": [round(st["q01"], 3), round(st["q99"], 3)],
            })
    evidence.sort(key=lambda e: -e["n"])
    return {"hit": bool(evidence), "n": int(total), "evidence": evidence[:25]}


# ══════════════════════════════════════════════════════════════════════
# C2  Parser encoding divergence
# ══════════════════════════════════════════════════════════════════════
def check_2_parser_encoding(df: pd.DataFrame, alt: pd.DataFrame | None) -> dict:
    """alt 为第二读取器（pyreadstat）读出的同一文件。"""
    evidence, total, divergent, total_meas = [], 0, 0, 0
    for c in numeric_cols(df):
        s = df[c]
        tiny = ((s.abs() < 1e-70) & (s != 0) & s.notna())
        n_tiny = int(tiny.sum())
        if n_tiny == 0:
            continue
        total += n_tiny
        if is_measurement(c):
            total_meas += n_tiny
        ev = {"col": str(c), "n_special_missing": n_tiny,
              "example": float(s[tiny].iloc[0]),
              "share_pct": round(100 * n_tiny / max(1, int(s.notna().sum())), 3)}
        if alt is not None and c in alt.columns:
            a = s.to_numpy(dtype=float, copy=False)
            b = alt[c].to_numpy(dtype=float, copy=False)
            if len(a) == len(b):
                # 同一 SAS 特殊缺失：pandas 给 5.4e-79，pyreadstat(float32) 下溢为 0
                d = int((tiny.to_numpy() & (b == 0)).sum())
                na = int((np.isnan(a) != np.isnan(b)).sum())
                ev["n_5e-79_vs_0"] = d
                ev["n_na_disagree"] = na
                divergent += d + na
        evidence.append(ev)
    evidence.sort(key=lambda e: -e["n_special_missing"])
    return {"hit": total > 0, "n": int(total),
            "n_measurement": int(total_meas),
            "evidence": evidence[:25], "n_divergent": int(divergent)}


# ══════════════════════════════════════════════════════════════════════
# C3  One-to-many inflation
# ══════════════════════════════════════════════════════════════════════
def check_3_row_inflation(df: pd.DataFrame, module: str = "") -> dict:
    key = next((c for c in df.columns if norm_col(c) == "SEQN"), None)
    if key is None:
        return {"hit": False, "n": 0, "evidence": ["no SEQN"]}
    n_rows, n_uniq = len(df), int(df[key].nunique(dropna=True))
    if n_uniq == 0:
        return {"hit": False, "n": 0, "evidence": ["SEQN all missing"]}
    ratio = n_rows / n_uniq
    event_level = ratio > 1.05
    documented = any(norm_col(module).startswith(norm_col(p)) for p in ONE_TO_MANY_OK)
    # 已知事件级文件：一人多药/多食物，本身不是缺陷，只记为结构性风险
    hit = bool(event_level and not documented)
    return {
        "hit": hit, "n": int(n_rows - n_uniq) if hit else 0,
        "evidence": [{"rows": n_rows, "unique_seqn": n_uniq,
                      "ratio": round(ratio, 3),
                      "tier": ("documented_event_level" if documented
                               else "structural_hazard" if event_level else "person_level"),
                      "documented_event_level": documented}],
    }


# ══════════════════════════════════════════════════════════════════════
# C4  Identifier key-space drift（单模块只采集，跨模块判定）
# ══════════════════════════════════════════════════════════════════════
def check_4_id_profile(df: pd.DataFrame) -> dict:
    key = next((c for c in df.columns if norm_col(c) == "SEQN"), None)
    if key is None:
        return {"hit": False, "n": 0,
                "evidence": [{"dtype": None, "has_seqn": False}]}
    s = df[key].dropna()
    if s.empty:
        return {"hit": False, "n": 0, "evidence": [{"dtype": str(df[key].dtype), "empty": True}]}
    numeric = pd.api.types.is_numeric_dtype(s)
    return {"hit": False, "n": 0, "evidence": [{
        "has_seqn": True,
        "dtype": str(df[key].dtype),
        "min": float(s.min()) if numeric else None,
        "max": float(s.max()) if numeric else None,
        "n_unique": int(s.nunique()),
        "has_fraction": bool(((s % 1) != 0).any()) if numeric else None,
        "width": (len(str(int(s.max()))) if numeric else None),
    }]}


# ══════════════════════════════════════════════════════════════════════
# C7  Format switch / 同一分析物多写法
# ══════════════════════════════════════════════════════════════════════
def _analyte_key(name: str) -> str | None:
    n = norm_col(name)
    n = re.sub(r"^(LBX|LBD|LB2|LBD2|LBXS|LBDS)", "", n)
    for suf in UNIT_SUFFIX:
        if n.endswith(suf) and len(n) > len(suf) + 1:
            n = n[: -len(suf)]
            break
    for suf in CODE_SUFFIX:
        if n.endswith(suf) and len(n) > len(suf) + 1:
            return None          # 解释码，不当成重复列
    return n or None


def check_7_analyte_duplication(df: pd.DataFrame) -> dict:
    """只做两件事：
    ① 登记本模块里每个 LBX/LBD 分析物的形态（conventional / si），供跨周期配对完整性判定；
    ② 统计本模块内成对齐全的分析物数（仅作描述，不计为缺陷）。
    非标准重名不再进入 hit —— 名称近似不等于同一变量，误报率太高。
    """
    inventory = []
    for c in df.columns:
        n = norm_col(c)
        if not n.startswith(("LBX", "LBD", "LB2", "LBD2")):
            continue
        k = _analyte_key(c)
        if not k or not df[c].notna().any():
            continue
        inventory.append({
            "analyte": k,
            "form": "si" if n.endswith("SI") else "conventional",
            "col": str(c),
        })
    forms_by_analyte: dict[str, set] = {}
    for it in inventory:
        forms_by_analyte.setdefault(it["analyte"], set()).add(it["form"])
    n_std = sum(1 for v in forms_by_analyte.values() if len(v) == 2)
    return {"hit": False, "n": 0, "evidence": [],
            "n_standard_dual_unit": n_std, "inventory": inventory}


# ══════════════════════════════════════════════════════════════════════
# C8  Subsample weight mismatch
# ══════════════════════════════════════════════════════════════════════
def check_8_weight_mismatch(df: pd.DataFrame) -> dict:
    cols = [norm_col(c) for c in df.columns]
    full = [c for c in cols if any(c.startswith(p) for p in FULL_W)]
    sub = [c for c in cols if any(c.startswith(p) for p in SUB_W)]
    wt_tiny = {}
    for c in df.columns:
        n = norm_col(c)
        if not n.startswith("WT"):
            continue
        s = df[c]
        k = int((((s.abs() < 1e-70) & (s != 0) & s.notna())).sum())
        if k:
            wt_tiny[str(c)] = k
    return {
        "hit": bool(full) and bool(sub),
        "n": len(sub) if (full and sub) else 0,
        "evidence": [{"full_weights": full, "subsample_weights": sub,
                      "weight_special_missing": wt_tiny}] if (full or sub) else [],
        "weight_special_missing_total": int(sum(wt_tiny.values())),
    }


# ══════════════════════════════════════════════════════════════════════
# C9  Short-name collision empties a column
# ══════════════════════════════════════════════════════════════════════
def check_9_shortname_collision(df: pd.DataFrame) -> dict:
    cols = [str(c) for c in df.columns]
    if len(df) == 0:
        return {"hit": False, "n": 0, "evidence": []}
    norm = {c: norm_col(c) for c in cols}
    confirmed, candidates = [], []
    for short in cols:
        s = norm[short]
        if len(s) < 3:
            continue
        fill_short = 1 - df[short].isna().mean()
        if fill_short >= 0.05:
            continue
        for long in cols:
            if long == short:
                continue
            lg = norm[long]
            if not (lg.startswith(s) and len(lg) > len(s)):
                continue
            fill_long = 1 - df[long].isna().mean()
            if fill_long <= 0.95:
                continue
            rec = {"short": short, "long": long,
                   "fill_short": round(float(fill_short), 4),
                   "fill_long": round(float(fill_long), 4)}
            # 确认是否同一变量：在两者都非空的交集上比较
            both = df[[short, long]].dropna()
            same = False
            if len(both) >= 5 and pd.api.types.is_numeric_dtype(df[short]) \
                    and pd.api.types.is_numeric_dtype(df[long]):
                x = both[short].to_numpy(float)
                y = both[long].to_numpy(float)
                same = bool(np.allclose(x, y, rtol=1e-6, atol=1e-9))
            rec["same_values"] = same
            if lg.endswith("IN") or lg.endswith("LC"):
                rec["adjudication"] = "interpretation_code_not_duplicate"
                candidates.append(rec)
            elif same:
                confirmed.append(rec)
            else:
                rec["adjudication"] = "prefix_only_not_duplicate"
                candidates.append(rec)
    return {"hit": bool(confirmed), "n": len(confirmed),
            "evidence": (confirmed + candidates)[:10],
            "n_candidates": len(candidates)}


# ══════════════════════════════════════════════════════════════════════
# C10  Column-name unit drift（跨周期，单模块只输出 profile）
# ══════════════════════════════════════════════════════════════════════
def check_10_unit_profile(df: pd.DataFrame, max_cols: int = 400) -> dict:
    prof = {}
    for c in numeric_cols(df)[:max_cols]:
        if not is_measurement(c):
            continue
        s = df[c].dropna()
        if len(s) < 30 or is_identifier(s, c):
            continue
        st = _non_sentinel_stats(s)
        if st is None:
            continue
        v = s.to_numpy(float)
        prof[str(c)] = {
            "q05": round(st["q05"], 6), "q50": round(st["q50"], 6),
            "q95": round(st["q95"], 6), "span": round(st["span"], 6),
            "decimal_share": round(float(np.mean((v % 1) != 0)), 4),
            "missing_rate": round(float(df[c].isna().mean()), 4),
        }
    return {"hit": False, "n": 0, "evidence": [prof]}


def run_all(df: pd.DataFrame, alt: pd.DataFrame | None = None, module: str = "") -> dict:
    out = {}
    fns = {
        "C1_missing_code": lambda: check_1_missing_code_collision(df),
        "C2_parser_encoding": lambda: check_2_parser_encoding(df, alt),
        "C3_row_inflation": lambda: check_3_row_inflation(df, module),
        "C4_id_drift": lambda: check_4_id_profile(df),
        "C7_analyte_dup": lambda: check_7_analyte_duplication(df),
        "C8_weight_mismatch": lambda: check_8_weight_mismatch(df),
        "C9_shortname_collision": lambda: check_9_shortname_collision(df),
        "C10_unit_drift": lambda: check_10_unit_profile(df),
    }
    for k, fn in fns.items():
        try:
            out[k] = fn()
        except Exception as e:
            out[k] = {"hit": False, "n": 0, "evidence": [f"ERROR: {type(e).__name__}: {e}"]}
    return out
