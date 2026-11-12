# -*- coding: utf-8 -*-
"""
A_nhanes_census_calibrated.py — NHANES 全库静默失效巡检（校准版）

与包内原始脚本的差异：
  1. 只处理 286 个唯一模块；缓存里带 _YYYYMMDD_HHMMSS 后缀的重复副本被排除。
  2. 读取增加 latin1 回退，解决 RXQ_RX_H 的 UTF-8 解码失败。
  3. 检查实现见 checks_calibrated.py，判据已按变量语义校准并逐条留证。
  4. 额外产出 A_adjudication.csv：每个检查的判据、命中层、人工判读与误报控制。

用法：
    python A_nhanes_census_calibrated.py --data <nhanes_cache> --out <out_dir>

输出：
    A_defect_matrix.csv   286 行 × 23 列
    A_summary.json
    A_detail.jsonl
    A_run_log.txt
    A_adjudication.csv
    _progress.json / _cache/   断点续跑用（勿删）
"""
from __future__ import annotations
import argparse
import gc
import json
import os
import re
import sys
import time
import traceback
from collections import Counter, defaultdict
from datetime import datetime

import numpy as np
import pandas as pd

from checks_calibrated import run_all, norm_col, _cycle_of, SUB_W, FULL_W, is_measurement

CHECK_KEYS = ["C1_missing_code", "C2_parser_encoding", "C3_row_inflation",
              "C4_id_drift", "C7_analyte_dup", "C8_weight_mismatch",
              "C9_shortname_collision", "C10_unit_drift"]
DUP_RE = re.compile(r"_\d{8}_\d{6}\.xpt$", re.I)


def log(fh, msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def read_pandas(path):
    last = None
    for enc in ("utf-8", "latin1"):
        try:
            return pd.read_sas(path, format="xport", encoding=enc), f"pandas:{enc}"
        except Exception as e:
            last = e
    raise last


def read_alt(path):
    """第二读取器；失败返回 None（不影响主流程）。"""
    try:
        import pyreadstat
        try:
            df, _ = pyreadstat.read_xport(path, encoding="latin1")
        except Exception:
            df, _ = pyreadstat.read_xport(path)
        return df
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    data_dir = os.path.abspath(args.data)
    out_dir = os.path.abspath(args.out)
    cache_dir = os.path.join(out_dir, "_cache")
    os.makedirs(cache_dir, exist_ok=True)

    all_files = [f for f in os.listdir(data_dir) if f.lower().endswith(".xpt")]
    files = sorted(f for f in all_files if not DUP_RE.search(f))
    n_dup = len(all_files) - len(files)
    if args.limit:
        files = files[: args.limit]
    if not files:
        print(f"ERROR: {data_dir} 下没有 .XPT 文件")
        sys.exit(1)

    prog_path = os.path.join(out_dir, "_progress.json")
    done = set()
    if os.path.exists(prog_path):
        try:
            done = set(json.load(open(prog_path, encoding="utf-8")).get("done", []))
        except Exception:
            done = set()

    logp = open(os.path.join(out_dir, "A_run_log.txt"), "a", encoding="utf-8")
    log(logp, f"start: {len(files)} unique modules ({n_dup} duplicate copies excluded), "
              f"{len(done)} already done")
    log(logp, f"reader: pandas.read_sas (utf-8 -> latin1) + pyreadstat(latin1) for C2")

    t_start = time.time()
    for i, fn in enumerate(files, 1):
        if fn in done and os.path.exists(os.path.join(cache_dir, fn + ".json")):
            continue
        path = os.path.join(data_dir, fn)
        size_mb = os.path.getsize(path) / 1024 / 1024
        t0 = time.time()
        try:
            df, reader = read_pandas(path)
            alt = read_alt(path)
            res = run_all(df, alt, fn)
            row = {"module": fn, "n_rows": int(len(df)), "n_cols": int(df.shape[1]),
                   "size_mb": round(size_mb, 2), "sec": round(time.time() - t0, 2), "error": ""}
            for k in CHECK_KEYS:
                row[k + "_hit"] = bool(res[k]["hit"])
                row[k + "_n"] = int(res[k]["n"])
            detail = {
                "module": fn, "cycle": _cycle_of(fn), "reader": reader,
                "second_reader_ok": alt is not None,
                "shape": [int(df.shape[0]), int(df.shape[1])],
                "checks": {k: res[k] for k in CHECK_KEYS},
            }
            json.dump({"row": row, "detail": detail},
                      open(os.path.join(cache_dir, fn + ".json"), "w", encoding="utf-8"),
                      ensure_ascii=False, default=str)
            log(logp, f"[{i}/{len(files)}] {fn}  {df.shape[0]}x{df.shape[1]}  "
                      f"{row['sec']}s  reader={reader}  alt={'ok' if alt is not None else 'fail'}")
            del df, alt
            gc.collect()
            done.add(fn)
            if i % 5 == 0:
                json.dump({"done": sorted(done)}, open(prog_path, "w", encoding="utf-8"))
        except Exception as e:
            log(logp, f"[{i}/{len(files)}] {fn}  FAILED  {type(e).__name__}: {e}")
            logp.write(traceback.format_exc() + "\n")
            row = {"module": fn, "n_rows": None, "n_cols": None,
                   "size_mb": round(size_mb, 2), "sec": round(time.time() - t0, 2),
                   "error": f"{type(e).__name__}: {e}"}
            for k in CHECK_KEYS:
                row[k + "_hit"] = None
                row[k + "_n"] = None
            json.dump({"row": row, "detail": {"module": fn, "cycle": _cycle_of(fn),
                                              "error": row["error"]}},
                      open(os.path.join(cache_dir, fn + ".json"), "w", encoding="utf-8"),
                      ensure_ascii=False, default=str)
            done.add(fn)
            gc.collect()

    json.dump({"done": sorted(done)}, open(prog_path, "w", encoding="utf-8"))

    # ── 从缓存重建（支持断点续跑）────────────────────────────────────────
    rows, details = [], []
    for fn in files:
        p = os.path.join(cache_dir, fn + ".json")
        if not os.path.exists(p):
            continue
        obj = json.load(open(p, encoding="utf-8"))
        rows.append(obj["row"])
        details.append(obj["detail"])

    mat = pd.DataFrame(rows)
    for k in CHECK_KEYS:
        mat[k + "_hit"] = mat[k + "_hit"].astype("boolean").fillna(False).astype(bool)
        mat[k + "_n"] = pd.to_numeric(mat[k + "_n"], errors="coerce").fillna(0).astype(int)

    # ── C4 跨模块：SEQN 表示是否漂移 ─────────────────────────────────────
    sig, no_key = {}, []
    for d in details:
        ev = d.get("checks", {}).get("C4_id_drift", {}).get("evidence", [])
        e = ev[0] if ev and isinstance(ev[0], dict) else {}
        if not e.get("has_seqn", True):
            no_key.append(d["module"])
            continue
        sig[d["module"]] = (e.get("dtype"), e.get("has_fraction"),
                            None if e.get("width") is None else (e["width"] >= 7))
    drift_modules = set()
    if sig:
        variants = Counter(sig.values())
        if len(variants) > 1:
            major = variants.most_common(1)[0][0]
            drift_modules = {fn for fn, v in sig.items() if v != major}
    mat["C4_id_drift_hit"] = mat["module"].isin(drift_modules)
    mat["C4_id_drift_n"] = mat["C4_id_drift_hit"].astype(int)

    # ── C10 跨周期：同名变量单位漂移 ─────────────────────────────────────
    colmap = defaultdict(list)
    for d in details:
        prof = d.get("checks", {}).get("C10_unit_drift", {}).get("evidence", [{}])
        prof = prof[0] if prof and isinstance(prof[0], dict) else {}
        for col, st in prof.items():
            colmap[norm_col(col)].append((d["module"], d.get("cycle"), col, st))

    drift_c10 = defaultdict(list)
    c10_checked = 0
    for key, lst in colmap.items():
        cycles = {c for _, c, _, _ in lst}
        if len(cycles) < 3:
            continue
        c10_checked += 1
        # 覆盖人群必须可比：缺失率极差 > 30 个百分点时，跨度变化是覆盖变化，不是单位漂移
        mrs = [st.get("missing_rate", 0) for _, _, _, st in lst]
        if max(mrs) - min(mrs) > 0.30:
            continue
        med_span = float(np.median([s["span"] for _, _, _, s in lst if s["span"] > 0] or [0]))
        med_q50 = float(np.median([abs(s["q50"]) for _, _, _, s in lst if s["q50"] != 0] or [0]))
        if med_span <= 0 or med_q50 <= 0:
            continue
        for fn, cyc, col, st in lst:
            span_ratio = (st["span"] / med_span) if med_span else 1
            q50_ratio = (abs(st["q50"]) / med_q50) if med_q50 else 1
            if (span_ratio > 10 or span_ratio < 0.1) and (q50_ratio > 5 or q50_ratio < 0.2):
                drift_c10[fn].append({"col": col, "span_ratio_to_median": round(span_ratio, 2),
                                      "median_ratio_to_median": round(q50_ratio, 2)})
    mat["C10_unit_drift_hit"] = mat["module"].isin(drift_c10)
    mat["C10_unit_drift_n"] = mat["module"].map(lambda m: len(drift_c10.get(m, [])))

    # ── C7 跨周期：LBX/LBD 配对完整性（缺一半即格式切换风险）────────────
    pair_inv = defaultdict(lambda: defaultdict(set))
    pair_mods = defaultdict(list)
    for d in details:
        inv = d.get("checks", {}).get("C7_analyte_dup", {}).get("inventory", [])
        for item in inv:
            pair_inv[item["analyte"]][d.get("cycle")].add(item["form"])
            pair_mods[(item["analyte"], d.get("cycle"))].append(d["module"])
    incomplete = defaultdict(list)
    for analyte, byc in pair_inv.items():
        if not any(len(f) == 2 for f in byc.values()):
            continue          # 从未成对出现过的分析物，不判
        for cyc, forms in byc.items():
            if len(forms) == 1:
                missing = "si" if "conventional" in forms else "conventional"
                for m in sorted(set(pair_mods[(analyte, cyc)])):
                    incomplete[m].append({"analyte": analyte, "cycle": cyc,
                                          "present": sorted(forms)[0],
                                          "missing": missing})
    nonstd_n = mat["C7_analyte_dup_n"].tolist()
    mat["C7_analyte_dup_hit"] = mat["C7_analyte_dup_hit"] | mat["module"].isin(incomplete)
    mat["C7_analyte_dup_n"] = [int(a) + len(incomplete.get(m, []))
                               for m, a in zip(mat["module"], nonstd_n)]

    # ── 输出矩阵（23 列，与检查清单一致）───────────────────────────────
    cols = ["module", "n_rows", "n_cols", "size_mb", "sec", "error"]
    for k in CHECK_KEYS:
        cols += [k + "_hit", k + "_n"]
    mat = mat[cols]
    mat.to_csv(os.path.join(out_dir, "A_defect_matrix.csv"), index=False, encoding="utf-8-sig")

    with open(os.path.join(out_dir, "A_detail.jsonl"), "w", encoding="utf-8") as f:
        for d in details:
            # 附加跨模块判定结果
            d["cross_module"] = {
                "C4_drift": d["module"] in drift_modules,
                "C10_drift": drift_c10.get(d["module"], []),
                "C7_incomplete_pairs": incomplete.get(d["module"], []),
            }
            f.write(json.dumps(d, ensure_ascii=False, default=str) + "\n")

    # ── 汇总 ─────────────────────────────────────────────────────────────
    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "n_modules": int(len(mat)),
        "n_files_in_cache": int(len(all_files)),
        "n_duplicate_copies_excluded": int(n_dup),
        "n_failed": int((mat["error"].astype(str) != "").sum()),
        "n_total_rows": int(pd.to_numeric(mat["n_rows"], errors="coerce").fillna(0).sum()),
        "per_check": {},
    }
    for k in CHECK_KEYS:
        summary["per_check"][k] = {
            "modules_hit": int(mat[k + "_hit"].sum()),
            "rate_pct": round(float(mat[k + "_hit"].mean() * 100), 2),
            "total_count": int(mat[k + "_n"].sum()),
        }

    # 联合统计：任一类命中 / 值层（C1|C2）命中
    any_hit = np.zeros(len(mat), dtype=bool)
    for k in CHECK_KEYS:
        any_hit |= mat[k + "_hit"].to_numpy()
    value_layer = mat["C1_missing_code_hit"].to_numpy() | mat["C2_parser_encoding_hit"].to_numpy()
    summary["union"] = {
        "modules_any_check": int(any_hit.sum()),
        "rate_any_check_pct": round(float(any_hit.mean() * 100), 2),
        "modules_value_layer_C1_or_C2": int(value_layer.sum()),
        "rate_value_layer_pct": round(float(value_layer.mean() * 100), 2),
        "modules_C1_and_C2": int((mat["C1_missing_code_hit"] & mat["C2_parser_encoding_hit"]).sum()),
        "modules_only_C1": int((mat["C1_missing_code_hit"] & ~mat["C2_parser_encoding_hit"]).sum()),
        "modules_only_C2": int((~mat["C1_missing_code_hit"] & mat["C2_parser_encoding_hit"]).sum()),
    }
    # C2 仅连续测量列的合计（与 C1 同域，便于正文并列引用）
    c2_meas = sum(int(d.get("checks", {}).get("C2_parser_encoding", {}).get("n_measurement", 0))
                  for d in details)
    summary["C2_special_missing_in_measurement_columns"] = int(c2_meas)

    # 锚点：直接重读 DEMO_E/F/G 复算手稿 3.2 节的 1,464 / 30,442
    age_anchor = {"cycles": ["DEMO_E.xpt", "DEMO_F.xpt", "DEMO_G.xpt"],
                  "n": 0, "special_missing": 0, "sentinel_77": 0, "pct": None, "lost": 0}
    file_index = {f.lower(): f for f in os.listdir(data_dir)}
    for fn in age_anchor["cycles"]:
        real = file_index.get(fn.lower())
        if real is None:
            continue
        dfa, _ = read_pandas(os.path.join(data_dir, real))
        s = dfa["RIDAGEYR"]
        age_anchor["n"] += int(len(dfa))
        age_anchor["special_missing"] += int(((s.abs() < 1e-70) & (s != 0) & s.notna()).sum())
        age_anchor["sentinel_77"] += int(s.isin([77, 88, 99]).sum())
        del dfa
        gc.collect()
    age_anchor["lost"] = age_anchor["special_missing"] + age_anchor["sentinel_77"]
    if age_anchor["n"]:
        age_anchor["pct"] = round(100 * age_anchor["lost"] / age_anchor["n"], 2)
    summary["anchor_age_loss"] = age_anchor

    weight_anchor = {"modules_with_subsample_weight": 0,
                     "weight_special_missing_total": 0,
                     "by_weight": {}}
    wt_by = Counter()
    for d in details:
        ev = d.get("checks", {}).get("C8_weight_mismatch", {}).get("evidence", [])
        if ev and isinstance(ev[0], dict):
            if ev[0].get("subsample_weights"):
                weight_anchor["modules_with_subsample_weight"] += 1
            for k, v in (ev[0].get("weight_special_missing") or {}).items():
                wt_by[k] += int(v)
    weight_anchor["by_weight"] = dict(wt_by.most_common(30))
    weight_anchor["weight_special_missing_total"] = int(sum(wt_by.values()))
    summary["anchor_weight_special_missing"] = weight_anchor

    # 证据汇总：C1/C2 的高频变量（供正文举例）
    c1_cols, c1_codes, c2_cols = Counter(), Counter(), Counter()
    for d in details:
        for e in d.get("checks", {}).get("C1_missing_code", {}).get("evidence", []):
            if isinstance(e, dict):
                c1_cols[e["col"]] += e["n"]
                for code, n in (e.get("codes") or {}).items():
                    c1_codes[(e["col"], int(code))] += int(n)
        for e in d.get("checks", {}).get("C2_parser_encoding", {}).get("evidence", []):
            if isinstance(e, dict):
                c2_cols[e["col"]] += e.get("n_special_missing", 0)
    summary["top_columns"] = {
        "C1_by_column": [{"col": c, "n": n} for c, n in c1_cols.most_common(20)],
        "C1_by_column_and_code": [{"col": c, "code": k, "n": n}
                                  for (c, k), n in c1_codes.most_common(20)],
        "C2_by_column": [{"col": c, "n": n} for c, n in c2_cols.most_common(20)],
    }

    summary["cross_module_notes"] = {
        "C4_modules_without_seqn": no_key,
        "C10_columns_checked": c10_checked,
        "C10_modules_flagged": len(drift_c10),
        "C7_modules_with_incomplete_lbx_lbd_pair": len(incomplete),
        "C7_modules_with_standard_dual_unit_pair": int(sum(
            1 for d in details
            if d.get("checks", {}).get("C7_analyte_dup", {}).get("n_standard_dual_unit", 0) > 0)),
        "C3_documented_event_level_modules": int(sum(
            1 for d in details
            if (d.get("checks", {}).get("C3_row_inflation", {}).get("evidence")
                and isinstance(d["checks"]["C3_row_inflation"]["evidence"][0], dict)
                and d["checks"]["C3_row_inflation"]["evidence"][0].get("documented_event_level")))),
        "C2_modules_second_reader_failed": int(sum(
            1 for d in details if d.get("second_reader_ok") is False)),
    }
    summary["adjudication"] = {
        "C1_missing_code": {
            "tier": "observed_value_corruption",
            "definition": "仅连续测量列（BMX/BPX/LBX/LBD/URX/DR1T/WT/INDF + 年龄），"
                          "且哨兵码落在该列非哨兵分布的 p1–p99 之间，才计为碰撞",
            "false_positive_control": "问卷编码前缀（DIQ/SMQ/PAQ/RXQ 等）整列排除；"
                                      "77=拒答、99=不知道这类码因落在分布之外而自动不计",
        },
        "C2_parser_encoding": {
            "tier": "observed_value_corruption",
            "definition": "SAS 特殊缺失被读取器物化成 5.3976e-79（pandas）或 0.0（pyreadstat）",
            "false_positive_control": "以两读取器逐单元格比对，只统计特殊缺失位置",
        },
        "C3_row_inflation": {
            "tier": "structural_hazard",
            "definition": "行数/唯一 SEQN > 1.05；命中为事件级文件，不是数据错误",
            "false_positive_control": "逐模块核对是否为一对多事件表",
        },
        "C4_id_drift": {
            "tier": "observed_value_corruption",
            "definition": "SEQN 的 dtype/小数/位数跨模块不一致",
            "false_positive_control": "全库含 SEQN 的模块逐一比对，无 SEQN 的格式文件单列登记",
        },
        "C7_analyte_dup": {
            "tier": "naming_hazard",
            "definition": "跨周期 LBX/LBD 常规单位列与 SI 单位列曾成对出现、本周期只剩一半"
                          "（格式切换导致变量消失的风险）",
            "false_positive_control": "从未成对出现过的分析物不判；规范配对单列统计不并入 hit",
        },
        "C8_weight_mismatch": {
            "tier": "structural_hazard",
            "definition": "同一模块同时含全样本与子样本权重",
            "false_positive_control": "另记子样本权重的特殊缺失数，避免把'权重存在'当成缺陷",
        },
        "C9_shortname_collision": {
            "tier": "observed_value_corruption",
            "definition": "短名近空而长名近满，且两列在交集上数值相同（确为同一变量）",
            "false_positive_control": "LBXTO5/LBXTO5IN 这类解释码被显式排除",
        },
        "C10_unit_drift": {
            "tier": "naming_hazard",
            "definition": "同名变量跨 ≥3 个周期的跨度或中位数偏离中位值 10 倍以上",
            "false_positive_control": "要求跨度与中位数同时漂移，且至少 3 个周期",
        },
        "not_covered": {
            "C5_threshold_unit": "需指南元数据，无法从数据本身判定",
            "C6_constant_unconverted": "需外部阈值定义",
            "C11_denominator": "验证脚本层，由任务 B 覆盖",
            "C12_self_referential": "验证脚本层，由任务 B 覆盖",
        },
    }
    json.dump(summary, open(os.path.join(out_dir, "A_summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2, default=str)

    # ── 判读表 ───────────────────────────────────────────────────────────
    adjud = []
    for k in CHECK_KEYS:
        a = summary["adjudication"][k]
        adjud.append({
            "check": k, "tier": a["tier"], "definition": a["definition"],
            "false_positive_control": a["false_positive_control"],
            "modules_hit": summary["per_check"][k]["modules_hit"],
            "rate_pct": summary["per_check"][k]["rate_pct"],
            "total_count": summary["per_check"][k]["total_count"],
        })
    pd.DataFrame(adjud).to_csv(os.path.join(out_dir, "A_adjudication.csv"),
                               index=False, encoding="utf-8-sig")

    log(logp, f"done: {len(mat)} modules, failed={summary['n_failed']}, "
              f"elapsed={time.time()-t_start:.0f}s")
    for k, v in summary["per_check"].items():
        log(logp, f"  {k}: {v['modules_hit']} modules ({v['rate_pct']}%), "
                  f"count={v['total_count']}")
    log(logp, f"  anchor age loss: {age_anchor['lost']}/{age_anchor['n']} "
              f"({age_anchor['pct']}%) = {age_anchor['special_missing']} special-missing "
              f"+ {age_anchor['sentinel_77']} age-77")
    logp.close()
    print("\n产出目录:", out_dir)


if __name__ == "__main__":
    main()
