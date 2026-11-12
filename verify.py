#!/usr/bin/env python3
"""verify.py - check that the shipped results carry every headline number in the manuscript.

Reads results/ only. No survey data and no institutional data required.
Exits non-zero if any check fails.
"""
import csv
import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "results")
rows, fails = [], []


def check(name, got, want):
    ok = got == want
    rows.append((name, want, got, "ok" if ok else "FAIL"))
    if not ok:
        fails.append(name)


def read_csv(name):
    with open(os.path.join(R, name), encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


# ---- 1. census of 286 survey modules (Table 4, Section 4.9) ----
j = json.load(open(os.path.join(R, "A_summary_rerun_20261112.json"), encoding="utf-8"))
check("modules analysed", j["n_modules"], 286)
check("rows analysed", j["n_total_rows"], 2279734)
check("class 1 modules", j["per_check"]["C1_missing_code"]["modules_hit"], 133)
check("class 1 legitimate values", j["per_check"]["C1_missing_code"]["total_count"], 45637)
check("class 2 modules", j["per_check"]["C2_parser_encoding"]["modules_hit"], 183)
check("class 2 cells", j["per_check"]["C2_parser_encoding"]["total_count"], 2077656)
check("modules with any hazard", j["union"]["modules_any_check"], 225)
check("share with any hazard (%)", round(j["union"]["rate_any_check_pct"], 1), 78.7)

a = j["anchor_age_loss"]
check("ages analysed", a["n"], 30442)
check("ages lost", a["lost"], 1464)
check("  of which SAS special missing", a["special_missing"], 1303)
check("  of which sentinel 77", a["sentinel_77"], 161)
check("percent lost", round(a["pct"], 1), 4.8)

# ---- 2. detection matrix under injection (Table 3, Section 4.10) ----
inj = read_csv("B_headtohead.csv")
check("injection rows", len(inj), 600)
c = Counter((r["defect"], r["method"], r["detected"]) for r in inj)
EXPECT = {}
for d in ("D1_sentinel", "D2_unit_rescale", "D3_row_dup", "D4_col_shift"):
    EXPECT[(d, "M1_double_programming")] = 0
    EXPECT[(d, "M2_reproducibility")] = 0
    EXPECT[(d, "M3_kahn_checklist")] = 30
    EXPECT[(d, "M4_four_rule_protocol")] = 30
EXPECT[("D5_round_int", "M1_double_programming")] = 0
EXPECT[("D5_round_int", "M2_reproducibility")] = 0
EXPECT[("D5_round_int", "M3_kahn_checklist")] = 7
EXPECT[("D5_round_int", "M4_four_rule_protocol")] = 30
for (d, m), want in EXPECT.items():
    check("%-20s %-24s detected" % (d, m), c[(d, m, "True")], want)

# ---- 3. boundary row: same-unit substitution, invisible to all four ----
b = read_csv("S3c_boundary.csv")
check("boundary replicates", len(b), 150)
check("boundary checks fired", sum(1 for x in b if x["fired"].strip().lower() == "true"), 0)

# ---- 4. 64-variable validation (Section 4.4) ----
v = read_csv("FULL64_validation_results.csv")
check("variables compared", len(v), 64)
check("variables with median difference > 0.01", sum(
    1 for x in v if x["Median_Diff"].strip() not in ("", "NA") and float(x["Median_Diff"]) > 0.01), 0)

# ---- 5. rounding denominators (Figure 3, Section 4.5) ----
s = read_csv("S1_A_rounding.csv")
cre = [x for x in s if x["label"] == "Creatinine"]
check("creatinine row present", len(cre), 1)
if cre:
    check("creatinine measurements", int(cre[0]["n"]), 19212)

# ---- report ----
w = max(len(r[0]) for r in rows)
print("=" * (w + 22))
print("%-*s  %-10s %-10s %s" % (w, "check", "expected", "found", ""))
print("=" * (w + 22))
for name, want, got, status in rows:
    print("%-*s  %-10s %-10s %s" % (w, name, want, got, status))
print("=" * (w + 22))
print("%d checks, %d failed" % (len(rows), len(fails)))
sys.exit(1 if fails else 0)
