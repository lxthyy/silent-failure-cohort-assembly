# -*- coding: utf-8 -*-
"""S4 synthetic demonstration.

Generates a schema-identical synthetic dataset (same column names, same types,
same rough sizes), then runs the full check suite of
`S4_retrospective_application.py` on it with zero manual steps.

Purpose: show that the checks are *executable* without access to the
institutional data, which cannot be shared. The numbers reported in the paper
come from the institutional data, not from this demo.

Usage:
    python S4_synthetic_demo.py
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "_synthetic")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261030)

# ---------------------------------------------------------------- D1
n1 = 4_000
ids = rng.integers(1_000_000, 9_999_999, n1)
names = np.array(["受试者%05d" % i for i in range(n1)])
ages = rng.integers(18, 90, n1)
sex = rng.choice(["男", "女"], n1)
d1 = pd.DataFrame({
    "住院号": ids,
    "姓名": names,
    "年龄": ages,
    "性别": sex,
    "采样日期": pd.to_datetime("2023-01-01") + pd.to_timedelta(rng.integers(0, 900, n1), unit="D"),
    "诊断": rng.choice(["抑郁状态", "焦虑状态", "其他"], n1, p=[0.5, 0.3, 0.2]),
    "血红蛋白": np.round(rng.normal(135, 15, n1), 1),
})
# injected defect: two diagnosis cells carry a timestamp instead of a label
d1.loc[[3, 17], "诊断"] = "2023/01/01 08:00:00"
d1.loc[5, "年龄"] = 131          # implausible age

# ---------------------------------------------------------------- D2
n2 = 900
d2 = pd.DataFrame({
    "样本号": ["126%09d" % i for i in range(n2)],
    "住院号": rng.integers(1_000_000, 9_999_999, n2),
    "年龄": rng.integers(18, 95, n2),
    "性别": rng.choice(["男", "女"], n2),
    "采样日期": pd.to_datetime("2022-06-01") + pd.to_timedelta(rng.integers(0, 700, n2), unit="D"),
    "诊断": rng.choice(["缺血性卒中", "出血性卒中", "对照"], n2, p=[0.45, 0.25, 0.30]),
    "白细胞计数": np.round(rng.gamma(2.0, 3.2, n2), 2),
})
# injected defect: three impossible white-cell counts
d2.loc[[10, 20, 30], "白细胞计数"] = [5_081.7, 2_400.0, 900.0]

# ---------------------------------------------------------------- D3
n3 = 3_000
d3 = pd.DataFrame({
    "样本编号": ["S%d" % i for i in range(n3)],
    "TSH": np.round(rng.lognormal(0.8, 0.5, n3), 2),
    "年龄": rng.integers(20, 80, n3),
    "性别": rng.choice(["男", "女"], n3),
    "采样日期": pd.to_datetime("2022-01-01") + pd.to_timedelta(rng.integers(0, 1_400, n3), unit="D"),
})

# ---------------------------------------------------------------- B (replication)
nb = int(n1 * 0.6)
idx = rng.choice(n1, nb, replace=False)
dfb = pd.DataFrame({"姓名": names[idx], "年龄": ages[idx], "性别": sex[idx]})

paths = {
    "S4_D1": os.path.join(OUT, "D1_synthetic.csv"),
    "S4_D2": os.path.join(OUT, "D2_synthetic.csv"),
    "S4_D3": os.path.join(OUT, "D3_synthetic.xlsx"),
    "S4_B": os.path.join(OUT, "B_replication_synthetic.xlsx"),
    "S4_OUT": os.path.join(OUT, "out"),
}
d1.to_csv(paths["S4_D1"], index=False, encoding="utf-8-sig")
d2.to_csv(paths["S4_D2"], index=False, encoding="utf-8-sig")
d3.to_excel(paths["S4_D3"], index=False)          # S4 reads D3 with read_excel
dfb.to_excel(paths["S4_B"], index=False)
os.makedirs(paths["S4_OUT"], exist_ok=True)

print("[synthetic] D1 %s rows, D2 %s, D3 %s, B %s" % (len(d1), len(d2), len(d3), len(dfb)))

env = dict(os.environ, **paths)
script = os.path.join(HERE, "S4_retrospective_application.py")
r = subprocess.run([sys.executable, script], env=env, capture_output=True, text=True)
print(r.stdout[-3000:])
if r.returncode != 0:
    print(r.stderr, file=sys.stderr)
    raise SystemExit("S4 failed on synthetic data")

csv = os.path.join(paths["S4_OUT"], "S4_retrospective_application.csv")
df = pd.read_csv(csv, encoding="utf-8-sig")
print("\n[synthetic] output rows:", len(df), "| datasets:",
      sorted(df["dataset"].unique()), "| checks:", sorted(df["check"].unique()))
print("[synthetic] OK — the check suite runs end to end with no manual steps")
