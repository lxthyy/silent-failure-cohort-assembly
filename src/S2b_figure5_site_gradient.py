# -*- coding: utf-8 -*-
"""S2b: Figure 5 - where the implausible glucose values come from (v2 layout).

Reads CHARLS 2011 venous blood, identifies collection sites from the identifier
prefix, and shows (A) the share of low-glucose records that originate in the two
affected sites, and (B) the glucose distribution of those sites against all others.

Reproduces the numbers quoted in Section 3.8:
  <70 mg/dL : 96 of 178 from the two sites (53.9%)
  <60 mg/dL : 73 of  88 (83.0%)
  <50 mg/dL : 38 of  43 (88.4%)
  two sites : n=136, median 58.5 vs 102.6 elsewhere, 1.2% of records

Set CHARLS_DIR to the folder holding Blood_20140429.dta.
Output: Figure_5_SiteGradient_v1.02.{svg,png,tif}
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.6,
    "axes.edgecolor": "#333333", "svg.fonttype": "none",
})
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CHARLS = os.environ.get("CHARLS_DIR",
                        os.path.join(os.path.dirname(os.path.abspath(__file__)), "charls"))
OUT = os.environ.get("FIGURE_OUT", os.path.dirname(os.path.abspath(__file__)))
os.makedirs(OUT, exist_ok=True)

INK, MID = "#1a1a1a", "#5a5a5a"
BLUE, BLUE_F = "#2E6DA4", "#D7E5F2"
ORANGE, ORANGE_F = "#C77F1A", "#FBE3C9"
PLUM = "#A85B86"

d = pd.read_stata(os.path.join(CHARLS, "Blood_20140429.dta"), convert_categoricals=False)
g = pd.to_numeric(d["newglu"], errors="coerce")
site = d["ID"].astype(str).str[:5]

impl = g.notna() & ((g < 40) | (g > 600))
two = site[impl].value_counts().head(2).index.tolist()
in2 = site.isin(two)
m = g.notna()

THR = [70, 60, 50]
shares, totals, hits = [], [], []
for t in THR:
    low = m & (g < t)
    h = low & in2
    shares.append(100 * h.sum() / max(low.sum(), 1))
    totals.append(int(low.sum()))
    hits.append(int(h.sum()))

share_sample = 100 * (m & in2).sum() / m.sum()

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw={"width_ratios": [1.05, 1]})

# ---- panel A
ax = axes[0]
x = np.arange(len(THR))
ax.bar(x, shares, width=0.55, color=BLUE_F, edgecolor=BLUE, linewidth=0.9, zorder=2)
for i, (s, n, h) in enumerate(zip(shares, totals, hits)):
    ax.text(i, s + 2.5, "%.1f%%" % s, ha="center", va="bottom",
            fontsize=7.4, color=INK)
    ax.text(i, s - 7, "%d/%d" % (h, n), ha="center", va="top",
            fontsize=6.4, color=INK, zorder=3)
ax.axhline(share_sample, color=PLUM, lw=0.9, ls=(0, (4, 2)), zorder=1)
ax.text(0.03, 0.97, "dashed line: the two sites hold %.1f%% of all records" % share_sample,
        transform=ax.transAxes, ha="left", va="top", fontsize=6.4, color=MID)
ax.set_xticks(x)
ax.set_xticklabels(["<%d mg/dL" % t for t in THR], fontsize=7.2)
ax.set_ylabel("Share of low values from the two sites (%)", fontsize=7.2)
ax.set_ylim(0, 108)
ax.tick_params(labelsize=6.8, length=2)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---- panel B
ax = axes[1]
a = g[m & in2].values
b = g[m & ~in2].values
bp = ax.boxplot([b, a], widths=0.45, showfliers=False, patch_artist=True,
                medianprops=dict(color=INK, lw=1.1),
                boxprops=dict(lw=0.8),
                whiskerprops=dict(color=MID, lw=0.7),
                capprops=dict(color=MID, lw=0.7))
bp["boxes"][0].set(facecolor=BLUE_F, edgecolor=BLUE)
bp["boxes"][1].set(facecolor=ORANGE_F, edgecolor=ORANGE)
ax.set_xticks([1, 2])
ax.set_xticklabels(["other %s\nsites" % format(len(b), ","),
                    "two affected\nsites %d" % len(a)], fontsize=7.0)
ax.set_ylabel("Fasting glucose (mg/dL)", fontsize=7.2)
ax.text(1, np.median(b) + 10, "%.1f" % np.median(b), ha="center",
        fontsize=6.8, color=INK)
ax.text(2, np.median(a) + 10, "%.1f" % np.median(a), ha="center",
        fontsize=6.8, color=INK)
ax.set_ylim(0, 220)
ax.tick_params(labelsize=6.8, length=2)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.tight_layout()
png = os.path.join(OUT, "Figure_5_SiteGradient_v1.02.png")
fig.savefig(png, dpi=400, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(OUT, "Figure_5_SiteGradient_v1.02.svg"), facecolor="white")
from PIL import Image  # noqa: E402
Image.open(png).convert("RGB").save(
    os.path.join(OUT, "Figure_5_SiteGradient_v1.02.tif"),
    compression="tiff_lzw", dpi=(400, 400))
plt.close(fig)

print("sites:", two)
print("share of sample in the two sites: %.2f%%" % share_sample)
for t, s, n, h in zip(THR, shares, totals, hits):
    print("  <%d: %d/%d = %.1f%%" % (t, h, n, s))
print("medians: other %.2f | two sites %.2f" % (np.median(b), np.median(a)))
print("[saved] Figure_5_SiteGradient_v1.02 (png/svg/tif)")
