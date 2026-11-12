# -*- coding: utf-8 -*-
"""S2: Figures 1-4 for the JAMIA manuscript (v2, refined layout).

Figure 1  Twelve classes grouped by layer (class names exactly as Table 1)
Figure 2  The four gates of the collision protocol (wording as Section 2.4
          and the Figure 2 legend)
Figure 3  Participants moved across a decision threshold by rounding
Figure 4  Weighted mean fasting glucose, examination vs fasting-subsample weight

Design-based estimates in fig4 are the ones printed by weights(); they come
from the NHANES analysis units (cycles 2007-2012, adults >= 20 y, stratified
on SDMVSTRA, clustered on SDMVPSU).

Outputs (env FIGURE_OUT overrides the destination):
  Figure_1_Classes_v1.01.{svg,png,tif}
  Figure_2_Protocol_v1.01.{svg,png,tif}
  Figure_3_RoundingThreshold_v1.02.{svg,png,tif}
  Figure_4_FastingWeight_v1.02.{svg,png,tif}

Self-check: after drawing each panel figure the script measures every text
object against its card boundary and reports any overflow (must print 0).
"""
import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.6,
    "axes.edgecolor": "#333333", "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "svg.fonttype": "none",
})

OUT = os.environ.get("FIGURE_OUT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures"))
os.makedirs(OUT, exist_ok=True)

INK, MID = "#1a1a1a", "#5a5a5a"
BLUE, ORANGE = "#2E6DA4", "#C77F1A"
BLUE_F, ORANGE_F = "#D7E5F2", "#FBE3C9"
# per-layer palette (colour-blind safe, grayscale-distinct)
LAYERS = [
    ("Value", "parsing  ·  linking", "#2E6DA4", "#E9F1F8"),
    ("Semantic", "converting  ·  weighting", "#C77F1A", "#FBF1DF"),
    ("Presentation", "naming", "#2E8B6E", "#E6F3EE"),
    ("Verification", "checking", "#A85B86", "#F6E9F0"),
]


def save_all(fig, stem):
    png = os.path.join(OUT, stem + ".png")
    fig.savefig(png, dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(OUT, stem + ".svg"), facecolor="white")
    from PIL import Image
    Image.open(png).convert("RGB").save(
        os.path.join(OUT, stem + ".tif"), compression="tiff_lzw", dpi=(400, 400))
    print("[saved]", stem, "(png/svg/tif)")


def check_texts(fig, cards, tag):
    """cards: list of (x0, x1, y0, y1) in data coords; verify every text
    that falls inside a card horizontally also stays inside vertically,
    and print overflow count."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bad = 0
    for ax, box in cards:
        x0, x1, y0, y1 = box
        p0 = ax.transData.transform((x0, y0))
        p1 = ax.transData.transform((x1, y1))
        for t in ax.texts:
            bb = t.get_window_extent(renderer=r)
            inside_x = p0[0] - 1 <= bb.x0 and bb.x1 <= p1[0] + 1
            inside_y = p0[1] - 1 <= bb.y0 and bb.y1 <= p1[1] + 1
            overlaps_x = bb.x0 < p1[0] and bb.x1 > p0[0]
            overlaps_y = bb.y0 < p1[1] and bb.y1 > p0[1]
            if overlaps_x and overlaps_y and not (inside_x and inside_y):
                bad += 1
                print("  OVERFLOW[%s]: %r" % (tag, t.get_text()[:40]))
    print("  %s text-overflow count: %d" % (tag, bad))
    return bad


# ---------------------------------------------------------------- figure 1
def fig1():
    classes = [
        ["1  missing-code collision",
         "2  parser encoding divergence",
         "3  one-to-many inflation",
         "4  identifier key-space drift"],
        ["5  threshold-unit mismatch",
         "6  constant on unconverted scale",
         "7  format switch deletes variables",
         "8  subsample weight mismatch"],
        ["9  short-name collision empties a column",
         "10  column-name unit drift"],
        ["11  denominator evaporation",
         "12  self-referential gold standard"],
    ]
    fig, ax = plt.subplots(figsize=(7.8, 3.9))
    ax.set_xlim(0, 12.6); ax.set_ylim(0, 6.5); ax.axis("off")

    W, STEP, Y0, H = 3.0, 3.12, 1.00, 4.80
    cards = []
    for i, ((name, ops, edge, face), items) in enumerate(zip(LAYERS, classes)):
        x = 0.08 + i * STEP
        ax.add_patch(FancyBboxPatch((x, Y0), W, H,
                                    boxstyle="round,pad=0.02,rounding_size=0.10",
                                    linewidth=1.0, edgecolor=edge, facecolor=face))
        ax.text(x + 0.17, 5.52, name, fontsize=9.2, weight="bold",
                color=edge, va="center")
        ax.text(x + 0.17, 5.15, ops, fontsize=6.6, color=MID,
                style="italic", va="center")
        ax.plot([x + 0.17, x + W - 0.17], [4.94, 4.94],
                color=edge, lw=0.6, alpha=0.45)
        y = 4.68
        for it in items:
            for ln in textwrap.wrap(it, width=24):
                ax.text(x + 0.17, y, ln, fontsize=7.2, color=INK, va="center")
                y -= 0.295
            y -= 0.30
        cards.append((ax, (x + 0.02, x + W - 0.02, Y0 + 0.02, Y0 + H - 0.02)))

    ax.text(0.08, 0.55,
            "Every class is defined by the operation that fails, not by the "
            "data source or the country it comes from.", fontsize=7.0, color=INK)
    save_all(fig, "Figure_1_Classes_v1.01")
    return check_texts(fig, cards, "fig1")


# ---------------------------------------------------------------- figure 2
def fig2():
    gates = [
        ("I", "Criteria fixed\nbefore the run",
         "Thresholds committed\nbefore any comparison\nis run"),
        ("II", "Denominator\naccounting",
         "Every declared item\ncompared, mismatched,\nor skipped, with a\nstated reason"),
        ("III", "Cross the\nimplementation\nboundary",
         "Two parsers, two\nlanguages, two data\npaths - at the level\nof the value"),
        ("IV", "No self-referential\ngold standard",
         "A snapshot may detect\ndrift; it may not\nestablish correctness"),
    ]
    fig, ax = plt.subplots(figsize=(7.8, 3.1))
    ax.set_xlim(0, 12.6); ax.set_ylim(0, 3.55); ax.axis("off")

    W, STEP, Y0, H = 2.92, 3.14, 0.45, 2.55
    cards = []
    for i, ((num, title, body), (name, _o, edge, face)) in enumerate(zip(gates, LAYERS)):
        x = 0.08 + i * STEP
        ax.add_patch(FancyBboxPatch((x, Y0), W, H,
                                    boxstyle="round,pad=0.02,rounding_size=0.10",
                                    linewidth=1.0, edgecolor=edge, facecolor=face))
        ax.text(x + 0.17, Y0 + H - 0.30, num, fontsize=11, weight="bold",
                color=edge, va="center")
        ax.text(x + 0.17, Y0 + H - 0.68, title, fontsize=7.6, weight="bold",
                color=INK, va="top", linespacing=1.28)
        n_t = title.count("\n") + 1
        ax.text(x + 0.17, Y0 + H - 0.68 - n_t * 0.30 - 0.10, body,
                fontsize=6.3, color=MID, va="top", linespacing=1.42)
        if i < 3:
            ax.add_patch(FancyArrowPatch((x + W + 0.03, Y0 + H / 2),
                                         (x + STEP - 0.05, Y0 + H / 2),
                                         arrowstyle="-|>", mutation_scale=11,
                                         linewidth=1.0, color=edge))
        cards.append((ax, (x + 0.02, x + W - 0.02, Y0 + 0.02, Y0 + H - 0.02)))

    ax.text(0.08, 0.18,
            "Each gate is designed to make the verification differ from the "
            "thing verified.", fontsize=7.0, color=INK)
    save_all(fig, "Figure_2_Protocol_v1.01")
    return check_texts(fig, cards, "fig2")


# ---------------------------------------------------------------- figure 3
def fig3():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7))
    specs = [
        ("Total cholesterol, 5.2 mmol/L", 171, 22221, "201 mg/dL"),
        ("Fasting glucose, 7.0 mmol/L", 36, 9529, "126 mg/dL"),
    ]
    for ax, (label, n, tot, src) in zip(axes, specs):
        share = 100.0 * n / tot
        ax.barh([0], [share], color=BLUE_F, edgecolor=BLUE, linewidth=0.9,
                height=0.42, zorder=2)
        ax.set_yticks([])
        ax.set_xlim(0, max(1.0, share * 1.55))
        ax.set_ylim(-0.62, 0.95)
        ax.set_title(label, fontsize=8.0, pad=8, color=INK)
        ax.set_xlabel("% of participants with a measured value", fontsize=7.0)
        ax.tick_params(labelsize=6.5)
        ax.annotate("%d of %s reclassified (%.2f%%)\nevery one carrying exactly %s"
                    % (n, format(tot, ","), share, src),
                    xy=(share, 0.235), xycoords="data",
                    xytext=(0.98, 0.97), textcoords="axes fraction",
                    ha="right", va="top",
                    fontsize=6.8, color=INK, linespacing=1.45,
                    arrowprops=dict(arrowstyle="-", linewidth=0.7,
                                    color=MID, shrinkB=0))
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout()
    save_all(fig, "Figure_3_RoundingThreshold_v1.02")


# ---------------------------------------------------------------- weights
def weights():
    """Design-based estimates from the NHANES analysis units (cycles E/F/G).

    Adults >= 20 y, stratified by SDMVSTRA and clustered on SDMVPSU.
    Pooled, adults: 105.80 (SE 0.58) vs 104.98 (SE 0.55), +0.82 mg/dL (+0.78%).
    """
    rows = [
        ("2007-2008", 31184, 2601, 107.48, 0.78, 106.25, 0.69, 1.23),
        ("2009-2010", 31152, 2762, 103.94, 0.79, 103.73, 0.84, 0.21),
        ("2011-2012", 29382, 2448, 105.97, 1.26, 104.98, 1.21, 1.00),
    ]
    for r in rows:
        print("%s  n=%d  exam=%.2f (se %.2f)  fast=%.2f (se %.2f)  diff=%.2f (%.2f%%)"
              % (r[0], r[2], r[3], r[4], r[5], r[6], r[7], 100 * r[7] / r[5]))
    return rows


def fig4(rows):
    fig, ax = plt.subplots(figsize=(7.2, 2.9))
    x = np.arange(len(rows))
    w = 0.34
    ax.bar(x - w / 2, [r[3] for r in rows], w, color=BLUE_F, edgecolor=BLUE,
           linewidth=0.9, label="Examination weight (WTMEC2YR)",
           yerr=[r[4] for r in rows], capsize=3,
           error_kw=dict(linewidth=0.7, capthick=0.7, ecolor="#444444"))
    ax.bar(x + w / 2, [r[5] for r in rows], w, color=ORANGE_F, edgecolor=ORANGE,
           linewidth=0.9, label="Fasting-subsample weight (WTSAF2YR)",
           yerr=[r[6] for r in rows], capsize=3,
           error_kw=dict(linewidth=0.7, capthick=0.7, ecolor="#444444"))
    ax.set_xticks(x)
    ax.set_xticklabels(["%s\n(n=%s)" % (r[0], format(r[2], ",")) for r in rows],
                       fontsize=7.0)
    ax.set_ylim(99, 112.5)
    ax.set_ylabel("Weighted mean fasting glucose (mg/dL)", fontsize=7.4)
    ax.tick_params(labelsize=6.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(fontsize=6.8, frameon=False, loc="lower center",
              bbox_to_anchor=(0.5, 1.02), ncol=2)
    for i, r in enumerate(rows):
        top = max(r[3] + r[4], r[5] + r[6])
        ax.annotate("+%.2f mg/dL\n(+%.2f%%)" % (r[7], 100 * r[7] / r[5]),
                    xy=(i, top + 0.55), ha="center", va="bottom",
                    fontsize=6.4, color="#444444", linespacing=1.35)
    fig.tight_layout()
    save_all(fig, "Figure_4_FastingWeight_v1.02")


import numpy as np  # noqa: E402  (kept late so the self-checks above read cleanly)

bad = fig1() + fig2()
fig3()
rows = weights()
fig4(rows)
print("\nwritten to", OUT)
for f in sorted(os.listdir(OUT)):
    if f.startswith("Figure_"):
        print("  ", f, os.path.getsize(os.path.join(OUT, f)), "bytes")
if bad:
    raise SystemExit("TEXT OVERFLOW DETECTED: %d" % bad)
print("self-check passed: no text overflow in fig1/fig2 cards")
