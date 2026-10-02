"""Deficiency panel: per model, the number of complexes n, the deficiency δ and the number of concordance classes.

δ = rank(M) is the dimension of all steady-state complex activity patterns, so the non-balanced complexes fall
into at least δ concordance classes (distinct directions of the nonzero rows of M, singletons included). The
panel shows that floor and how close concordance comes to it, on the 33 elementary-step-split benchmark models
and on the 12 networks of Küken et al. (2022) as they preprocessed them (blocked reactions removed).

    python figures/make_deficiency_panel.py [out_dir]     # run from benchmarks/; default figures/main
"""
import json, os, sys
import matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np

out = sys.argv[1] if len(sys.argv) > 1 else "figures/main"
SPLIT = json.load(open("deficiency/deficiency_33.json"))
KUKEN = json.load(open("deficiency/deficiency_kuken12.json"))

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "svg.fonttype": "none"})
INK, MUTE = "#1b2430", "#9aa3ad"
CLS, DEF = "#2f6db3", "#e6862b"                      # concordance classes, deficiency
fig, ax = plt.subplots(figsize=(8.5 / 2.54, 6.4 / 2.54))
lo, hi = 30, 2.5e4
ax.plot([lo, hi], [lo, hi], color=MUTE, lw=0.8, ls="--", zorder=1)
ax.text(1.6e4, 1.05e4, "n", color=MUTE, ha="right", va="bottom")
for rows, mk, ms in ((SPLIT, "o", 14), (KUKEN, "^", 18)):
    n = np.array([r["n"] for r in rows]); d = np.array([r["delta"] for r in rows]); c = np.array([r["directions"] for r in rows])
    for i in range(len(n)):                          # each model's class count stands on its own floor
        ax.plot([n[i], n[i]], [d[i], c[i]], color="#c9d3de", lw=0.6, zorder=2)
    ax.scatter(n, c, s=ms, marker=mk, facecolor=CLS, edgecolor="white", linewidth=0.4, zorder=3)
    ax.scatter(n, d, s=ms, marker=mk, facecolor=DEF, edgecolor="white", linewidth=0.4, zorder=3)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(500, hi); ax.set_ylim(lo, hi)
ax.set_xlabel("complexes n"); ax.set_ylabel("count")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.tick_params(colors=INK, length=2.5)
from matplotlib.lines import Line2D
h = [Line2D([], [], ls="", marker="s", ms=4, mfc=CLS, mec=CLS, label="concordance classes"),
     Line2D([], [], ls="", marker="s", ms=4, mfc=DEF, mec=DEF, label="deficiency δ = rank M"),
     Line2D([], [], ls="", marker="o", ms=3.5, mfc="white", mec=INK, mew=0.6, label="33 benchmark models (split)"),
     Line2D([], [], ls="", marker="^", ms=4, mfc="white", mec=INK, mew=0.6, label="12 models of Küken et al.")]
ax.legend(handles=h, loc="upper left", frameon=False, handletextpad=0.3, borderaxespad=0.2, labelspacing=0.35)
fig.tight_layout(pad=0.3)
os.makedirs(out, exist_ok=True)
for ext in ("svg", "png"):
    fig.savefig(f"{out}/deficiency_panel.{ext}", dpi=300)
