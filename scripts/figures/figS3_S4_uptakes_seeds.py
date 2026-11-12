"""Figs. S3 and S4 of the Supplement, in the figure conventions of the paper (svg_finish, svg_panels).

fig. S3 (FigS3_uptakes) attached against detached uptakes, per network (Kimonu, iterated route): kinetic modules (A), ACR
                metabolites (B) and independent ACRR ratios (C, k - 1 per group of k) among free metabolites. figdata.json, arms "kimonu" (attached,
                timing arm att) and "kimonu_det" (detached, timing arm det).
fig. S4 (FigS4_seeds) relations that change with the random seed, per network and arm: coupled complex pairs (A) and ACR
                metabolites plus independent ACRR ratios among free metabolites (B) that differ between the seeds (the
                rank of what some seed finds minus the rank of what all seeds share, seed_spread.py).
                seed_spread.json (seed_spread.py); Kimonu's r* seeds 0, 1, 2 change nothing.
Both JSON files are read from figures/source_data/ of the records (records.py).

    python figures/figS3_S4_uptakes_seeds.py [<out_dir>]      # from scripts/; default drafts/supp of KIMONU_OUT
"""
import io, json, os, sys, tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
sys.path[:0] = [os.path.dirname(os.path.abspath(__file__)), os.path.dirname(os.path.dirname(os.path.abspath(__file__)))]
import svg_finish, svg_panels
import records

plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
                     "ytick.major.size": 2.5, **svg_finish.mathtext_rc()})
INK, MUTE, RULE = "#1b2430", "#6b7480", "#d5d9de"
GREEN, BLUE, RED = "#4DAF4A", "#377EB8", "#E41A1C"
LETTER = dict(fontsize=10, fontweight="bold", family="Arial", color=INK, va="top", ha="left")   # panel letters 10 pt bold (Science)
W = 7.2
D = json.load(open(records.read(records.SOURCE_DATA, "figdata.json"))); nC = D["n_complexes"]
MODELS = sorted(D["kimonu_det"], key=lambda m: nC[m])
from names import network as SHORT                 # one display name per network (figures/names.py)


def save(fig, stem, panels):
    buf = io.StringIO(); fig.savefig(buf, format="svg"); raw = buf.getvalue()
    n, box = svg_finish.check(raw, raw, tempfile.mkdtemp())
    assert n == 0, f"{stem}: removing the clipping changes {n} pixels in {box}"
    fg = svg_panels.Figure(raw); miss = fg.assign(panels); assert not miss, len(miss)
    open(stem + ".svg", "w", encoding="utf-8").write(svg_finish.finish(fg.svg())[0])
    fig.savefig(stem + ".png", dpi=200)                      # for checking only; the vault keeps the SVG


def strip(ax, ylabel, title, names=True):
    ax.set_xlim(-0.6, len(MODELS) - 0.4); ax.set_xticks(range(len(MODELS)))
    ax.set_xticklabels([SHORT(m) for m in MODELS] if names else [], rotation=90, fontsize=6)
    ax.tick_params(axis="x", pad=2)                     # a tick at every network, the major size of the y axis
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.yaxis.grid(True, color=RULE, lw=0.4, zorder=0); ax.set_axisbelow(True)
    ax.set_ylabel(ylabel, labelpad=2); ax.set_title(title, fontsize=8, color=INK, loc="center", pad=4)


def uptakes(out):
    H = 5.3; fig = plt.figure(figsize=(W, H))
    keys = [("n_kinetic", "kinetic modules", (1, 3000)), ("acr_free", "ACR metabolites (free)", (0.5, 300)),
            ("acrr_ratios_free", "ACRR ratios (free)", (0.5, 3000))]
    hA, gap, top = 1.15, 0.4, 0.35
    for k, (key, title, ylim) in enumerate(keys):
        y0 = H - top - (k + 1) * hA - k * gap
        ax = fig.add_axes([0.55 / W, y0 / H, (W - 0.65) / W, hA / H])
        for i, m in enumerate(MODELS):
            a, d = D["kimonu"][m][key], D["kimonu_det"][m][key]
            ax.plot([i, i], [min(a, d), max(a, d)], color="#aab2bb", lw=0.9, solid_capstyle="butt", zorder=2)
            ax.scatter(i, d, s=14, facecolor=GREEN, edgecolor="white", linewidths=0.5, zorder=3)
            ax.scatter(i, a, s=14, facecolor="white", edgecolor=GREEN, linewidths=0.9, zorder=3)
        ax.set_yscale("symlog", linthresh=1); ax.set_ylim(-0.3, ylim[1]); ax.minorticks_off()
        ax.set_yticks([0, 1, 10, 100, 1000][: 5 if ylim[1] >= 1000 else 4]); ax.set_yticklabels(["0", "1", "10", "100", "1000"][: 5 if ylim[1] >= 1000 else 4])
        strip(ax, "count", title, names=k == 2)
        fig.text(0.0, (y0 + hA + 0.3) / H, "ABC"[k], **LETTER)
    fig.legend(handles=[Line2D([], [], ls="", marker="o", ms=4, mfc=GREEN, mec="white", mew=0.5, label="detached uptakes (default)"),
                        Line2D([], [], ls="", marker="o", ms=4, mfc="white", mec=GREEN, mew=0.9, label="attached uptakes")],
               loc="upper right", bbox_to_anchor=(1.0, 1.0), frameon=False, ncol=2, handletextpad=0.2, columnspacing=1.2)
    Hp = H * 72; b = [(top + (k + 1) * hA + k * gap + 0.1) * 72 for k in range(3)]
    save(fig, os.path.join(out, "FigS3_uptakes"), [("legend", [(300, 0, 518.4, 14)]), ("A", [(0, 0, 518.4, b[0])]),
                                                  ("B", [(0, b[0], 518.4, b[1])]), ("C", [(0, b[1], 518.4, Hp)])])


def seeds(out):
    S = json.load(open(records.read(records.SOURCE_DATA, "seed_spread.json")))
    ARMS = [("cocoa_false", "COCOA.jl, efficient=false", RED, True), ("cocoa_true", "COCOA.jl, efficient=true", RED, False),
            ("ref_sampled", "Reference, sampled (Langary et al., 2025)", BLUE, False), ("kimonu_det", "Kimonu (iterated route)", GREEN, True)]
    H = 3.9; fig = plt.figure(figsize=(W, H)); hA, gap, top = 1.15, 0.4, 0.35
    L = 0.65                                                 # left edge of the axes: room for the "1 000 000" ticks and "count"
    axs = []
    for k, (key, title) in enumerate([("pairs", "coupled complex pairs that change with the seed"),
                                      ("fixed_ratios", "ACR metabolites plus ACRR ratios (free) that change with the seed")]):
        y0 = H - top - (k + 1) * hA - k * gap
        ax = fig.add_axes([L / W, y0 / H, (W - 0.1 - L) / W, hA / H]); axs.append(ax)
        for j, (arm, lab, col, filled) in enumerate(ARMS):
            for i, m in enumerate(MODELS):
                r = S.get(arm, {}).get(m)
                if r is None: continue
                v = r[key][0]
                ax.scatter(i + (j - 1.5) * 0.17, v, s=12, marker="o" if arm.startswith("kimonu") else "^" if arm.startswith("cocoa") else "s",
                           facecolor=col if filled else "white", edgecolor="white" if filled else col, linewidths=0.5 if filled else 0.8, zorder=3)
        top_v = max(r[key][0] for a, *_ in ARMS for r in S.get(a, {}).values())
        ticks = [t for t in (0, 1, 10, 100, 1000, 10000, 100000, 1000000) if t <= top_v * 10]
        ax.set_yscale("symlog", linthresh=1); ax.set_ylim(-0.3, top_v * 3); ax.minorticks_off()
        ax.set_yticks(ticks); ax.set_yticklabels([f"{t:,}".replace(",", " ") for t in ticks])
        strip(ax, "count", title, names=k == 1)
        fig.text(0.0, (y0 + hA + 0.3) / H, "AB"[k], **LETTER)
    fig.canvas.draw()                                        # one x for both "count" labels, that of A (the wider ticks)
    xA = fig.transFigure.inverted().transform((axs[0].yaxis.label.get_position()[0], 0))[0]
    for ax in axs:
        ax.yaxis.set_label_coords(xA, 0.5, transform=blended_transform_factory(fig.transFigure, ax.transAxes))
    hd = [Line2D([], [], ls="", marker="o" if a.startswith("kimonu") else "^" if a.startswith("cocoa") else "s", ms=4,
                 mfc=c if f else "white", mec="white" if f else c, mew=0.5 if f else 0.8, label=l) for a, l, c, f in ARMS]
    fig.legend(handles=hd, loc="upper right", bbox_to_anchor=(1.0, 1.0), frameon=False, ncol=4, handletextpad=0.2, columnspacing=1.0)
    Hp = H * 72; b = [(top + hA + 0.1) * 72]
    save(fig, os.path.join(out, "FigS4_seeds"), [("legend", [(100, 0, 518.4, 14)]), ("A", [(0, 0, 518.4, b[0])]), ("B", [(0, b[0], 518.4, Hp)])])


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else records.write(records.DRAFTS, "supp", folder=True); os.makedirs(out, exist_ok=True)
    uptakes(out); seeds(out)
