"""Figs. S5 and S6, and the draft of Figure 3E, F, from kinetic_variants/float_four.py and float_exactcmp.py: floating point
against exact arithmetic, for both routes and four ways of computing and comparing the nullspace rows.

    exactcmp   exact nullspace (Kimonu's RREF bases), rows kept as fractions, element-wise tolerance, decided exactly
    exact      exact nullspace rounded to an orthonormal float64 basis, row-relative Euclidean tolerance
    svdabs     SVD nullspace (orthonormal), element-wise tolerance
    svd        SVD nullspace (orthonormal), row-relative Euclidean tolerance

Figure3EF_draft     draft of Figure 3E, F (drawn into the figure by fig3EF_float.py): per tolerance, what the comparison of the two routes certifies, exact nullspace
                    (exactcmp) against SVD nullspace (svdabs), one element-wise rule. Classes: the routes agree and are
                    right (the complete set, certified), the routes differ as they do exactly (the gap models), the
                    routes differ only because of floating point, the routes agree on a wrong answer (a false
                    certificate). The first bar is tolerance 0 (bitwise equal rows); the dashed line is what exact
                    arithmetic certifies. Answers are compared as ACR, ACRR and kinetic modules.
FigS5_certificate    fig. S5, the same for all four ways.
FigS6_float_errors   fig. S6, per route and way: the models whose ACR and ACRR come out right, only made up, only missed, or both.

Figs. S5 and S6 have the colour key in one row below the plots, so the plots take the full page width, and one SVG group
per plot (panel_A, panel_B, ... row by row) plus panel_legend.

    python figures/figS5_S6_float.py <out_dir>        # from benchmarks/
"""
import glob, json, os, sys
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import io, tempfile
import svg_finish, svg_panels

KF = os.environ.get("KF_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_float_four")
plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "font.family": "Arial",
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "svg.fonttype": "none",
                     "text.color": "#1b2430", "axes.labelcolor": "#1b2430", **svg_finish.mathtext_rc()})   # INK, as the labels set it
INK, MUTE = "#1b2430", "#6b7480"
GREEN, RED, BLUE_, BOTH, GREY = "#4DAF4A", "#d6604d", "#4393c3", "#9e7bb5", "#b8bec6"
WAYS = ["exactcmp", "exact", "svdabs", "svd"]
WAYNAME = {"exactcmp": "exact nullspace, exact comparison", "exact": "exact nullspace, rounded to float64",
           "svdabs": "SVD nullspace, element-wise", "svd": "SVD nullspace, row-relative"}
KEYS = ("h_acr", "h_acrr", "h_mod")
FILL_UNFINISHED = bool(os.environ.get("KF_FILL_UNFINISHED"))   # draft of the figure while the last parts run

# ---- merge the parts per model
runs = defaultdict(dict)
for f in sorted(glob.glob(f"{KF}/parts/*.json")):
    d = json.load(open(f))
    r = runs[d["model"]]
    r.setdefault("taus", d["taus"])
    for k, v in d.items():
        if k.endswith("_relations"):
            r[k] = v
        elif isinstance(v, dict) and "sweep" in v:
            a = r.setdefault(k, {"sweep": {}})
            for t, s in zip(v.get("taus", d["taus"]), v["sweep"]):
                a["sweep"][t] = s
            if "check" in v:
                a["check"] = v["check"]
TAUS = next(iter(runs.values()))["taus"]
complete = lambda m, arm: arm in runs[m] and len(runs[m][arm]["sweep"]) == len(TAUS) and "check" in runs[m][arm]
models_for = lambda arms: sorted(m for m in runs if all(complete(m, a) for a in arms)
                                 and "iterated_exact_relations" in runs[m] and "tangent_exact_relations" in runs[m])
X = np.arange(len(TAUS)); XE, WE = -2.5, 2.3                     # bar positions: tolerance 0 (bitwise equal rows), then the grid
dec = [k for k, t in enumerate(TAUS) if abs(np.log10(t) - round(np.log10(t))) < 1e-9 and round(np.log10(t)) % 3 == 1]


def save_svg(fig, stem, panels=None):
    """SVG with live text, embedded fonts and no clipping, plus a PNG for checking. panels: the groups, as for
    svg_panels.Figure.assign (default one group for everything)"""
    buf = io.StringIO(); fig.savefig(buf, format="svg"); raw = buf.getvalue()
    n, box = svg_finish.check(raw, raw, tempfile.mkdtemp())
    assert n == 0, f"{os.path.basename(stem)}: removing the clipping changes {n} pixels in {box} (pt)"
    fg = svg_panels.Figure(raw)
    miss = fg.assign(panels or [("all", [(0, 0, 1e4, 1e4)])])
    assert not miss, f"{os.path.basename(stem)}: {len(miss)} elements outside the panels"
    open(stem + ".svg", "w", encoding="utf-8").write(svg_finish.finish(fg.svg())[0])
    fig.savefig(stem + ".png", dpi=200)


# ---- figs. S5 and S6: plots over the full page width, the colour key in one row below them (all in pt)
W = 518.4                                            # page width (7.2 in)
LEFT, RIGHT, WGAP = 30.0, 4.0, 9.0                   # y label and ticks, right edge, between the plots
KEY = 21.0                                           # band at the bottom that holds a key of two-line labels


def place(fig, nrow, ncol, top, bottom, hgap=0.0):
    """Places the axes of plt.subplots in pt from the page edges and returns the plot width and height."""
    H = fig.get_figheight() * 72
    w = (W - LEFT - RIGHT - (ncol - 1) * WGAP) / ncol
    h = (H - top - bottom - (nrow - 1) * hgap) / nrow
    fig.subplots_adjust(left=LEFT / W, right=1 - RIGHT / W, top=1 - top / H, bottom=bottom / H, wspace=WGAP / w, hspace=hgap / h)
    return w, h


def groups(fig, ncol, w, row_edges, key=KEY):
    """One group per plot (A, B, ... row by row, the columns split halfway between the axes) and one for the key.
    row_edges: the y in pt from the top between the rows, from 0 to the top of the key band (key pt high)."""
    H = fig.get_figheight() * 72
    xs = [0.0] + [LEFT + j * (w + WGAP) - WGAP / 2 for j in range(1, ncol)] + [W]
    names = iter("ABCDEFGH")
    out = [("legend", [(0, H - key, W, H)])]
    for y0, y1 in zip(row_edges, row_edges[1:]):
        out += [(next(names), [(xs[j], y0, xs[j + 1], y1)]) for j in range(ncol)]
    return out


def key_below(fig, handles, labels):
    """The colour key in one row under the plots, centred on them."""
    H = fig.get_figheight() * 72
    fig.legend(handles, labels, frameon=False, loc="lower center", ncol=len(labels), fontsize=6, borderaxespad=0,
               bbox_to_anchor=((LEFT + W - RIGHT) / 2 / W, 3 / H), handlelength=1, handletextpad=.4, columnspacing=1.5)


def sparse_ticks(ax):
    """0, then every third decade from 1e-11 (a label at 1e-14 would touch the 0 beside it)."""
    ax.set_xticks([XE] + [X[k] for k in dec[1:]])
    ax.set_xticklabels(["0"] + [f"$\\mathregular{{10^{{{round(np.log10(TAUS[k]))}}}}}$" for k in dec[1:]])


def axes_style(ax, n):
    ax.set_xticks([XE] + [X[k] for k in dec])
    ax.set_xticklabels(["0"] + [f"$\\mathregular{{10^{{{round(np.log10(TAUS[k]))}}}}}$" for k in dec])    # Arial, as every tick label
    ax.set_xlim(XE - WE / 2 - 0.2, X[-1] + 0.6); ax.set_ylim(0, max(n, 1))
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.tick_params(axis="both", length=2, colors=INK)


def bars(ax, cls, zero, colors):
    bottom = np.zeros(len(TAUS)); bz = 0
    for c, col in enumerate(colors):
        ax.bar(X, cls[c], bottom=bottom, width=0.86, color=col, lw=0)
        ax.bar([XE], [zero[c]], bottom=[bz], width=WE, color=col, lw=0)
        bottom += cls[c]; bz += zero[c]


same = lambda a, b: all(a[k] == b[k] for k in KEYS)


LABELS = ["routes agree on a wrong answer\n(false certificate)", "routes differ, unlike in\nexact arithmetic",
          "routes differ as in exact\narithmetic (the gap)", "routes agree and are right\n(certified)"]


def draw_certificate(axs, ways, titles=None):
    """One stacked-bar panel per way on the given axes. While a way is unfinished a tolerance counts only the models it has."""
    out = {}
    for ax, way in zip(axs, ways):
        models = sorted(m for m in runs if f"tangent_{way}" in runs[m] and f"iterated_{way}" in runs[m] and "check" in runs[m][f"iterated_{way}"]
                        and "iterated_exact_relations" in runs[m] and "tangent_exact_relations" in runs[m])
        cls = np.zeros((4, len(TAUS))); zero = np.zeros(4); ex = np.zeros(4); filled = 0
        for m in models:
            xi, xt = runs[m]["iterated_exact_relations"], runs[m]["tangent_exact_relations"]
            gap = not same(xi, xt); ex[1 if gap else 0] += 1

            def klass(fi, ft):
                agree, ok = same(fi, ft), same(fi, xi) and same(ft, xt)
                return (0 if ok else 3) if agree else (1 if (gap and ok) else 2)

            zero[klass(runs[m][f"iterated_{way}"]["check"], runs[m][f"tangent_{way}"]["check"])] += 1
            for k, t in enumerate(TAUS):
                si, st = runs[m][f"iterated_{way}"]["sweep"].get(t), runs[m][f"tangent_{way}"]["sweep"].get(t)
                if si is not None and st is not None:
                    cls[klass(si, st), k] += 1
                elif FILL_UNFINISHED:                      # draft only: a part still running counts as "routes differ"
                    cls[2, k] += 1; filled += 1
        bars(ax, cls, zero, (GREEN, BLUE_, GREY, RED))
        axes_style(ax, len(models))
        ax.axhline(ex[0], color=INK, lw=0.6, ls=(0, (3, 2)))
        ax.set_title(titles[way] if titles else WAYNAME[way], fontsize=8, color=INK, pad=3)
        # the bar at tolerance 0 is the result: wider than the others, and labelled in place
        bottom = 0
        for c, word in ((0, "certified"), (1, "gap"), (2, "differ"), (3, "false certificates")):
            if zero[c] >= 8:
                ax.text(XE, bottom + zero[c] / 2, f"{zero[c]:.0f} {word}", rotation=90, ha="center", va="center", fontsize=7,
                        fontweight="bold", color="white", zorder=5)
            elif zero[c] > 0:
                ax.text(XE, bottom + zero[c] / 2, f"{zero[c]:.0f}", ha="center", va="center", fontsize=6, fontweight="bold",
                        color=INK if c == 2 else "white", zorder=5)
            bottom += zero[c]
        out[way] = dict(models=models, certified=cls[0].astype(int).tolist(), gap=cls[1].astype(int).tolist(),
                        float_differ=cls[2].astype(int).tolist(), false_certificate=cls[3].astype(int).tolist(),
                        per_tau=cls.sum(axis=0).astype(int).tolist(), at_zero=zero.astype(int).tolist(), exact=ex.astype(int).tolist(),
                        filled=filled)
    return out


def certificate_legend(fig, x, y, exact):
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (RED, GREY, BLUE_, GREEN)] + \
              [plt.Line2D([], [], color=INK, lw=0.6, ls=(0, (3, 2)))]
    fig.legend(handles, LABELS + [f"certified in exact arithmetic\n({exact[0]:.0f} of {sum(exact):.0f})"], frameon=False, loc="upper left",
               fontsize=6, bbox_to_anchor=(x, y), handlelength=1, handletextpad=.4, borderaxespad=0)


def certificate(ways, name):
    """Figure3EF_draft (two ways, the key at the right as in Figure 3) or, for more ways, fig. S5 (the key below)."""
    if len(ways) > 2:
        return certificate_supp(ways, name)
    fig, axs = plt.subplots(1, len(ways), figsize=(7.2, 2.2), sharey=True, squeeze=False)
    axs = axs[0]
    fig.subplots_adjust(left=0.07, right=0.76, bottom=0.2, top=0.86, wspace=0.08)
    out = draw_certificate(axs, ways, {"exactcmp": "exact nullspace", "svdabs": "SVD nullspace"})
    for ax in axs:
        ax.set_xlabel("tolerance (0 = rows exactly equal)", color=INK)
    axs[0].set_ylabel("networks", color=INK)
    certificate_legend(fig, 0.77, 0.86, out[ways[0]]["exact"])
    save_svg(fig, f"{OUT}/{name}"); plt.close(fig)
    return dict(n=max(len(v["models"]) for v in out.values()), ways=out)


def certificate_supp(ways, name):
    """fig. S5: one plot per way over the full width, the key in one row below."""
    fig, axs = plt.subplots(1, len(ways), figsize=(W / 72, 2.6), sharey=True, squeeze=False)
    axs = axs[0]
    w, h = place(fig, 1, len(ways), top=26, bottom=52)
    out = draw_certificate(axs, ways, {w_: WAYNAME[w_].replace(", ", ",\n") for w_ in ways})
    for ax in axs:
        sparse_ticks(ax)
        ax.set_xlabel("tolerance (0 = rows exactly equal)", color=INK)
    axs[0].set_ylabel("networks", color=INK)
    ex = out[ways[0]]["exact"]
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (GREEN, BLUE_, GREY, RED)] + \
              [plt.Line2D([], [], color=INK, lw=0.6, ls=(0, (3, 2)))]
    key_below(fig, handles, LABELS[::-1] + [f"certified in exact arithmetic\n({ex[0]:.0f} of {sum(ex):.0f})"])
    save_svg(fig, f"{OUT}/{name}", groups(fig, len(ways), w, [0, fig.get_figheight() * 72 - KEY])); plt.close(fig)
    return dict(n=max(len(v["models"]) for v in out.values()), ways=out)


def errors(ways, name):
    arms = [f"{r}_{w}" for r in ("iterated", "tangent") for w in ways]
    models = models_for(arms); n = len(models)
    fig, axs = plt.subplots(2, len(ways), figsize=(W / 72, 4.4), sharey=True, squeeze=False)
    w, h = place(fig, 2, len(ways), top=26, bottom=47, hgap=40)
    out = {}
    klass = lambda r: (lambda miss, fake: 0 if not (miss or fake) else 2 if not fake else 1 if not miss else 3)(
        r["acr_fn"] + r["acrr_fn"] > 0, r["acr_fp"] + r["acrr_fp"] > 0)
    for i, route in enumerate(("iterated", "tangent")):
        for j, way in enumerate(ways):
            ax, arm = axs[i, j], f"{route}_{way}"
            cls = np.zeros((4, len(TAUS))); zero = np.zeros(4)
            for m in models:
                zero[klass(runs[m][arm]["check"])] += 1
                for k, t in enumerate(TAUS):
                    cls[klass(runs[m][arm]["sweep"][t]), k] += 1
            bars(ax, cls, zero, (GREEN, RED, BLUE_, BOTH)); axes_style(ax, n); sparse_ticks(ax)
            ax.set_title(f"{route} route, " + WAYNAME[way].replace(", ", "\n", 1), fontsize=8, color=INK, pad=3)   # fits the plot width
            out[arm] = dict(right=cls[0].astype(int).tolist(), at_zero=zero.astype(int).tolist())
    for ax in axs[1]: ax.set_xlabel("tolerance (0 = rows exactly equal)", color=INK)
    for ax in axs[:, 0]: ax.set_ylabel("networks")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (GREEN, RED, BLUE_, BOTH)]
    key_below(fig, handles, ["all relations right", "relations made up", "relations missed", "missed and made up"])
    save_svg(fig, f"{OUT}/{name}", groups(fig, len(ways), w, [0, 26 + h + 16, fig.get_figheight() * 72 - 14], key=14)); plt.close(fig)
    return dict(models=models, n=n, arms=out)


if __name__ == "__main__":
    OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
    summary = {}
    avail = [w for w in WAYS if any(f"iterated_{w}" in runs[m] for m in runs)]
    if {"exactcmp", "svdabs"} <= set(avail):
        summary["Figure3EF_draft"] = certificate(["exactcmp", "svdabs"], "Figure3EF_draft")
    if avail:
        summary["FigS5_certificate"] = certificate(avail, "FigS5_certificate")
        summary["FigS6_float_errors"] = errors(avail, "FigS6_float_errors")
    json.dump(dict(taus=TAUS, **summary), open(f"{OUT}/float_four_summary.json", "w"), indent=1)
    for k, v in summary.items():
        print(k, f"{v['n']} models")
        for way, c in v.get("ways", {}).items():
            print(f"  {way:9s} at 0: {c['at_zero']}  best certified {max(c['certified'])}  max false certificates "
                  f"{max(c['false_certificate'] + [c['at_zero'][3]])}  exact {c['exact']}")
