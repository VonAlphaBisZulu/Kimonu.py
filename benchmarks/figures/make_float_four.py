"""Drafts from kinetic_variants/float_four.py and float_exactcmp.py: floating point against exact arithmetic, for
both routes and four ways of computing and comparing the nullspace rows.

    exactcmp   exact kernel (Kimonu's RREF bases), rows kept as fractions, element-wise tolerance, decided exactly
    exact      exact kernel rounded to an orthonormal float64 basis, row-relative Euclidean tolerance
    svdabs     SVD kernel (orthonormal), element-wise tolerance
    svd        SVD kernel (orthonormal), row-relative Euclidean tolerance

Fig_certificate     Figure 2I draft: per tolerance, what the comparison of the two routes certifies, exact kernel
                    (exactcmp) against SVD kernel (svdabs), one element-wise rule. Classes: the routes agree and are
                    right (the complete set, certified), the routes differ as they do exactly (the gap models), the
                    routes differ only because of floating point, the routes agree on a wrong answer (a false
                    certificate). The first bar is tolerance 0 (bitwise equal rows); the dashed line is what exact
                    arithmetic certifies. Answers are compared as ACR, ACRR and kinetic modules.
FigS_certificate    the same for all four ways.
FigS_float_errors   per route and way: the models whose ACR and ACRR come out right, only made up, only missed, or both.

    python figures/make_float_four.py <out_dir>        # from benchmarks/
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
                     **svg_finish.mathtext_rc()})
INK, MUTE = "#1b2430", "#6b7480"
GREEN, RED, BLUE_, BOTH, GREY = "#4DAF4A", "#d6604d", "#4393c3", "#9e7bb5", "#b8bec6"
WAYS = ["exactcmp", "exact", "svdabs", "svd"]
WAYNAME = {"exactcmp": "exact kernel, exact comparison", "exact": "exact kernel, rounded to float64",
           "svdabs": "SVD kernel, element-wise", "svd": "SVD kernel, row-relative"}
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


def save_svg(fig, stem):
    """SVG with live text, embedded fonts and no clipping, plus a PNG for checking"""
    buf = io.StringIO(); fig.savefig(buf, format="svg"); raw = buf.getvalue()
    n, box = svg_finish.check(raw, raw, tempfile.mkdtemp())
    assert n == 0, f"{os.path.basename(stem)}: removing the clipping changes {n} pixels in {box} (pt)"
    fg = svg_panels.Figure(raw)
    fg.assign([("all", [(0, 0, 1e4, 1e4)])])
    open(stem + ".svg", "w", encoding="utf-8").write(svg_finish.finish(fg.svg())[0])
    fig.savefig(stem + ".png", dpi=200)


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
    fig, axs = plt.subplots(1, len(ways), figsize=(7.2, 2.2), sharey=True, squeeze=False)
    axs = axs[0]
    right = 0.76 if len(ways) == 2 else 0.8
    fig.subplots_adjust(left=0.07, right=right, bottom=0.2, top=0.8 if len(ways) > 2 else 0.86, wspace=0.08)
    titles = {"exactcmp": "exact kernel", "svdabs": "SVD kernel"} if len(ways) == 2 else {w: WAYNAME[w].replace(", ", ",\n") for w in ways}
    out = draw_certificate(axs, ways, titles)
    if len(ways) > 2:                                    # narrow panels: fewer tick labels, so that 0 stays clear
        for ax in axs:
            ax.set_xticks([XE] + [X[k] for k in dec[1::2]])
            ax.set_xticklabels(["0"] + [f"$\\mathregular{{10^{{{round(np.log10(TAUS[k]))}}}}}$" for k in dec[1::2]])
            ax.title.set_fontsize(7)
    for ax in axs:
        ax.set_xlabel("tolerance (0: rows exactly equal)", color=INK)
    axs[0].set_ylabel("models", color=INK)
    certificate_legend(fig, right + 0.01, 0.8 if len(ways) > 2 else 0.86, out[ways[0]]["exact"])
    save_svg(fig, f"{OUT}/{name}"); plt.close(fig)
    return dict(n=max(len(v["models"]) for v in out.values()), ways=out)


def errors(ways, name):
    arms = [f"{r}_{w}" for r in ("iterated", "tangent") for w in ways]
    models = models_for(arms); n = len(models)
    fig, axs = plt.subplots(2, len(ways), figsize=(7.2, 4.0), sharey=True, squeeze=False)
    fig.subplots_adjust(left=0.06, right=0.84, bottom=0.11, top=0.9, wspace=0.08, hspace=0.5)
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
            bars(ax, cls, zero, (GREEN, RED, BLUE_, BOTH)); axes_style(ax, n)
            ax.set_title(f"{route} route\n{WAYNAME[way]}", fontsize=7, color=INK, pad=3)
            out[arm] = dict(right=cls[0].astype(int).tolist(), at_zero=zero.astype(int).tolist())
    for ax in axs[1]: ax.set_xlabel("tolerance (0: rows exactly equal)", color=INK)
    for ax in axs[:, 0]: ax.set_ylabel("models")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (BOTH, BLUE_, RED, GREEN)]
    fig.legend(handles, ["missed and made up", "relations missed", "relations made up", "all relations right"],
               frameon=False, loc="upper left", fontsize=6, bbox_to_anchor=(0.85, 0.9), handlelength=1, handletextpad=.4)
    save_svg(fig, f"{OUT}/{name}"); plt.close(fig)
    return dict(models=models, n=n, arms=out)


if __name__ == "__main__":
    OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
    summary = {}
    avail = [w for w in WAYS if any(f"iterated_{w}" in runs[m] for m in runs)]
    if {"exactcmp", "svdabs"} <= set(avail):
        summary["Fig_certificate"] = certificate(["exactcmp", "svdabs"], "Fig_certificate")
    if avail:
        summary["FigS_certificate"] = certificate(avail, "FigS_certificate")
        summary["FigS_float_errors"] = errors(avail, "FigS_float_errors")
    json.dump(dict(taus=TAUS, **summary), open(f"{OUT}/float_four_summary.json", "w"), indent=1)
    for k, v in summary.items():
        print(k, f"{v['n']} models")
        for way, c in v.get("ways", {}).items():
            print(f"  {way:9s} at 0: {c['at_zero']}  best certified {max(c['certified'])}  max false certificates "
                  f"{max(c['false_certificate'] + [c['at_zero'][3]])}  exact {c['exact']}")
