"""Figure 3A, B. Proven <= true <= complete on the 33 benchmark models: the iterated route against the tangent route.

Per model (ordered by number of complexes) an interval from the iterated route's count (lower bound,
proven for all rate constants) to the tangent route's count (upper bound, complete but with possible false
couplings); the true count lies in between. Where both coincide the model is decided and
gets one marker. Left: pairs of complexes with a fixed concentration ratio (kinetic coupling). Right:
fixed metabolite relations, ACR metabolites plus independent ACRR ratios (k - 1 per group of k free metabolites),
free metabolites as in Figure 4. Detached
uptakes. Data: HIERARCHY_DIR, default kinetic_variants/hierarchy_detached/<model>.json (hierarchy_detached.py from
uptake_modes.py), the input of Figure 3A, B. HIERARCHY_DIR=kinetic_variants/hierarchy draws the attached panels of
October (hierarchy_data.py, removed, see HISTORY.md).

    python figures/fig3AB_hierarchy.py [out_dir]   # from scripts/; default figures/main

Writes <out_dir>/Figure3AB_draft.svg. The published Figure 3 carries hand edits and was updated in place by
fig3B_update_ratios.py, fig3AB_update_legend.py, fig3CD_coverage.py and fig3_fig4_update_terms.py.
"""
import glob, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

out = sys.argv[1] if len(sys.argv) > 1 else "figures/main"
# the combined Figure 2 of October used FIG_HEIGHT_IN=2.3 FIG_LETTERS=0 FIG_LEGEND_BOTTOM=1 (make_figure2.py, see HISTORY.md)
HEIGHT_IN = float(os.environ.get("FIG_HEIGHT_IN", 2.9)); LETTERS = os.environ.get("FIG_LETTERS", "1") == "1"
NAME = os.environ.get("FIG_NAME", "Figure3AB_draft")
LEGEND_BOTTOM = os.environ.get("FIG_LEGEND_BOTTOM", "0") == "1"   # combined figure: the top band holds the panel letters
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "svg.fonttype": "none", "pdf.fonttype": 42})
INK, MUTE, RULE = "#1b2430", "#6b7480", "#d5d9de"
GREY, DARK = "#6b7480", "#3a4250"              # no colours: the two bounds differ by marker shape only
HIERARCHY_DIR = os.environ.get("HIERARCHY_DIR", "kinetic_variants/hierarchy_detached")
D = sorted((json.load(open(f)) for f in glob.glob(f"{HIERARCHY_DIR}/*.json")), key=lambda d: d["n_complexes"])
from names import NAME as SHORT                    # one display name per network (figures/names.py)
LABEL_BELOW = {("A", "iCN718")}
bad = [d["model"] for d in D if any(d["violations"].values())]
assert not bad, f"iterated relations missing from the tangent route: {bad}"

fig, axs = plt.subplots(1, 2, figsize=(17.5 / 2.54, HEIGHT_IN), sharex=True)
panels = [("pairs", "coupled complex pairs", "A"),
          (lambda r: r["acr_free"] + r["acrr_ratios_free"], "ACR metabolites + ACRR ratios (free)", "B")]
for ax, (key, ylabel, letter) in zip(axs, panels):
    get = (lambda r: r[key]) if isinstance(key, str) else key
    for x, d in enumerate(D):
        lo, hi = get(d["it"]), get(d["tan"])
        lo_, hi_ = max(lo, 0.7), max(hi, 0.7)            # zero counts sit on the floor of the log axis
        if lo == hi:
            ax.scatter(x, lo_, s=10, color=GREY, zorder=3, linewidths=0)
        else:
            # the bracket as an error bar: iterated route (proven) at the bottom, tangent route (local) at the top
            ax.errorbar(x, hi_, yerr=[[hi_ - lo_], [0]], fmt="none", ecolor=DARK, elinewidth=0.9, capsize=2.2,
                        capthick=0.9, zorder=3)
            below = (letter, d["model"]) in LABEL_BELOW       # placed by hand where a neighbour sits to the upper right
            ax.annotate(f"+{hi - lo:,}", (x, hi_), xytext=(3, -3 if below else 2), textcoords="offset points",
                        ha="left", va="top" if below else "bottom", fontsize=6, color=DARK)
    ax.set_yscale("log"); ax.set_ylabel(ylabel, color=INK)
    ax.yaxis.set_major_locator(LogLocator(base=10)); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else ""))
    ax.set_xticks(range(len(D))); ax.set_xticklabels([SHORT.get(d["model"], d["model"]) for d in D], rotation=90)
    ax.tick_params(axis="both", colors=INK, length=2)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(MUTE)
    ax.grid(axis="y", color=RULE, lw=0.4); ax.set_axisbelow(True)
    if LETTERS:
        ax.text(-0.12, 1.04, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="bottom")
axs[0].set_xlabel(""); axs[1].set_xlabel("")
if not LEGEND_BOTTOM:
    fig.text(0.5, 0.005, "models, by number of complexes", ha="center", va="bottom", color=INK)
bar = axs[0].errorbar([np.nan], [np.nan], yerr=[[np.nan], [np.nan]], fmt="none", ecolor=DARK, elinewidth=0.9, capsize=2.2,
                      capthick=0.9, label="undecided, from iterated route (proven, bottom) to tangent route (complete, top)       +n candidates that require solving polynomials")
handles = [Line2D([], [], marker="o", ls="", color=GREY, ms=3.5, label="decided, iterated = tangent"), bar]
if LEGEND_BOTTOM:
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
else:
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0.03, 1, 0.92))
os.makedirs(out, exist_ok=True)
for ext in ("svg", "png"):
    fig.savefig(f"{out}/{NAME}.{ext}", dpi=300)
print("gap models:", [d["model"] for d in D if d["gap"]])
