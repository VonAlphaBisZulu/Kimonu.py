"""Proven <= true <= local on the 33 benchmark models: the iterated route against the tangent route.

Per model (ordered by number of complexes) an interval from the iterated route's count (lower bound,
proven for all rate constants) to the tangent route's count (upper bound, local ACR in the sense of
Pascual-Escudero & Feliu); the true count lies in between. Where both coincide the model is decided and
gets one marker. Left: pairs of complexes with a fixed concentration ratio (kinetic coupling). Right:
fixed metabolite relations, ACR metabolites plus ACRR pairs, free metabolites as in Figure 3. Attached
uptakes. Data: kinetic_variants/hierarchy/<model>.json (hierarchy_data.py).

    python figures/make_hierarchy.py [out_dir]   # from benchmarks/; default figures/main
"""
import glob, json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

out = sys.argv[1] if len(sys.argv) > 1 else "figures/main"
# composed into Figure 2 by make_figure2.py: a shorter panel, no letters (the composition adds them)
HEIGHT_IN = float(os.environ.get("FIG_HEIGHT_IN", 2.9)); LETTERS = os.environ.get("FIG_LETTERS", "1") == "1"
NAME = os.environ.get("FIG_NAME", "hierarchy")
LEGEND_BOTTOM = os.environ.get("FIG_LEGEND_BOTTOM", "0") == "1"   # Figure 2: the top band holds the panel letters
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "svg.fonttype": "none", "pdf.fonttype": 42})
INK, MUTE, RULE = "#1b2430", "#6b7480", "#d5d9de"
GREY, DARK = "#6b7480", "#3a4250"              # no colours: the two bounds differ by marker shape only
D = sorted((json.load(open(f)) for f in glob.glob("kinetic_variants/hierarchy/*.json")), key=lambda d: d["n_complexes"])
SHORT = {"Chlorella_variabilis_iAJ526": "iAJ526", "Saha2011iRS1597": "iRS1597", "ArabidopsisCoreModel": "AraCore",
         "iIS312_Trypomastigote": "iIS312 T", "iIS312_Amastigote": "iIS312 A", "iIS312_Epimastigote": "iIS312 E"}
LABEL_BELOW = {("A", "iCN718")}
bad = [d["model"] for d in D if any(d["violations"].values())]
assert not bad, f"iterated relations missing from the tangent route: {bad}"

fig, axs = plt.subplots(1, 2, figsize=(17.5 / 2.54, HEIGHT_IN), sharex=True)
panels = [("pairs", "coupled complex pairs", "A"),
          (lambda r: r["acr_free"] + r["acrr_pairs_free"], "ACR metabolites + ACRR pairs (free)", "B")]
for ax, (key, ylabel, letter) in zip(axs, panels):
    get = (lambda r: r[key]) if isinstance(key, str) else key
    for x, d in enumerate(D):
        lo, hi = get(d["it"]), get(d["tan"])
        lo_, hi_ = max(lo, 0.7), max(hi, 0.7)            # zero counts sit on the floor of the log axis
        if lo == hi:
            ax.scatter(x, lo_, s=10, color=GREY, zorder=3, linewidths=0)
        else:
            # the two bounds side by side, so a gap that is small on the log axis stays visible
            ax.plot([x - 0.22, x + 0.22], [lo_, hi_], color=GREY, lw=1.0, zorder=2)
            ax.scatter(x - 0.22, lo_, s=20, facecolor="white", edgecolor=DARK, linewidths=0.9, zorder=3)
            ax.scatter(x + 0.22, hi_, s=18, marker="D", facecolor="white", edgecolor=DARK, linewidths=0.9, zorder=3)
            below = (letter, d["model"]) in LABEL_BELOW       # placed by hand where a neighbour sits to the upper right
            ax.annotate(f"+{hi - lo:,}", (x + 0.22, hi_), xytext=(4, -3 if below else 2), textcoords="offset points",
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
handles = [Line2D([], [], marker="o", ls="", color=GREY, ms=3.5, label="decided: iterated = tangent"),
           Line2D([], [], marker="o", ls="", mfc="white", mec=DARK, mew=0.9, ms=5, label="iterated route (proven)"),
           Line2D([], [], marker="D", ls="", mfc="white", mec=DARK, mew=0.9, ms=4, label="tangent route (local ACR); +n = undecided"),]
if LEGEND_BOTTOM:
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
else:
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0.03, 1, 0.92))
os.makedirs(out, exist_ok=True)
for ext in ("svg", "png"):
    fig.savefig(f"{out}/{NAME}.{ext}", dpi=300)
print("gap models:", [d["model"] for d in D if d["gap"]])
