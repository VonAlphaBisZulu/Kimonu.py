"""Proven <= true <= local on the 33 benchmark models: the iterated route against the tangent route.

Per model (ordered by number of complexes) an interval from the iterated route's count (lower bound,
proven for all rate constants) to the tangent route's count (upper bound, local ACR in the sense of
Pascual-Escudero & Feliu); the true count lies in between. Where both coincide the model is decided and
gets one marker. Left: pairs of complexes with a fixed concentration ratio (kinetic coupling). Right:
fixed metabolite relations, ACR metabolites plus ACRR pairs, free metabolites as in Figure 2. Attached
uptakes. Data: kinetic_variants/hierarchy/<model>.json (hierarchy_data.py).

    python figures/make_hierarchy.py [out_dir]   # from benchmarks/; default figures/main
"""
import glob, json, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

out = sys.argv[1] if len(sys.argv) > 1 else "figures/main"
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "svg.fonttype": "none", "pdf.fonttype": 42})
INK, MUTE, RULE = "#1b2430", "#6b7480", "#d5d9de"
GREEN, PURPLE = "#4DAF4A", "#984EA3"          # Kimonu green as in Figure 2; local ACR in Set1 purple
D = sorted((json.load(open(f)) for f in glob.glob("kinetic_variants/hierarchy/*.json")), key=lambda d: d["n_complexes"])
SHORT = {"Chlorella_variabilis_iAJ526": "iAJ526", "Saha2011iRS1597": "iRS1597", "ArabidopsisCoreModel": "AraCore",
         "iIS312_Trypomastigote": "iIS312 T", "iIS312_Amastigote": "iIS312 A", "iIS312_Epimastigote": "iIS312 E"}
LABEL_BELOW = {("A", "iCN718")}
bad = [d["model"] for d in D if any(d["violations"].values())]
assert not bad, f"iterated relations missing from the tangent route: {bad}"

fig, axs = plt.subplots(1, 2, figsize=(17.5 / 2.54, 2.9), sharex=True)
panels = [("pairs", "coupled complex pairs", "A"),
          (lambda r: r["acr_free"] + r["acrr_pairs_free"], "ACR metabolites + ACRR pairs (free)", "B")]
for ax, (key, ylabel, letter) in zip(axs, panels):
    get = (lambda r: r[key]) if isinstance(key, str) else key
    for x, d in enumerate(D):
        lo, hi = get(d["it"]), get(d["tan"])
        lo_, hi_ = max(lo, 0.7), max(hi, 0.7)            # zero counts sit on the floor of the log axis
        if lo == hi:
            ax.scatter(x, lo_, s=10, color=INK, zorder=3, linewidths=0)
        else:
            # the two bounds side by side, so a gap that is small on the log axis stays visible
            ax.plot([x - 0.22, x + 0.22], [lo_, hi_], color=MUTE, lw=1.0, zorder=2)
            ax.scatter(x - 0.22, lo_, s=20, color=GREEN, edgecolor="white", linewidths=0.5, zorder=3)
            ax.scatter(x + 0.22, hi_, s=18, marker="D", facecolor="white", edgecolor=PURPLE, linewidths=1.0, zorder=3)
            below = (letter, d["model"]) in LABEL_BELOW       # placed by hand where a neighbour sits to the upper right
            ax.annotate(f"+{hi - lo:,}", (x + 0.22, hi_), xytext=(4, -3 if below else 2), textcoords="offset points",
                        ha="left", va="top" if below else "bottom", fontsize=6, color=PURPLE)
    ax.set_yscale("log"); ax.set_ylabel(ylabel, color=INK)
    ax.yaxis.set_major_locator(LogLocator(base=10)); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else ""))
    ax.set_xticks(range(len(D))); ax.set_xticklabels([SHORT.get(d["model"], d["model"]) for d in D], rotation=90)
    ax.tick_params(axis="both", colors=INK, length=2)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(MUTE)
    ax.grid(axis="y", color=RULE, lw=0.4); ax.set_axisbelow(True)
    ax.text(-0.12, 1.04, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK, va="bottom")
axs[0].set_xlabel(""); axs[1].set_xlabel("")
fig.text(0.5, 0.005, "models, by number of complexes", ha="center", va="bottom", color=INK)
handles = [Line2D([], [], marker="o", ls="", color=INK, ms=3.5, label="decided: iterated = tangent"),
           Line2D([], [], marker="o", ls="", color=GREEN, mec="white", ms=5.5, label="iterated route (proven)"),
           Line2D([], [], marker="D", ls="", mfc="white", mec=PURPLE, ms=4, label="tangent route (local ACR); +n = undecided"),]
fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0.03, 1, 0.92))
os.makedirs(out, exist_ok=True)
for ext in ("svg", "png"):
    fig.savefig(f"{out}/hierarchy.{ext}", dpi=300)
print("gap models:", [d["model"] for d in D if d["gap"]])
