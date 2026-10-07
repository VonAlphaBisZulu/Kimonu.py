"""Figure 3 of the paper: Kimonu in its reference-equivalent mode reproduces the published analysis exactly, in seconds.

A  the key quantities on the eleven models where the exhaustive reference finished. Every value is identical
   in both tools (asserted below), so each cell shows the shared value.
B  wall time of that same computation, row-aligned with A, and the speedup.
C  wall time vs network size on all 33 models for every arm, with a log-log least-squares line per arm.
Supplementary figure (figures/supp/FigS_stage_timing): where the time goes on those eleven models: the concordance stage against everything after it (kinetic
   modules, ACR/ACRR), per tool, as the range over the models with the median marked. Measured vs derived:
   Kimonu's stage times are measured (median of three replicates); the reference's concordance stage is its
   own measured time, and its R stage the wall time between the chain log's start stamp and the end of its
   ACR log (the R stage reads the stoichiometric coupling from a file, whose computation is not included);
   COCOA's concordance stage is its own logged "Elapsed time", which ends before kinetic_analysis is called,
   and the rest is the measured total minus that.

Kimonu here is run_kimonu_refmode.py (triplicate_refmode/, the median of three replicates), which follows the
reference's definitions: concordance from ker(N), upstream sets merged by overlap only, ACR/ACRR by Corollary
S2-10. Stochastic arms show the median over three seeds. Layout is in inches at 17.5 cm (double column) width;
the table columns are sized from their measured text, so the axes of A and B are placed by hand.

    python figures/make_figure3.py [out_dir]    # run from benchmarks/ after gather_figdata.py; default figures/main
"""
import datetime as dt, glob, json, os, re, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

out = sys.argv[1] if len(sys.argv) > 1 else "figures/main"
D = json.load(open("figdata.json")); nC = D["n_complexes"]
REF = sorted(D["ref_exh"], key=lambda m: nC[m])

# the agreement the figure claims, checked on every quantity the reference reports (free and all metabolites)
KEYS = ["balanced", "conc", "n_kinetic", "giant", "acr_free", "acrr_free", "acr", "acrr"]
for m in REF:
    for k in KEYS:
        assert D["kimonu_ref"][m][k] == D["ref_exh"][m][k], (m, k, D["kimonu_ref"][m][k], D["ref_exh"][m][k])

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                     "svg.fonttype": "none", "pdf.fonttype": 42})
INK, MUTE, RULE, STRIPE = "#1b2430", "#6b7480", "#d5d9de", "#f1f3f5"
GREEN, BLUE, RED = "#4DAF4A", "#377EB8", "#E41A1C"          # Set1, as in the earlier drafts; shapes carry identity too
# arm: label, colour, marker, filled, fit line style
ARMS = {"kimonu_ref":  ("Kimonu, reference-equivalent mode", GREEN, "o", True, "-"),
        "ref_exh":     ("Reference, exhaustive",             BLUE,  "s", True, "-"),
        "ref_sampled": ("Reference, sampled",                BLUE,  "s", False, "--"),
        "cocoa_false": ("COCOA.jl, efficient=false",         RED,   "^", True, "-"),
        "cocoa_true":  ("COCOA.jl, efficient=true",          RED,   "^", False, "--")}
SHORT = {"iIS312_Trypomastigote": "iIS312 Trypomastigote", "iIS312_Amastigote": "iIS312 Amastigote",
         "iIS312_Epimastigote": "iIS312 Epimastigote", "Saha2011iRS1597": "iRS1597"}
COLS = [("nC", "Complexes"), ("balanced", "Balanced\ncomplexes"), ("conc", "Concordance\nmodules"),
        ("n_kinetic", "Kinetic\nmodules"), ("giant", "Largest\nmodule"), ("acr_free", "ACR\nspecies"), ("acrr_free", "ACRR\npairs")]
T = {1: "1 s", 60: "1 min", 3600: "1 h", 86400: "1 day"}

def seconds(arm, m):
    v = D[arm].get(m)
    if v is None: return None
    return float(np.median([r["seconds"] for r in v])) if isinstance(v, list) else v["seconds"]

def stages(arm, m):
    """(concordance stage, everything after it) in seconds for one model, or None; see the docstring for
    which part is measured and which derived."""
    if arm == "kimonu_ref":
        R = [json.load(open(f)) for f in sorted(glob.glob(f"triplicate_refmode/kimonu_{m}_rep*.json"))]
        return (float(np.median([r["seconds_concordance"] for r in R])),
                float(np.median([r["seconds_kinetic"] + r["seconds_acr"] for r in R]))) if R else None
    if arm == "ref_exh":
        r = D["ref_exh"][m]
        stem = f"{m}_concordant_fixed_retest" if "retest" in r else m if os.path.isfile(f"reference_final/acr/{m}.acr.log") else f"{m}_concordant_fixed"
        start = {}
        for f in ("reference_final/acr/chain.log", "reference_final/acr/chain2.log", "reference_final/chain_553338.log",
                  "reference_final/acr/chain_retest_iIT341b.log"):
            start.update({n: dt.datetime.fromisoformat(t) for n, t in re.findall(r"^== kinetic (\S+) (\S+)$", open(f).read(), re.M)})
        end = dt.datetime.fromtimestamp(os.path.getmtime(f"reference_final/acr/{stem}.acr.log")).astimezone()
        return r["seconds"], (end - start[stem]).total_seconds()
    d, pre = {"cocoa_true": ("cocoa_eff", "cocoa_eff"), "cocoa_false": ("cocoa_v110", "cocoa_v110")}[arm]
    conc, rest = [], []
    for seed in (1234, 2026, 31337):
        j, lg = f"sweep/{d}/{pre}_{m}_HiGHS_p64_seed{seed}.json", f"sweep/{d}/cell_{m}_seed{seed}.log"
        if not (os.path.isfile(j) and os.path.isfile(lg)): continue
        el = re.findall(r"Elapsed time: (\d+):(\d+):([\d.]+)", open(lg).read())
        if not el: continue
        c = int(el[-1][0]) * 3600 + int(el[-1][1]) * 60 + float(el[-1][2]); t = json.load(open(j))["seconds"]
        conc.append(c); rest.append(t - c)
    return (float(np.median(conc)), float(np.median(rest))) if conc else None

def mark(ax, x, y, arm, s=16, z=3):
    _, col, mk, filled, _ = ARMS[arm]
    ax.scatter(x, y, s=s, marker=mk, facecolor=col if filled else "white", edgecolor=col if not filled else "white",
               linewidths=0.9 if not filled else 0.5, zorder=z)

# --- layout, in inches
W, RH, HEAD, PAD = 6.89, 0.2, 2.6, 0.1                  # width, table row height, header rows, cell padding
n = len(REF)
hA = (n + HEAD) * RH
H = 0.06 + hA + 0.5 + 0.18 + 2.3 + 0.42
fig = plt.figure(figsize=(W, H)); rend = fig.canvas.get_renderer()
def width_in(s, size):
    t = fig.text(0, 0, s, fontsize=size); w = t.get_window_extent(rend).width / fig.dpi; t.remove(); return w
labw = max(width_in(SHORT.get(m, m), 7) for m in REF) + PAD
colw = [max(width_in(h, 6.5), *(width_in(str(nC[m] if k == "nC" else D["kimonu_ref"][m][k]), 7) for m in REF)) + PAD for k, h in COLS]
wA = labw + sum(colw); x0A = 0.02
SPD = 0.62                                               # speedup column right of B
x0B = x0A + wA + 0.12; wB = W - x0B - SPD - 0.02
yA = H - 0.06 - hA
axA = fig.add_axes([x0A / W, yA / H, wA / W, hA / H]); axB = fig.add_axes([x0B / W, yA / H, wB / W, hA / H])
for ax in (axA, axB):
    ax.set_ylim(n, -HEAD)
    for i in range(0, n, 2): ax.axhspan(i, i + 1, color=STRIPE, lw=0, zorder=0)

# --- A: the shared values
axA.set_xlim(0, wA); axA.axis("off")
right = labw + np.cumsum(colw) - PAD / 2                # right edge of each column's text
for j, (k, h) in enumerate(COLS):
    axA.text(right[j], -0.2, h, ha="right", va="bottom", color=MUTE if k == "nC" else INK, fontsize=6.5, linespacing=1.05)
    for i, m in enumerate(REF):
        v = nC[m] if k == "nC" else D["kimonu_ref"][m][k]
        axA.text(right[j], i + 0.5, str(v), ha="right", va="center", color=MUTE if k == "nC" else INK, fontsize=7)
for i, m in enumerate(REF):
    axA.text(PAD / 2, i + 0.5, SHORT.get(m, m), ha="left", va="center", fontsize=7, color=INK)
x1 = labw + colw[0] + PAD / 2                            # the spanner covers the six compared quantities
axA.plot([x1, wA - PAD / 4], [-1.95, -1.95], color=INK, lw=0.6, clip_on=False)
axA.text((x1 + wA) / 2, -2.05, "identical in Kimonu and the exhaustive reference", ha="center", va="bottom", fontsize=6.5, color=INK)
axA.plot([0, wA], [n, n], color=RULE, lw=0.6, clip_on=False)

# --- B: the same computation, timed
axB.set_xscale("log"); axB.set_xlim(0.15, 2.5e6)
axB.vlines(list(T), 0, n, color=RULE, lw=0.5, zorder=1)      # gridlines stop at the header, clear of the legend
ratios = []
for i, m in enumerate(REF):
    k, r = seconds("kimonu_ref", m), seconds("ref_exh", m); ratios.append(r / k)
    axB.plot([k, r], [i + 0.5] * 2, color="#aab2bb", lw=0.9, zorder=2)
    mark(axB, k, i + 0.5, "kimonu_ref", s=18); mark(axB, r, i + 0.5, "ref_exh", s=16)
    e = int(np.floor(np.log10(r / k)))
    axB.text(1.0 + 0.08 / wB, i + 0.5, rf"${r / k / 10 ** e:.1f}\times10^{{{e}}}$", transform=axB.get_yaxis_transform(),
             ha="left", va="center", fontsize=7, color=INK, clip_on=False)
axB.text(1.0 + 0.08 / wB, -0.2, "Speedup", transform=axB.get_yaxis_transform(), ha="left", va="bottom", fontsize=6.5, color=INK)
axB.set_xticks(list(T)); axB.set_xticklabels(list(T.values())); axB.minorticks_off()
axB.set_yticks([]); axB.tick_params(axis="x", colors=INK, length=2)
for s in ("top", "right", "left"): axB.spines[s].set_visible(False)
axB.spines["bottom"].set_color(MUTE)
axB.set_xlabel("wall time", labelpad=2)
axB.legend(handles=[Line2D([], [], ls="", marker="o", ms=4.5, mfc=GREEN, mec="white", mew=0.5, label="Kimonu, reference-equivalent mode"),
                    Line2D([], [], ls="", marker="s", ms=4, mfc=BLUE, mec="white", mew=0.5, label="Reference, exhaustive")],
           loc="upper left", bbox_to_anchor=(0.06, 1.0), frameon=False, handletextpad=0.2, borderaxespad=0.05, labelspacing=0.25)

# --- C: every model, every arm
hC, yC = 2.3, 0.42
axC = fig.add_axes([0.5 / W, yC / H, 3.75 / W, hC / H])
for y, lab in {60: "1 min", 3600: "1 h", 86400: "1 day", 604800: "1 week"}.items():
    axC.axhline(y, color=RULE, lw=0.5, zorder=1); axC.text(1.9e4, y * 1.18, lab, fontsize=6, color=MUTE, ha="right", va="bottom")
handles = []
for arm, (lab, col, mk, filled, ls) in ARMS.items():
    X = np.array([nC[m] for m in D["models"] if seconds(arm, m) is not None], float)
    Y = np.array([seconds(arm, m) for m in D["models"] if seconds(arm, m) is not None], float)
    b, a0 = np.polyfit(np.log10(X), np.log10(Y), 1); gx = np.linspace(np.log10(X.min()), np.log10(X.max()), 50)
    axC.plot(10 ** gx, 10 ** (a0 + b * gx), color=col, lw=1.1, ls=ls, zorder=2)
    mark(axC, X, Y, arm, s=14 if mk != "^" else 16)
    handles.append(Line2D([], [], color=col, lw=1.1, ls=ls, marker=mk, ms=4 if mk != "^" else 4.5, mfc=col if filled else "white",
                          mec="white" if filled else col, mew=0.5 if filled else 0.9, label=f"{lab}\n{len(X)} models, slope {b:.2f}"))
axC.set_xscale("log"); axC.set_yscale("log"); axC.set_xlim(900, 2e4); axC.set_ylim(0.15, 2e6)
axC.set_xticks([1000, 2000, 5000, 10000, 20000]); axC.set_xticklabels(["1000", "2000", "5000", "10 000", "20 000"]); axC.minorticks_off()
axC.set_yticks([1, 10, 100, 1e3, 1e4, 1e5, 1e6]); axC.set_yticklabels(["1", "10", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$", r"$10^6$"])
axC.set_xlabel("complexes", labelpad=2); axC.set_ylabel("wall time (s)", labelpad=2)
for s in ("top", "right"): axC.spines[s].set_visible(False)
axC.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.04, 1.0), frameon=False, handlelength=2.2, labelspacing=0.7)
fig.text(4.5 / W, (yC + 0.02) / H, "Lines: least-squares fit of log wall time on\nlog size per arm. Stochastic arms: median of\nthree seeds; Kimonu: median of three replicates.\nCOCOA.jl did not finish iS_1188 and iSDY_1059\n(both settings) or iAF1260b (efficient=false)\nwithin 7 days.",
         fontsize=6, color=MUTE, va="bottom", linespacing=1.3)

for L, x, y in (("A", x0A, H - 0.06), ("B", x0B, H - 0.06), ("C", x0A, yC + hC + 0.2)):
    fig.text(x / W, y / H, L, fontsize=10, fontweight="bold", va="top", ha="left", color=INK)

os.makedirs(out, exist_ok=True)
for ext in ("svg", "pdf", "png"):
    fig.savefig(os.path.join(out, f"Figure3_draft.{ext}"), dpi=300)
# --- supplementary figure: concordance stage against the rest, per tool, over the eleven models of A
hD, yD, HS = 1.0, 0.36, 1.5
figS = plt.figure(figsize=(W, HS))
axD = figS.add_axes([1.75 / W, yD / HS, 3.6 / W, hD / HS])
ROWS = ["ref_exh", "cocoa_false", "cocoa_true", "kimonu_ref"]
axD.set_xscale("log"); axD.set_xlim(0.05, 2.5e6); axD.set_ylim(len(ROWS) - 0.4, -0.6)
for y, lab in {1: "1 s", 60: "1 min", 3600: "1 h", 86400: "1 day", 604800: "1 week"}.items():
    axD.axvline(y, color=RULE, lw=0.5, zorder=1)
for i, arm in enumerate(ROWS):
    lab, col, mk, filled, _ = ARMS[arm]
    S = np.array([stages(arm, m) for m in REF if stages(arm, m) is not None])
    for k, dy, fill in ((0, -0.14, True), (1, 0.14, False)):
        lo, mid, hi = S[:, k].min(), np.median(S[:, k]), S[:, k].max()
        axD.plot([lo, hi], [i + dy] * 2, color=col, lw=1.0, alpha=0.6, solid_capstyle="round", zorder=2)
        axD.scatter(mid, i + dy, s=18 if mk != "^" else 20, marker=mk, facecolor=col if fill else "white", edgecolor="white" if fill else col,
                    linewidths=0.5 if fill else 0.9, zorder=3)
    axD.text(0.035, i, lab.replace(", reference-equivalent mode", ""),
             ha="right", va="center", fontsize=6.5, color=INK)
axD.set_yticks([]); axD.minorticks_off()
axD.set_xticks([1, 60, 3600, 86400, 604800]); axD.set_xticklabels(["1 s", "1 min", "1 h", "1 day", "1 week"])
for sp in ("top", "right", "left"): axD.spines[sp].set_visible(False)
axD.spines["bottom"].set_color(MUTE); axD.tick_params(axis="x", colors=INK, length=2)
axD.set_xlabel("wall time per stage (range over the 11 models of Fig. 3A, median marked)", labelpad=2)
axD.legend(handles=[Line2D([], [], ls="", marker="o", ms=4.5, mfc=INK, mec="white", mew=0.5, label="concordance"),
                    Line2D([], [], ls="", marker="o", ms=4.5, mfc="white", mec=INK, mew=0.9, label="kinetic modules, ACR/ACRR")],
           loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, handletextpad=0.2, labelspacing=0.4)

supp = os.path.join(os.path.dirname(os.path.abspath(out)), "supp")
os.makedirs(supp, exist_ok=True)
for ext in ("svg", "png"):
    figS.savefig(os.path.join(supp, f"FigS_stage_timing.{ext}"), dpi=300)
print(f"speedup over the exhaustive reference: {min(ratios):.3g} to {max(ratios):.3g}, median {np.median(ratios):.3g} (11 models)")
