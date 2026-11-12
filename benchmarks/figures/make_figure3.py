"""Figure 4 of the paper (A to D, E and F through make_fig3_methods.py) and figs. S1, S2: Kimonu's relations against the
published tools, and the reference reproduced in seconds. The file name dates from when this was Figure 3.

A  wall time vs network size on all 33 models for every arm, with a log-log least-squares line per arm.
B  the tangent route at scale: the published symbolic local-ACR criterion against one exact nullspace, on the small
   networks and the size ladder (symbolic and exact) and the 33 models (exact, detached uptakes). Symbolic times from
   kinetic_variants/local_acr_symbolic.py ($KIMONU_LACR), exact times from kinetic_variants/timing_exclusive.sbatch
   ($KIMONU_TIMING, median of the replicates).
C  coupled complex pairs, D independent ACR and ACRR relations (k - 1 per ACRR group) among free metabolites, per arm on the eleven
   models where the exhaustive reference finished (bars: value or median over seeds, whiskers: seed range);
   centred titles, "count" on each y axis and one legend below both, as Figure 3A, B.
E, F  the map of relations and methods, and what each tool computes on iLJ478 (make_fig3_methods.py).
fig. S2 (figures/supp/FigS_reference_agreement): the key quantities of Kimonu's reference-equivalent
   mode and the exhaustive reference on the eleven models (identical, asserted below) and the wall time of that
   same computation, with the speedup.
fig. S1 (figures/supp/FigS_stage_timing): where the time goes on those eleven models: the concordance stage against everything after it (kinetic
   modules, ACR/ACRR), per tool, as the range over the models with the median marked. Measured vs derived:
   Kimonu's stage times are measured (median of three replicates); the reference's concordance stage is its
   own measured time, and its R stage the wall time between the chain log's start stamp and the end of its
   ACR log (the R stage reads the stoichiometric coupling from a file, whose computation is not included);
   COCOA's concordance stage is its own logged "Elapsed time", which ends before kinetic_analysis is called,
   and the rest is the measured total minus that.

Kimonu here is run_kimonu_refmode.py ($KIMONU_REF_DIR through figdata.json, the median of three runs), which follows the
reference's definitions: concordance from ker(N), upstream sets merged by overlap only, ACR/ACRR by the
reference's pairwise rule (the Shinar-Feinberg criterion within kinetic modules). Stochastic arms show the median over three seeds. Layout is in inches at 7.2 in (518.4 pt) width; the
supplementary table columns are sized from their measured text, so its axes are placed by hand.

    python figures/make_figure3.py [out_dir]    # run from benchmarks/ after gather_figdata.py and pair_sets.py; default figures/main

KIMONU_TIMING (default the paper's run, /scratch/sadreyev_lab/pschneider/kimonu_rerun_1111/timing) holds the tangent-route
times of B, KIMONU_LACR the symbolic baselines, KIMONU_REF_DIR the runs of the reference-equivalent mode for fig. S1.
"""
import datetime as dt, glob, json, os, re, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_finish, svg_panels                             # live text, embedded fonts, no clipping; one group per panel

out = sys.argv[1] if len(sys.argv) > 1 else "figures/main"
D = json.load(open("figdata.json")); nC = D["n_complexes"]
REF = sorted(D["ref_exh"], key=lambda m: nC[m])

# the agreement the figure claims, checked on every quantity the reference reports (free and all metabolites)
KEYS = ["balanced", "conc", "n_kinetic", "giant", "acr_free", "acrr_free", "acr", "acrr", "acrr_ratios_free", "acrr_ratios"]
for m in REF:
    for k in KEYS:
        assert D["kimonu_ref"][m][k] == D["ref_exh"][m][k], (m, k, D["kimonu_ref"][m][k], D["ref_exh"][m][k])

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "legend.fontsize": 6, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                     **svg_finish.mathtext_rc()})               # words Arial, formulas Latin Modern (LaTeX), both live text
CITE = {"ref_exh": "Langary et al., 2025", "ref_sampled": "Langary et al., 2025"}   # the publication behind each reference arm
LETTER = dict(fontsize=18, fontweight="bold", family="Arial")                 # panel letters as in Figures 1 and 2 (18 pt bold)
INK, MUTE, RULE, STRIPE = "#1b2430", "#6b7480", "#d5d9de", "#f1f3f5"
GREEN, BLUE, RED = "#4DAF4A", "#377EB8", "#E41A1C"          # Set1, as in the earlier drafts; shapes carry identity too
# arm: label, colour, marker, filled, fit line style
ARMS = {"kimonu_det":  ("Kimonu (iterated route)",           GREEN, "o", True, "-"),
        "kimonu_ref":  ("Kimonu, reference-equivalent mode", GREEN, "o", False, "--"),
        "ref_exh":     ("Reference, exhaustive",             BLUE,  "s", True, "-"),
        "ref_sampled": ("Reference, sampled",                BLUE,  "s", False, "--"),
        "cocoa_false": ("COCOA.jl, efficient=false",         RED,   "^", True, "-"),
        "cocoa_true":  ("COCOA.jl, efficient=true",          RED,   "^", False, "--")}
from names import NAME as SHORT                    # one display name per network (figures/names.py)
COLS = [("nC", "Complexes"), ("balanced", "Balanced\ncomplexes"), ("conc", "Concordance\nmodules"),
        ("n_kinetic", "Kinetic\nmodules"), ("giant", "Largest\nmodule"), ("acr_free", "ACR\nmetabolites"), ("acrr_ratios_free", "ACRR\nratios")]
T = {1: "1 s", 60: "1 min", 3600: "1 h", 86400: "1 day"}

def seconds(arm, m):
    v = D[arm].get(m)
    if v is None: return None
    return float(np.median([r["seconds"] for r in v])) if isinstance(v, list) else v["seconds"]

def stages(arm, m):
    """(concordance stage, everything after it) in seconds for one model, or None; see the docstring for
    which part is measured and which derived."""
    if arm == "kimonu_ref":
        R = [json.load(open(f)) for f in sorted(glob.glob(f"{os.environ.get('KIMONU_REF_DIR', '/scratch/sadreyev_lab/pschneider/kimonu_rerun_1111/timing/ref')}/kimonu_{m}_rep*.json"))]
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

def save_finished(f, stem, pdf=True, panels=None):
    """SVG with one group per panel (svg_panels), live text, fonts embedded and no clipping (svg_finish), plus PDF
    (fonts embedded) and PNG for checking; the vault keeps only the SVG. Fails if removing the clipping would
    change the drawing, i.e. if a mark reaches outside its axes. panels: [(name, [(x0, y0, x1, y1) in pt])]."""
    import io, tempfile
    buf = io.StringIO(); f.savefig(buf, format="svg"); raw = buf.getvalue()
    if os.environ.get("RAW_SVG"): open(os.environ["RAW_SVG"], "w").write(raw)
    n, box = svg_finish.check(raw, raw, tempfile.mkdtemp())
    assert n == 0, f"{os.path.basename(stem)}: removing the clipping changes {n} pixels in {box} (pt)"
    fg = svg_panels.Figure(raw)
    miss = fg.assign(panels or [("all", [(0, 0, 1e4, 1e4)])])
    assert not miss, f"{os.path.basename(stem)}: {len(miss)} elements outside every panel"
    svg, used = svg_finish.finish(fg.svg())
    open(stem + ".svg", "w", encoding="utf-8").write(svg)
    if pdf: f.savefig(stem + ".pdf")
    f.savefig(stem + ".png", dpi=300)
    return used

# --- layout, in inches: A scaling and B tangent route (top), C and D results per method (middle), E map and F computations (bottom, make_fig3_methods)
W, RH, HEAD, PAD = 7.2, 0.2, 3.2, 0.1                  # width, table row height, header rows, cell padding
n = len(REF)
hTab = (n + HEAD) * RH
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_fig3_methods import draw as draw_methods, HCD
hS, hBar, gapSB, gapBC = 2.05, 1.05, 1.26, 0.8      # A, B: axis label and the legend of A below; C, D: titles above, tick labels and the shared legend below
xA, wS, xB, wB = 0.5, 3.0, 4.2, 2.9                # the two scaling panels
H = 0.3 + hS + gapSB + hBar + gapBC + HCD
fig = plt.figure(figsize=(W, H)); rend = fig.canvas.get_renderer()
def width_in(s, size):
    t = fig.text(0, 0, s, fontsize=size); w = t.get_window_extent(rend).width / fig.dpi; t.remove(); return w

# --- A: every model, every arm
yS = H - 0.3 - hS
axS = fig.add_axes([xA / W, yS / H, wS / W, hS / H])
for y, lab in {60: "1 min", 3600: "1 h", 86400: "1 day", 604800: "1 week"}.items():
    axS.axhline(y, color=RULE, lw=0.5, zorder=1); axS.text(2.6e4, y * 1.18, lab, fontsize=6, color=MUTE, ha="right", va="bottom")
handles = []
A_ORDER = ["ref_exh", "ref_sampled", "cocoa_false", "cocoa_true", "kimonu_det", "kimonu_ref"]   # legend top to bottom, as the lines run
for arm in A_ORDER:
    lab, col, mk, filled, ls = ARMS[arm]
    X = np.array([nC[m] for m in D["models"] if seconds(arm, m) is not None], float)
    Y = np.array([seconds(arm, m) for m in D["models"] if seconds(arm, m) is not None], float)
    b, a0 = np.polyfit(np.log10(X), np.log10(Y), 1); gx = np.linspace(np.log10(X.min()), np.log10(X.max()), 50)
    axS.plot(10 ** gx, 10 ** (a0 + b * gx), color=col, lw=1.1, ls=ls, zorder=2)
    mark(axS, X, Y, arm, s=14 if mk != "^" else 16)
    handles.append(Line2D([], [], color=col, lw=1.1, ls=ls, marker=mk, ms=4 if mk != "^" else 4.5, mfc=col if filled else "white",
                          mec="white" if filled else col, mew=0.5 if filled else 0.9, label=f"{lab}\n" + (f"({CITE[arm]})\n" if arm in CITE else "") + f"{len(X)} networks, slope {b:.2f}"))
axS.set_xscale("log"); axS.set_yscale("log"); axS.set_xlim(900, 2.7e4); axS.set_ylim(0.05, 2e6)
axS.set_xticks([1000, 2000, 5000, 10000, 20000]); axS.set_xticklabels(["1000", "2000", "5000", "10 000", "20 000"]); axS.minorticks_off()
axS.set_yticks([0.1, 1, 10, 100, 1e3, 1e4, 1e5, 1e6]); axS.set_yticklabels(["0.1", "1", "10", r"$10^2$", r"$10^3$", r"$10^4$", r"$10^5$", r"$10^6$"])
axS.set_xlabel("complexes", labelpad=2); axS.set_ylabel("wall time (s)", labelpad=2)
for s in ("top", "right"): axS.spines[s].set_visible(False)
# one column per tool: filled arm above, open arm below (handles are in that order)
for h_ in handles:
    h_.set_label(h_.get_label().replace("Kimonu, reference-equivalent mode", "Kimonu, reference-\nequivalent mode"))
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.06 / W, (yS - 0.36) / H), ncol=3, frameon=False, handlelength=2.1,
           labelspacing=0.45, columnspacing=0.55, borderaxespad=0, handletextpad=0.4)

# --- B: the tangent route at scale, symbolic criterion against one exact nullspace
LACR = os.environ.get("KIMONU_LACR", "/scratch/sadreyev_lab/pschneider/kimonu_local_acr")
TIMING = os.environ.get("KIMONU_TIMING", "/scratch/sadreyev_lab/pschneider/kimonu_rerun_1111/timing")
LIMIT = 3600
PURPLE, ORANGE = "#984EA3", "#FF7F00"
def med_times(arm, stem=""):
    runs = {}
    for f in glob.glob(f"{TIMING}/{arm}/tangent_*_rep*.json"):
        r = json.load(open(f)); runs.setdefault(r["network"], []).append(r)
    return {k: (v[0]["complexes"], float(np.median([x["seconds"] for x in v]))) for k, v in runs.items()}
small, genome = med_times("small"), med_times("tan")
rows = json.load(open(f"{LACR}/small.json")) + [json.load(open(f)) for f in sorted(glob.glob(f"{LACR}/ladder_*.json"))]
axB = fig.add_axes([xB / W, yS / H, wB / W, hS / H])
for y, lab in ((1, "1 s"), (60, "1 min")):
    axB.axhline(y, color=RULE, lw=0.5, zorder=1); axB.text(1.5, y * 1.25, lab, fontsize=6, color=MUTE, ha="left", va="bottom")
axB.axhline(LIMIT, color=RULE, lw=0.5, zorder=1); axB.text(1.5, LIMIT * 1.25, "1 h limit (open markers did not finish)", fontsize=6, color=MUTE, ha="left", va="bottom")
b_handles = []
for key, col, mk, lab in (("minors", PURPLE, "s", "Symbolic minors (Pascual-Escudero & Feliu, 2022)"),
                          ("rref", ORANGE, "^", "Symbolic elimination (same rank condition, this work)")):
    done = [(r["complexes"], r[key]["seconds"]) for r in rows if r.get(key)]
    out_ = [(r["complexes"], LIMIT) for r in rows if not r.get(key) and r.get(f"{key}_timeout") is not None]
    x, y = zip(*done); axB.scatter(x, np.maximum(y, 1e-3), s=14, color=col, marker=mk, lw=0, zorder=3)
    if out_:
        x, y = zip(*out_); axB.scatter(x, y, s=16, facecolors="white", edgecolors=col, marker=mk, lw=0.8, zorder=3)
    b_handles.append(Line2D([], [], ls="", marker=mk, ms=4, color=col, mew=0, label=lab))
exact = [small[r["network"]] for r in rows if r["network"] in small] + list(genome.values())
x, y = zip(*exact); axB.scatter(x, np.maximum(y, 1e-3), s=14, color=GREEN, marker="o", edgecolor="white", linewidths=0.3, zorder=4)
b_handles.append(Line2D([], [], ls="", marker="o", ms=4, color=GREEN, mew=0, label=f"One exact nullspace (Kimonu, tangent route)\nat most {max(t for _, t in genome.values()):.0f} s on {len(genome)} networks, one core"))
axB.set_xscale("log"); axB.set_yscale("log"); axB.set_xlim(1.4, 2.7e4); axB.set_ylim(6e-4, LIMIT * 6); axB.minorticks_off()
axB.set_xticks([10, 100, 1000, 10000]); axB.set_xticklabels(["10", r"$10^2$", r"$10^3$", r"$10^4$"])
axB.set_yticks([1e-3, 1e-2, 0.1, 1, 10, 100, 1e3]); axB.set_yticklabels(["0.001", "0.01", "0.1", "1", "10", r"$10^2$", r"$10^3$"])
axB.set_xlabel("complexes", labelpad=2); axB.set_ylabel("wall time (s)", labelpad=2)
for s_ in ("top", "right"): axB.spines[s_].set_visible(False)
fig.legend(handles=b_handles, loc="upper left", bbox_to_anchor=((xB - 0.2) / W, (yS - 0.36) / H), frameon=False, handletextpad=0.3,
           labelspacing=0.45, borderaxespad=0)

# --- C, D: what each method reports, on the eleven models where the exhaustive reference finished
BARS = ["cocoa_false", "cocoa_true", "kimonu_det", "kimonu_ref", "ref_exh", "ref_sampled"]   # reference-equivalent next to the reference it matches
PC = json.load(open("figures/pair_counts.json"))          # figures/pair_sets.py
def counts(arm, m, k):
    """(bar, low, high): the value, or for stochastic arms the median and range over seeds"""
    v = PC[m].get(arm)
    if v is None or v == []: return None
    if isinstance(v, list):
        x = sorted(r[k] for r in v); return x[len(x) // 2], x[0], x[-1]
    return v[k], v[k], v[k]
yBar = yS - gapSB - hBar
wGap = 0.75                                               # room for C's own y label
wBar = (W - 0.5 - wGap - 0.05) / 2
axK = fig.add_axes([0.5 / W, yBar / H, wBar / W, hBar / H])
axR = fig.add_axes([(0.5 + wBar + wGap) / W, yBar / H, wBar / W, hBar / H])
bw = 0.8 / len(BARS)
from matplotlib.colors import to_rgb
tint = lambda c, t=0.6: tuple(t + (1 - t) * v for v in to_rgb(c))   # open arms: a light fill of the arm colour
YLIM = {"pairs": (40, 40000), "fixed": (0.5, 200)}
for ax, k, title in ((axK, "pairs", "coupled complex pairs"), (axR, "fixed", "fixed relations, ACR and ACRR (free metabolites)")):
    for j, arm in enumerate(BARS):
        lab, col, mk, filled, _ = ARMS[arm]
        for i, m in enumerate(REF):
            c = counts(arm, m, k)
            if c is None or c[0] <= 0: continue
            x = i - 0.4 + (j + 0.5) * bw
            lo, hi = YLIM[k]                       # bars start at the axis floor and whiskers stop at its top, so no
            b0 = lo * 1.03                         # half an edge width above the floor: the edge stays on the spine
            ax.bar(x, c[0] - b0, bottom=b0, width=bw * 0.92, color=col if filled else tint(col), edgecolor=col,   # clipping is needed
                   lw=0.0 if filled else 0.6, zorder=2)
            if c[2] > c[1]:
                ax.plot([x, x], [max(c[1], 0.8, lo), min(c[2], hi)], color=INK, lw=0.5, zorder=3, solid_capstyle="butt")
                if c[2] > hi: print(f"whisker cut at the axis top: {arm} {m} {k} {c}")
    ax.set_yscale("log"); ax.set_xlim(-0.55, n - 0.45); ax.minorticks_off()
    ax.set_xticks(range(n)); ax.set_xticklabels([SHORT.get(m, m) for m in REF],
                                                rotation=40, ha="right", rotation_mode="anchor", fontsize=6)
    ax.tick_params(axis="x", length=0, pad=2); ax.tick_params(axis="y", colors=INK, length=2)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.yaxis.grid(True, color=RULE, lw=0.4, zorder=0); ax.set_axisbelow(True)
    ax.set_title(title, fontsize=8, color=INK, loc="center", pad=4)
    ax.set_ylabel("count", labelpad=2)
axK.set_ylim(*YLIM["pairs"]); axR.set_ylim(*YLIM["fixed"])
axK.set_yticks([100, 1000, 10000]); axK.set_yticklabels(["100", "1000", "10 000"])
axR.set_yticks([1, 10, 100]); axR.set_yticklabels(["1", "10", "100"])
for ax in (axK, axR): ax.tick_params(axis="y", labelsize=6)
from matplotlib.patches import Patch
bar_handles = [Patch(facecolor=ARMS[a][1] if ARMS[a][3] else tint(ARMS[a][1]), edgecolor=ARMS[a][1], lw=0 if ARMS[a][3] else 0.6, label=ARMS[a][0] + (f" ({CITE[a]})" if a in CITE else "")) for a in BARS]
yLeg = yBar - 0.45                                       # the legend of C and D sits below their tick labels, as in Figure 2
fig.legend(handles=bar_handles, loc="upper left", bbox_to_anchor=(0.5 / W, yLeg / H), ncol=3, frameon=False,
           fontsize=6, handlelength=1.2, handleheight=0.8, columnspacing=1.2, labelspacing=0.25, borderaxespad=0)

# --- E, F: the map and the computations per tool
axCD = fig.add_axes([0, 0, 1, HCD / H]); axCD.set_xlim(0, W); axCD.set_ylim(0, HCD); axCD.axis("off")
math_texts = draw_methods(axCD, W, HCD)

# --- supplementary figure (figures/supp/FigS_reference_agreement): the shared values and the same computation, timed
HT = hTab + 0.55
figT = plt.figure(figsize=(W, HT))
labw = max(width_in(SHORT.get(m, m), 7) for m in REF) + PAD
colw = [max(width_in(h, 6.5), *(width_in(str(nC[m] if k == "nC" else D["kimonu_ref"][m][k]), 7) for m in REF)) + PAD for k, h in COLS]
wA = labw + sum(colw); x0A = 0.02
SPD = 0.62                                               # speedup column right of D
x0B = x0A + wA + 0.12; wB = W - x0B - SPD - 0.02
yA = 0.42
axA = figT.add_axes([x0A / W, yA / HT, wA / W, hTab / HT]); axB = figT.add_axes([x0B / W, yA / HT, wB / W, hTab / HT])
for ax in (axA, axB):
    ax.set_ylim(n, -HEAD)
    for i in range(0, n, 2): ax.axhspan(i, i + 1, color=STRIPE, lw=0, zorder=0)
axA.set_xlim(0, wA); axA.axis("off")
right = labw + np.cumsum(colw) - PAD / 2                # right edge of each column's text
for j, (k, h) in enumerate(COLS):
    axA.text(right[j], -0.2, h, ha="right", va="bottom", color=MUTE if k == "nC" else INK, fontsize=6, linespacing=1.05)
    for i, m in enumerate(REF):
        v = nC[m] if k == "nC" else D["kimonu_ref"][m][k]
        axA.text(right[j], i + 0.5, str(v), ha="right", va="center", color=MUTE if k == "nC" else INK, fontsize=7)
for i, m in enumerate(REF):
    axA.text(PAD / 2, i + 0.5, SHORT.get(m, m), ha="left", va="center", fontsize=7, color=INK)
x1 = labw + colw[0] + PAD / 2                            # the spanner covers the six compared quantities
axA.plot([x1, wA - PAD / 4], [-1.95, -1.95], color=INK, lw=0.6, clip_on=False)
axA.text((x1 + wA) / 2, -2.05, "identical in Kimonu's reference-equivalent mode and the exhaustive\nreference, on the 11 networks where the exhaustive reference finished", ha="center", va="bottom", fontsize=6, color=INK, linespacing=1.1)
axA.plot([0, wA], [n, n], color=RULE, lw=0.6, clip_on=False)

axB.set_xscale("log"); axB.set_xlim(0.05, 2.5e6)
axB.vlines(list(T), 0, n, color=RULE, lw=0.5, zorder=1)      # gridlines stop at the header, clear of the legend
ratios = []
for i, m in enumerate(REF):
    k, r = seconds("kimonu_ref", m), seconds("ref_exh", m); ratios.append(r / k)
    axB.plot([k, r], [i + 0.5] * 2, color="#aab2bb", lw=0.9, zorder=2)
    mark(axB, k, i + 0.5, "kimonu_ref", s=18); mark(axB, r, i + 0.5, "ref_exh", s=16)
    e = int(np.floor(np.log10(r / k)))
    axB.text(1.0 + 0.08 / wB, i + 0.5, rf"${r / k / 10 ** e:.1f}\times10^{{{e}}}$", transform=axB.get_yaxis_transform(),
             ha="left", va="center", fontsize=7, color=INK, clip_on=False)
axB.text(1.0 + 0.08 / wB, -0.2, "Speed-up", transform=axB.get_yaxis_transform(), ha="left", va="bottom", fontsize=6, color=INK)
axB.set_xticks(list(T)); axB.set_xticklabels(list(T.values())); axB.minorticks_off()
axB.set_yticks([]); axB.tick_params(axis="x", colors=INK, length=2)
for s in ("top", "right", "left"): axB.spines[s].set_visible(False)
axB.spines["bottom"].set_color(MUTE)
axB.set_xlabel("wall time", labelpad=2)
axB.legend(handles=[Line2D([], [], ls="", marker="o", ms=4.5, mfc="white", mec=GREEN, mew=0.9, label="Kimonu, reference-equivalent mode"),
                    Line2D([], [], ls="", marker="s", ms=4, mfc=BLUE, mec="white", mew=0.5, label=f"Reference, exhaustive ({CITE['ref_exh']})")],
           loc="upper left", bbox_to_anchor=(0.06, 1.0), frameon=False, handletextpad=0.2, borderaxespad=0.05, labelspacing=0.25)

for L, x, y in (("A", 0.0, H - 0.02), ("B", xB - 0.5, H - 0.02), ("C", 0.0, yBar + hBar + 0.3), ("D", wBar + wGap, yBar + hBar + 0.3)):
    fig.text(x / W, y / H, L, va="top", ha="left", color=INK, **LETTER)
os.makedirs(out, exist_ok=True)
Hpt, yAB, yCD = H * 72, (H - yBar - hBar - 0.3) * 72, (H - HCD) * 72  # panel bounds in pt from the top: below the legends of A, B; above E, F
xAB = (xB - 0.5) * 72                                                     # between A and B
xBC, yL = (0.5 + wBar + 0.15) * 72, (H - yLeg + 0.03) * 72              # between the two bar charts; top of their legend
save_finished(fig, os.path.join(out, "Figure3_draft"), panels=[
    ("A", [(0, 0, xAB, yAB)]), ("B", [(xAB, 0, 518.4, yAB)]), ("C", [(0, yAB, xBC, yL)]), ("D", [(xBC, yAB, 518.4, yL)]),
    ("legend_CD", [(0, yL, 518.4, yCD)]), ("E", [(0, yCD, 225, Hpt)]), ("F", [(225, yCD, 518.4, Hpt)])])
supp = os.path.join(os.path.dirname(os.path.abspath(out)), "supp"); os.makedirs(supp, exist_ok=True)
save_finished(figT, os.path.join(supp, "FigS_reference_agreement"))
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
             ha="right", va="center", fontsize=6, color=INK)
axD.set_yticks([]); axD.minorticks_off()
axD.set_xticks([1, 60, 3600, 86400, 604800]); axD.set_xticklabels(["1 s", "1 min", "1 h", "1 day", "1 week"])
for sp in ("top", "right", "left"): axD.spines[sp].set_visible(False)
axD.spines["bottom"].set_color(MUTE); axD.tick_params(axis="x", colors=INK, length=2)
axD.set_xlabel("wall time per stage (range over the 11 networks of Fig. 4C, median marked)", labelpad=2)
axD.legend(handles=[Line2D([], [], ls="", marker="o", ms=4.5, mfc=INK, mec="white", mew=0.5, label="concordance"),
                    Line2D([], [], ls="", marker="o", ms=4.5, mfc="white", mec=INK, mew=0.9, label="kinetic modules, ACR/ACRR")],
           loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, handletextpad=0.2, labelspacing=0.4)

save_finished(figS, os.path.join(supp, "FigS_stage_timing"), pdf=False)
print(f"speedup over the exhaustive reference: {min(ratios):.3g} to {max(ratios):.3g}, median {np.median(ratios):.3g} (11 models)")
