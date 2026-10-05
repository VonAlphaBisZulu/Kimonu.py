"""Tangent-route benchmark: the published local-ACR criterion (symbolic) against one exact nullspace.

A  runtime against network size (species). Symbolic minors = Pascual-Escudero & Feliu's Theorem 5.3(3) as
   stated, symbolic elimination = the same rank condition by Gaussian elimination over Q(v); both are our
   implementations (sympy). Exact nullspace = the tangent route of Kimonu on the small networks, the size
   ladder (k disjoint copies of the Figure 1 network, e_coli_core) and the 33 genome-scale models. Runs that
   reached the 1 h limit are drawn at the limit as open markers.
B  float64 instead of exact arithmetic on the 33 models (ordered by size): wrong ACR and ACRR calls against the
   exact readout, at numpy's default rank and a row tolerance of 1e-8, and at the best of nine tolerance
   pairs per model ("oracle", chosen knowing the exact answer).

    python figures/make_tangent_benchmark.py      # from benchmarks/
"""
import glob, json, os
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

LACR = "/scratch/sadreyev_lab/pschneider/kimonu_local_acr"
TF = "/scratch/sadreyev_lab/pschneider/kimonu_tangent_float"
OUT = os.environ.get("FIG_OUT", "figures/main"); LIMIT = 3600
# composed into Figure 3 by make_figure3.py: a shorter panel, no letters (the composition adds them)
HEIGHT_CM = float(os.environ.get("FIG_HEIGHT_CM", 6.2)); LETTERS = os.environ.get("FIG_LETTERS", "1") == "1"
NAME = os.environ.get("FIG_NAME", "tangent_benchmark")
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "svg.fonttype": "none", "pdf.fonttype": 42})
INK, MUTE, RULE = "#1b2430", "#6b7480", "#d5d9de"
GREEN, PURPLE, ORANGE = "#4DAF4A", "#984EA3", "#FF7F00"

# ---- data: small networks and the ladder (all three methods), genome scale (exact only)
rows = json.load(open(f"{LACR}/small.json")) + [json.load(open(f)) for f in sorted(glob.glob(f"{LACR}/ladder_*.json"))]
def series(key):
    done = [(r["species"], r[key]["seconds"]) for r in rows if r.get(key)]
    out = [(r["species"], LIMIT) for r in rows if not r.get(key) and r.get(f"{key}_timeout") is not None]
    return done, out
genome = []
for f in glob.glob("kinetic_variants/fast_*.json"):
    M = json.load(open(f))["model"]
    crn = json.load(open(f"unblocked/{M}_fixed.crn.json"))
    genome.append((crn["n_metabolites"], json.load(open(f))["attached"]["seconds"]))
floats = [json.load(open(f)) for f in sorted(glob.glob(f"{TF}/*.json"))]

fig = plt.figure(figsize=(17.5 / 2.54, HEIGHT_CM / 2.54))
ax = fig.add_axes([0.07, 0.17, 0.42, 0.74]); bx = fig.add_axes([0.6, 0.17, 0.38, 0.74])
# ---- A
for key, col, mk, lab in (("minors", PURPLE, "s", "symbolic minors (published criterion)"),
                          ("rref", ORANGE, "^", "symbolic elimination"),
                          ("kimonu", GREEN, "o", "exact nullspace (this work)")):
    done, out = series(key)
    if key == "kimonu": done = done + genome
    if done:
        x, y = zip(*done); ax.scatter(x, np.maximum(y, 1e-3), s=14, color=col, marker=mk, lw=0, zorder=3, label=lab)
    if out:
        x, y = zip(*out); ax.scatter(x, y, s=16, facecolors="none", edgecolors=col, marker=mk, lw=0.8, zorder=3)
ax.axhline(LIMIT, color=RULE, lw=0.6, zorder=1); ax.text(1.05, LIMIT * 1.45, "1 h limit (open: not finished)", color=MUTE, fontsize=5.5)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-3 * 0.6, LIMIT * 6)
ax.set_xlabel("species"); ax.set_ylabel("wall time (s)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(frameon=False, loc="lower right")
if LETTERS: ax.text(-0.13, 1.02, "A", transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
# ---- B: wrong calls against the dynamic range of the interior flux, the source of float trouble
if floats:
    w = lambda v: v["acr_fp"] + v["acr_fn"] + v["acrr_fp"] + v["acrr_fn"]
    span = [d["rate_span_log10"] for d in floats]
    for get, col, mk, lab, dx in ((lambda d: d["float"]["numpy_default|1e-08"], MUTE, "o", "float64, default tolerances", -0.04),
                                  (lambda d: d["oracle"], INK, "D", "float64, best of 9 tolerances", 0.04)):
        bx.scatter(np.array(span) + dx, [w(get(d)) for d in floats], s=12, color=col, marker=mk, lw=0, label=lab, zorder=3)
    bx.set_yscale("symlog", linthresh=1); bx.set_ylim(-0.3, 1e5)
    bx.set_xlabel("dynamic range of interior flux r* (log10 max/min)"); bx.set_ylabel("wrong ACR + ACRR calls per model")
    n_def = sum(w(d["float"]["numpy_default|1e-08"]) > 0 for d in floats); n_or = sum(w(d["oracle"]) > 0 for d in floats)
    bx.text(0.02, 0.97, f"exact nullspace: 0 errors (reference)\nmodels with errors: {n_def}/{len(floats)} default, {n_or}/{len(floats)} best tolerance",
            transform=bx.transAxes, fontsize=6, color=INK, va="top")
    for s_ in ("top", "right"): bx.spines[s_].set_visible(False)
    bx.legend(frameon=False, loc="lower right")
if LETTERS: bx.text(-0.12, 1.02, "B", transform=bx.transAxes, fontsize=10, fontweight="bold", va="bottom")
os.makedirs(OUT, exist_ok=True)
for ext in ("svg", "png"):
    fig.savefig(f"{OUT}/{NAME}.{ext}", dpi=300)
print("wrote", OUT)
