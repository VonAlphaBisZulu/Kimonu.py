"""Tangent-route benchmark: the published local-ACR criterion (symbolic) against one exact nullspace.

A  runtime against network size (species). Symbolic minors = Pascual-Escudero & Feliu's Theorem 5.3(3) as
   stated, symbolic elimination = the same rank condition by Gaussian elimination over Q(v); both are our
   implementations (sympy). Exact nullspace = the tangent route of Kimonu on the small networks, the size
   ladder (k disjoint copies of the Figure 1 network, e_coli_core) and the 33 genome-scale models. Runs that
   reached the 1 h limit are drawn at the limit as open markers.
B  float64 instead of exact arithmetic on the 33 models (detached uptakes): per tolerance for equal rows, how
   many models come out exactly right, only miss relations, only make some up, or both
   (kinetic_variants/float_tradeoff.py).

    python figures/make_tangent_benchmark.py      # from benchmarks/
"""
import glob, json, os
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

LACR = "/scratch/sadreyev_lab/pschneider/kimonu_local_acr"
TF = "/scratch/sadreyev_lab/pschneider/kimonu_tangent_float"
OUT = os.environ.get("FIG_OUT", "figures/main"); LIMIT = 3600
# composed into Figure 2 by make_figure2.py: a shorter panel, no letters (the composition adds them)
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
ax = fig.add_axes([0.07, 0.17, 0.42, 0.74]); bx = fig.add_axes([0.6, 0.17, 0.235, 0.74])
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
ax.axhline(LIMIT, color=RULE, lw=0.6, zorder=1); ax.text(1.05, LIMIT / 3.2, "1 h limit (open: not finished)", color=MUTE, fontsize=5.5)
for sec, lab in ((1, "1 s"), (60, "1 min")):                 # reference lines
    ax.axhline(sec, color=RULE, lw=0.6, ls=(0, (3, 2)), zorder=1); ax.text(1.05, sec * 1.25, lab, color=MUTE, fontsize=5.5)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-3 * 0.6, LIMIT * 6)
ax.set_xlabel("species"); ax.set_ylabel("wall time (s)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(frameon=False, loc="lower right")
if LETTERS: ax.text(-0.13, 1.02, "A", transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
# ---- B: float64 instead of exact needs a tolerance to call two rows equal (or a row zero). Too tight misses
#      relations, too loose makes them up. Per row tolerance (numpy's default rank), the models split into
#      exactly right, only missing, only made up, and both; the zero-error window is narrow where it exists at
#      all, and no single tolerance gets every model right. Data: kinetic_variants/float_tradeoff.py.
FT = os.environ.get("FT_OUT", "/scratch/sadreyev_lab/pschneider/kimonu_float_tradeoff")
trade = [json.load(open(f)) for f in sorted(glob.glob(f"{FT}/*.json"))]
if trade:
    taus = np.array(trade[0]["row_tols"]); n = len(trade)
    cls = np.zeros((4, len(taus)))                     # right, missed only, made up only, both
    for d in trade:
        for k, r in enumerate(d["rank"]["numpy"]["sweep"]):
            miss, fake = r["acr_fn"] + r["acrr_fn"] > 0, r["acr_fp"] + r["acrr_fp"] > 0
            cls[0 if not (miss or fake) else 1 if not fake else 2 if not miss else 3, k] += 1
    RED, BLUE_, BOTH = "#d6604d", "#4393c3", "#9e7bb5"
    bx.stackplot(taus, cls[0], cls[2], cls[1], cls[3], colors=(GREEN, RED, BLUE_, BOTH), alpha=.85, lw=0,   # right at the bottom
                 labels=("all relations right", "relations made up", "relations missed", "missed and made up"))
    bx.set_xscale("log"); bx.set_xlim(taus[0], taus[-1]); bx.set_ylim(0, n)
    best = cls[0].max(); kb = int(np.argmax(cls[0]))
    bx.text(taus[kb] * 1.6, best + 0.4, f"at best {best:.0f} of {n} right", fontsize=6, color=INK, va="bottom")
    bx.set_xlabel("tolerance for equal rows (float64, default rank)"); bx.set_ylabel("models")
    for s_ in ("top", "right"): bx.spines[s_].set_visible(False)
    bx.legend(frameon=False, loc="upper left", fontsize=6, handlelength=1, handletextpad=.4, borderaxespad=0,
              bbox_to_anchor=(1.03, 1.0), reverse=True)
    bx.text(1.03, 0.0, "too tight: relations missed\ntoo loose: relations made up\nsome models: both at every\ntolerance",
            transform=bx.transAxes, fontsize=5.5, color=MUTE, va="bottom")
if LETTERS: bx.text(-0.12, 1.02, "B", transform=bx.transAxes, fontsize=10, fontweight="bold", va="bottom")
os.makedirs(OUT, exist_ok=True)
for ext in ("svg", "png"):
    fig.savefig(f"{OUT}/{NAME}.{ext}", dpi=300)
print("wrote", OUT)
