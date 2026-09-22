"""Main-text candidate figures M1-M3 from figdata.json, appended to the figures page.

M1 kinetic modules and M2 ACRR (free) on the ten reference models in order of size; M3 runtime vs
network size on every model, with a log-log least-squares line per arm and its 95% confidence band
for the mean. Run from benchmarks/:  python figures/make_main.py
"""
import json, os, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
sys.path.insert(0, "figures")
import make_prelim as P            # renders F1-F5 on import; reuses D, ARMS, val, grouped

D, ARMS, nC = P.D, P.ARMS, P.nC
MAIN = [m for m in sorted(D["ref_exh"] | D["ref_sampled"], key=lambda m: nC[m]) if m != "iIS312_Epimastigote"][:10]
P.models = MAIN
os.makedirs("figures/main", exist_ok=True)

import textwrap
def save(fig, name, note):
    fig.text(0.01, 0.005, "\n".join(textwrap.wrap(note, 190)), fontsize=7, color="#555", va="bottom"); fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(f"figures/main/{name}.png", dpi=160); fig.savefig(f"figures/main/{name}.svg"); plt.close(fig)

fig, ax = plt.subplots(figsize=(11, 4.6))
P.grouped(ax, "n_kinetic", ["kimonu", "cocoa_false", "cocoa_true", "ref_exh"], "log", "M1  Kinetic modules (non-singleton)", "count (log)",
          labels={"ref_exh": "Reference exhaustive"})
save(fig, "M1_kinetic_modules", P.BASE + "  Reference on the landed exhaustive runs; the three exact arms agree wherever the reference has landed.")

fig, ax = plt.subplots(figsize=(11, 4.6))
P.grouped(ax, "acrr_free", ["kimonu", "kimonu_s210", "cocoa_false", "cocoa_true", "ref_exh"], "symlog", "M2  ACRR metabolite pairs (free)", "pairs (symlog)",
          labels={"ref_exh": "Reference exhaustive †"})
save(fig, "M2_acrr", P.BASE + "  † reference criterion (Corollary S2-10) is a sufficient condition; * Kimonu applying that criterion to its own modules reproduces the reference's pairs exactly on every landed model.")

# M3 runtime vs size with fitted lines and 95% confidence bands
fig, ax = plt.subplots(figsize=(11, 5.4))
arms = [("kimonu", "seconds"), ("cocoa_false", "seconds"), ("cocoa_true", "seconds"), ("ref_exh", "seconds"), ("ref_sampled", "seconds")]
xs_all = np.array([nC[m] for m in D["models"]], float)
for a, key in arms:
    lab, col, _ = ARMS[a]; X, Y = [], []
    for m in D["models"]:
        v = P.val(a, m, key)
        if v is None: continue
        X.append(nC[m]); Y.append(v[0])
        if len(v) > 3: ax.scatter([nC[m]] * len(v[3]), v[3], s=8, color=col, alpha=0.35, linewidths=0, zorder=2)
    if len(X) < 3: continue
    X, Y = np.array(X, float), np.array(Y, float); lx, ly = np.log10(X), np.log10(Y)
    ax.scatter(X, Y, s=26, color=col, edgecolor="white", linewidths=0.4, zorder=3, label=f"{lab} (n={len(X)}, slope {np.polyfit(lx, ly, 1)[0]:.2f})")
    b, a0 = np.polyfit(lx, ly, 1); n = len(lx); resid = ly - (a0 + b * lx); s2 = resid @ resid / max(n - 2, 1)
    gx = np.linspace(lx.min(), lx.max(), 80); se = np.sqrt(s2 * (1 / n + (gx - lx.mean()) ** 2 / ((lx - lx.mean()) ** 2).sum()))
    t = stats.t.ppf(0.975, max(n - 2, 1)); fit = a0 + b * gx
    ax.plot(10 ** gx, 10 ** fit, color=col, lw=1.4, zorder=2); ax.fill_between(10 ** gx, 10 ** (fit - t * se), 10 ** (fit + t * se), color=col, alpha=0.12, linewidths=0)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("complexes"); ax.set_ylabel("seconds")
for y, t_ in [(60, "1 min"), (3600, "1 h"), (86400, "1 day"), (604800, "1 week")]:
    ax.axhline(y, color="#ddd", lw=0.6); ax.text(xs_all.min() * 1.02, y * 1.15, t_, fontsize=6.5, color="#888")
ax.set_title("M3  Runtime vs network size, all models", loc="left", fontsize=11)
ax.legend(fontsize=7, frameon=False, loc="lower right"); ax.grid(color="#eee", lw=0.5); [ax.spines[s].set_visible(False) for s in ("top", "right")]
save(fig, "M3_runtime", "Lines: least-squares fit of log seconds on log complexes per arm, on per-model medians; bands: 95% confidence for the fitted mean. Stochastic arms show seeds as faint points. Reference times are from the fixed bridge; the exhaustive arm is concordance-dominated.")
print("main figures:", sorted(os.listdir("figures/main")), "models:", MAIN)
