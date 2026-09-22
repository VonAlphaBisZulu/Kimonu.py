"""Preliminary figures F1-F5 from figdata.json.

Run from benchmarks/:  python figures/make_prelim.py
Every arm's provenance and trust level is in FINDINGS.md. Error bars are min-max over seeds and
appear only on stochastic arms; exact arms have none by construction.

figdata.json schema:
  models, n_complexes[m]
  kimonu[m]: balanced, conc, giant, giant_overlap, acr, acrr, acr_free, acrr_free, seconds, sec_conc, sec_kin, sec_acr
  cocoa_false[m] / cocoa_true[m]: [{seed, balanced, conc, giant, acr, acrr, acr_free, acrr_free, seconds}, ...]
  ref_exh[m]: balanced, conc, seconds, sec_conc (+ giant, acr, acrr, acr_free, acrr_free where the R stage ran)
  ref_sampled[m]: [{seed, balanced, conc, seconds}, ...]
Free counts: metabolite index < index of "E_1" in the model's own metabolite list.
"""
import json, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch

D = json.load(open("figdata.json")); nC = D["n_complexes"]
models = sorted(D["models"], key=lambda m: nC[m]); os.makedirs("figures/prelim", exist_ok=True)
ARMS = {"kimonu":      ("Kimonu (Upstream, S2-1)",      "#1C6B58", ""),
        "kimonu4":     ("Kimonu, overlap merge only*",   "#7FB8A6", "//"),
        "kimonu_s210": ("Kimonu, S2-10 criterion*",     "#7FB8A6", "//"),
        "cocoa_false": ("COCOA efficient=false",        "#B8641E", ""),
        "cocoa_true":  ("COCOA efficient=true",         "#E0A26A", ""),
        "ref_exh":     ("Reference exhaustive",         "#4A5A88", ""),
        "ref_sampled": ("Reference sampled*",           "#9AA6C8", "//")}
BASE = "Stochastic arms: bar at the median over seeds, each seed a point; exact arms have no spread by construction.  * comparability-only arm."

def val(arm, m, key):
    if arm == "kimonu":  v = D["kimonu"].get(m, {}).get(key); return None if v is None else (v, None, None)
    if arm == "kimonu4": v = D["kimonu"].get(m, {}).get("giant_overlap") if key == "giant" else None; return None if v is None else (v, None, None)
    if arm == "kimonu_s210":
        v = D["kimonu"].get(m, {}).get({"acr_free": "acr_s210_free", "acrr_free": "acrr_s210_free", "acr": "acr_s210", "acrr": "acrr_s210"}.get(key, "-")); return None if v is None else (v, None, None)
    if arm == "ref_exh": v = D["ref_exh"].get(m, {}).get(key); return None if v is None else (v, None, None)
    runs = [r[key] for r in D[arm].get(m, []) if r.get(key) is not None]
    return None if not runs else (float(np.median(runs)), float(min(runs)), float(max(runs)), runs)

def grouped(ax, key, arms, scale="linear", title="", ylabel="count", labels=None):
    labels = labels or {}; n = len(arms); w = 0.8 / n; x = np.arange(len(models)); drawn = {a: 0 for a in arms}
    for i, a in enumerate(arms):
        lab, col, hatch = ARMS[a]
        for j, m in enumerate(models):
            v = val(a, m, key)
            if v is None: continue
            med, lo, hi = v[0], v[1], v[2]; runs = v[3] if len(v) > 3 else None; drawn[a] += 1; xx = x[j] + (i - n / 2 + 0.5) * w
            ax.bar(xx, med, w, color=col, hatch=hatch, edgecolor=col if hatch else "white", linewidth=0.4)
            if runs and len(runs) > 1:
                ax.scatter([xx] * len(runs), runs, s=6, color="#222", zorder=3, linewidths=0)
    ax.set_xticks(x); ax.set_xticklabels(models, rotation=65, ha="right", fontsize=7)
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left", fontsize=11)
    if scale == "log": ax.set_yscale("log")
    if scale == "symlog": ax.set_yscale("symlog", linthresh=1); ax.set_ylim(bottom=0)
    ax.grid(axis="y", color="#e5e5e5", linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(handles=[Patch(facecolor=ARMS[a][1], hatch=ARMS[a][2], edgecolor=ARMS[a][1] if ARMS[a][2] else "white",
                             label=f"{labels.get(a, ARMS[a][0])}  (n={drawn[a]})") for a in arms],
              fontsize=7, frameon=False, ncol=2, loc="upper left")
    return drawn

def save(fig, name, note):
    fig.text(0.01, 0.005, note, fontsize=7, color="#555"); fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(f"figures/prelim/{name}.png", dpi=160); fig.savefig(f"figures/prelim/{name}.svg"); plt.close(fig)

# F1, F2, F3
fig, ax = plt.subplots(figsize=(14, 5)); grouped(ax, "balanced", ["kimonu", "cocoa_false", "cocoa_true", "ref_exh", "ref_sampled"], "log", "F1  Balanced complexes", "count (log)")
save(fig, "F1_balanced", BASE + "  Reference values are from the final rerun (fixed algorithm and bridge); balanced matches Kimonu exactly on every model, exhaustive and sampled.")
fig, ax = plt.subplots(figsize=(14, 5)); grouped(ax, "conc", ["kimonu", "cocoa_false", "cocoa_true", "ref_exh", "ref_sampled"], "log", "F2  Concordance modules (non-singleton)", "count (log)",
        labels={"ref_exh": "Reference exhaustive"})
save(fig, "F2_concordance", BASE + "  Reference concordance is the final rerun: exhaustive is identical to Kimonu on every landed model (4 of 11 so far); sampled (11 models, 3 seeds) is a strict subset of Kimonu's pairs on every cell, as its prefilter can only remove candidates.")
fig, ax = plt.subplots(figsize=(14, 5)); grouped(ax, "giant", ["kimonu", "kimonu4", "cocoa_false", "cocoa_true", "ref_exh"], "log", "F3  Largest kinetic module", "count (log)",
        labels={"ref_exh": "Reference (R stage, corrected NA)*"})
save(fig, "F3_giant", BASE + "  * Kimonu without the Proposition S3-4 coupling merge = the reference's merge rule; the four-phase and two-phase upstream sets are identical on all 33 models. Reference giants: corrected R stage on the final exhaustive concordance (4 models so far), equal to Kimonu on each.")
# F4: ACR linear (mostly zeros), ACRR symlog
fig, axs = plt.subplots(2, 1, figsize=(14, 8.5)); a4 = ["kimonu", "cocoa_false", "cocoa_true", "kimonu_s210", "ref_exh"]; lab4 = {"ref_exh": "Reference exhaustive †"}
grouped(axs[0], "acr_free", a4, "linear", "F4a  ACR metabolites (free)", labels=lab4)
grouped(axs[1], "acrr_free", a4, "symlog", "F4b  ACRR metabolite pairs (free), symlog", labels=lab4)
save(fig, "F4_acr_acrr", BASE + "  † reference ACR/ACRR uses the combinatorial criterion (Corollary S2-10), a subset of Kimonu's kernel test.  * Kimonu applying that criterion to its own modules; equals the reference's ACR/ACRR name sets on every landed model.")
# F5
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [2, 1.2]})
for a in ["kimonu", "cocoa_false", "cocoa_true", "ref_exh", "ref_sampled"]:
    lab, col, _ = ARMS[a]; xs, ys, lo, hi = [], [], [], []
    for m in models:
        v = val(a, m, "sec_conc" if a == "ref_exh" else "seconds")
        if v is None: continue
        xs.append(nC[m]); ys.append(v[0])
        if len(v) > 3: a1.scatter([nC[m]] * len(v[3]), v[3], s=7, color=col, alpha=0.45, linewidths=0, zorder=2)
    if xs: a1.scatter(xs, ys, s=22, color=col, label=f"{lab} (n={len(xs)})", zorder=3, edgecolor="white", linewidths=0.4)
a1.set_xscale("log"); a1.set_yscale("log"); a1.set_xlabel("complexes"); a1.set_ylabel("seconds"); a1.set_title("F5a  Runtime vs network size", loc="left", fontsize=11)
for y, t in [(60, "1 min"), (3600, "1 h"), (86400, "1 day"), (604800, "1 week")]:
    a1.axhline(y, color="#ddd", lw=0.6); a1.text(1050, y * 1.15, t, fontsize=6.5, color="#888")
a1.legend(fontsize=7, frameon=False); a1.grid(color="#eee", lw=0.5); [a1.spines[s].set_visible(False) for s in ("top", "right")]
km = [(D["kimonu"][m]["sec_conc"], D["kimonu"][m]["sec_kin"], D["kimonu"][m]["sec_acr"]) for m in models]; tot = [sum(t) for t in km]
x = np.arange(len(models)); b = np.zeros(len(models))
for i, (lab, col) in enumerate([("concordance", "#1C6B58"), ("kinetic", "#7FB8A6"), ("ACR/ACRR", "#C9DED6")]):
    frac = np.array([t[i] / tt for t, tt in zip(km, tot)]); a2.bar(x, frac, 0.8, bottom=b, color=col, label=lab, edgecolor="white", lw=0.3); b += frac
a2.set_xticks(x); a2.set_xticklabels(models, rotation=65, ha="right", fontsize=6); a2.set_ylabel("fraction of Kimonu runtime"); a2.set_title("F5b  Where Kimonu spends its time", loc="left", fontsize=11)
a2.legend(fontsize=7, frameon=False, loc="lower left"); [a2.spines[s].set_visible(False) for s in ("top", "right")]
rf = [D["ref_exh"][m]["sec_conc"] / D["ref_exh"][m]["seconds"] for m in D["ref_exh"]]
save(fig, "F5_runtime", f"Reference spends {100*min(rf):.1f}–{100*max(rf):.1f}% of its time in concordance ({len(rf)} models landed); reference times are from the fixed bridge (cold re-solve on any non-optimal LP status).  Stochastic arms: median marker, seeds as faint points.")
print("rendered:", sorted(os.listdir("figures/prelim")))
# the discrepancies F4 makes visible, for the record
for m in ["iIT341", "iHN637", "iJN678", "iLJ478"]:
    k = D["kimonu"][m]; cf = D["cocoa_false"].get(m, [])
    print(f"{m:10s} ACR free K={k['acr_free']} Cfalse={[r['acr_free'] for r in cf]} | ACRR free K={k['acrr_free']} Cfalse={[r['acrr_free'] for r in cf]} | giant K={k['giant']} Cfalse={[r['giant'] for r in cf]}")
