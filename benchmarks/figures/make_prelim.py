"""Preliminary figures F1-F5 from figdata.json.

Run from benchmarks/:  python figures/make_prelim.py
Every arm's provenance and trust level is in FINDINGS.md. Error bars are min-max over seeds and
appear only on stochastic arms; exact arms have none by construction.

figdata.json schema:
  models, n_complexes[m]
  kimonu[m]: balanced, conc, giant, giant_4phase, acr, acrr, acr_free, acrr_free, seconds, sec_conc, sec_kin, sec_acr
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
        "kimonu4":     ("Kimonu (tentative §S.2.2)*",   "#7FB8A6", "//"),
        "cocoa_false": ("COCOA efficient=false",        "#B8641E", ""),
        "cocoa_true":  ("COCOA efficient=true",         "#E0A26A", ""),
        "ref_exh":     ("Reference exhaustive",         "#4A5A88", ""),
        "ref_sampled": ("Reference sampled*",           "#9AA6C8", "//")}
BASE = "Error bars: min–max over seeds; exact arms have none by construction.  * comparability-only arm."

def val(arm, m, key):
    if arm == "kimonu":  v = D["kimonu"].get(m, {}).get(key); return None if v is None else (v, None, None)
    if arm == "kimonu4": v = D["kimonu"].get(m, {}).get("giant_4phase") if key == "giant" else None; return None if v is None else (v, None, None)
    if arm == "ref_exh": v = D["ref_exh"].get(m, {}).get(key); return None if v is None else (v, None, None)
    runs = [r[key] for r in D[arm].get(m, []) if r.get(key) is not None]
    return None if not runs else (float(np.mean(runs)), float(min(runs)), float(max(runs)))

def grouped(ax, key, arms, scale="linear", title="", ylabel="count", labels=None):
    labels = labels or {}; n = len(arms); w = 0.8 / n; x = np.arange(len(models)); drawn = {a: 0 for a in arms}
    for i, a in enumerate(arms):
        lab, col, hatch = ARMS[a]
        for j, m in enumerate(models):
            v = val(a, m, key)
            if v is None: continue
            mean, lo, hi = v; drawn[a] += 1; xx = x[j] + (i - n / 2 + 0.5) * w
            ax.bar(xx, mean, w, color=col, hatch=hatch, edgecolor=col if hatch else "white", linewidth=0.4)
            if lo is not None and hi > lo:
                ax.errorbar(xx, mean, yerr=[[mean - lo], [hi - mean]], fmt="none", ecolor="#222", elinewidth=0.7, capsize=1.5)
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
save(fig, "F1_balanced", BASE + "  Reference values are from the pre-fix bridge; balanced uses no pair LP and matched Kimonu exactly on 9/9.")
fig, ax = plt.subplots(figsize=(14, 5)); grouped(ax, "conc", ["kimonu", "cocoa_false", "cocoa_true", "ref_exh", "ref_sampled"], "log", "F2  Concordance modules (non-singleton)", "count (log)",
        labels={"ref_exh": "Reference exhaustive (pre-rerun)*"})
save(fig, "F2_concordance", BASE + "  Reference concordance predates the bridge fix and the three algorithm fixes; it is sound on 6 of 9 models by adjudication and will be rerun.")
fig, ax = plt.subplots(figsize=(14, 5)); grouped(ax, "giant", ["kimonu", "kimonu4", "cocoa_false", "cocoa_true", "ref_exh"], "log", "F3  Largest kinetic module", "count (log)",
        labels={"ref_exh": "Reference (R stage, corrected NA)*"})
save(fig, "F3_giant", BASE + "  Reference giants are from the corrected R stage on the reference's own pre-rerun concordance (6 models).")
# F4: ACR linear (mostly zeros), ACRR symlog
fig, axs = plt.subplots(2, 1, figsize=(14, 8.5)); a4 = ["kimonu", "cocoa_false", "cocoa_true", "ref_exh"]; lab4 = {"ref_exh": "Reference exhaustive †"}
grouped(axs[0], "acr_free", a4, "linear", "F4a  ACR metabolites (free)", labels=lab4)
grouped(axs[1], "acrr_free", a4, "symlog", "F4b  ACRR metabolite pairs (free), symlog", labels=lab4)
save(fig, "F4_acr_acrr", BASE + "  † reference ACR/ACRR uses the combinatorial criterion (Corollary S2-10), a subset of Kimonu's kernel test; every reference pair is in Kimonu's set on all 6 models.")
# F5
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [2, 1.2]})
for a in ["kimonu", "cocoa_false", "cocoa_true", "ref_exh", "ref_sampled"]:
    lab, col, _ = ARMS[a]; xs, ys, lo, hi = [], [], [], []
    for m in models:
        v = val(a, m, "sec_conc" if a == "ref_exh" else "seconds")
        if v is None: continue
        xs.append(nC[m]); ys.append(v[0]); lo.append(v[0] - (v[1] if v[1] is not None else v[0])); hi.append((v[2] if v[2] is not None else v[0]) - v[0])
    if xs: a1.errorbar(xs, ys, yerr=[lo, hi], fmt="o", ms=4.5, color=col, ecolor=col, elinewidth=0.7, capsize=1.5, label=f"{lab} (n={len(xs)})")
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
save(fig, "F5_runtime", f"Reference spends {100*min(rf):.1f}–{100*max(rf):.1f}% of its time in concordance (9 models); reference times are from the pre-fix bridge.  Error bars: min–max over seeds.")
print("rendered:", sorted(os.listdir("figures/prelim")))
# the discrepancies F4 makes visible, for the record
for m in ["iIT341", "iHN637", "iJN678", "iLJ478"]:
    k = D["kimonu"][m]; cf = D["cocoa_false"].get(m, [])
    print(f"{m:10s} ACR free K={k['acr_free']} Cfalse={[r['acr_free'] for r in cf]} | ACRR free K={k['acrr_free']} Cfalse={[r['acrr_free'] for r in cf]} | giant K={k['giant']} Cfalse={[r['giant'] for r in cf]}")
