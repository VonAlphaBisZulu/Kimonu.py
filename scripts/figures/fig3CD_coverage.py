"""Figure 3C, D: how much of each network the relations cover, drawn into the vault Figure 3.

The vault Figure 3 (A, B the counts per network, C, D the tolerance sweep) gets a new middle row:

C  share of the complexes that lie in a kinetic module (two or more complexes), against the number of complexes
D  share of the free metabolites with ACR or in an ACRR group, against the number of free metabolites
   (SHARE_OF=all: all metabolites, enzyme forms included)

The old C, D and their legend move down and become E, F. Marks as in A, B: one dot where the two routes agree, so the
share is exact, and a bar from the iterated route (proven) to the tangent route (complete) on the networks with a
gap. Data: kinetic_variants/member_counts/<model>.json (member_counts.py, detached uptakes). Complexes are counted as
in table S2, the zero complex included and the complexes that occur in no reaction (51 in iJB785, 29 in iSynCJ816)
left out.

    python figures/fig3CD_coverage.py <vault Figure3.svg> <out.svg> [--all-rows]     # from benchmarks/

Everything that was in the input must render as before, checked at the end (rows compared pixel by pixel, the two
relabelled letters left out). An input that already has the panels C and D (groups panel_C, panel_D) gets them
redrawn in place instead, and everything outside their row must render as before. With --all-rows the complexes are
counted as before 2026-11-11 (every complex row of the crn.json), which must redraw the vault figure of that time as it
was (the check of the redraw).
"""
import glob, hashlib, io, json, os, sys, tempfile
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_finish, svg_panels

src, dest = sys.argv[1], sys.argv[2]
SHARE_OF = os.environ.get("SHARE_OF", "free")
text = open(src, encoding="utf-8").read()
print("input", hashlib.sha1(text.encode()).hexdigest()[:12])
W = 518.4
H_OLD = float(ET.fromstring(text).get("viewBox").split()[3])
Y0 = 163.5                                                # pt, where the old row C, D starts
ROW = 119.0                                               # height of the new row
REPLACE = 'id="panel_C"' in text                          # the panels exist: redraw them in place
ALL_ROWS = "--all-rows" in sys.argv[3:]                   # every complex row of the crn.json, as before 2026-11-11 (the check)
H = H_OLD if REPLACE else round(H_OLD + ROW, 1)
XA, XB, WAX = 48.29, 299.77, 209.95                       # the axes of A and B, which the new panels line up with
AX_TOP, AX_H = Y0 + 19.0, 68.0
INK, MUTE, RULE, GREY, DARK = "#1b2430", "#6b7480", "#d5d9de", "#6b7480", "#3a4250"
from names import NAME as SHORT                    # one display name per network (figures/names.py)

plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "svg.fonttype": "none",
                     **svg_finish.mathtext_rc()})

D = [json.load(open(f)) for f in sorted(glob.glob("kinetic_variants/member_counts/*.json")) if not f.endswith(".attached.json")]
assert len(D) == 33 and all(d["uptakes"] == "detached" for d in D), len(D)
assert sorted(d["model"] for d in D if d["gap"]) == ["iCN718", "iHN637", "iNJ661"]

def n_complexes(m):
    """complexes of the network with the zero complex, without those that occur in no reaction (the rows of A in use)"""
    return len({i for i, _, _ in json.load(open(f"unblocked/{m}_fixed.crn.json"))["A"]})
def complexes(d, r):
    n = d["n_complexes"] if ALL_ROWS else n_complexes(d["model"])
    return n, 100 * r["complexes_in_modules"] / n
def metabolites(d, r):
    if SHARE_OF == "all": return d["n_metabolites"], 100 * r["fixed_members"] / d["n_metabolites"]
    return d["n_free"], 100 * r["fixed_members_free"] / d["n_free"]

# title, x label, y label, x limits, x ticks, y limits, y ticks, labelled networks: (dx, dy in pt, ha)
PANELS = [
    (complexes, "complexes in kinetic modules", "complexes", "% of complexes", (900, 17500), [1000, 2000, 5000, 10000],
     (0, 72), [0, 20, 40, 60], {"iJB785": (4, 0, "left"), "iSynCJ816": (4, 0, "left"), "iCN718": (0, -8, "center")}),
    (metabolites, "free metabolites with ACR or ACRR" if SHARE_OF == "free" else "metabolites with ACR or ACRR",
     "free metabolites" if SHARE_OF == "free" else "metabolites", "% of free metabolites" if SHARE_OF == "free" else "% of metabolites",
     (150, 1250) if SHARE_OF == "free" else (600, 11000), [200, 500, 1000] if SHARE_OF == "free" else [1000, 2000, 5000, 10000],
     (0, 17) if SHARE_OF == "free" else (0, 13), [0, 5, 10, 15] if SHARE_OF == "free" else [0, 4, 8, 12],
     {"iCN718": (-4, 0, "right"), "iHN637": (-4, 0, "right"), "iNJ661": (4, 0, "left")}),
]

fig = plt.figure(figsize=(W / 72, H / 72))
fig.patch.set_visible(False)
ay = lambda pt: 1 - pt / H
summary = {}
for (get, title, xlabel, ylabel, xlim, xticks, ylim, yticks, names), x0, letter, lx, tx in zip(PANELS, (XA, XB), "CD", (7.87, 259.07), (60.77, 318.49)):
    ax = fig.add_axes([x0 / W, ay(AX_TOP + AX_H), WAX / W, AX_H / H])
    ax.patch.set_visible(False)
    vals = []
    for d in D:
        x, lo = get(d, d["it"])
        hi = get(d, d["tan"])[1] if d["gap"] else lo
        vals.append((d["model"], x, lo, hi))
        if d["gap"]:                                     # iterated route (proven) at the bottom, tangent route (complete) at the top
            ax.errorbar(x, hi, yerr=[[hi - lo], [0]], fmt="none", ecolor=DARK, elinewidth=0.9, capsize=2.2, capthick=0.9, zorder=4)
        else:
            ax.scatter(x, lo, s=10, color=GREY, zorder=3, linewidths=0)
        if d["model"] in names:
            dx, dy, ha = names[d["model"]]
            ax.annotate(SHORT.get(d["model"], d["model"]), (x, (lo + hi) / 2), xytext=(dx, dy), textcoords="offset points",
                        ha=ha, va="center", fontsize=6, color=DARK)
    ax.set_xscale("log"); ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    ax.set_xticks(xticks); ax.set_xticklabels([f"{t:,}" for t in xticks]); ax.minorticks_off()
    ax.set_yticks(yticks)
    ax.set_xlabel(xlabel, labelpad=2, color=INK); ax.set_ylabel(ylabel, labelpad=2, color=INK)
    ax.tick_params(axis="both", colors=INK, length=2)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(MUTE)
    ax.grid(axis="y", color=RULE, lw=0.4); ax.set_axisbelow(True)
    fig.text(tx / W, ay(Y0 + 12.99), title, fontsize=8, fontweight="bold", color=INK, ha="left", va="baseline")
    fig.text(lx / W, ay(Y0 + 16.87), letter, fontsize=18, fontweight="bold", color=INK, ha="left", va="baseline")
    exact = sorted(v[2] for v in vals if v[2] == v[3])
    summary[title] = dict(n_exact=len(exact), min=round(exact[0], 1), median=round(float(np.median(exact)), 1), max=round(exact[-1], 1),
                          gap={m: (round(lo, 1), round(hi, 1)) for m, _, lo, hi in vals if lo != hi})
buf = io.StringIO(); fig.savefig(buf, format="svg"); raw = buf.getvalue(); plt.close(fig)
n, box = svg_finish.check(raw, raw, tempfile.mkdtemp())
assert n == 0, f"removing the clipping changes {n} pixels in {box} (pt)"

# --- the old figure: the row below Y0 moves down, its letters C, D become E, F; the new leaves join it
if REPLACE:                                               # or the old C, D go and the new ones take their place
    root = ET.fromstring(text)
    for g in [g for g in root if g.get("id") in ("panel_C", "panel_D")]:
        root.remove(g)
    old = svg_panels.Figure(ET.tostring(root, encoding="unicode"))
else:
    old = svg_panels.Figure(text)
new = svg_panels.Figure(raw)
T = np.array([[1, 0, 0], [0, 1, ROW], [0, 0, 1.0]])
relabelled = 0
for lf in old.leaves if not REPLACE else []:
    b = lf.bbox()
    if b is not None and (b[1] + b[3]) / 2 > Y0:
        lf.M = T @ lf.M
        if lf.tag() == "text" and lf.el.get("font-size") == "18":
            for e in lf.el.iter():
                if e.text and e.text.strip() in ("C", "D"):
                    e.text = {"C": "E", "D": "F"}[e.text.strip()]; relabelled += 1
assert relabelled == (0 if REPLACE else 2), relabelled
for lf in new.leaves:
    lf.order += len(old.leaves)
old.leaves += new.leaves
have = {d.get("id") for d in old.nested_defs} | {e.get("id") for d in old.defs for e in d.iter()}
old.nested_defs += [d for d in new.nested_defs if d.get("id") not in have]     # a marker both drew is defined once
XD = 259.0                                                 # between the two new panels, as between A and B
missing = old.assign([("A", [(0, 0, 259, 151.5)]), ("B", [(259, 0, W, 151.5)]), ("legend_AB", [(0, 151.5, W, Y0)]),
                      ("C", [(0, Y0, XD - 6, Y0 + ROW)]), ("D", [(XD - 6, Y0, W, Y0 + ROW)]),
                      ("E", [(0, Y0 + ROW, 206, H)]), ("F", [(206, Y0 + ROW, 404, H)]), ("legend_EF", [(404, Y0 + ROW, W, H)])])
assert not missing, [ET.tostring(m.el)[:100] for m in missing]
old.root.set("viewBox", f"0 0 {W} {H}")
merged, used = svg_finish.finish(old.svg())
open(dest, "w", encoding="utf-8").write(merged)

# --- the old drawing must be unchanged: top row in place, bottom row moved by ROW (letters C, D excepted)
tmp = tempfile.mkdtemp()
a, b = os.path.join(tmp, "a.png"), os.path.join(tmp, "b.png")
svg_finish._render(text, a); svg_finish._render(merged, b)
A = np.asarray(Image.open(a).convert("L"), dtype=int); B = np.asarray(Image.open(b).convert("L"), dtype=int)
k = A.shape[1] / W
cut, row = int(round(Y0 * k)), int(round(ROW * k))
if REPLACE:                                               # outside the row of C, D nothing may change
    assert A.shape == B.shape, (A.shape, B.shape)
    ch = np.abs(A - B) > 40
    inside, end = int(ch[cut:cut + row].sum()), cut + row
    outside = int(ch[:cut].sum() + ch[end:].sum())
    ys, xs = np.nonzero(ch[cut:end])
    print(f"pixels changed: outside the row {outside}, inside {inside}"
          + (f" (x {xs.min() / k:.1f} to {xs.max() / k:.1f} pt, y {Y0 + ys.min() / k:.1f} to {Y0 + ys.max() / k:.1f} pt)" if inside else "")
          + f" | {W} x {H} pt | fonts {sorted(used)}")
    assert outside == 0
    print(json.dumps(summary, indent=1))
    sys.exit(0)
top = int((np.abs(A[:cut] - B[:cut]) > 40).sum())
lo_a, lo_b = A[cut:].copy(), B[cut + row:cut + row + A.shape[0] - cut].copy()
for x in (0, 206):                                        # the two relabelled letters
    for M in (lo_a, lo_b):
        M[:int(22 * k), int(x * k):int((x + 18) * k)] = 255
bottom = int((np.abs(lo_a - lo_b) > 40).sum())
print(f"pixels changed: top row {top}, moved row {bottom} | {W} x {H} pt | fonts {sorted(used)}")
assert top == 0 and bottom == 0
print(json.dumps(summary, indent=1))
