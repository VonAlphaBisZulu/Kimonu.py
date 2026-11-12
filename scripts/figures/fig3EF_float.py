"""Figure 3E, F, drawn into the combined Figure 2 of early November (Phil's Illustrator file): panels A-G stay as they
are, the bottom row becomes the certificate of the two routes, exact nullspace (H) against SVD nullspace (I), both with one
element-wise tolerance rule. fig2_fig3_split.py then made H, I the panels E, F of Figure 3 (via C, D and fig3CD_coverage.py).

    python figures/fig3EF_float.py <vault Figure2_draft.svg> <out.svg>     # from scripts/

Drawn on a canvas of the full figure height, so that the new panels land where the old ones were; the render check at the
end compares everything above the bottom row with the input.
"""
import copy, hashlib, io, os, sys, tempfile
import xml.etree.ElementTree as ET
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_finish, svg_panels
import figS5_S6_float as ff
import matplotlib.pyplot as plt
from PIL import Image

SVG = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG)
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
src, dest = sys.argv[1], sys.argv[2]
text = open(src, encoding="utf-8").read()
print("input", hashlib.sha1(text.encode()).hexdigest()[:12])
root = ET.fromstring(text)
vb = [float(v) for v in root.get("viewBox").split()]
Wpt, Hpt = vb[2], vb[3]
TOP = 536.0                                               # start of the bottom row, pt from the top
AX0, AX1 = 549.2, 665.4                                   # plot height of the old H, I
XH, XI, WP = 36.3, 228.0, 170.0                           # left edges and width of the two plots, pt
W, H = Wpt / 72, Hpt / 72


def ay(pt):
    return 1 - pt / Hpt


fig = plt.figure(figsize=(W, H))
fig.patch.set_visible(False)
axH = fig.add_axes([XH / Wpt, ay(AX1), WP / Wpt, (AX1 - AX0) / Hpt])
axI = fig.add_axes([XI / Wpt, ay(AX1), WP / Wpt, (AX1 - AX0) / Hpt], sharey=axH)
out = ff.draw_certificate([axH, axI], ["exactcmp", "svdabs"], {"exactcmp": "exact nullspace", "svdabs": "SVD nullspace"})
plt.setp(axI.get_yticklabels(), visible=False)
for t in axI.texts:                                       # Phil's wording for the SVD bar at tolerance 0
    if "false certificates" in t.get_text():
        t.set_text("no rows ever exactly equal")
axH.set_ylabel("networks", labelpad=2, color=ff.INK)
for ax in (axH, axI):
    ax.set_xlabel("tolerance (0 = rows exactly equal)", labelpad=2, color=ff.INK)
    ax.tick_params(axis="x", pad=2)
    ax.tick_params(axis="y", length=2)
LEG_X = XI + WP + 8
ff.certificate_legend(fig, LEG_X / Wpt, ay(AX0 - 12), out["exactcmp"]["exact"])
x, s = out["exactcmp"], out["svdabs"]
fig.text(LEG_X / Wpt, ay(622), f"Exact rows at tolerance 0 settle\n{x['at_zero'][0]} of {sum(x['at_zero'])} networks, and the other {x['at_zero'][1]} are\n"
         f"the gap. Float rows are never\nexactly equal, so both routes find\nnothing and agree. No tolerance\n"
         f"repairs this. The SVD nullspace\ncertifies at most {max(s['certified'])} of {sum(s['exact'])}.",
         fontsize=6, color=ff.INK, va="top", ha="left", linespacing=1.25)
LET = dict(fontsize=10, fontweight="bold", family="Arial", color=ff.INK, va="top", ha="left")   # panel letters 10 pt bold (Science)
fig.text(0.0, ay(TOP + 1), "H", **LET)
fig.text((XI - 22) / Wpt, ay(TOP + 1), "I", **LET)
buf = io.StringIO()
fig.savefig(buf, format="svg")
raw = buf.getvalue()
plt.close(fig)
n, box = svg_finish.check(raw, raw, tempfile.mkdtemp())
assert n == 0, f"removing the clipping changes {n} pixels in {box} (pt)"
fg = svg_panels.Figure(raw)
miss = fg.assign([("H", [(0, TOP, XI - 23, Hpt)]), ("I", [(XI - 23, TOP, LEG_X - 4, Hpt)]), ("legend_HI", [(LEG_X - 4, TOP, Wpt, Hpt)])])
assert not miss, f"{len(miss)} elements outside every panel"
new = ET.fromstring(fg.svg())
groups = {g.get("id"): g for g in new.iter(f"{{{SVG}}}g") if (g.get("id") or "").startswith("panel_")}
for name in ("panel_H", "panel_I", "panel_legend_HI"):
    for old in [g for g in root.iter(f"{{{SVG}}}g") if g.get("id") == name]:
        next(p for p in root.iter() if old in list(p)).remove(old)
parent = next((p for p in root.iter() if any(g.get("id") == "panel_G" for g in list(p))), root)   # Illustrator may ungroup F, G
for k in ("panel_H", "panel_I", "panel_legend_HI"):
    parent.append(copy.deepcopy(groups[k]))
for d in new.iter(f"{{{SVG}}}defs"):                      # definitions the new panels refer to; fonts are redone below
    d2 = copy.deepcopy(d)
    for st in list(d2):
        if st.tag == f"{{{SVG}}}style":
            d2.remove(st)
    if len(d2):
        root.append(d2)
merged, used = svg_finish.finish('<?xml version="1.0" encoding="UTF-8"?>\n' + ET.tostring(root, encoding="unicode"))
open(dest, "w", encoding="utf-8").write(merged)
# everything above the bottom row must render as before
tmp = tempfile.mkdtemp()
a, b = os.path.join(tmp, "a.png"), os.path.join(tmp, "b.png")
svg_finish._render(text, a)
svg_finish._render(merged, b)
A = np.asarray(Image.open(a).convert("L"), dtype=int)
B = np.asarray(Image.open(b).convert("L"), dtype=int)
cut = int(TOP * A.shape[1] / Wpt)
print("cells counted as gray without a result:", {w: v["filled"] for w, v in out.items()})
print("pixels changed above the bottom row:", int((np.abs(A[:cut] - B[:cut]) > 40).sum()), "| size", A.shape, B.shape)
for way, v in out.items():
    print(way, len(v["models"]), "models | at 0:", v["at_zero"], "| exact:", v["exact"], "| best certified", max(v["certified"]),
          "| per tolerance", min(v["per_tau"]), "-", max(v["per_tau"]))
