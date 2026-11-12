"""Figure 1, panel G as in the paper (a test variant, made the vault Figure1.svg on 2026-11-11): W above both routes,
u and the feedback rows in the iterated route.

In the vault figure the definition W = N diag(r*) S^T sits inside the first box of the iterated route, although the
tangent route (W Y^T h = 0) uses the same matrix, and the unknown u of the iterated route appears only in the
feedback label. Here

  - W = N diag(r*) S^T and the note on r* move up to the line of the "mass action" pill, above both routes,
  - the first box reads ker(W), W u = 0 and what u is (u_c = psi_c / psi*_c),
  - the feedback label says what the feedback does: it adds rows u_c = u_d to W.

Nothing else moves. Everything outside the three edited places must render as before, checked pixel by pixel.

    python fig1G_update_shared_W.py <vault Figure1.svg> <out.svg>
"""
import hashlib, os, re, sys, tempfile
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "figures"))
import svg_finish
from svg_panels import NS, text_width, typeset

src, out = sys.argv[1:3]
text = open(src, encoding="utf-8").read()
print("input", hashlib.sha1(text.encode()).hexdigest()[:12])
INK, MUTE = "#1b2430", "#9aa3ad"
ET.register_namespace("", "http://www.w3.org/2000/svg")
root = ET.fromstring(text)
assert not list(root.iter(f"{NS}metadata")), "Illustrator's private copy is still in the file"
G = [g for g in root.iter(f"{NS}g") if g.get("id") == "panel_G"][0]


def pos(el):
    m = re.search(r"translate\(([-\d.]+)[ ,]([-\d.]+)\)", el.get("transform") or "")
    return (float(m.group(1)), float(m.group(2))) if m else None


def shift(el, dx, dy):
    """Move a <text> (its own translate) or an outline <path> (a translate in front)."""
    p = pos(el)
    if p and el.get("transform").startswith("translate"):
        el.set("transform", re.sub(r"translate\([^)]*\)", f"translate({p[0] + dx:.3f} {p[1] + dy:.3f})", el.get("transform"), count=1))
    else:
        el.set("transform", f"translate({dx:.3f} {dy:.3f})" + (" " + el.get("transform") if el.get("transform") else ""))


def said(el): return re.sub(r"\s+", "", "".join(el.itertext()))
def texts(el): return list(el.iter(f"{NS}text"))
def child(test):
    hits = [c for c in G if test(c)]
    assert len(hits) == 1, len(hits)
    return hits[0]


def tex(x, y, s, size, fill, anchor="left"):
    el, box = typeset(s, size, fill)
    if anchor == "center": x -= (box[0] + box[2]) / 2
    el.set("transform", f"translate({x:.3f} {y:.3f})")
    G.append(el)
    return x + box[2]


def words(x, y, t, size=6, fill=MUTE):
    el = ET.SubElement(G, f"{NS}text", {"transform": f"translate({x:.3f} {y:.3f})", "fill": fill, "font-size": f"{size:g}", "font-family": "ArialMT, Arial"})
    ET.SubElement(el, f"{NS}tspan", {"x": "0", "y": "0"}).text = t
    return x + text_width(t, "ArialMT", size)


# ---- the pieces of the vault panel that move or go
W_def = child(lambda c: said(c).startswith("W=Ndiag"))
W_T = child(lambda c: c.tag == f"{NS}path" and (c.get("d") or "").startswith("M277.51,458.66"))      # its transpose sign, an outline
r_note = child(lambda c: said(c) == ":asteady-stateflux")
r_star = child(lambda c: said(c) == "r*" and texts(c)[0].get("fill") == MUTE and abs(pos(texts(c)[0])[1] - 475.73) < 0.5)
feedback = child(lambda c: said(c) == "feedback")
feedback_u = child(lambda c: said(c) == "uc=ud")
box1 = child(lambda c: c.tag == f"{NS}rect" and c.get("x") == "206.21" and c.get("width") == "86")
box2 = child(lambda c: c.tag == f"{NS}rect" and c.get("width") == "76")
pill = child(lambda c: c.tag == f"{NS}path" and c.get("fill") == "#2f6f9f")
BX0, BW = float(box1.get("x")), float(box1.get("width"))
B2X0 = float(box2.get("x"))

# ---- W and the note on r* go up to the line of the pill, which both routes share
HEAD = 428.4                                                # baseline, the 7 pt capitals centred on the pill
x_W = pos(texts(W_def)[0])[0]
dx, dy = 253.0 - x_W, HEAD - pos(texts(W_def)[0])[1]
for t in texts(W_def): shift(t, dx, dy)
shift(W_T, dx, dy)
w_end = 277.51 + dx + 3.5                                   # right end of S^T
dx, dy = w_end + 7.0 - pos(texts(r_star)[0])[0], HEAD - pos(texts(r_note)[0])[1]
for t in texts(r_star) + texts(r_note): shift(t, dx, dy)

# ---- the first box: W u = 0, and what u is
mid = BX0 + BW / 2
tex(mid, 464.05, r"\mathbf{W}\mathbf{u} = \mathbf{0}", 7, INK, "center")
tex(mid, 475.0, r"u_c = \psi_c\,/\,\psi^*_c", 6, MUTE, "center")

# ---- the feedback label: what the feedback does to W
G.remove(feedback); G.remove(feedback_u)
parts = [("w", "feedback adds rows"), ("t", r"u_c = u_d"), ("w", "to"), ("t", r"\mathbf{W}")]
GAP = text_width(" ", "ArialMT", 6)                         # between the parts (a space at the end of a <text> would collapse)
width = sum(text_width(s, "ArialMT", 6) if k == "w" else typeset(s, 6, MUTE)[1][2] for k, s in parts) + 3 * GAP
x, y = (BX0 + BW + B2X0) / 2 - width / 2, 491.3
for k, s in parts:
    x = (words(x, y, s) if k == "w" else tex(x, y, s, 6, MUTE)) + GAP

merged, used = svg_finish.finish(ET.tostring(root, encoding="unicode"))
open(out, "w", encoding="utf-8").write(merged)

# ---- everything outside the three edited places must render as before
tmp = tempfile.mkdtemp()
a, b = os.path.join(tmp, "a.png"), os.path.join(tmp, "b.png")
svg_finish._render(text, a); svg_finish._render(merged, b)
A = np.asarray(Image.open(a).convert("L"), dtype=int); B = np.asarray(Image.open(b).convert("L"), dtype=int)
assert A.shape == B.shape, (A.shape, B.shape)
k = A.shape[1] / 518.4
D = np.abs(A - B) > 40
for x0, y0, x1, y1 in ((248, 420, 420, 432),                 # the line of the pill
                       (BX0 + 2, 457, BX0 + BW - 2, 480),    # the first box below ker(W)
                       (262, 484, 364, 495)):                # the feedback label
    D[int(y0 * k):int(y1 * k), int(x0 * k):int(x1 * k)] = False
print(f"pixels changed outside the three edited places: {int(D.sum())} | fonts {sorted(used)}")
if D.any():
    ys, xs = np.nonzero(D)
    print("changed box (pt):", xs.min() / k, ys.min() / k, xs.max() / k, ys.max() / k)
assert not D.any()
