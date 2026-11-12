"""Figure 2, panel C, lower plot: the nullclines of C and P on log-log axes.

On linear axes the nullcline dP/dt = 0 is the hyperbola P = 4/C. On log-log axes it is the straight line through
both steady states, A at (1, 4) and B at (4, 1), which is the line ln x_A + B of the text (B the differences of
ln x between steady states, here spanned by ln x_B - ln x_A). The nullcline dC/dt = 0 crosses it at the two
steady states only.

Only the drawing of the lower plot is replaced (axes, ticks, the two curves, the two dots and their labels).
The title, the axis labels and the formulas on the right stay where they are, and so does everything else in
the figure, checked pixel by pixel at the end.

    python fig2C_update_loglog.py <vault Figure2.svg> <out.svg>
"""
import hashlib, math, os, sys, tempfile
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "figures"))
import svg_finish
from svg_panels import NS, text_width

src, out = sys.argv[1:3]
text = open(src, encoding="utf-8").read()
print("input", hashlib.sha1(text.encode()).hexdigest()[:12])
ET.register_namespace("", "http://www.w3.org/2000/svg")
root = ET.fromstring(text)
panel = {g.get("id")[6:]: g for g in root.iter(f"{NS}g") if (g.get("id") or "").startswith("panel_")}
C = panel["C"]
INK, MUTE, BLUE, ORANGE = "#1b2430", "#9aa3ad", "#5a9bd5", "#e6862b"
X0, X1, Y0, Y1 = 358.15, 484.03, 229.38, 161.41               # the axes of the lower plot in pt (Y0 the bottom)
LO, HI = 0.5, 8.0                                             # both axes, log scale
TICKS = [0.5, 1, 2, 4]                                        # as before, the upper end carries no tick


def px(c): return X0 + (X1 - X0) * math.log(c / LO) / math.log(HI / LO)
def py(p): return Y0 + (Y1 - Y0) * math.log(p / LO) / math.log(HI / LO)


# ---- what goes: the drawing of the old lower plot
def old_plot(el):
    tag = el.tag.replace(NS, "")
    if tag == "line":
        return float(el.get("y1")) > 150
    if tag == "polyline":
        return float(el.get("points").split()[1]) > 150
    if tag == "circle":
        return float(el.get("cy")) > 150
    if tag == "text" and el.get("transform"):
        x, y = (float(v) for v in el.get("transform")[10:-1].split())
        txt = "".join(el.itertext()).strip()
        return 150 < y < 245 and x < 445 and (el.get("font-size") == "6" or (el.get("font-size") == "7" and txt in ("A", "B")))
    return False


gone = [el for el in list(C) if old_plot(el)]
assert len(gone) == 20, len(gone)                             # 8 lines, 6 tick labels, 2 curves, 2 dots, A and B
at = list(C).index(gone[0])
for el in gone:
    C.remove(el)

# ---- what comes
new = []
def line(x1, y1, x2, y2):
    new.append(ET.Element(f"{NS}line", {"x1": f"{x1:.2f}", "y1": f"{y1:.2f}", "x2": f"{x2:.2f}", "y2": f"{y2:.2f}", "fill": "none",
                                        "stroke": INK, "stroke-linejoin": "round", "stroke-width": ".6"}))
def words(x, y, t, size=6, fill=INK):
    el = ET.Element(f"{NS}text", {"fill": fill, "font-family": "ArialMT, Arial", "font-size": f"{size:g}", "transform": f"translate({x:.3f} {y:.3f})"})
    ET.SubElement(el, f"{NS}tspan", {"x": "0", "y": "0"}).text = t
    new.append(el)
def curve(pts, colour):
    new.append(ET.Element(f"{NS}polyline", {"points": " ".join(f"{x:.2f} {y:.2f}" for x, y in pts), "fill": "none", "stroke": colour,
                                            "stroke-linejoin": "round", "stroke-width": "1.2"}))
def dot(c, p):
    new.append(ET.Element(f"{NS}circle", {"cx": f"{px(c):.2f}", "cy": f"{py(p):.2f}", "r": "3.02", "fill": INK, "stroke": "#fff",
                                          "stroke-linejoin": "round", "stroke-width": ".8"}))

line(X0, Y0, X1, Y0); line(X0, Y0, X0, Y1)
for v in TICKS:
    s = f"{v:g}"
    line(px(v), Y0, px(v), Y0 + 3.15)
    words(px(v) - text_width(s, "ArialMT", 6) / 2, Y0 + 12.59, s)
    line(X0, py(v), X0 - 3.15, py(v))
    words(X0 - 7.22 - text_width(s, "ArialMT", 6), py(v) + 3.15, s)
# dC/dt = 0: P = 2 sqrt((5 - C)/C), that is C = 20/(P^2 + 4), from the left edge down to the lower edge
c_left = LO
p_top = 2 * math.sqrt((5 - c_left) / c_left)
P = np.exp(np.linspace(math.log(p_top), math.log(LO), 160))
curve([(px(20 / (p * p + 4)), py(p)) for p in P], BLUE)
# dP/dt = 0: P = 4/C, a straight line on these axes
curve([(px(c), py(4 / c)) for c in (LO, 2.0, HI)], ORANGE)
dot(1, 4); words(px(1) + 4.5, py(4) - 3.0, "A", 7)
dot(4, 1); words(px(4) - 9.5, py(1) + 8.5, "B", 7)
words(X0 + 5.0, Y0 - 4.5, "log–log axes", 6, MUTE)
for k, el in enumerate(new):
    C.insert(at + k, el)

# both steady states lie on both curves
for c, p in ((1, 4), (4, 1)):
    assert abs(5 - c - c * p * p / 4) < 1e-12 and abs(c * p * p / 4 - p) < 1e-12

merged, used = svg_finish.finish(ET.tostring(root, encoding="unicode"))
open(out, "w", encoding="utf-8").write(merged)

# ---- everything outside the lower plot must render as before
tmp = tempfile.mkdtemp()
a, b = os.path.join(tmp, "a.png"), os.path.join(tmp, "b.png")
svg_finish._render(text, a); svg_finish._render(merged, b)
A = np.asarray(Image.open(a).convert("L"), dtype=int); B = np.asarray(Image.open(b).convert("L"), dtype=int)
assert A.shape == B.shape, (A.shape, B.shape)
k = A.shape[1] / 518.4
D = np.abs(A - B) > 40
D[int(156 * k):int(245 * k), int(340 * k):int(488 * k)] = False
print(f"pixels changed outside the lower plot of C: {int(D.sum())} | fonts {sorted(used)}")
assert not D.any()
