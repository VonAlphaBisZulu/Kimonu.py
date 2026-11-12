"""Figure 2, panel E: the third row of the kinetic-module table and the route headers of both tables.

The kinase has seven complexes (∅, X, 2X, C, C+2P, 3P, P) and two kinetic modules, {∅, X, 2X} and {P, C+2P}. The
tangent route reports one module of all seven complexes, which contains both true modules and is false as a whole.
The third row of the right table used to be labelled "17 other pairs", as if the table listed pairs. It now names
the two complexes that complete the seven, "{C}, {3P}", and the tangent-route column gets one box across all three
rows. Its upper two thirds keep the two greens of the true modules, its lower third is red, and the outline of the
whole box is the red of a false positive, because the module as a whole is a false coupling. One "coupled" label sits
at the centre of the box.

The headers "iterated route" and "tangent route" of both tables lose their underline and bold and become regular
text like the "ground truth" headers next to them. The regular face is narrower, so each word is re-centred on the
centre it had before. Nothing else in the figure changes, checked pixel by pixel at the end.

    python fig2E_update_table.py <vault Figure2.svg> <out.svg> [<panelE_before.png> <panelE_after.png>]
"""
import hashlib, os, sys, tempfile
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
E = panel["E"]
OLD_LABEL, NEW_LABEL = "17 other pairs", "{C}, {3P}"
GROUND_TRUTH_X = 330.54                                       # left edge of the ground-truth column of the right table
HEADERS = {(196.760, 285.790): "iterated", (200.650, 294.190): "route", (229.800, 285.790): "tangent", (233.690, 294.190): "route",
           (367.390, 285.870): "iterated", (371.280, 294.270): "route", (400.430, 285.870): "tangent", (404.320, 294.270): "route"}


def xy(el):
    return tuple(round(float(v), 3) for v in el.get("transform")[10:-1].split())


def tspan(el):
    return el.find(f"{NS}tspan")


# ---- 1. the third row of the right table names the one module of all seven complexes
labels = [el for el in E.iter(f"{NS}text") if "".join(el.itertext()).strip() == OLD_LABEL]
assert len(labels) == 1, len(labels)
label = labels[0]
assert label.get("font-family") == "Arial-BoldMT, Arial" and label.get("font-size") == "8", label.attrib
x_label = xy(label)[0]
w_old = text_width(OLD_LABEL, "Arial-BoldMT", 8)
w_new = text_width(NEW_LABEL, "Arial-BoldMT", 8)
gap = GROUND_TRUTH_X - (x_label + w_new)
assert gap >= 3, gap
tspan(label).text = NEW_LABEL
print(f"row label {OLD_LABEL!r} ({w_old:.2f} pt) -> {NEW_LABEL!r} ({w_new:.2f} pt), gap to the ground-truth column {gap:.2f} pt")

# ---- 2. the tangent-route column of the right table: one box across the three rows
X, W, R = 398.54, 30.0, 2.0                                   # the column's boxes as drawn
Y0, Y1, Y2, Y3 = 299.5, 310.8, 322.1, 335.1                   # top, the two row boundaries, bottom (row 3 ends at 335.1)
LIGHT, DARK, RED_FILL, RED, INK = "#d3efd4", "#a8dcb0", "#f6c9c2", "#c8453b", "#1b2430"
TEXT_X, TEXT_DY = 401.28, 2.5                                 # "coupled" left edge, baseline below a box centre


def in_column(el):
    tag = el.tag.replace(NS, "")
    if tag == "rect":
        return float(el.get("x")) == X and Y0 <= float(el.get("y")) < Y3
    if tag == "text" and el.get("transform"):
        x, y = xy(el)
        return x == TEXT_X and Y0 < y < Y3 and "".join(el.itertext()).strip() == "coupled"
    return False


old = [el for el in list(E) if in_column(el)]
assert len(old) == 6, len(old)                                # 2 fills + outline + "coupled", red box + "coupled"
at = list(E).index(old[0])
for el in old:
    E.remove(el)
new = [
    ET.Element(f"{NS}path", {"d": f"M{X + R},{Y0}H{X + W - R}A{R},{R},0,0,1,{X + W},{Y0 + R}V{Y1}H{X}V{Y0 + R}A{R},{R},0,0,1,{X + R},{Y0}Z",
                             "fill": LIGHT}),
    ET.Element(f"{NS}rect", {"x": f"{X}", "y": f"{Y1}", "width": f"{W:g}", "height": f"{Y2 - Y1:.1f}", "fill": DARK}),
    ET.Element(f"{NS}path", {"d": f"M{X},{Y2}H{X + W}V{Y3 - R}A{R},{R},0,0,1,{X + W - R},{Y3}H{X + R}A{R},{R},0,0,1,{X},{Y3 - R}Z",
                             "fill": RED_FILL}),
    ET.Element(f"{NS}rect", {"x": f"{X}", "y": f"{Y0}", "width": f"{W:g}", "height": f"{Y3 - Y0:.1f}", "rx": f"{R:g}", "ry": f"{R:g}",
                             "fill": "none", "stroke": RED, "stroke-linejoin": "round", "stroke-width": ".8"}),
]
label_el = ET.Element(f"{NS}text", {"fill": INK, "font-family": "ArialMT, Arial", "font-size": "7", "isolation": "isolate",
                                    "transform": f"translate({TEXT_X:.3f} {(Y0 + Y3) / 2 + TEXT_DY:.3f})"})
ET.SubElement(label_el, f"{NS}tspan", {"x": "0", "y": "0"}).text = "coupled"
new.append(label_el)
assert abs(TEXT_X - (X + (W - text_width("coupled", "ArialMT", 7)) / 2)) < 0.05
for k, el in enumerate(new):
    E.insert(at + k, el)
print(f"tangent column: one box {X} to {X + W} x {Y0} to {Y3} pt, greens to {Y2}, red below, red outline, label baseline {(Y0 + Y3) / 2 + TEXT_DY}")

# ---- 3. the route headers: regular face, no underline, same centre
done = 0
for el in E.iter(f"{NS}text"):
    if el.get("transform") is None or xy(el) not in HEADERS:
        continue
    word = "".join(el.itertext()).strip()
    assert word == HEADERS[xy(el)], (xy(el), word)
    assert el.get("font-weight") == "700" and el.get("text-decoration") == "underline" and el.get("font-size") == "7", el.attrib
    x, y = xy(el)
    centre = x + text_width(word, "Arial-BoldMT", 7) / 2
    x_new = centre - text_width(word, "ArialMT", 7) / 2
    del el.attrib["font-weight"], el.attrib["text-decoration"]
    el.set("font-family", "ArialMT, Arial")
    el.set("transform", f"translate({x_new:.3f} {y:.3f})")
    done += 1
assert done == len(HEADERS), done

merged, used = svg_finish.finish(ET.tostring(root, encoding="unicode"))
open(out, "w", encoding="utf-8").write(merged)
vb = [float(v) for v in ET.fromstring(merged).get("viewBox").replace(",", " ").split()]
print(f"output {hashlib.sha1(merged.encode()).hexdigest()[:12]} | {vb[2]:g} x {vb[3]:g} pt | fonts {sorted(used)}")

# ---- everything outside the edited cells of panel E must render as before
tmp = tempfile.mkdtemp()
a, b = os.path.join(tmp, "a.png"), os.path.join(tmp, "b.png")
svg_finish._render(text, a); svg_finish._render(merged, b)
A = np.asarray(Image.open(a).convert("L"), dtype=int); B = np.asarray(Image.open(b).convert("L"), dtype=int)
assert A.shape == B.shape, (A.shape, B.shape)
k = A.shape[1] / 518.4
D = np.abs(A - B) > 40
for x0, y0, x1, y1 in ((194, 279, 257, 297), (365, 279, 428, 297), (278, 325, 333, 336), (397, 298, 430, 337)):
    # headers left, headers right, row-3 label, tangent column of the right table
    D[int(y0 * k):int(y1 * k), int(x0 * k):int(x1 * k)] = False
print(f"pixels changed outside the panel E table area: {int(D.sum())}")
assert not D.any()
if len(sys.argv) > 4:                                          # optional crops of panel E before and after
    box = tuple(int(v * k) for v in (120, 268, 440, 372))
    Image.open(a).crop(box).save(sys.argv[3]); Image.open(b).crop(box).save(sys.argv[4])
