"""On Phil's Illustrator save of Figure 1: matrix names K and M at the top left of their matrices (as B and Y^T B in F),
the letters C and E with their pills a little higher, B = ker(V^T) moved up and left with what V and B are on its right,
then fonts embedded and the size set in points (svg_finish).

Illustrator's own copy of the document (metadata with i:aipgf, written with "Preserve Illustrator Editing Capabilities")
is dropped: it would no longer match the SVG, and Illustrator refuses such a file as invalid.

    python fig1_update_matrix_labels.py <Figure1_draft2.svg from Illustrator> <out.svg>
"""
import os, re, sys
import xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "figures"))
import svg_finish
from svg_panels import NS, text_width, typeset

src, out = sys.argv[1:3]
INK, MUTE = "#1b2430", "#9aa3ad"
FB = "#c0632b"
ET.register_namespace("", "http://www.w3.org/2000/svg")
ET.register_namespace("i", "http://ns.adobe.com/AdobeIllustrator/10.0/")
root = ET.parse(src).getroot()
for md in list(root.iter(f"{NS}metadata")):                 # Illustrator's private copy of the document
    for parent in root.iter():
        if md in list(parent):
            parent.remove(md)
panel = {g.get("id")[6:]: g for g in root.iter(f"{NS}g") if (g.get("id") or "").startswith("panel_")}


def pos(el):
    m = re.search(r"translate\(([-\d.]+)[ ,]([-\d.]+)\)", el.get("transform") or "")
    return (float(m.group(1)), float(m.group(2))) if m else None


def shift(el, dx, dy):
    p = pos(el)
    if p and (el.get("transform") or "").startswith("translate"):
        el.set("transform", re.sub(r"translate\([^)]*\)", f"translate({p[0] + dx:.3f} {p[1] + dy:.3f})", el.get("transform"), count=1))
    else:
        el.set("transform", f"translate({dx:.3f} {dy:.3f})" + (" " + el.get("transform") if el.get("transform") else ""))


def tex(g, x, y, s, size, fill):
    el, box = typeset(s, size, fill)
    el.set("transform", f"translate({x:.3f} {y:.3f})")
    g.append(el)
    return box[2]


def words(g, x, y, t, size=6, fill=MUTE):
    el = ET.SubElement(g, f"{NS}text", {"transform": f"translate({x:.3f} {y:.3f})", "fill": fill, "font-size": f"{size:g}", "font-family": "ArialMT, Arial"})
    ET.SubElement(el, f"{NS}tspan", {"x": "0", "y": "0"}).text = t


# ---- C and E: panel letter and pill 3 pt higher
UP = -3.0
for k, letter in (("C", "C"), ("E", "E")):
    for el in panel[k].iter():
        tag = el.tag.replace(NS, "")
        txt = "".join(el.itertext()).strip() if tag == "text" else ""
        if tag == "text" and ((el.get("font-size") in ("18", "10") and txt == letter) or (el.get("fill") == "#fff" and txt == "flux balance")):
            shift(el, 0, UP)
        elif tag == "path" and el.get("fill") == FB:
            shift(el, 0, UP)

# ---- K and M at the top left of their matrices, on the column of the row labels (8 pt bold, as B and Y^T B in F)
tex(panel["C"], 38.6, 251.0, r"\mathbf{K}", 8, INK)
tex(panel["E"], 191.5, 251.0, r"\mathbf{M}", 8, INK)

# ---- F: B = ker(V^T) up and left, what V and B are on its right
DX, DY = -48.8, -5.0
for el in list(panel["F"].iter()):
    tag = el.tag.replace(NS, "")
    if tag == "text" and "LMRoman" in (el.get("font-family") or ""):
        p = pos(el)
        if p and 405 < p[0] < 475 and abs(p[1] - 245.1) < 0.5:
            shift(el, DX, DY)
    elif tag == "path" and (el.get("d") or "").startswith("M465.59,237.52"):   # the transpose sign, an outline
        shift(el, DX, DY)
X = 430.0
for y, sym, rest in ((236.5, r"\mathbf{V}", ": fixed relationships"), (244.0, r"\mathbf{B}", r": free directions for ln x")):
    w = tex(panel["F"], X, y, sym, 6, MUTE)
    words(panel["F"], X + w + 0.4, y, rest)

text = ET.tostring(root, encoding="unicode")
open(out, "w", encoding="utf-8").write(svg_finish.finish(text)[0])
print("written", out)
