"""On Phil's Illustrator save of Figure 2: Figure 1's colours for ACR and kinetic modules, panel E's legend as true,
missed (false negative) and false (false positive), "local" for the tangent route in the legend of F and G, the
formulas of panel C in Latin Modern, and C's tick labels in the same colour as the other plots. Illustrator's private
copy of the document is dropped (it would not match). The bottom row (H, I) is redrawn afterwards by
figures/fig3EF_float.py, which also embeds the fonts.

    python fig2_update_conventions.py <Figure2_draft.svg from Illustrator> <out.svg>
"""
import os, re, sys
import xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "figures"))
from svg_panels import NS, text_width, typeset

src, out = sys.argv[1:3]
INK, BLUE, ORANGE = "#1b2430", "#5a9bd5", "#e6862b"
COLOURS = {"#f2c14e": "#5b8fd1",       # ACR, blue as in Figure 1
           "#cdf3df": "#d3efd4",       # kinetic module {0, X, 2X}, light green
           "#c9ecf4": "#a8dcb0"}       # kinetic module {P, C+2P}, green
ET.register_namespace("", "http://www.w3.org/2000/svg")
ET.register_namespace("i", "http://ns.adobe.com/AdobeIllustrator/10.0/")
root = ET.parse(src).getroot()
parent = {c: p for p in root.iter() for c in p}
for md in list(root.iter(f"{NS}metadata")):
    parent[md].remove(md)

for el in root.iter():
    for attr in ("fill", "stroke"):
        v = (el.get(attr) or "").lower()
        if v in COLOURS:
            el.set(attr, COLOURS[v])


def xy(el):
    m = re.search(r"translate\(([-\d.]+)[ ,]([-\d.]+)\)", el.get("transform") or "")
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)


texts = [el for el in root.iter(f"{NS}text")]
content = lambda el: "".join(el.itertext())


def one(pred):
    hits = [el for el in texts if pred(el)]
    assert len(hits) == 1, (len(hits), [content(h) for h in hits])
    return hits[0]


def set_text(el, new):
    spans = list(el.iter(f"{NS}tspan"))
    if spans:
        spans[0].text = new
        for s in spans[1:]:
            parent[s].remove(s)
    else:
        el.text = new


# ---- E legend: what the boxes mean
set_text(one(lambda e: content(e).strip() == "reported"), "true (TP)")
set_text(one(lambda e: content(e).strip() == "missed"), "missed (FN)")
set_text(one(lambda e: content(e).strip() == "false"), "false (FP)")
# ---- E legend 12 pt to the left, so that the longer labels stay inside the page
def anchor(el):
    tag = el.tag.replace(NS, "")
    if tag == "text":
        return xy(el)
    if tag == "rect":
        return float(el.get("x", 0)), float(el.get("y", 0))
    if tag in ("path", "polygon", "polyline"):
        n = re.findall(r"-?\d+\.?\d*", el.get("d") or el.get("points") or "")
        return (float(n[0]), float(n[1])) if len(n) > 1 else (None, None)
    if tag == "line":
        return float(el.get("x1", 0)), float(el.get("y1", 0))
    return None, None
for el in list(root.iter()):
    if el.tag.replace(NS, "") in ("text", "rect", "path", "polygon", "polyline", "line") and not any(a.tag == f"{NS}text" for a in [parent.get(el)] if a is not None):
        x, y = anchor(el)
        if x is not None and 448 < x < 518 and 285 < y < 350:
            el.set("transform", "translate(-12 0)" + (" " + el.get("transform") if el.get("transform") else ""))

# ---- legend of F and G
el = one(lambda e: "tangent route (local ACR, top)" in content(e))
for node in el.iter():
    for attr in ("text", "tail"):
        v = getattr(node, attr)
        if v and "local ACR" in v:
            setattr(node, attr, v.replace("(local ACR, top)", "(local, top)"))

# ---- C: tick labels in ink, as in F to I
for el in texts:
    x, y = xy(el)
    if x is not None and x > 330 and y < 262 and el.get("fill") == "#9aa3ad" and content(el).strip().isdigit():
        el.set("fill", INK)

# ---- C: the formulas in Latin Modern, right-aligned where the Arial ones ended
FORMULAS = [  # (start of the first Arial run, its baseline, pieces to remove, TeX, colour)
    ("dC/dt = 5 − C − CP", 164.0, ["dC/dt = 5 − C − CP", "/4 = 0"], r"dC/dt = 5 - C - CP^2/4 = 0", BLUE),
    ("dP/dt = CP", 173.5, ["dP/dt = CP", "/4 − P = 0"], r"dP/dt = CP^2/4 - P = 0", ORANGE),
    ("P = 4/C:", 186.0, ["P = 4/C:"], r"P = 4/C\colon", INK),
    ("C", 195.5, ["C", " − 5C + 4 = 0"], r"C^2 - 5C + 4 = 0", INK),
    ("(C − 1)(C − 4) = 0", 205.0, ["(C − 1)(C − 4) = 0"], r"(C - 1)(C - 4) = 0", INK),
]
for first, base, pieces, tex, colour in FORMULAS:
    runs = [el for el in texts if content(el) in pieces and xy(el)[1] is not None and abs(xy(el)[1] - base) < 0.3 and xy(el)[0] > 430]
    sups = [el for el in texts if el.get("font-size") == "4.5" and xy(el)[0] is not None and xy(el)[0] > 430 and 0 < base - xy(el)[1] < 4]
    assert len(runs) == len(pieces), (first, [content(r) for r in runs])
    last = max(runs, key=lambda e: xy(e)[0])
    right = xy(last)[0] + text_width(content(last), "ArialMT", 6)
    g = parent[runs[0]]
    for el in runs + sups:
        parent[el].remove(el)
    new, box = typeset(tex, 6, colour)
    new.set("transform", f"translate({right - box[2]:.3f} {base:.3f})")
    g.append(new)

open(out, "w", encoding="utf-8").write(ET.tostring(root, encoding="unicode"))
print("written", out)
