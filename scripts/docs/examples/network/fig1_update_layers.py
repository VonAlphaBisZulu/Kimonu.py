"""Figure 1 with flux balance and mass action told apart: warm colours and a "flux balance" pill for everything that
follows from N r = 0 alone, cold colours and a "mass action" pill for everything that needs mass-action kinetics.

    python fig1_update_layers.py <Figure1_draft.svg> <out.svg>

Flux balance (warm): concordance modules, balanced complexes, the kernel columns k1 to k3 and the coupled reactions
of panel C. Mass action (cold): ACR and ACRR groups, kinetic modules. Each equation that introduces a layer and each
panel carries the pill of its layer. The fill of the ker(W) box in G goes, since nothing else shares its meaning.
"""
import os, re, sys
import xml.etree.ElementTree as ET
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "figures"))
import svg_finish
from svg_panels import NS, text_width, typeset

src, out = sys.argv[1:3]
FB, MA, INK = "#c0632b", "#2f6f9f", "#1b2430"
WARM = {"#5a9bd5": "#f9dd6f",       # balanced complexes, yellow (B nodes, legend)
        "#a7cff5": "#fcefb0",       # balanced complexes, lighter yellow (E rows)
        "#fff2cd": "#f7b267",       # concordance module 1, orange
        "#ffe0ff": "#f08a7e",       # concordance module 2, coral
        "#e5d9ff": "#d2a679",       # concordance module 3, tan
        "#3a9d5d": "#c8414b",       # k2, red
        "#56708c": "#7a4a22"}       # k3, brown (k1 stays orange)
COLD = {"#f2c14e": "#5b8fd1",       # ACR, blue (A halo, legend)
        "#fce3a6": "#9dbfe9",       # ACR row in F
        "#d9e7c4": "#c9dcf4",       # ACRR groups, light blue (A halos, F rows, legend)
        "#c9ecf4": "#d3efd4",       # kinetic module 1, light green
        "#cdf3df": "#a8dcb0"}       # kinetic module 2, green

ET.register_namespace("", "http://www.w3.org/2000/svg")
s = re.sub(r"<style.*?</style>", "", open(src, encoding="utf-8").read(), flags=re.S)
root = ET.fromstring(s)
panel = {g.get("id")[6:]: g for g in root.iter(f"{NS}g") if (g.get("id") or "").startswith("panel_")}

for el in root.iter():
    for attr in ("fill", "stroke"):
        v = (el.get(attr) or "").lower()
        if v in WARM:
            el.set(attr, WARM[v])
        elif v in COLD:
            el.set(attr, COLD[v])
for el in panel["E"]:                                   # E's rows take lighter tints, so the grey zeros stay legible
    if el.tag == f"{NS}rect" and el.get("fill") in ("#f7b267", "#f08a7e", "#d2a679"):
        el.set("fill", {"#f7b267": "#fbd3a8", "#f08a7e": "#f8c2bb", "#d2a679": "#e9d2ba"}[el.get("fill")])
for el in panel["G"]:                                   # the ker(W) box loses its kinetic-module fill
    if el.tag == f"{NS}rect" and el.get("fill") == "#d3efd4":
        el.set("fill", "#fff")


def words(g, x, y, t, size, fill, bold=False):
    a = {"transform": f"translate({x:.3f} {y:.3f})", "fill": fill, "font-size": f"{size:g}",
         "font-family": "Arial-BoldMT, Arial" if bold else "ArialMT, Arial"}
    if bold:
        a["font-weight"] = "700"
    el = ET.SubElement(g, f"{NS}text", a)
    ET.SubElement(el, f"{NS}tspan", {"x": "0", "y": "0"}).text = t


def pill(g, x, y, kind):
    """A pill with its top-left corner at (x, y); returns its right edge."""
    t, c = ("flux balance", FB) if kind == "FB" else ("mass action", MA)
    w, h = text_width(t, "Arial-BoldMT", 6) + 6, 7.5
    ET.SubElement(g, f"{NS}rect", {"x": f"{x:.3f}", "y": f"{y:.3f}", "width": f"{w:.3f}", "height": f"{h}", "rx": f"{h / 2}", "ry": f"{h / 2}", "fill": c})
    words(g, x + 3, y + 5.6, t, 6, "#fff", True)
    return x + w


def by_text(g, t):
    return next(c for c in g if c.tag == f"{NS}text" and "".join(c.itertext()).strip() == t)


def move(el, dy):
    x, y = map(float, re.match(r"translate\(([-\d.]+) ([-\d.]+)\)", el.get("transform")).groups())
    el.set("transform", f"translate({x:.3f} {y + dy:.3f})")


# ---- A and B: the first equation is flux balance; mass action takes all lines, because a mass-action model also
#      satisfies flux balance. A gets the mass-action rate in metabolite form, B has it through the complexes.
A, B = panel["A"], panel["B"]
el, _ = typeset(r"\ln\mathbf{r} = \ln\mathbf{k} + \mathbf{S}^\top\mathbf{Y}^\top\ln\mathbf{x}", 9, INK)
el.set("transform", "translate(26.373 36.664)"); A.append(el)
for g, x_fb, x_ma in ((A, 96, None), (B, 412, None)):
    right = pill(g, x_fb, 6.5, "FB") + 4
    ET.SubElement(g, f"{NS}path", {"d": f"M{right},4 h2 v35.5 h-2", "stroke": MA, "stroke-width": "0.75", "fill": "none"})
    pill(g, right + 5, 18, "MA")
pill(A, 233, 171.5, "MA")
for t in ("balanced complexes", "concordance module 1", "concordance module 2", "concordance module 3", "kinetic module 1", "kinetic module 2"):
    B.remove(by_text(B, t))
sw = {c.get("fill"): c for c in B if c.tag == f"{NS}rect" and float(c.get("y")) > 172}
for f, (x, y) in {"#f7b267": (302.083, 195.8), "#f08a7e": (311.5, 195.8), "#d2a679": (320.9, 195.8), "#f9dd6f": (302.083, 209.8),
                  "#d3efd4": (416.462, 195.8), "#a8dcb0": (416.462, 209.8)}.items():
    sw[f].set("x", f"{x:g}"); sw[f].set("y", f"{y:g}")
pill(B, 302.083, 181.5, "FB"); pill(B, 416.462, 181.5, "MA")
words(B, 332, 201.2, "concordance modules", 8, INK)
words(B, 314.1, 215.2, "balanced complexes", 8, INK)
words(B, 428.4, 201.2, "kinetic module 1", 8, INK)
words(B, 428.4, 215.2, "kinetic module 2", 8, INK)
# ---- C to G: the pill of the layer under the panel letter (D's note moves down to make room)
for t in ("examples of", "steady-state fluxes", "that are vectors of the", "nullspace basis"):
    move(by_text(panel["D"], t), 4.5)
for c in panel["D"]:                                    # the K that ends the note moves with it
    if c.tag == f"{NS}text" and "LMRoman" in (c.get("font-family") or "") + "".join(t.get("font-family") or "" for t in c.iter()):
        move(c, 4.5)
for k, (x, y, kind) in {"C": (5.5, 235.5, "FB"), "E": (172, 237.5, "FB"), "F": (313, 245.5, "MA"), "D": (6.5, 421.5, "FB"), "G": (206.5, 417.5, "MA")}.items():
    pill(panel[k], x, y, kind)

text = ET.tostring(root, encoding="unicode")
open(out, "w", encoding="utf-8").write(svg_finish.finish(text)[0])
print("written", out)
