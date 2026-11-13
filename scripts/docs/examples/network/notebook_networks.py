"""Network drawings for the example notebook (examples/kimonu_examples.ipynb), cut from Figures 1 and 2.

Each drawing holds panels A (metabolites and reactions) and B (complexes and reactions) of one figure, side by side,
without the panel letters and the equations above the panels, in two versions.
    fig<n>_network.svg   the plain network, every highlight of a relation and the legend removed
    fig<n>_truth.svg     the ground truth as the paper shows it (ACR, ACRR groups, concordance modules, balanced
                         complexes, kinetic modules) with its legend
The fonts stay embedded as in the published figures (figures/svg_finish.py). The boxes of the elements come from
Inkscape (inkscape --query-all), in px, which are 4/3 of the figures' pt.

    python notebook_networks.py <Figure1.svg> <Figure2.svg> <out dir>   # from scripts/; out dir examples/figures
"""
import copy, csv, io, os, subprocess, sys, tempfile
import xml.etree.ElementTree as ET

NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", NS)
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
PX = 0.75                                   # Inkscape px to pt
MARGIN = 4.0                                # pt around the drawing

# Per figure and panel, in pt. "header": elements ending above this line are the equations and their labels.
# "legend": the box (x0, y0, x1, y1) of the legend. "letter": the panel letter.
RULES = {
    1: {"panel_A": dict(letter="A", header=32, legend=(225, 172, 300, 220)),
        "panel_B": dict(letter="B", header=45, legend=(285, 175, 520, 225))},
    2: {"panel_A": dict(letter="A", header=0, legend=(0, 100, 145, 140)),
        "panel_B": dict(letter="B", header=0, legend=None)},
}
# The plain network: the module boxes go and the complexes become white. The coloured circles of the metabolites
# (ACR, ACRR groups) go as well, so that every metabolite is drawn as D and E are in Figure 1 and C and P in Figure 2.
MODULE_FILLS = {"#d3efd4", "#a8dcb0"}
COMPLEX_FILLS = {"#f9dd6f", "#f7b267", "#d2a679", "#f08a7e"}
METABOLITE_FILLS = {"#c9dcf4", "#5b8fd1", "#92b5e2"}


def boxes(svg_bytes):
    with tempfile.NamedTemporaryFile(suffix=".svg") as f:
        f.write(svg_bytes); f.flush()
        out = subprocess.run(["inkscape", "--query-all", f.name], capture_output=True, text=True, check=True).stdout
    return {r[0]: [float(v) * PX for v in r[1:]] for r in csv.reader(io.StringIO(out)) if len(r) == 5}


def inside(b, box):
    x, y, w, h = b
    return box is not None and x >= box[0] and y >= box[1] and x + w <= box[2] and y + h <= box[3]


def cut(path, fig, plain):
    tree = ET.parse(path)
    root = tree.getroot()
    keep = [g for g in root if g.get("id") in RULES[fig]]
    for e in list(root):                                     # only the defs (fonts) and panels A and B stay
        if e.tag != f"{{{NS}}}defs" and e not in keep:
            root.remove(e)
    tag = {}
    for g in keep:
        for i, e in enumerate(g):
            tag[e] = f"{g.get('id')}__{i}"
            e.set("id", tag[e])
    bb = boxes(ET.tostring(root))
    used = []
    for g in keep:
        rule = RULES[fig][g.get("id")]
        for e in list(g):
            b = bb.get(e.get("id"), [0, 0, 0, 0])
            text = "".join(e.itertext()).strip()
            fills = {x.get("fill") for x in e.iter()}
            gone = (text == rule["letter"] and b[1] < 12) or (b[1] + b[3] < rule["header"])
            if plain:
                gone = gone or inside(b, rule["legend"]) or bool(fills & MODULE_FILLS)
                gone = gone or (e.tag.endswith("circle") and e.get("fill") in METABOLITE_FILLS)
                if e.tag.endswith("rect") and e.get("fill") in COMPLEX_FILLS:
                    e.set("fill", "#fff")
            if gone:
                g.remove(e)
            else:
                if b[2] > 0 and b[3] > 0:
                    used.append(b)
            del e.attrib["id"]
    x0 = min(b[0] for b in used) - MARGIN; y0 = min(b[1] for b in used) - MARGIN
    x1 = max(b[0] + b[2] for b in used) + MARGIN; y1 = max(b[1] + b[3] for b in used) + MARGIN
    root.set("viewBox", f"{x0:.2f} {y0:.2f} {x1 - x0:.2f} {y1 - y0:.2f}")
    root.set("width", f"{x1 - x0:.2f}pt"); root.set("height", f"{y1 - y0:.2f}pt")
    return tree


if __name__ == "__main__":
    fig1, fig2, out = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    for n, path in ((1, fig1), (2, fig2)):
        for plain, name in ((True, "network"), (False, "truth")):
            dst = os.path.join(out, f"fig{n}_{name}.svg")
            cut(path, n, plain).write(dst, xml_declaration=True, encoding="utf-8")
            print("wrote", dst)
