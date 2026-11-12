"""Figures 1 and 2 (Phil's Illustrator files): one group per panel, font sizes on the figure scale, formulas as
live Latin Modern text, then svg_finish (fonts embedded, no clipping).

Size scale, shared with Figure 3 (figures/fig4_figS1_S2_benchmark.py):
    10 bold   panel letters (Science: part labels 10 pt bold; 18 pt before 2026-11-12)
    14 bold   network nodes; 11 bold for the small copies in Figure 2D
    8 or 9    panel titles, diagram and table labels
    7         axis labels, point labels
    6         tick labels, legends, notes (the grey annotations)
    scripts   70 % of their base, rounded up to half points (4.5 under 6 pt)
Formulas take the size of the text around them. Figure 2's formulas typed in Arial stay Arial.

The formula outlines (pdflatex glyphs from figure1_latex.py and fig2_kinase.py) are found by their boxes below, in
the coordinates of the current files; each is replaced by the same formula set by matplotlib (svg_panels.typeset).

    python fig1_fig2_panels.py <in Figure1.svg> <in Figure2.svg> <out dir>
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "figures"))
import numpy as np
import xml.etree.ElementTree as ET
import svg_finish
from svg_panels import Figure, Leaf, NS, text_width

W = 518.4
INK, MUTE, PSI, BLUE = "#1b2430", "#9aa3ad", "#2f6fae", "#5a9bd5"
GREY = {"#9aa3ad", "#6b7480", "#8a929b"}

PANELS = {
    1: [("A", [(0, 0, 285, 216)]), ("B", [(285, 0, W, 216)]), ("C", [(0, 216, 170, 400)]), ("E", [(170, 216, 310, 400)]),
        ("F", [(310, 216, W, 400), (425, 400, W, 440)]), ("D", [(0, 400, 205, 538)]), ("G", [(205, 400, W, 538)])],
    2: [("A", [(0, 0, 145, 140)]), ("B", [(145, 0, 330, 140)]), ("C", [(330, 0, W, 268)]), ("D", [(0, 140, 330, 268)]),
        ("E", [(0, 268, W, 362)]), ("legend_FG", [(0, 520, W, 538)]), ("F", [(0, 362, 248, 535)]), ("G", [(248, 362, W, 535)]),
        ("H", [(0, 535, 260, 694.1)]), ("I", [(260, 535, W, 694.1)])],
}

# (box of the outlines, TeX, size, colour, aligned on the outlines' left / center / right, label it follows)
FORMULAS = {
    1: [((26, 5, 91, 15), r"d\mathbf{x}/dt = \mathbf{N}\mathbf{r} = \mathbf{0}", 9, INK, "left", None),
        ((311, 5, 408, 15), r"d\mathbf{x}/dt = \mathbf{N}\mathbf{r} = \mathbf{Y}\mathbf{A}\mathbf{r} = \mathbf{0}", 9, INK, "left", None),
        ((27, 17.5, 47, 25), r"\mathbf{r} \geq \mathbf{0}", 9, INK, "left", None),
        ((313, 15.5, 376, 26), r"\mathbf{r} = \mathrm{diag}(\mathbf{k})\,\mathbf{S}^\top\boldsymbol{\psi}", 9, INK, "left", None),
        ((311, 28.5, 370, 39.5), r"\ln\boldsymbol{\psi} = \mathbf{Y}^\top\ln\mathbf{x}", 9, INK, "left", None),
        ((68, 237, 120, 248), r"\mathbf{K} = \ker(\mathbf{N})", 10, INK, "center", None),
        ((230, 242, 272, 250), r"\mathbf{M} = \mathbf{A}\mathbf{K}", 10, INK, "center", None),
        ((394, 238.5, 451, 250.5), r"\mathbf{B} = \ker(\mathbf{V}^\top)", 10, INK, "center", None),
        ((358, 266, 396, 275), r"\ker(\mathbf{W}\mathbf{Y}^\top)", 7, INK, "center", None),
        ((404.5, 267.5, 410, 275), r"\subseteq", 7, INK, "center", None),
        ((424.5, 267.5, 431, 275), r"\mathbf{B}", 7, INK, "center", None),
        ((446, 267.5, 452, 275), r"\subseteq", 7, INK, "center", None),
        ((457, 266, 501.5, 276), r"\ker((\mathbf{Y}\Delta)^\top)", 7, INK, "center", None),
        ((364, 284, 371, 291), r"\mathbf{B}", 8, INK, "center", None),
        ((459, 281, 480, 289.5), r"\mathbf{Y}^\top\mathbf{B}", 8, INK, "center", None),
        ((367, 408, 374, 414.5), r"\mathbf{B}", 8, INK, "center", None),
        ((235, 438.5, 262.5, 448), r"\ker(\mathbf{W})", 8, INK, "center", None),
        ((430.5, 443.5, 495.5, 454.5), r"\ker((\mathbf{Y}\Delta)^\top) \supseteq \mathbf{B}", 8, INK, "left", None),
        ((76, 450.5, 83, 457), r"\mathbf{K}", 6, MUTE, "left", "nullspace basis"),
        ((218, 449, 281, 458), r"\mathbf{W} = \mathbf{N}\,\mathrm{diag}(\mathbf{r}^*)\,\mathbf{S}^\top", 7, INK, "center", None),
        ((347, 449, 398, 458), r"\ln\boldsymbol{\psi} = \mathbf{Y}^\top\ln\mathbf{x}", 7, INK, "center", None),
        ((314.5, 478, 341, 483), r"u_c = u_d", 6, MUTE, "left", "feedback"),
        ((430, 501.5, 489.5, 511.5), r"\ker(\mathbf{W}\mathbf{Y}^\top) \subseteq \mathbf{B}", 8, INK, "left", None),
        ((308, 507, 315.5, 513.5), r"\mathbf{r}^*", 8, INK, "left", None),
        ((324, 505.5, 370, 514), r"\mathbf{W}\mathbf{Y}^\top\mathbf{h} = \mathbf{0}", 8, INK, "left", None)],
    2: [((251, 1.5, 260.5, 7.5), r"x_{\mathrm{C}}", 8, PSI, "center", None),
        ((169.5, 11.5, 179.5, 21.5), r"x_{\mathrm{X}}^2", 8, PSI, "right", None),
        ((169.5, 53.5, 179.5, 59), r"x_{\mathrm{X}}", 8, PSI, "right", None),
        ((313.5, 52, 317, 57), r"1", 8, PSI, "left", None),
        ((169.5, 89, 179, 94.5), r"x_{\mathrm{P}}", 8, PSI, "right", None),
        ((158, 110.5, 178.5, 120.5), r"x_{\mathrm{C}}\,x_{\mathrm{P}}^2", 8, PSI, "right", None),
        ((299, 110.5, 308.5, 120.5), r"x_{\mathrm{P}}^3", 8, PSI, "left", None),
        ((432.5, 37.5, 467, 46.5), r"x + x^2/10", 6, BLUE, "left", "loss"),
        ((44.5, 216.5, 50, 224), r"\psi", 8, PSI, "left", None)],
}


def sizes(n):
    def rule(lf, panel, text, size, face):
        fill = (lf.el.get("fill") or "").lower()
        if size == 18 or size == 14:
            return size
        if n == 1:
            if size == 6.5 or (size == 7 and fill in GREY):
                return 6                                        # notes
            return size
        table = {4.6: 4.5, 4.755: 4.5, 5.747: 6, 6.27: 6, 6.324: 6, 6.792: 6, 7.315: 7, 7.5: 8, 11.2: 11}
        if size in table:
            return table[size]
        if panel == "C" and size == 8:                          # the two plots of C, as the plots of F-I
            t = text.strip()
            if t.replace(".", "").isdigit():
                return 6                                        # ticks
            if t in ("X", "rate", "P", "C", "A, B", "A", "B"):
                return 7                                        # axis and point labels
            if t == "supply 20":
                return 6                                        # curve label, as "loss"
        return size
    return rule


# Figure 2F, G as Figure 3B, C: "count" on the y axis, the old axis label as a centred title above the plot
# (plot frame x0, x1 and top in pt; frames from the gridlines and spines of the current file)
BAR_PLOTS = [("F", "coupled complex pairs", (48.3, 258.2, 383.0, 473.0)),
             ("G", "ACR metabolites + ACRR pairs (free)", (299.8, 509.7, 383.0, 473.0))]
TITLE_SIZE = 8                                              # plot titles, as C's "X: equal at both steady states"


def bar_titles(F):
    texts = {(F.panel_of[id(l)], "".join(r[3] for r in l.runs()).strip()): l for l in F.leaves if l.tag() == "text"}
    ref = texts[("C", "X: equal at both steady states")]
    for panel, label, (x0, x1, top, bottom) in BAR_PLOTS:
        lab = texts[(panel, label)]
        # the y label becomes "count", centred on the frame as before (rotated text reads upwards from its origin)
        run = next(e for e in lab.el.iter() if e.text and e.text.strip())
        size = float(lab.el.get("font-size")); run.text = "count"
        w = text_width("count", "ArialMT", size)
        O = lab.M[:2, 2]; d = lab.M[:2, 0] / np.hypot(*lab.M[:2, 0])
        mid = np.array([O[0], (top + bottom) / 2]) if abs(d[0]) < 1e-6 else O
        start = mid - d * w / 2
        lab.M = lab.M.copy(); lab.M[0, 2], lab.M[1, 2] = start
        # the old label as a title above the frame, in the size and colour of the figure's other plot titles
        tw = text_width(label, "ArialMT", TITLE_SIZE)
        el = ET.Element(f"{NS}text", {"fill": ref.el.get("fill") or "#1b2430", "font-family": "ArialMT, Arial",
                                      "font-size": f"{TITLE_SIZE:g}"})
        ET.SubElement(el, f"{NS}tspan", {"x": "0", "y": "0"}).text = label
        t = Leaf(el, np.array([[1, 0, (x0 + x1) / 2 - tw / 2], [0, 1, top - 4.0], [0, 0, 1]]), lab.order)
        t.inherited = {}
        F.leaves.append(t); F.panel_of[id(t)] = panel


def build(src, n):
    F = Figure(open(src, encoding="utf-8").read())
    miss = F.assign(PANELS[n])
    assert not miss, [m.bbox() for m in miss]
    shift = F.resize(sizes(n))
    done = F.replace_formulas(FORMULAS[n], shift)
    if n == 2:
        bar_titles(F)
    out, used = svg_finish.finish(F.svg())
    return out, done, used


if __name__ == "__main__":
    f1, f2, outdir = sys.argv[1:4]
    for n, src in ((1, f1), (2, f2)):
        out, done, used = build(src, n)
        open(os.path.join(outdir, f"Figure{n}_draft.svg"), "w", encoding="utf-8").write(out)
        print(f"Figure {n}: {len(done)} formulas", {k: len(v) for k, v in used.items()})
        for d in done:
            print("   ", d)
