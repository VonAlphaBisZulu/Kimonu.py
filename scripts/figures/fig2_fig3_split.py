"""Split the vault Figure 2 into Figure 2 (the kinase example, panels A-E) and Figure 3 (the 33 genome-scale models and
the certificate, old panels F-I, now A-D), and fix the label typo in the old Figure 3, which becomes Figure 4.

    python figures/fig2_fig3_split.py <vault figures dir> <out dir>      # from scripts/

Every drawn element goes to the panel whose rectangle holds the centre of its bounding box (svg_panels), the bottom
half moves up to the top of its own artboard, and its letters F-I become A-D. Fonts are subset again (svg_finish).
"""
import os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_finish, svg_panels

src, out = sys.argv[1:3]
W = 518.4
CUT = 366.0          # pt, between panel E and panels F, G
SHIFT = 369.5        # pt, moves the letter F (top at 373.5 pt) to the height of the letter A (3.9 pt)
H_OLD = 694.1

TOP = [("A", [(0, 0, 145, 135)]), ("B", [(145, 0, 320, 135)]), ("C", [(320, 0, W, 262)]),
       ("D", [(0, 135, 320, 262)]), ("E", [(0, 262, W, CUT)])]
BOTTOM = [("A", [(0, CUT, 259, 521)]), ("B", [(259, CUT, W, 521)]), ("legend_AB", [(0, 521, W, 533)]),
          ("C", [(0, 533, 206, H_OLD)]), ("D", [(206, 533, 404, H_OLD)]), ("legend_CD", [(404, 533, W, H_OLD)])]


def part(panels, height, shift, relabel, fig_id):
    text = open(os.path.join(src, "Figure2.svg"), encoding="utf-8").read()
    fig = svg_panels.Figure(text)
    keep = []
    for lf in fig.leaves:
        b = lf.bbox()
        cy = None if b is None else (b[1] + b[3]) / 2
        if cy is not None and (cy < CUT) == (shift == 0):
            keep.append(lf)
    fig.leaves = keep
    if shift:
        T = np.array([[1, 0, 0], [0, 1, -shift], [0, 0, 1.0]])
        for lf in fig.leaves:
            lf.M = T @ lf.M
        panels = [(n, [(x0, y0 - shift, x1, y1 - shift) for x0, y0, x1, y1 in r]) for n, r in panels]
    missing = fig.assign(panels)
    assert not missing, [svg_panels.ET.tostring(m.el)[:100] for m in missing]
    for lf in fig.leaves:                                   # panel letters
        if lf.tag() == "text" and (lf.el.get("font-size") == "18"):
            for e in lf.el.iter():
                if e.text and e.text.strip() in relabel:
                    e.text = relabel[e.text.strip()]
    fig.root.set("viewBox", f"0 0 {W} {height}")
    fig.root.set("id", fig_id)
    svg, _ = svg_finish.finish(fig.svg())
    return svg


f2 = part(TOP, CUT, 0, {}, "Figure2")
f3 = part(BOTTOM, round(H_OLD - SHIFT, 1), SHIFT, dict(zip("FGHI", "ABCD")), "Figure3")
for name, svg in (("Figure2.svg", f2), ("Figure3.svg", f3)):
    open(os.path.join(out, name), "w", encoding="utf-8").write(svg)

old3 = open(os.path.join(src, "Figure3.svg"), encoding="utf-8").read()
assert old3.count("ACCR") == 1
f4 = old3.replace("ACCR", "ACRR").replace('id="Figure3"', 'id="Figure4"', 1)
open(os.path.join(out, "Figure4.svg"), "w", encoding="utf-8").write(f4)
for n in ("Figure2.svg", "Figure3.svg", "Figure4.svg"):
    print(n, os.path.getsize(os.path.join(out, n)))
