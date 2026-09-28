"""Complex network of the Langary et al. example, in the style of the hand-drawn metabolic panel.

Same canvas (468 x 324), fonts (Arial 24 for nodes, 13 for reactions), 2-unit strokes, full heads for
irreversible reactions, harpoons for reversible pairs with the reverse direction in light gray,
offset by the same distance. One row per linkage class.

    python complex_network.py <out.svg>
"""
import sys

GRAY, GAP, FS_NODE, FS_RXN = "#bbb", 3.6, 24, 13
W = {"A": .667, "B": .667, "C": .722, "D": .722, "E": .667, "F": .611, "G": .778, "H": .722, "I": .278,
     "+": .584, "R": .722, "/": .278, **{str(d): .556 for d in range(10)}}          # Arial advance widths
width = lambda t, fs: fs * sum(W[c] for c in t)
CAP = 0.716                                                                        # Arial cap height

# linkage classes: nodes left to right, and (forward, reverse or None) between neighbours
ROWS = [(["A", "B", "C", "D"], [("R1", "R2"), ("R3", "R4"), ("R5", None)]),
        (["D+E", "F", "B+G"], [("R6", "R7"), ("R8", None)]),
        (["C+G", "H", "C+E"], [("R9", "R10"), ("R11", None)]),
        (["A+G", "I", "A+E"], [("R12", "R13"), ("R14", None)])]
LEFT, RIGHT = 48, 426                # outermost node centres; each row spreads evenly between them
Y = [44, 134, 224, 314]              # node baselines

def full_head(x, y, d):
    """full arrowhead, tip at (x, y), pointing along d = +1 (right) or -1 (left)"""
    return f'<polygon points="{x:.3f} {y:.3f} {x - d * 7:.3f} {y - 4:.3f} {x - d * 7:.3f} {y + 4:.3f}"/>'

def harpoon(x, y, d, barb, fill):
    """one-sided arrowhead: tip at (x, y) pointing along d, barb on side barb (-1 up, +1 down)"""
    return f'<polygon points="{x:.3f} {y:.3f} {x - d * 8:.3f} {y:.3f} {x - d * 8:.3f} {y + barb * 5.5:.3f}" fill="{fill}"/>'

def text(x, y, t, fs, fill=None):
    f = f' fill="{fill}"' if fill else ""
    return f'<text transform="translate({x:.3f} {y:.3f})" font-family="ArialMT, Arial" font-size="{fs}"{f}><tspan x="0" y="0">{t}</tspan></text>'

nodes, fwd, rev, labels = [], [], [], []
for r, (cx, rxns) in enumerate(ROWS):
    y = Y[r]; mid = y - CAP * FS_NODE / 2                       # arrow height: middle of the capitals
    X = [LEFT + i * (RIGHT - LEFT) / (len(cx) - 1) for i in range(len(cx))]
    for i, c in enumerate(cx):
        nodes.append(text(X[i] - width(c, FS_NODE) / 2, y, c, FS_NODE))
    for i, (f, b) in enumerate(rxns):
        x0 = X[i] + width(cx[i], FS_NODE) / 2 + 7                # clear of the node labels
        x1 = X[i + 1] - width(cx[i + 1], FS_NODE) / 2 - 7
        g = [f'<g id="{f}">']
        if b is None:
            g += [f'<line x1="{x0:.3f}" y1="{mid:.3f}" x2="{x1 - 6:.3f}" y2="{mid:.3f}" fill="none" stroke="#000" stroke-width="2"/>',
                  full_head(x1, mid, 1)]
        else:
            yf, yb = mid - GAP / 2, mid + GAP / 2                # forward above, reverse below: the ⇌ glyph
            g += [f'<line x1="{x0:.3f}" y1="{yf:.3f}" x2="{x1 - 6:.3f}" y2="{yf:.3f}" fill="none" stroke="#000" stroke-width="2"/>',
                  harpoon(x1, yf, 1, -1, "#000")]
            rev.append(f'<g id="{b}_reverse"><line x1="{x0 + 6:.3f}" y1="{yb:.3f}" x2="{x1:.3f}" y2="{yb:.3f}" fill="none" stroke="{GRAY}" stroke-width="2"/>'
                       + harpoon(x0, yb, -1, 1, GRAY) + '</g>')
        fwd.append("\n    ".join(g) + "\n  </g>")
        lab = f + ("/" + b if b else ""); lx = (x0 + x1) / 2 - width(lab, FS_RXN) / 2; ly = mid - 9
        t = f'<text transform="translate({lx:.3f} {ly:.3f})" font-family="ArialMT, Arial" font-size="{FS_RXN}"><tspan x="0" y="0">{f}</tspan>'
        if b: t += f'<tspan x="{width(f, FS_RXN):.3f}" y="0" fill="{GRAY}">/{b}</tspan>'
        labels.append(t + "</text>")

svg = ['<?xml version="1.0" encoding="UTF-8"?>',
       '<svg id="complex_network" xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 468 324">',
       '  <g id="reverse_directions">\n    ' + "\n    ".join(rev) + '\n  </g>',
       '  <g id="reactions">\n  ' + "\n  ".join(fwd) + '\n  </g>',
       '  <g id="complexes">\n    ' + "\n    ".join(nodes) + '\n  </g>',
       '  <g id="reaction_labels">\n    ' + "\n    ".join(labels) + '\n  </g>', '</svg>']
open(sys.argv[1], "w").write("\n".join(svg) + "\n")
