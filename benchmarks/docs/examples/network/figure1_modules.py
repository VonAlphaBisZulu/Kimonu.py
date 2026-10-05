"""Figure 1 after the terminal-complex fix: B+G joins the kinetic module, A+E and C+E form a second one,
their Y^T.T rows get module chips, and a small panel G sketches the iterated route with T, Y^T.T as its
one-step local readout."""
import re, sys
src, out = sys.argv[1:3]
s = open(src).read()
INK, MUTE, KIN = "#1b2430", "#9aa3ad", "#c9ecf4"
DASH = f'fill="{KIN}" stroke="#000" stroke-dasharray="6 2" stroke-miterlimit="10" stroke-width="1.5"'
# --- B: module outline now also encloses B+G (right column, row 2); D above it stays outside
old = re.search(r'<path d="M451\.147,48\.911h-140\.451[^"]*"[^>]*stroke-dasharray="6 2"[^>]*/>', s)
assert old
outline = ('M310.696,48.911 L451.146,48.911 A12,12 0 0 1 463.146,60.911 L463.146,75.9 A4,4 0 0 0 467.146,79.9 '
           'L509,79.9 A4,4 0 0 1 513,83.9 L513,104 A4,4 0 0 1 509,108 L438.183,108 A4,4 0 0 0 434.183,112 '
           'L434.183,162.887 A12,12 0 0 1 422.181,174.889 L310.696,174.889 A12,12 0 0 1 298.696,162.889 '
           'L298.696,60.911 A12,12 0 0 1 310.696,48.911 Z')
# second module: C+E and A+E (pills x 469.65-508.45, y 115.34-167.08)
second = '<rect x="465.6" y="112.6" width="47" height="58.4" rx="8" ry="8" ' + DASH + '/>'
s = s[:old.start()] + f'<path d="{outline}" {DASH}/>\n  ' + second + s[old.end():]
# --- F: Y^T.T chips for the terminal complexes that join (label baseline - 7.2 = chip top)
chip = lambda y: f'<rect x="364.346" y="{y:.3f}" width="19.624" height="8.6" rx="2" ry="2" fill="{KIN}"/>'
anchor = '<rect x="364.346" y="283.176" width="19.624" height="8.6" rx="2" ry="2" fill="#c9ecf4"/>'
assert s.count(anchor) == 1
s = s.replace(anchor, anchor + "\n    " + "\n    ".join(chip(y) for y in (343.176, 373.176, 403.176)))
# --- G: the iterated route as a loop, and where the one-step readout of F sits
def t(x, y, txt, size=8, fill=INK, bold=False, anchor=None):
    fam = 'font-family="Arial-BoldMT, Arial" font-weight="700"' if bold else 'font-family="ArialMT, Arial"'
    a = f' text-anchor="{anchor}"' if anchor else ""
    return f'<text x="{x}" y="{y}" {fam} font-size="{size}" fill="{fill}"{a} xml:space="preserve">{txt}</text>'
def box(x, y, w, h, fill="#fff", stroke=INK):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" ry="4" fill="{fill}" stroke="{stroke}" stroke-width="1"/>'
def arrow(x1, y1, x2, y2, c=INK, curve=None):
    head = f'<polygon points="{x2},{y2} {x2 - 4.5 * (1 if x2 > x1 else -1)},{y2 - 2.2} {x2 - 4.5 * (1 if x2 > x1 else -1)},{y2 + 2.2}" fill="{c}"/>'
    d = f'M{x1},{y1} Q{curve[0]},{curve[1]} {x2},{y2}' if curve else f'M{x1},{y1} L{x2},{y2}'
    return f'<path d="{d}" fill="none" stroke="{c}" stroke-width="1"/>' + head
g = [f'<text transform="translate(208 438)" fill="{INK}" font-family="Helvetica-Bold, Helvetica" font-size="18" font-weight="700">G</text>',
     t(228, 437, "iterated route (Fig. 3)", 8, INK, True),
     box(212, 452, 74, 30, KIN), t(249, 465, "ker(", 8, INK, anchor="middle"),
     box(330, 452, 74, 30, "#fff"),
     t(249, 465, "", 8),
     ]
g[3] = t(249, 464.5, 'ker(<tspan font-family="Arial-BoldMT, Arial" font-weight="700">W</tspan>)', 8, INK, anchor="middle")
g += [t(249, 475, "linear in ψ", 7, MUTE, anchor="middle"),
      t(367, 464.5, "log closure", 8, INK, anchor="middle"), t(367, 475, "via ln ψ = Yᵀ·ln x", 7, MUTE, anchor="middle"),
      arrow(286, 459, 329, 459), t(307.5, 456, "coupled sets", 6.5, MUTE, anchor="middle"),
      arrow(330, 477, 287, 477, curve=(308.5, 490)), t(308.5, 492.5, "feedback u_c = u_d", 6.5, MUTE, anchor="middle"),
      arrow(404, 467, 422, 467), t(426, 463, "modules,", 8), t(426, 473, "ACR, ACRR", 8), t(426, 482, "(proven)", 7, MUTE),
      t(212, 512, 'one-step local readout (F): <tspan font-family="Arial-BoldMT, Arial" font-weight="700">T</tspan>, '
                  '<tspan font-family="Arial-BoldMT, Arial" font-weight="700">Y</tspan>ᵀ·<tspan font-family="Arial-BoldMT, Arial" font-weight="700">T</tspan>'
                  ' (Pascual-Escudero &amp; Feliu)', 7, MUTE),
      t(212, 521, "proven ⊆ true ⊆ local; equal on 30 of 33 models", 7, MUTE)]
g = [x for x in g if 'y="465"' not in x or "ker(" not in x]
s = s.replace("</svg>", '  <g id="panel_G_iterated_route">' + "".join(g) + "</g>\n</svg>")
open(out, "w").write(s)
