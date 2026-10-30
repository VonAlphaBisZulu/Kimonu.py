"""Update Phil's Illustrator file of Figure 2 from the generator, keeping his layout.

Phil edits Figure 2 in Illustrator, so the composed SVG in the vault is the master for panels A-C and E-I.
This script takes that file and
  1. replaces the two network copies of panel D (D_steady_state_1/2) by those of figure2.py, which draws B's
     network with Phil's shapes and the fluxes and concentrations of each steady state;
  2. removes the kinetic-module highlight from B (box and legend line): the figure marks only ground truth,
     and the iterated route's one module {P, C+2P} is not all of it ({0, X, 2X} is a module too, since X
     has ACR);
  3. typesets the formula runs of panels C and D with LaTeX (latex_text.py), so the figure needs no math
     glyph that Arial lacks.

    python figure2_update.py <phil.svg> <toy.svg from figure2.py> <out.svg>
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex_text import replace

BLUE, ORANGE, PSI = "#5a9bd5", "#e6862b", "#2f6fae"
F = [
    (40.8, (355, 366), 7.5, "left", [("text", "loss ", BLUE), ("tex", r"x + x^2/10", BLUE)]),
    (139, (380, 395), 7.5, "left", [("tex", r"d\mathrm{C}/dt = 0", BLUE)]),
    (148, (380, 395), 7.5, "left", [("tex", r"d\mathrm{P}/dt = 0", ORANGE)]),
    (156.1, (380, 395), 7.5, "left", [("tex", r"\mathrm{C}\cdot\mathrm{P} = 4", ORANGE)]),
    (225, (280, 290), 8, "left", [("text", "at complexes: concentration ", PSI), ("tex", r"\psi", PSI)]),
]


def group(s, gid):
    """(start, end) of the element <g id="gid"> ... </g>, nesting respected"""
    i = s.index(f'<g id="{gid}"')
    depth = 0
    for m in re.finditer(r"<g\b|</g>", s[i:]):
        depth += 1 if m.group(0) == "<g" else -1
        if depth == 0:
            return i, i + m.end()


phil, toy, out = sys.argv[1:4]
s, t = open(phil, encoding="utf-8").read(), open(toy, encoding="utf-8").read()
for k in (1, 2):                                       # 1) panel D from the generator
    a, b = group(s, f"D_steady_state_{k}"); c, d = group(t, f"D_steady_state_{k}")
    s = s[:a] + t[c:d] + s[b:]
a, b = group(s, "B_complex_view")                       # 2) no module highlight in B
B = s[a:b]
B, n_box = re.subn(r'\s*<rect\b[^>]*stroke-dasharray[^>]*/>', "", B)
B, n_leg = re.subn(r'\s*<g isolation="isolate">\s*<text\b[^>]*>(?:(?!</text>).)*kinetic module \(iterated route\)(?:(?!</g>).)*</g>', "", B, flags=re.S)
assert (n_box, n_leg) == (2, 1), (n_box, n_leg)
s = s[:a] + B + s[b:]
s, n = replace(s, F)                                    # 3) formulas
open(out, "w", encoding="utf-8").write(s)
print(f"D replaced, module box and legend removed, {len(F)} formulas ({n} runs)")
