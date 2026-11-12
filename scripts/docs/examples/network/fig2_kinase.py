"""Figure 2, panels A-E: a network whose answer is known, where the iterated and the tangent route disagree.

A protein X is made, degraded and assembled into an inactive kinase (2X -> C); the active kinase P activates
C cooperatively (C + 2P -> 3P), and C and P are degraded.
A: metabolite view, as Figure 1A. In C + 2P -> 3P, P is consumed twice and produced three times, so its net
   change is +1 and the net stoichiometry (C -> P) hides the autocatalysis. Halos mark true ACR (X).
B: complex view, as Figure 1B, with the two true kinetic modules boxed, {0, X, 2X} and {P, C+2P}, and the
   concentration of every complex as a formula. The iterated route finds only {P, C+2P}.
C: steady states at k = (20, 1, 1/20, 1/4, 1, 1) for R1..R6, chosen so that every steady state is an integer
   (test_integer_example_for_the_completeness_figure). X has one (supply 20 = loss x + x^2/10 at x = 10); in
   the (C, P) plane at x = 10 the nullclines cross twice, at (1, 4) and (4, 1).
D: the two steady states on the network of B: exact flux r* on each arrow, shaft width 2·√(r*/10), and the
   concentration ψ of every complex. X and the fluxes through it are the same in both; C is routed mostly
   into C + 2P -> 3P in one and mostly into decay in the other. ψ_C and ψ_P swap, while the module's ratio
   ψ_{C+2P} / ψ_P = 4 is the same in both.
E: ACR per metabolite, truth against the iterated route, the tangent route (local ACR) and the sign certificate.
A and B follow Phil's layout of 2026-10-06 (fig2AB_layout_phil_2026-10-06.svg) in Figure 1's units: bold 14 pt
nodes, 8 pt labels, 1.2 pt strokes, Figure 1's heads, and harpoons with the reverse direction in #bbb whose
heads sit half a stroke out, flush with the outer edge of the shaft. Every reaction is its own group, so the
drawing can be rearranged in Illustrator. All panels are drawn at final size, 518.4 pt wide. Arial has no
glyph for the empty set or for check marks, so those are paths.

    python fig2_kinase.py <out.svg>
"""
import math, os, sys
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex_text import place                           # formulas as LaTeX outlines, no math font needed

INK, MUTE, GRAY, FONT = "#1b2430", "#9aa3ad", "#bbb", 'font-family="ArialMT, Arial"'
BOLD = 'font-family="Arial-BoldMT, Arial" font-weight="700"'
ACR_FILL, KIN, BLUE, ORANGE = "#f2c14e", "#c9ecf4", "#5a9bd5", "#e6862b"
FALSE_FILL, FALSE_LINE, MISS_LINE, OK = "#f6c9c2", "#c8453b", "#c9a227", "#3a9d5d"
REG = {**dict(zip("ABCDEFGHIJKLMNOPQRSTUVWXYZ", [.667, .667, .722, .722, .667, .611, .778, .722, .278, .5, .667, .556, .833,
                                                  .722, .778, .667, .778, .722, .667, .611, .722, .667, .944, .667, .667, .611])),
       **dict(zip("abcdefghijklmnopqrstuvwxyz", [.556, .556, .5, .556, .556, .278, .556, .556, .222, .222, .5, .222, .833,
                                                  .556, .556, .556, .556, .333, .5, .278, .556, .5, .722, .5, .5, .5])),
       **{str(d): .556 for d in range(10)}, " ": .278, ".": .278, ",": .278, "=": .584, "/": .278, "(": .333, ")": .333,
       "+": .584, "-": .333, ":": .278, "·": .278, "–": .556}                        # Arial advance widths
BLD = {**REG, "A": .722, "B": .722, "P": .667, "C": .722, "X": .667}                             # Arial Bold, for the complexes
width = lambda s, fs, w=REG: fs * sum(w[c] for c in s)

# ---- the model: reactions with the rate constants of the test, and its positive steady states ----
from fractions import Fraction as F
RX = [("0", "X", 20), ("X", "0", 1), ("2X", "C", F(1, 20)), ("C", "0", 1), ("C+2P", "3P", F(1, 4)), ("P", "0", 1)]
# mass action: dX = 20 - X - X^2/10, dC = X^2/20 - C - C P^2/4, dP = C P^2/4 - P
SS = [dict(X=F(10), C=F(1), P=F(4)), dict(X=F(10), C=F(4), P=F(1))]
for v in SS:
    x, c, p = v["X"], v["C"], v["P"]
    assert 20 - x - x * x / 10 == 0 and x * x / 20 - c - c * p * p / 4 == 0 and c * p * p / 4 - p == 0

def t(x, y, s, size=8, fill=INK, bold=False, anchor=None):
    """text at baseline (x, y); anchor "middle" or "end" shifts by the Arial advance width"""
    if anchor == "middle": x -= width(s, size, BLD if bold else REG) / 2
    if anchor == "end": x -= width(s, size, BLD if bold else REG)
    return f'<text transform="translate({x:.2f} {y:.2f})" {BOLD if bold else FONT} font-size="{size}" fill="{fill}"><tspan x="0" y="0">{s}</tspan></text>'
def letter(x, y, L):
    return f'<text transform="translate({x} {y})" fill="{INK}" font-family="Helvetica-Bold, Helvetica" font-size="18" font-weight="700"><tspan x="0" y="0">{L}</tspan></text>'
def rect(x, y, w, h, fill, rx=2, stroke=None, sw=0.8, dash=None):
    st = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    st += f' stroke-dasharray="{dash}"' if dash else ""
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" ry="{rx}" fill="{fill}"{st}/>'
def line(x1, y1, x2, y2, stroke=INK, sw=1.2):
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" fill="none" stroke="{stroke}" stroke-width="{sw}"/>'
def poly(pts, fill=INK):
    return f'<polygon points="{" ".join(f"{x:.2f} {y:.2f}" for x, y in pts)}" fill="{fill}"/>'
def check(x, y):
    return f'<path d="M{x - 3.2:.2f},{y - .2:.2f} L{x - 1:.2f},{y + 2.3:.2f} L{x + 3.4:.2f},{y - 3:.2f}" fill="none" stroke="{OK}" stroke-width="1.3" stroke-linecap="round" stroke-linejoin="round"/>'
def cross(x, y):
    return (f'<path d="M{x - 2.6:.2f},{y - 2.6:.2f} L{x + 2.6:.2f},{y + 2.6:.2f} M{x + 2.6:.2f},{y - 2.6:.2f} L{x - 2.6:.2f},{y + 2.6:.2f}" '
            f'fill="none" stroke="{FALSE_LINE}" stroke-width="1.3" stroke-linecap="round"/>')

# ---- shapes shared by A, B and D, in Figure 1's sizes ----
FULL = {"A": (4.07, 2.35), "B": (4.63, 2.64)}           # full heads: length, half base (equilateral, side 4.7 / 5.29)
HARP = {"A": (4.7, 3.3), "B": (5.29, 3.63)}             # harpoons: length along the shaft, barb
def unit(p, q):
    L = math.hypot(q[0] - p[0], q[1] - p[1]); return ((q[0] - p[0]) / L, (q[1] - p[1]) / L)
def full_head(tip, v, size, col="#000"):
    L, h = FULL[size]; b = (tip[0] - L * v[0], tip[1] - L * v[1]); n = (v[1], -v[0])
    return poly([tip, (b[0] + h * n[0], b[1] + h * n[1]), (b[0] - h * n[0], b[1] - h * n[1])], col)
def arrow(p, q, size, col="#000", w=1.2):
    """straight arrow p -> q with its tip at q; the shaft runs 0.3 pt into the head"""
    v = unit(p, q); L = FULL[size][0]
    return [line(p[0], p[1], q[0] - (L - .3) * v[0], q[1] - (L - .3) * v[1], col, w), full_head(q, v, size, col)]
def pair(p, q, size, w=(1.2, 1.2)):
    """reversible p <-> q: forward harpoon black, reverse #bbb, each on the left of its travel with the barb
    outward. Each shaft runs to its tip, and the head sits half a stroke out from the shaft's centre line,
    so its flat edge closes out with the outer edge of the shaft."""
    o = []
    for a, b, col, wi in ((p, q, "#000", w[0]), (q, p, GRAY, w[1])):
        v = unit(a, b); n = (v[1], -v[0]); c = .59 + wi / 2           # shaft centre: inner edges 1.18 pt apart
        back, barb = HARP[size]
        s0, s1 = (a[0] + c * n[0], a[1] + c * n[1]), (b[0] + c * n[0], b[1] + c * n[1])
        tip = (s1[0] + wi / 2 * n[0], s1[1] + wi / 2 * n[1])
        base = (tip[0] - back * v[0], tip[1] - back * v[1])
        o += [line(s0[0], s0[1], s1[0], s1[1], col, wi), poly([tip, base, (base[0] + barb * n[0], base[1] + barb * n[1])], col)]
    return o
def turn(p, v0, tip, v1, size, col="#000", w=1.2):
    """smooth quarter turn: leaves p along v0, arrives at a full head at tip along v1 (axis-parallel v0, v1)"""
    L = FULL[size][0]; e = (tip[0] - (L - .3) * v1[0], tip[1] - (L - .3) * v1[1])
    d0 = abs((e[0] - p[0]) * v0[0] + (e[1] - p[1]) * v0[1]); d1 = abs((e[0] - p[0]) * v1[0] + (e[1] - p[1]) * v1[1])
    c0 = (p[0] + .552 * d0 * v0[0], p[1] + .552 * d0 * v0[1]); c1 = (e[0] - .552 * d1 * v1[0], e[1] - .552 * d1 * v1[1])
    return [f'<path d="M{p[0]:.2f},{p[1]:.2f} C{c0[0]:.2f},{c0[1]:.2f} {c1[0]:.2f},{c1[1]:.2f} {e[0]:.2f},{e[1]:.2f}" '
            f'fill="none" stroke="{col}" stroke-width="{w}"/>', full_head(tip, v1, size, col)]
def group(name, items): return [f'<g id="{name}">', *items, '</g>']

# ---- B: complex view, as Figure 1B, on Phil's layout of 2026-10-06 ----
PH, PAD = 20.26, 5.8                                   # pill height and side padding of Figure 1B
def pill_w(c): return max(22.46, width(c, 14, BLD) + 2 * PAD) if c != "0" else 22.46
POS = {"2X": (184.1, 19.21), "C": (240.1, 19.21), "X": (180.97, 54.0), "0": (280.2, 54.0), "P": (199.5, 89.63),
       "C+2P": (199.5, 114.43), "3P": (264.8, 114.43)}  # 3P 6 pt right of Phil's, room for R4; rows moved by Phil 2026-10-30; X flush left under 2X
L_ = lambda c, gap=3: POS[c][0] - pill_w(c) / 2 - gap  # just outside a pill's left, right, top and bottom
R_ = lambda c, gap=3: POS[c][0] + pill_w(c) / 2 + gap
T_ = lambda c, gap=2.5: POS[c][1] - PH / 2 - gap
B__ = lambda c, gap=2.5: POS[c][1] + PH / 2 + gap

def pill(c, fill="#fff"):
    x, y = POS[c]; w = pill_w(c)
    o = [rect(x - w / 2, y - PH / 2, w, PH, fill, rx=9.98, stroke="#000", sw=1)]
    if c == "0":                                        # the zero complex: an empty-set sign
        o += [f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.6" fill="none" stroke="#000" stroke-width="1.5"/>',
              line(x - 4.6, y + 6.2, x + 4.6, y - 6.2, "#000", 1.5)]
    else:
        o.append(t(x, y + 5.1, c, 14, "#000", bold=True, anchor="middle"))
    return o

# the true kinetic modules, boxed as in Figure 1 (dashed, filled, one colour per module): {P, C+2P}
# (psi_{C+2P} / psi_P = C·P = k_R6 / k_R4) in blue, and {0, X, 2X} (X has ACR, so psi_X and psi_2X are fixed
# like psi_0 = 1) in green, an L around 2X, X and the zero complex that leaves C out.
MOD2, MOD1 = KIN, "#cdf3df"
mx0, mx1 = POS["C+2P"][0] - pill_w("C+2P") / 2 - 4, POS["C+2P"][0] + pill_w("C+2P") / 2 + 4
my0, my1 = POS["P"][1] - PH / 2 - 4, POS["C+2P"][1] + PH / 2 + 4
PSI = "#2f6fae"                                        # complex concentrations in B and D
def rounded(pts, radii):
    """closed path through pts, corner k rounded by radii[k] (a quadratic through the corner)"""
    d, n = [], len(pts)
    for k in range(n):
        p0, p1, p2 = pts[k - 1], pts[k], pts[(k + 1) % n]
        a = unit(p1, p0); b = unit(p1, p2); r = radii[k]
        ra = min(r, math.dist(p0, p1) / 2); rb = min(r, math.dist(p1, p2) / 2)
        s0 = (p1[0] + ra * a[0], p1[1] + ra * a[1]); s1 = (p1[0] + rb * b[0], p1[1] + rb * b[1])
        d.append(f"{'M' if k == 0 else 'L'}{s0[0]:.2f},{s0[1]:.2f} Q{p1[0]:.2f},{p1[1]:.2f} {s1[0]:.2f},{s1[1]:.2f}")
    return " ".join(d) + " Z"
def module_boxes():
    g = 2.5                                            # pad around 2X; the band around X and the zero complex
    l, t, r2, b2 = POS["2X"][0] - pill_w("2X") / 2 - g, POS["2X"][1] - PH / 2 - g, POS["2X"][0] + pill_w("2X") / 2 + g, 37.5
    r0, bb = POS["0"][0] + pill_w("0") / 2 + g, POS["X"][1] + PH / 2 + g   # band from 37.5: room for R1/R2 inside
    L = [(l, t), (r2, t), (r2, b2), (r0, b2), (r0, bb), (l, bb)]
    st = 'stroke="#000" stroke-width="1.5" stroke-dasharray="6 2"'
    return [f'<path d="{rounded(L, (10, 8, 4, 8, 10, 10))}" fill="{MOD1}" {st}/>',
            rect(mx0, my0, mx1 - mx0, my1 - my0, MOD2, rx=10, stroke="#000", sw=1.5, dash="6 2")]

# Phil's Illustrator edits of 2026-10-30: R5 and R6 run straight and turn down/up just before the zero complex,
# R4 leaves 0.4 pt earlier, A's R4 loop is a rounded rectangle, C moves 66.32 pt left, and D's flux key is gone.
R5_PATH = "M254.33,19.21h4.357c12.847,0,21.513,5.71,21.513,17.83"
R6_PATH = "M213.73,89.63h46.95c16.88,0,19.52-7.71,19.52-18.67"
A_LOOP = "M100.685,67.03c0-17.7-29.195-17.7-29.195,0v13.988c0,10.006,6.375,14.252,14.597,14.252h6.587"
def shaft(d, w): return f'<path d="{d}" fill="none" stroke="#000" stroke-width="{round(w, 3)}"/>'

# where a complex's concentration goes: left of the column, above C, right of 3P and of the zero complex
PSI_AT = lambda: {"2X": (164, POS["2X"][1] + 3.5, "right"), "X": (164, POS["X"][1] + 3.5, "right"),
                  "P": (164, POS["P"][1] + 3.5, "right"), "C+2P": (164, POS["C+2P"][1] + 3.5, "right"),
                  "C": (POS["C"][0], T_("C") - 2, "center"), "3P": (R_("3P") + 1, POS["3P"][1] + 3.5, "left"),
                  "0": (R_("0") + 2, POS["0"][1] + 3.5, "left")}
PSI_TEX = {"2X": r"x_{\mathrm{X}}^2", "X": r"x_{\mathrm{X}}", "C": r"x_{\mathrm{C}}", "P": r"x_{\mathrm{P}}",
           "C+2P": r"x_{\mathrm{C}}\,x_{\mathrm{P}}^2", "3P": r"x_{\mathrm{P}}^3", "0": "1"}

def complex_view(flux=None, psi=None, fs=8, pre="B", boxes=False):
    """the network of B. With flux (one steady state), shafts are 2·√(r*/10) wide and carry r* instead of the
    reaction names, and psi labels every complex with its concentration (numbers; B gets formulas instead)."""
    W = lambda r: 1.2 if flux is None else WIDTH(flux[r])
    name = lambda r: r if flux is None else str(flux[r])
    o = module_boxes() if boxes else []
    yX, y0 = POS["X"][1], POS["0"][1]
    o += group(f"{pre}_R1_R2_zero_X", pair((L_("0"), y0), (R_("X"), yX), "B", (W("R1"), W("R2")))   # R1 0 -> X travels left: black below
               + [t((R_("X") + L_("0")) / 2 - 2.3, yX - 7.6, name("R1"), fs, INK, anchor="end"),
                  t((R_("X") + L_("0")) / 2 - 2.3, yX - 7.6, "/" + name("R2"), fs, GRAY)])
    o += group(f"{pre}_R3_2X_C", arrow((R_("2X", 3.5), POS["2X"][1]), (L_("C"), POS["C"][1]), "B", w=W("R3"))
               + [t((R_("2X", 3.5) + L_("C")) / 2, POS["C"][1] - 4, name("R3"), fs, INK, anchor="middle")])
    o += group(f"{pre}_R5_C_zero", [shaft(R5_PATH, W("R5")), full_head((POS["0"][0], T_("0")), (0, 1), "B")]
               + [t(R_("C") + 13, POS["C"][1] - 4, name("R5"), fs, INK, anchor="middle")])
    o += group(f"{pre}_R6_P_zero", [shaft(R6_PATH, W("R6")), full_head((POS["0"][0], B__("0")), (0, -1), "B")]
               + [t(250, POS["P"][1] - 4, name("R6"), fs, INK, anchor="middle")])
    o += group(f"{pre}_R4_C2P_3P", arrow((mx1 + .75, POS["C+2P"][1]), (L_("3P"), POS["3P"][1]), "B", w=W("R4"))
               + [t((mx1 + .75 + L_("3P")) / 2, POS["3P"][1] - 4, name("R4"), fs, INK, anchor="middle")])
    for c in POS:
        o += group(f"{pre}_complex_{c.replace('+', '_').replace('0', 'zero')}", pill(c))
    if psi is not None:
        anchor = {"right": "end", "center": "middle", "left": None}
        o += [t(x, y, str(psi[c]), fs, PSI, anchor=anchor[a]) for c, (x, y, a) in PSI_AT().items()]
    return o

# B: the true modules boxed, every complex with its concentration as a formula (D gives the numbers)
A_ = complex_view(boxes=True) + place([(x, y, 8, a, [("tex", PSI_TEX[c], PSI)]) for c, (x, y, a) in PSI_AT().items()]) + [
    rect(mx0 + 4.5, 136.2, 7.85, 6.43, MOD1, rx=3.17, stroke="#000", sw=1.5, dash="6 2"),
    rect(mx0 + 15.5, 136.2, 7.85, 6.43, MOD2, rx=3.17, stroke="#000", sw=1.5, dash="6 2"), t(mx0 + 29.4, 141.5, "kinetic modules", 8)]

# ---- A: metabolite view, as Figure 1A, on Phil's layout; also the two copies of D ----
NODE = {"X": (51.455, 34.2), "C": (100.685, 34.2), "P": (100.685, 95.27)}
HALF = {"X": (9.5, 9.5), "C": (5.05, 5.0), "P": (4.67, 5.0)}       # halo radius for X; half glyph for C and P
TRUE_ACR = {"X"}                                       # halo: the true answer, as Figure 1A marks ACR
# Phil's numbering: R1 0 -> X, R2 X -> 0, R3 2X -> C, R4 C + 2P -> 3P, R5 C -> 0, R6 P -> 0
RID = {"R1": 0, "R2": 1, "R3": 2, "R4": 4, "R5": 3, "R6": 5}
def rate(i, v):
    """mass-action rate of RX[i] at concentrations v, exact"""
    s, _, k = RX[i]; r = F(k)
    for term in s.split("+"):
        if term != "0":
            n, m = (int(term[:-1]), term[-1]) if term[0].isdigit() else (1, term)
            r *= v[m] ** n
    return r
FLUX = [{r: rate(i, v) for r, i in RID.items()} for v in SS]
assert [FLUX[0][r] for r in RID] == [20, 10, 5, 4, 1, 4]
assert [FLUX[1][r] for r in RID] == [20, 10, 5, 1, 4, 1]
WIDTH = lambda r: 2 * math.sqrt(r / 10)                # shaft width in pt: 20 -> 2.8, 5 -> 1.4, 1 -> 0.6

def metabolite_view(flux=None, conc=None, fs=8, pre="A"):
    """the network of A. With flux (one steady state), shafts are 2·√(r*/10) wide and carry r* instead of the
    reaction names, and conc labels the species; the stoichiometric numbers are left out."""
    (xX, yR), (xC, _), (_, yP) = NODE["X"], NODE["C"], NODE["P"]
    W = lambda r: 1.2 if flux is None else WIDTH(flux[r])
    name = lambda r: r if flux is None else str(flux[r])
    o = [f'<circle cx="{xX}" cy="{yR}" r="9.5" fill="{ACR_FILL}"/>']
    # R1/R2: uptake and decay of X, a reversible pair from open space (there is no zero complex here)
    o += group(f"{pre}_R1_R2", pair((13.0, yR), (xX - 9.5 - 2, yR), "A", (W("R1"), W("R2"))) + (
        [t(15.2, 47.0, "R1", fs, INK), t(15.2 + width("R1", fs), 47.0, "/R2", fs, GRAY)] if flux is None else
        [t(26.5, 27.6, name("R1"), fs, INK, anchor="middle"), t(26.5, 47.0, name("R2"), fs, GRAY, anchor="middle")]))
    x0, x1 = xX + 9.5 + 2.5, xC - HALF["C"][0] - 2.5
    o += group(f"{pre}_R3", arrow((x0, yR), (x1, yR), "A", w=W("R3")) + [t((x0 + x1) / 2, 47.0, name("R3"), fs, INK, anchor="middle")]
               + ([t(x0 + 1.5, 29.6, "2", fs, INK), t(x1 - 6.5, 29.6, "1", fs, INK)] if flux is None else []))
    xo0, xo1 = xC + HALF["C"][0] + 3, 137.2
    o += group(f"{pre}_R5", arrow((xo0, yR), (xo1, yR), "A", w=W("R5"))
               + [t((xo0 + xo1) / 2, 29.6 if flux is None else 47.0, name("R5"), fs, INK, anchor="middle")])
    o += group(f"{pre}_R6", arrow((xo0, yP), (xo1, yP), "A", w=W("R6")) + [t((xo0 + xo1) / 2, 90.7, name("R6"), fs, INK, anchor="middle")])
    # R4, C + 2P -> 3P: C comes down the vertical; the loop takes 2P from P's left into the junction J, and
    # 3P arrive at P. Two cubics with a common vertical tangent at the loop's far side.
    J, xl, ye = (xC, 67.03), 71.49, (xC - HALF["P"][0] - 3.3, yP)
    loop = (f"M{J[0]},{J[1]} C{J[0]},{J[1] - 17.7} {xl},{J[1] - 17.7} {xl},{J[1]} "
            f"C{xl},{J[1] + .552 * (yP - J[1]):.2f} {ye[0] - .552 * (ye[0] - xl):.2f},{yP} {ye[0]:.2f},{yP}")
    o += group(f"{pre}_R4", arrow((xC, yR + HALF["C"][1] + 3), (xC, yP - HALF["P"][1] - 3), "A", w=W("R4"))
               + [f'<path d="{A_LOOP if flux is None else loop}" fill="none" stroke="#000" stroke-width="{W("R4"):.2f}"/>',
                  t(66.5, 70.0, name("R4"), fs, INK, anchor="end")]
               + ([t(104.4, 52.5, "1", fs, INK), t(104.4, 82.5, "3", fs, INK), t(85.5, 90.0, "2", fs, INK)] if flux is None else []))
    for m, (x, y) in NODE.items():
        o.append(t(x, y + 5.0, m, 14, "#000", bold=True, anchor="middle"))
    if conc is not None:
        o += [t(xX, 20.8, f"x = {conc['X']}", fs, INK, anchor="middle"), t(xC, 23.6, f"c = {conc['C']}", fs, INK, anchor="middle"),
              t(xC, 108.6, f"p = {conc['P']}", fs, INK, anchor="middle")]
    return o

M_ = metabolite_view() + [f'<circle cx="19.0" cy="88.5" r="4.5" fill="{ACR_FILL}"/>', t(26.5, 91.4, "ACR (true)", 8)]

# ---- D: the two steady states on the network of B ----
def complex_conc(c, v):
    """concentration of complex c at metabolite concentrations v, exact"""
    out = F(1)
    for term in c.split("+"):
        if term != "0":
            n, m = (int(term[:-1]), term[-1]) if term[0].isdigit() else (1, term)
            out *= v[m] ** n
    return out
PSIS = [{c: complex_conc(c, v) for c in POS} for v in SS]
assert all(p["C+2P"] / p["P"] == 4 for p in PSIS)       # the module's ratio is the same in both
DS, DY, DX0 = 0.8, 150, 146                            # copies of B at 80 %; 10 pt labels print at 8 pt
D_ = [letter(4.302, 166, "D")]
for k, ox in enumerate((18, 150)):
    D_ += [f'<g id="D_steady_state_{k + 1}" transform="translate({ox - DS * DX0:.2f} {DY}) scale({DS})">',
           *complex_view(FLUX[k], PSIS[k], fs=10, pre=f"D{k + 1}"), t(DX0, 6, f"steady state {k + 1}", 10, INK, bold=True), '</g>']
# key: on arrows the exact flux r* (width grows with its square root, see the caption), at complexes psi
KX, KY = 285, 166
D_ += place([(KX, KY + 59, 8, "left", [("text", "at complexes: concentration ", PSI), ("tex", r"\psi", PSI)])])
D_ += [t(KX, KY + 70, "ratio C+2P / P = 4 in both", 8, MUTE)]


# ---- C: nullclines, stacked in the right column ----
def plot(x0, y0, w, h, xmax, ymax, ticks, curves, points, xlab, ylab, title, yticks=None, ylab_dx=17):
    """data box of w x h pt with its lower left corner at (x0, y0); curves as (points, colour, label, (x, y))"""
    X = lambda v: x0 + v / xmax * w
    Yp = lambda v: y0 - v / ymax * h
    o = [t(x0 + w / 2, y0 - h - 8, title, anchor="middle"),
         line(x0, y0, x0 + w, y0, INK, .6), line(x0, y0, x0, y0 - h, INK, .6)]
    for v in ticks:
        o += [line(X(v), y0, X(v), y0 + 2.5, INK, .6), t(X(v), y0 + 10, f"{v:g}", 7, MUTE, anchor="middle")]
    for v in (ticks if yticks is None else yticks):
        o += [line(x0, Yp(v), x0 - 2.5, Yp(v), INK, .6), t(x0 - 4.5, Yp(v) + 2.5, f"{v:g}", 7, MUTE, anchor="end")]
    o += [t(x0 + w / 2, y0 + 20, xlab, 8, anchor="middle"), t(x0 - ylab_dx, y0 - h / 2 + 3, ylab, 8, anchor="middle")]
    for pts, col, lab, (lx, ly) in curves:
        seg = [(X(a), Yp(b)) for a, b in pts if 0 <= a <= xmax and 0 <= b <= ymax]
        o.append(f'<polyline points="{" ".join(f"{a:.2f},{b:.2f}" for a, b in seg)}" fill="none" stroke="{col}" stroke-width="1.2"/>')
        if lab:
            o.append(t(X(lx), Yp(ly), lab, 7.5, col))
    for (a, b), lab, (dx, dy), anchor in points:
        o += [f'<circle cx="{X(a):.2f}" cy="{Yp(b):.2f}" r="2.4" fill="{INK}" stroke="#fff" stroke-width=".8"/>',
              t(X(a) + dx, Yp(b) + dy, lab, 7, INK, anchor=anchor)]
    return o

grid = lambda lo, hi, n: [lo + (hi - lo) * i / n for i in range(n + 1)]
NX, NW, NH = 337.68, 100, 54                           # data boxes of the two stacked plots (Phil: 66.32 pt left)
B_ = []
# X: supply 20 against loss x + x^2/10, one crossing at x = 10 (both steady states)
B_ += plot(NX, 84, NW, NH, 16, 40, [0, 10], [([(x, 20) for x in grid(0, 16, 2)], ORANGE, "supply 20", (.5, 23.5)),
                                            ([(x, x + x * x / 10) for x in grid(0, 16, 200)], BLUE, "", (3.6, 32))],
           [((10, 20), "1, 2", (4, 9), None)], "X", "rate", "X: one steady state", yticks=[0, 20, 40], ylab_dx=23)
# C, P at x = 10: dC/dt = 0 is c (1 + p^2/4) = 5, dP/dt = 0 is c p = 4; they cross twice
B_ += plot(NX, 184, NW, NH, 6, 6, [0, 2, 4],
           [([(5 / (1 + p * p / 4), p) for p in grid(0, 6, 300)], BLUE, "", (3, 5)),
            ([(4 / p, p) for p in grid(2 / 3, 6, 300)], ORANGE, "", (3, 4))],
           [((1, 4), "1", (4, 2.5), None), ((4, 1), "2", (0, -5), "middle")], "C", "P", "C, P: two steady states")
# the nullclines as formulas: dC/dt = 0 is C (4 + P^2) = 20, dP/dt = 0 is C·P = 4 (the module {P, C+2P})
LXN = NX + NW + 6                                      # right of the plot, where Phil's shift left made room
B_ += place([(NX + 3.6 / 16 * NW, 84 - 32 / 40 * NH, 7.5, "left", [("text", "loss ", BLUE), ("tex", r"x + x^2/10", BLUE)]),
             (LXN, 141, 7.5, "left", [("tex", r"d\mathrm{C}/dt = 0", BLUE)]),
             (LXN, 151, 7.5, "left", [("tex", r"\mathrm{P} = \sqrt{20/\mathrm{C} - 4}", BLUE)]),
             (LXN, 164, 7.5, "left", [("tex", r"d\mathrm{P}/dt = 0", ORANGE)]),
             (LXN, 174, 7.5, "left", [("tex", r"\mathrm{P} = 4/\mathrm{C}", ORANGE)])])

# ---- E: what each route reports, ACR per metabolite (left) and kinetic modules (right) ----
Y0 = 272                                # below D
RH = 13
def chip(x, y, label, fill, w=23, stroke=None, dash=None, text=INK):
    return [rect(x - w / 2, y - 7.3, w, 9.6, fill, rx=2, stroke=stroke, sw=.8, dash=dash), t(x, y, label, 7, text, anchor="middle")]
def cell(x, y, claimed, true, label="ACR", fill=ACR_FILL, w=23):
    """claimed and true: chip; claimed but false: red chip; missed: dashed outline around a dash"""
    if claimed:
        return chip(x, y, label, fill if true else FALSE_FILL, w, stroke=None if true else FALSE_LINE)
    if true:
        return chip(x, y, "–", "none", w, stroke=MISS_LINE, dash="2 1.2", text=MUTE)
    return [t(x, y, "–", 7, MUTE, anchor="middle")]
ROUTES = ("iterated", "tangent", "certificate")
HEAD = {"truth": ("", "truth"), "iterated": ("proven", "iterated"), "tangent": ("local", "tangent"),
        "certificate": ("", "certificate")}
def table(x_label, cols, rows, title, label, fill, w, span=None):
    """rows: (name, truth, {route: claimed}, extra cells[, fill]); cols: {column: x}; span: {route: (first, last,
    fills)} draws one box over those rows (the tangent route merging the two modules); returns svg and the y of
    the check rows"""
    o = [t(x_label - 16, Y0 - 7, title)]
    for c, x in cols.items():
        top, bottom = HEAD.get(c, ("", c))
        o += [t(x, Y0 + 9, bottom, 7, MUTE, anchor="middle")] + ([t(x, Y0 + 1.5, top, 6.5, MUTE, anchor="middle")] if top else [])
    span = span or {}
    for i, (name, truth, calls, extra, *rf) in enumerate(rows):
        y = Y0 + 10 + (i + 1) * RH; f = rf[0] if rf else fill
        o += [t(x_label, y, name, 8, INK, bold=True)] + [t(cols[k], y, v, 7.5, anchor="middle") for k, v in extra.items()]
        o += cell(cols["truth"], y, truth, truth, label, f, w)
        o += sum((cell(cols[r], y, calls[r], truth, label, f, w) for r in ROUTES
                  if not (r in span and span[r][0] <= i <= span[r][1])), [])
    for r, (i0, i1, fills) in span.items():            # one box over several rows, one colour band per row
        x, ya, yb = cols[r], Y0 + 10 + (i0 + 1) * RH - 7.3, Y0 + 10 + (i1 + 1) * RH + 2.3
        band = (yb - ya) / len(fills)
        o += [rect(x - w / 2, ya + k * band, w, band, f, rx=0) for k, f in enumerate(fills)]
        o += [rect(x - w / 2, ya, w, yb - ya, "none", rx=2, stroke=MUTE, sw=.6), t(x, (ya + yb) / 2 + 2.5, label, 7, INK, anchor="middle")]
    ys = Y0 + 10 + (len(rows) + 1) * RH + 4
    o.append(line(x_label - 4, ys - 9.5, max(cols.values()) + 16, ys - 9.5, "#dde1e6", .6))
    for j, prop in enumerate(("correct", "complete")):     # no false claim / nothing missed
        y = ys + j * 12
        for r in ROUTES:
            wrong = [row[0] for row in rows if row[2][r] != row[1] and row[2][r] == (prop == "correct")]
            o.append(cross(cols[r], y - 2.6) if wrong else check(cols[r], y - 2.6))
    return o, ys
C_ = [letter(4.302, Y0 - 5, "E")]
# ACR: per metabolite, with its concentration in both steady states
VAL = {m: (str(SS[0][m]), str(SS[1][m])) for m in "XCP"}      # exact, here integers: 10, 1, 4
ACR_CALLS = {"iterated": set(), "tangent": {"X", "C", "P"}, "certificate": {"X"}}
ACOL = {"1": 62, "2": 94, "truth": 128, "iterated": 160, "tangent": 192, "certificate": 226}
acr_rows = [(m, m == "X", {r: m in ACR_CALLS[r] for r in ROUTES}, {"1": VAL[m][0], "2": VAL[m][1]}) for m in "XCP"]
o, ys = table(36, ACOL, acr_rows, "ACR, per metabolite", "ACR", ACR_FILL, 23)
C_ += o + [t(78, Y0 + 1.5, "steady state", 6.5, MUTE, anchor="middle")]
C_ += [t(122, ys + j * 12, prop, 7.5, INK, anchor="end") for j, prop in enumerate(("correct", "complete"))]
# kinetic modules: the two true ones, and the 17 other pairs of the 7 complexes that the tangent route couples
# too (T = {0} at both steady states, so every complex looks fixed); the certificate proves X, and with it
# the module {0, X, 2X}
MCOL = {"truth": 316, "iterated": 350, "tangent": 384, "certificate": 418}
mod_rows = [("{Ø, X, 2X}", True, {"iterated": False, "tangent": True, "certificate": True}, {}, MOD1),
            ("{P, C+2P}", True, {"iterated": True, "tangent": True, "certificate": True}, {}, MOD2),
            ("17 other pairs", False, {"iterated": False, "tangent": True, "certificate": False}, {})]
o, _ = table(250, MCOL, mod_rows, "Kinetic modules", "coupled", KIN, 30, span={"tangent": (0, 1, (MOD1, MOD2))})
C_ += o
# legend: reported (ACR, coupled), missed, false
LX = 440
C_ += chip(LX + 15, Y0 + 23, "ACR", ACR_FILL) + chip(LX + 15, Y0 + 36, "coupled", KIN, 30) + [t(LX + 34, Y0 + 32.7, "reported", 7.5)]
C_ += chip(LX + 15, Y0 + 49, "–", "none", 23, stroke=MISS_LINE, dash="2 1.2", text=MUTE) + [t(LX + 34, Y0 + 49.2, "missed", 7.5)]
C_ += chip(LX + 15, Y0 + 62, "ACR", FALSE_FILL, stroke=FALSE_LINE) + [t(LX + 34, Y0 + 62.2, "false", 7.5)]

H = Y0 + 95 - RH
svg = ['<?xml version="1.0" encoding="UTF-8"?>',
       f'<svg id="Figure3" xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 518.4 {H}" width="518.4pt" height="{H}pt">',
       '<rect x="0" y="0" width="518.4" height="%d" fill="#fff"/>' % H,
       letter(4.302, 16.991, "A"), letter(151.109, 16.991, "B"), letter(307.68, 16.991, "C"),
       '<g id="A_metabolite_view">', *M_, '</g>', '<g id="B_complex_view">', *A_, '</g>',
       '<g id="C_nullclines">', *B_, '</g>', '<g id="D_steady_states">', *D_, '</g>',
       '<g id="E_acr_table">', *C_, '</g>', '</svg>']
open(sys.argv[1], "w").write("\n".join(svg) + "\n")
