"""Figure 3, panels A-E: a network whose answer is known, where the iterated and the tangent route disagree.

A protein X is made, degraded and assembled into an inactive kinase (2X -> C); the active kinase P activates
C cooperatively (C + 2P -> 3P), and C and P are degraded.
A: metabolite view, as Figure 1A. In C + 2P -> 3P, P is consumed twice and produced three times, so its net
   change is +1 and the net stoichiometry (C -> P) hides the autocatalysis. Halos mark true ACR (X).
B: complex view, as Figure 1B. The iterated route's one kinetic module, {P, C+2P}, is outlined; it fixes
   C·P, not C or P.
C: steady states at the rate constants of test_natural_candidates_for_the_completeness_figure. X has one
   (supply 2 = loss x + x^2 at x = 1); in the (C, P) plane at x = 1 the nullclines cross twice, at
   (1/10, 2/5) and (2/5, 1/10).
D: the two steady states on the network of A: exact flux r* on each arrow, shaft width 2·√r*, and the
   concentrations. X and the fluxes through it are the same in both; C is routed mostly into C + 2P -> 3P
   in one and mostly into decay in the other.
E: ACR per metabolite, truth against the iterated route, the tangent route (local ACR) and the sign certificate.
A and B follow Phil's layout of 2026-10-06 (figure3_phil_2026-10-06.svg) in Figure 1's units: bold 14 pt
nodes, 8 pt labels, 1.2 pt strokes, Figure 1's heads, and harpoons with the reverse direction in #bbb whose
heads sit half a stroke out, flush with the outer edge of the shaft. Every reaction is its own group, so the
drawing can be rearranged in Illustrator. All panels are drawn at final size, 518.4 pt wide. Arial has no
glyph for the empty set or for check marks, so those are paths.

    python figure3.py <out.svg>
"""
import math, sys
from scipy.optimize import brentq

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
RX = [("0", "X", 2), ("X", "0", 1), ("2X", "C", F(1, 2)), ("C", "0", 1), ("C+2P", "3P", 25), ("P", "0", 1)]
# mass action: dX = 2 - X - X^2, dC = X^2/2 - C - 25 C P^2, dP = 25 C P^2 - P
SS = [dict(X=F(1), C=F(1, 10), P=F(2, 5)), dict(X=F(1), C=F(2, 5), P=F(1, 10))]
for v in SS:
    x, c, p = v["X"], v["C"], v["P"]
    assert 2 - x - x * x == 0 and x * x / 2 - c - 25 * c * p * p == 0 and 25 * c * p * p - p == 0

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
POS = {"2X": (184.1, 26.7), "C": (240.1, 26.7), "X": (199.5, 54.0), "0": (280.2, 54.0), "P": (199.5, 83.5),
       "C+2P": (199.5, 108.3), "3P": (264.8, 108.3)}     # 3P 6 pt right of Phil's, room for R4
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

A_ = []
# the iterated route's kinetic module {P, C+2P}: psi_{C+2P} / psi_P = C·P is fixed, so C·P = k_R4 / k_R6
mx0, mx1 = POS["C+2P"][0] - pill_w("C+2P") / 2 - 6, POS["C+2P"][0] + pill_w("C+2P") / 2 + 6
my0, my1 = POS["P"][1] - PH / 2 - 6, POS["C+2P"][1] + PH / 2 + 6
A_.append(rect(mx0, my0, mx1 - mx0, my1 - my0, KIN, rx=12, stroke="#000", sw=1.5, dash="6 2"))
yX, y0 = POS["X"][1], POS["0"][1]
A_ += group("B_R1_R2_zero_X", pair((L_("0"), y0), (R_("X"), yX), "B")      # R1 0 -> X travels left: black below
            + [t((R_("X") + L_("0")) / 2 - 2.3, yX - 7.6, "R1", 8, INK, anchor="end"), t((R_("X") + L_("0")) / 2 - 2.3, yX - 7.6, "/R2", 8, GRAY)])
A_ += group("B_R3_2X_C", arrow((R_("2X", 3.5), POS["2X"][1]), (L_("C"), POS["C"][1]), "B")
            + [t((R_("2X", 3.5) + L_("C")) / 2, POS["C"][1] - 4, "R3", 8, INK, anchor="middle")])
A_ += group("B_R5_C_zero", turn((R_("C"), POS["C"][1]), (1, 0), (POS["0"][0], T_("0")), (0, 1), "B")
            + [t(R_("C") + 13, POS["C"][1] - 4, "R5", 8, INK, anchor="middle")])
A_ += group("B_R6_P_zero", turn((R_("P"), POS["P"][1]), (1, 0), (POS["0"][0], B__("0")), (0, -1), "B")
            + [t(250, POS["P"][1] - 4, "R6", 8, INK, anchor="middle")])
A_ += group("B_R4_C2P_3P", arrow((R_("C+2P"), POS["C+2P"][1]), (L_("3P"), POS["3P"][1]), "B")
            + [t((mx1 + .75 + L_("3P")) / 2, POS["3P"][1] - 4, "R4", 8, INK, anchor="middle")])
for c in POS:
    A_ += group(f"B_complex_{c.replace('+', '_').replace('0', 'zero')}", pill(c))
A_ += [rect(mx0 + 4.5, 136.2, 7.85, 6.43, KIN, rx=3.17, stroke="#000", sw=1.5, dash="6 2"),
       t(mx0 + 18.4, 141.5, "kinetic module (iterated route)", 8)]

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
assert [FLUX[0][r] for r in RID] == [2, 1, F(1, 2), F(2, 5), F(1, 10), F(2, 5)]
assert [FLUX[1][r] for r in RID] == [2, 1, F(1, 2), F(1, 10), F(2, 5), F(1, 10)]

def metabolite_view(flux=None, conc=None, fs=8, pre="A"):
    """the network of A. With flux (one steady state), shafts are 2·√r* wide and carry r* instead of the
    reaction names, and conc labels the species; the stoichiometric numbers are left out."""
    (xX, yR), (xC, _), (_, yP) = NODE["X"], NODE["C"], NODE["P"]
    W = lambda r: 1.2 if flux is None else 2 * math.sqrt(flux[r])
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
               + [f'<path d="{loop}" fill="none" stroke="#000" stroke-width="{W("R4"):.2f}"/>',
                  t(66.5, 70.0, name("R4"), fs, INK, anchor="end")]
               + ([t(104.4, 52.5, "1", fs, INK), t(104.4, 82.5, "3", fs, INK), t(85.5, 90.0, "2", fs, INK)] if flux is None else []))
    for m, (x, y) in NODE.items():
        o.append(t(x, y + 5.0, m, 14, "#000", bold=True, anchor="middle"))
    if conc is not None:
        o += [t(xX, 20.8, f"x = {conc['X']}", fs, INK, anchor="middle"), t(xC, 23.6, f"c = {conc['C']}", fs, INK, anchor="middle"),
              t(xC, 108.6, f"p = {conc['P']}", fs, INK, anchor="middle")]
    return o

M_ = metabolite_view() + [f'<circle cx="19.0" cy="88.5" r="4.5" fill="{ACR_FILL}"/>', t(26.5, 91.4, "ACR (true)", 8)]

# ---- D: the two steady states on the network of A ----
DS, DY = 0.8, 140.8                                    # copies of A at 80 %; 10 pt labels print at 8 pt
D_ = [letter(4.302, 166, "D")]
for k, ox in enumerate((18, 150)):
    D_ += [f'<g id="D_steady_state_{k + 1}" transform="translate({ox} {DY}) scale({DS})">',
           *metabolite_view(FLUX[k], SS[k], fs=10, pre=f"D{k + 1}"), t(4, 99.5, f"steady state {k + 1}", 10, INK, bold=True), '</g>']
# key: numbers on arrows are the steady-state flux r*, exact; shaft width grows with its square root
KX, KY = 285, 166
D_ += [t(KX, KY, "on arrows: flux r*", 8), t(KX, KY + 11, "width ∝ √r*", 8, MUTE)]
for j, r in enumerate((F(2), F(1, 2), F(1, 10))):
    yy = KY + 23 + j * 9
    D_ += [line(KX, yy, KX + 20, yy, "#000", DS * 2 * math.sqrt(r)), t(KX + 25, yy + 2.6, str(r), 8)]
D_ += [t(KX, KY + 59, "x, c, p: concentrations", 8, MUTE)]


# ---- C: nullclines, stacked in the right column ----
def plot(x0, y0, w, h, xmax, ymax, ticks, curves, points, xlab, ylab, title, yticks=None):
    """data box of w x h pt with its lower left corner at (x0, y0); curves as (points, colour, label, (x, y))"""
    X = lambda v: x0 + v / xmax * w
    Yp = lambda v: y0 - v / ymax * h
    o = [t(x0 + w / 2, y0 - h - 8, title, anchor="middle"),
         line(x0, y0, x0 + w, y0, INK, .6), line(x0, y0, x0, y0 - h, INK, .6)]
    for v in ticks:
        o += [line(X(v), y0, X(v), y0 + 2.5, INK, .6), t(X(v), y0 + 10, f"{v:g}", 7, MUTE, anchor="middle")]
    for v in (ticks if yticks is None else yticks):
        o += [line(x0, Yp(v), x0 - 2.5, Yp(v), INK, .6), t(x0 - 4.5, Yp(v) + 2.5, f"{v:g}", 7, MUTE, anchor="end")]
    o += [t(x0 + w / 2, y0 + 20, xlab, 8, anchor="middle"), t(x0 - 17, y0 - h / 2 + 3, ylab, 8, anchor="middle")]
    for pts, col, lab, (lx, ly) in curves:
        seg = [(X(a), Yp(b)) for a, b in pts if 0 <= a <= xmax and 0 <= b <= ymax]
        o.append(f'<polyline points="{" ".join(f"{a:.2f},{b:.2f}" for a, b in seg)}" fill="none" stroke="{col}" stroke-width="1.2"/>')
        o.append(t(X(lx), Yp(ly), lab, 7.5, col))
    for (a, b), lab, (dx, dy), anchor in points:
        o += [f'<circle cx="{X(a):.2f}" cy="{Yp(b):.2f}" r="2.4" fill="{INK}" stroke="#fff" stroke-width=".8"/>',
              t(X(a) + dx, Yp(b) + dy, lab, 7, INK, anchor=anchor)]
    return o

grid = lambda lo, hi, n: [lo + (hi - lo) * i / n for i in range(n + 1)]
NX, NW, NH = 404, 100, 54                              # data boxes of the two stacked plots
B_ = []
# X: supply 2 against loss x + x^2, one crossing at x = 1 (both steady states)
B_ += plot(NX, 84, NW, NH, 1.6, 4, [0, 1], [([(x, 2) for x in grid(0, 1.6, 2)], ORANGE, "supply 2", (.05, 2.35)),
                                            ([(x, x + x * x) for x in grid(0, 1.6, 200)], BLUE, "loss x + x²", (.45, 3.2))],
           [((1, 2), "1, 2", (4, 9), None)], "X", "rate", "X: one steady state", yticks=[0, 2, 4])
# C, P at x = 1: dC/dt = 0 is c (1 + 25 p^2) = 1/2, dP/dt = 0 is 25 c p = 1; they cross twice
B_ += plot(NX, 184, NW, NH, .6, .6, [0, .2, .4],
           [([(.5 / (1 + 25 * p * p), p) for p in grid(0, .6, 300)], BLUE, "dC/dt = 0", (.3, .5)),
            ([(1 / (25 * p), p) for p in grid(1 / 15, .6, 300)], ORANGE, "dP/dt = 0", (.3, .4))],
           [((.1, .4), "1", (4, 2.5), None), ((.4, .1), "2", (0, -5), "middle")], "C", "P", "C, P: two steady states")
# the orange curve is where the iterated route's kinetic module holds: C·P = 1/25 at every steady state
B_.append(t(NX + .3 / .6 * NW, 184 - .31 / .6 * NH, "C·P = 1/25", 7.5, ORANGE))

# ---- E: ACR per metabolite ----
Y0 = 255
C_ = [letter(4.302, Y0 - 5, "E"), t(24, Y0 - 7, "ACR, per metabolite")]
COLS = [("steady state 1", 84), ("steady state 2", 140), ("truth", 196), ("iterated", 250), ("tangent", 304),
        ("certificate", 358)]
for name, x in COLS:
    C_.append(t(x, Y0 + 6, name, 7, MUTE, anchor="middle"))
C_.append(t(304, Y0 - 2, "local ACR", 6.5, MUTE, anchor="middle"))   # the tangent route decides local ACR
C_.append(t(250, Y0 - 2, "proven", 6.5, MUTE, anchor="middle"))
ROWS = "XCP"
VAL = {m: (str(SS[0][m]), str(SS[1][m])) for m in ROWS}      # exact fractions: 1, 1/10, 2/5
TRUTH = {"X": True, "C": False, "P": False}
CALLS = {"iterated": set(), "tangent": {"X", "C", "P"}, "certificate": {"X"}}
RH = 13
def cell(x, y, claimed, true):
    """ACR claimed: amber chip; claimed but false: red chip; missed: dashed outline around a dash"""
    if claimed:
        o = [rect(x - 11.5, y - 7.3, 23, 9.6, ACR_FILL if true else FALSE_FILL, rx=2,
                  stroke=None if true else FALSE_LINE)]
        return o + [t(x, y, "ACR", 7, INK, anchor="middle")]
    if true:
        return [rect(x - 11.5, y - 7.3, 23, 9.6, "none", rx=2, stroke=MISS_LINE, sw=.8, dash="2 1.2"),
                t(x, y, "–", 7, MUTE, anchor="middle")]
    return [t(x, y, "–", 7, MUTE, anchor="middle")]
for i, m in enumerate(ROWS):
    y = Y0 + 7 + (i + 1) * RH
    C_ += [t(40, y, m, 8, INK, bold=True, anchor="middle"),
           t(84, y, VAL[m][0], 7.5, anchor="middle"), t(140, y, VAL[m][1], 7.5, anchor="middle")]
    C_ += cell(196, y, TRUTH[m], TRUTH[m])
    for name, x in COLS[3:]:
        C_ += cell(x, y, m in CALLS[name], TRUTH[m])
# the two properties the figure is about: no false ACR claims (correct), no missed ACR (complete)
ys = Y0 + 7 + (len(ROWS) + 1) * RH + 4
C_.append(line(30, ys - 9.5, 374, ys - 9.5, "#dde1e6", .6))
for j, prop in enumerate(("correct", "complete")):
    y = ys + j * 12
    C_.append(t(222, y, prop, 7.5, INK, anchor="end"))
    for name, x in COLS[3:]:
        wrong = [m for m in ROWS if (m in CALLS[name]) != TRUTH[m] and (m in CALLS[name]) == (prop == "correct")]
        C_.append(cross(x, y - 2.6) if wrong else check(x, y - 2.6))
# legend
LX = 408
C_ += [rect(LX, Y0 + 12.7, 23, 9.6, ACR_FILL, rx=2), t(LX + 11.5, Y0 + 20, "ACR", 7, INK, anchor="middle"),
       t(LX + 28, Y0 + 20.2, "reported", 8),
       rect(LX, Y0 + 25.7, 23, 9.6, "none", rx=2, stroke=MISS_LINE, sw=.8, dash="2 1.2"),
       t(LX + 11.5, Y0 + 33, "–", 7, MUTE, anchor="middle"), t(LX + 28, Y0 + 33.2, "missed", 8),
       rect(LX, Y0 + 38.7, 23, 9.6, FALSE_FILL, rx=2, stroke=FALSE_LINE), t(LX + 11.5, Y0 + 46, "ACR", 7, INK, anchor="middle"),
       t(LX + 28, Y0 + 46.2, "false", 8)]

H = Y0 + 92 - RH
svg = ['<?xml version="1.0" encoding="UTF-8"?>',
       f'<svg id="Figure3" xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 518.4 {H}" width="518.4pt" height="{H}pt">',
       '<rect x="0" y="0" width="518.4" height="%d" fill="#fff"/>' % H,
       letter(4.302, 16.991, "A"), letter(151.109, 16.991, "B"), letter(374, 16.991, "C"),
       '<g id="A_metabolite_view">', *M_, '</g>', '<g id="B_complex_view">', *A_, '</g>',
       '<g id="C_nullclines">', *B_, '</g>', '<g id="D_steady_states">', *D_, '</g>',
       '<g id="E_acr_table">', *C_, '</g>', '</svg>']
open(sys.argv[1], "w").write("\n".join(svg) + "\n")
