"""Figure 3, panels A-D: a network whose answer is known, where the iterated and the tangent route disagree.

A protein X is made, degraded and assembled into an inactive kinase (2X -> C); the active kinase P activates
C cooperatively (C + 2P -> 3P), and C and P are degraded.
A: metabolite view, as Figure 1A. In C + 2P -> 3P, P is consumed twice and produced three times, so its net
   change is +1 and the net stoichiometry (C -> P) hides the autocatalysis. Halos mark true ACR (X).
B: complex view, as Figure 1B. The iterated route's one kinetic module, {P, C+2P}, is outlined; it fixes
   C·P, not C or P.
C: steady states at the rate constants of test_natural_candidates_for_the_completeness_figure. X has one
   (supply 2 = loss x + x^2 at x = 1); in the (C, P) plane at x = 1 the nullclines cross twice, at
   (1/10, 2/5) and (2/5, 1/10).
D: ACR per metabolite, truth against the iterated route, the tangent route (local ACR) and the sign certificate.
A and B are drawn in Figure 1's units (bold 14 pt nodes, 1.2 pt strokes, harpoons with the reverse direction
in #bbb) and shrunk by S to sit side by side; every reaction is its own group, so the drawing can be
rearranged in Illustrator. C and D are drawn at final size, 518.4 pt wide in all. Arial has no glyph for the
empty set or for check marks, so those are paths.

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

# ---- A: complex network, pills as in Figure 1B ----
PH, PAD = 20.26, 5.8                                   # pill height and side padding of Figure 1B
def pill_w(c): return max(22.46, width(c, 14, BLD) + 2 * PAD) if c != "0" else 22.46
POS = {"0": (118, 118), "X": (38, 72), "2X": (112, 30), "C": (196, 30), "P": (196, 130),
       "C+2P": (196, 182), "3P": (256, 182)}

def pill(c, fill="#fff"):
    x, y = POS[c]; w = pill_w(c)
    o = [rect(x - w / 2, y - PH / 2, w, PH, fill, rx=9.98, stroke="#000", sw=1)]
    if c == "0":                                        # the zero complex: an empty-set sign
        o += [f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.6" fill="none" stroke="#000" stroke-width="1.5"/>',
              line(x - 4.6, y + 6.2, x + 4.6, y - 6.2, "#000", 1.5)]
    else:
        o.append(t(x, y + 5.1, c, 14, "#000", bold=True, anchor="middle"))
    return o

def rim(c, u, gap=2.5):
    """point where the ray from the centre of pill c along unit vector u leaves the pill, plus a gap"""
    x, y = POS[c]; hw, hh = pill_w(c) / 2, PH / 2
    s = min(hw / abs(u[0]) if u[0] else math.inf, hh / abs(u[1]) if u[1] else math.inf)
    return x + (s + gap) * u[0], y + (s + gap) * u[1]

def arrow(src, dst, reversible):
    """irreversible: one line and a full head; reversible: the forward harpoon in black on the left of
    its travel direction and the reverse one in gray, as in Figure 1B"""
    (x0, y0), (x1, y1) = POS[src], POS[dst]
    L = math.hypot(x1 - x0, y1 - y0); u = ((x1 - x0) / L, (y1 - y0) / L)
    o = []
    for k, (a, b, col) in enumerate([(src, dst, "#000"), (dst, src, GRAY)][:1 + reversible]):
        v = u if k == 0 else (-u[0], -u[1]); n = (v[1], -v[0])        # left of travel, in SVG coordinates
        off = 1.19 if reversible else 0
        p = rim(a, v); q = rim(b, (-v[0], -v[1]))
        p = (p[0] + off * n[0], p[1] + off * n[1]); q = (q[0] + off * n[0], q[1] + off * n[1])
        if reversible:
            tip = (q[0] - .6 * n[0], q[1] - .6 * n[1]); base = (tip[0] - 5.29 * v[0], tip[1] - 5.29 * v[1])
            head = [tip, base, (base[0] + 3.63 * n[0], base[1] + 3.63 * n[1])]
        else:
            tip = q; base = (q[0] - 4.63 * v[0], q[1] - 4.63 * v[1])
            head = [tip, (base[0] + 2.64 * n[0], base[1] + 2.64 * n[1]), (base[0] - 2.64 * n[0], base[1] - 2.64 * n[1])]
        end = (base[0] + .3 * v[0], base[1] + .3 * v[1])
        o += [line(p[0], p[1], end[0], end[1], col), poly(head, col)]
    return o

S = 0.7                                             # A and B are drawn in Figure 1 units and shrunk by S
A_ = []
# the iterated route's kinetic module {P, C+2P}: psi_{C+2P} / psi_P = C·P is fixed, so C·P = k9 / k10
mx0 = POS["C+2P"][0] - pill_w("C+2P") / 2 - 6; mx1 = POS["C+2P"][0] + pill_w("C+2P") / 2 + 6
A_.append(rect(mx0, POS["P"][1] - PH / 2 - 6, mx1 - mx0, POS["C+2P"][1] - POS["P"][1] + PH + 12, KIN, rx=12,
               stroke="#000", sw=1.5, dash="6 2"))
for s_, p_, _ in RX:
    if (p_, s_) in {(a, b) for a, b, _ in RX} and (p_, s_) < (s_, p_):
        continue                                        # reverse of a pair: drawn with its forward reaction
    rid = f"{s_}_to_{p_}".replace("+", "_").replace("0", "zero")
    A_ += [f'<g id="B_{rid}">', *arrow(s_, p_, (p_, s_) in {(a, b) for a, b, _ in RX}), '</g>']
for c in POS:
    A_ += pill(c)
A_ += [rect(12, 213, 10.5, 8.6, KIN, rx=4.2, stroke="#000", sw=1.5, dash="6 2"),
       t(28, 221, "kinetic module (iterated route)", 10.5)]

# ---- metabolite view, as Figure 1A: species with halos, reactions as arrows between them ----
SP = {"X": (62, 112), "C": (150, 50), "P": (150, 172)}
NODE_R = 12                                            # halo radius 9.5 plus a gap, Figure 1 units
def zero(x, y, r=4.2):
    """where an exchange reaction meets the environment: nothing is drawn. The zero complex belongs to
    the complex view; in the metabolite view uptakes and exports simply start and end in open space."""
    return []
def head(tip, v, col="#000", harpoon_side=None):
    """full arrowhead at tip pointing along unit v, or a harpoon on side n (Figure 1 sizes)"""
    n = (v[1], -v[0])
    if harpoon_side is not None:
        base = (tip[0] - 5.29 * v[0], tip[1] - 5.29 * v[1])
        return poly([tip, base, (base[0] + 3.63 * n[0], base[1] + 3.63 * n[1])], col)
    base = (tip[0] - 4.63 * v[0], tip[1] - 4.63 * v[1])
    return poly([tip, (base[0] + 2.64 * n[0], base[1] + 2.64 * n[1]), (base[0] - 2.64 * n[0], base[1] - 2.64 * n[1])], col)
def unit(p, q):
    L = math.hypot(q[0] - p[0], q[1] - p[1]); return ((q[0] - p[0]) / L, (q[1] - p[1]) / L)
def shrink(p, q, rp, rq):
    """segment p->q cut back by rp at p and rq at q"""
    u = unit(p, q); return (p[0] + rp * u[0], p[1] + rp * u[1]), (q[0] - rq * u[0], q[1] - rq * u[1]), u
def pair(p, q, rp, rq):
    """reversible pair p <-> q: black forward harpoon, gray reverse harpoon, offset as in Figure 1"""
    o = []
    for k, (a, b, ra, rb, col) in enumerate([(p, q, rp, rq, "#000"), (q, p, rq, rp, GRAY)]):
        s0, s1, v = shrink(a, b, ra, rb); n = (v[1], -v[0])
        s0 = (s0[0] + 1.19 * n[0], s0[1] + 1.19 * n[1]); s1 = (s1[0] + 1.19 * n[0], s1[1] + 1.19 * n[1])
        tip = (s1[0] - .6 * n[0], s1[1] - .6 * n[1])
        o += [line(s0[0], s0[1], tip[0] - 5 * v[0], tip[1] - 5 * v[1], col), head(tip, v, col, harpoon_side=1)]
    return o
def single(p, q, rp, rq, col="#000"):
    s0, s1, v = shrink(p, q, rp, rq)
    return [line(s0[0], s0[1], s1[0] - 4.3 * v[0], s1[1] - 4.3 * v[1], col), head(s1, v, col)]
def bez(p, c, q, col="#000", arrow=True):
    """quadratic curve p -> q with control c; arrowhead at q along the end tangent"""
    o = [f'<path d="M{p[0]:.2f},{p[1]:.2f} Q{c[0]:.2f},{c[1]:.2f} {q[0]:.2f},{q[1]:.2f}" fill="none" stroke="{col}" stroke-width="1.2"/>']
    if arrow: o.append(head(q, unit(c, q), col))
    return o
def on_circle(centre, towards, r):
    u = unit(centre, towards); return (centre[0] + r * u[0], centre[1] + r * u[1])

M_ = []
TRUE_ACR = {"X"}                                       # halos: the true answer, as Figure 1A marks ACR
g = lambda name, items: [f'<g id="A_{name}">', *items, '</g>']
# supply and decay of X: a reversible pair with the environment (open space, no zero complex here)
M_ += g("X_supply_decay", pair((12, 112), SP["X"], 0, NODE_R))
# 2X -> C, labelled with the two molecules it consumes
M_ += g("2X_to_C", single(SP["X"], SP["C"], NODE_R, NODE_R) + [t(88, 76, "2", 10.5)])
M_ += g("C_decay", single(SP["C"], (150, 6), NODE_R, 0))
M_ += g("P_decay", single(SP["P"], (150, 218), NODE_R, 0))
# C + 2P -> 3P: C and P merge into the reaction, which returns to P; numbers give consumed and produced
J2 = (212, 112)
c_ = on_circle(SP["C"], J2, NODE_R); p_out = on_circle(SP["P"], (180, 140), NODE_R); p_in = on_circle(SP["P"], (190, 184), NODE_R)
M_ += g("C_2P_to_3P", bez(c_, (210, 62), J2, arrow=False) + bez(p_out, (196, 140), J2, arrow=False)
        + bez(J2, (240, 182), p_in) + [t(170, 140, "2", 10.5), t(214, 186, "3", 10.5)])
for m, (x, y) in SP.items():
    M_ += g(f"species_{m}", [f'<circle cx="{x}" cy="{y}" r="9.5" fill="{ACR_FILL if m in TRUE_ACR else "#fff"}"'
                             f'{"" if m in TRUE_ACR else " stroke=\"#c9ced4\" stroke-width=\"1\""}/>',
                             t(x, y + 5.1, m, 14, "#000", bold=True, anchor="middle")])
M_ += [f'<circle cx="18" cy="225" r="5.2" fill="{ACR_FILL}"/>', t(28, 229, "ACR (true)", 10.5)]

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

# ---- D: ACR per metabolite ----
Y0 = 226
C_ = [letter(4.302, Y0 - 5, "D"), t(24, Y0 - 7, "ACR, per metabolite")]
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
       letter(4.302, 16.991, "A"), letter(184, 16.991, "B"), letter(374, 16.991, "C"),
       f'<g id="A_metabolite_view" transform="translate(4 14) scale({S})">', *M_, '</g>',
       f'<g id="B_complex_view" transform="translate(186 14) scale({S})">', *A_, '</g>',
       '<g id="C_nullclines">', *B_, '</g>', '<g id="D_acr_table">', *C_, '</g>', '</svg>']
open(sys.argv[1], "w").write("\n".join(svg) + "\n")
