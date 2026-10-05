"""Figure 3: a network whose answer is known, where the iterated and the tangent route disagree.

A: complex network. Dimerising proteins A, B (supply, decay, homo- and heterodimer removal) next to an
   autocatalytic kinase switch C, P (C + 2P -> 3P). The iterated route's one kinetic module, {P, C+2P},
   is outlined; it fixes C·P, not C or P.
B: nullclines at the rate constants of test_example_network_for_the_completeness_figure. A and B cross
   once, C and P twice: two positive steady states that differ in C and P.
C: ACR per metabolite, truth against the iterated route, the tangent route (local ACR) and the sign certificate.
Style as Figure 1: 518.4 pt wide, Arial 8 pt text, bold 14 pt complexes, 1.2 pt strokes, harpoons with
the reverse direction in #bbb. Arial has no glyph for the empty set or for check marks, so those are paths.

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
BLD = {**REG, "A": .722, "B": .722, "P": .667, "C": .722}                             # Arial Bold, for the complexes
width = lambda s, fs, w=REG: fs * sum(w[c] for c in s)

# ---- the model: reactions with the rate constants of the test, and its positive steady states ----
RX = [("0", "A", 2), ("A", "0", 1), ("2A", "0", .5), ("0", "B", 1), ("B", "0", 1), ("A+B", "0", 1),
      ("0", "C", 4), ("C", "0", 1), ("P", "0", 1), ("C+2P", "3P", 1)]
# mass action: dA = 2 - A - 2*.5*A^2 - A*B, dB = 1 - B - A*B, dC = 4 - C - C*P^2, dP = -P + C*P^2
a_ss = brentq(lambda a: 2 - a - a * a - a / (1 + a), .01, 2); b_ss = 1 / (1 + a_ss)
cp_ss = [(2 - 3 ** .5, 2 + 3 ** .5), (2 + 3 ** .5, 2 - 3 ** .5)]      # C*P = 1 and C + P = 4
assert abs(2 - a_ss - a_ss ** 2 - a_ss * b_ss) < 1e-12 and abs(1 - b_ss - a_ss * b_ss) < 1e-12
for c, p in cp_ss:
    assert abs(4 - c - c * p * p) < 1e-12 and abs(-p + c * p * p) < 1e-12

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
POS = {"0": (112, 106), "A": (38, 44), "2A": (26, 92), "A+B": (32, 136), "B": (48, 180),
       "C": (184, 58), "P": (176, 136), "C+2P": (176, 180), "3P": (240, 180)}

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

A_ = [letter(4.302, 16.991, "A")]
# the iterated route's kinetic module {P, C+2P}: psi_{C+2P} / psi_P = C·P is fixed, so C·P = k9 / k10
mx0 = POS["C+2P"][0] - pill_w("C+2P") / 2 - 6; mx1 = POS["C+2P"][0] + pill_w("C+2P") / 2 + 6
A_.append(rect(mx0, POS["P"][1] - PH / 2 - 6, mx1 - mx0, POS["C+2P"][1] - POS["P"][1] + PH + 12, KIN, rx=12,
               stroke="#000", sw=1.5, dash="6 2"))
arrows = []
for s, p, _ in RX:
    if (p, s) in {(a, b) for a, b, _ in RX} and (p, s) < (s, p):
        continue                                        # reverse of a pair: drawn with its forward reaction
    arrows += arrow(s, p, (p, s) in {(a, b) for a, b, _ in RX})
A_ += arrows
for c in POS:
    A_ += pill(c)
A_ += [rect(12, 199, 7.85, 6.43, KIN, rx=3.17, stroke="#000", sw=1.5, dash="6 2"),
       t(24, 205, "kinetic module (iterated route)")]

# ---- B: nullclines ----
def plot(x0, y0, w, h, xmax, ymax, ticks, curves, points, xlab, ylab, title):
    """data box of w x h pt with its lower left corner at (x0, y0); curves as (points, colour, label, (x, y))"""
    X = lambda v: x0 + v / xmax * w
    Yp = lambda v: y0 - v / ymax * h
    o = [t(x0 + w / 2, y0 - h - 8, title, anchor="middle"),
         line(x0, y0, x0 + w, y0, INK, .6), line(x0, y0, x0, y0 - h, INK, .6)]
    for v in ticks:
        o += [line(X(v), y0, X(v), y0 + 2.5, INK, .6), t(X(v), y0 + 10, f"{v:g}", 7, MUTE, anchor="middle"),
              line(x0, Yp(v), x0 - 2.5, Yp(v), INK, .6), t(x0 - 4.5, Yp(v) + 2.5, f"{v:g}", 7, MUTE, anchor="end")]
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
B_ = [letter(262, 16.991, "B")]
B_ += plot(290, 160, 92, 110, 1.5, 1.5, [0, 1],
           [([(a, (2 - a - a * a) / a) for a in grid(.62, 1, 200)], BLUE, "dA/dt = 0", (1.0, 1.3)),
            ([(a, 1 / (1 + a)) for a in grid(0, 1.5, 200)], ORANGE, "dB/dt = 0", (.05, 1.1))],
           [((a_ss, b_ss), "1, 2", (4, -4), None)], "A", "B", "A, B: one steady state")
B_ += plot(420, 160, 92, 110, 4.5, 4.5, [0, 2, 4],
           [([(4 / (1 + p * p), p) for p in grid(0, 4.5, 300)], BLUE, "dC/dt = 0", (2.0, 2.75)),
            ([(c, 1 / c) for c in grid(.2, 4.5, 300)], ORANGE, "dP/dt = 0", (2.0, 2.15))],
           [(cp_ss[0], "1", (4, 2.5), None), (cp_ss[1], "2", (0, -5), "middle")], "C", "P", "C, P: two steady states")
# the orange curve is the iterated route's kinetic module: C·P = k9/k10 on every steady state
B_.append(t(420 + 2.0 / 4.5 * 92, 160 - 1.7 / 4.5 * 110, "C·P = 1", 7.5, ORANGE))

# ---- C: ACR per metabolite ----
Y0 = 236
C_ = [letter(4.302, Y0 - 5, "C"), t(24, Y0 - 7, "ACR, per metabolite")]
COLS = [("steady state 1", 84), ("steady state 2", 140), ("truth", 196), ("iterated", 250), ("tangent", 304),
        ("certificate", 358)]
for name, x in COLS:
    C_.append(t(x, Y0 + 6, name, 7, MUTE, anchor="middle"))
C_.append(t(304, Y0 - 2, "local ACR", 6.5, MUTE, anchor="middle"))   # the tangent route decides local ACR
C_.append(t(250, Y0 - 2, "proven", 6.5, MUTE, anchor="middle"))
SS = {"A": (a_ss, a_ss), "B": (b_ss, b_ss), "C": (cp_ss[0][0], cp_ss[1][0]), "P": (cp_ss[0][1], cp_ss[1][1])}
TRUTH = {"A": True, "B": True, "C": False, "P": False}
CALLS = {"iterated": set(), "tangent": {"A", "B", "C", "P"}, "certificate": {"A", "B"}}
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
for i, m in enumerate("ABCP"):
    y = Y0 + 7 + (i + 1) * RH
    C_ += [t(40, y, m, 8, INK, bold=True, anchor="middle"),
           t(84, y, f"{SS[m][0]:.3f}", 7.5, anchor="middle"), t(140, y, f"{SS[m][1]:.3f}", 7.5, anchor="middle")]
    C_ += cell(196, y, TRUTH[m], TRUTH[m])
    for name, x in COLS[3:]:
        C_ += cell(x, y, m in CALLS[name], TRUTH[m])
# the two properties the figure is about: no false ACR claims (correct), no missed ACR (complete)
ys = Y0 + 7 + 5 * RH + 4
C_.append(line(30, ys - 9.5, 374, ys - 9.5, "#dde1e6", .6))
for j, prop in enumerate(("correct", "complete")):
    y = ys + j * 12
    C_.append(t(222, y, prop, 7.5, INK, anchor="end"))
    for name, x in COLS[3:]:
        wrong = [m for m in "ABCP" if (m in CALLS[name]) != TRUTH[m] and (m in CALLS[name]) == (prop == "correct")]
        C_.append(cross(x, y - 2.6) if wrong else check(x, y - 2.6))
# legend
LX = 408
C_ += [rect(LX, Y0 + 12.7, 23, 9.6, ACR_FILL, rx=2), t(LX + 11.5, Y0 + 20, "ACR", 7, INK, anchor="middle"),
       t(LX + 28, Y0 + 20.2, "reported", 8),
       rect(LX, Y0 + 25.7, 23, 9.6, "none", rx=2, stroke=MISS_LINE, sw=.8, dash="2 1.2"),
       t(LX + 11.5, Y0 + 33, "–", 7, MUTE, anchor="middle"), t(LX + 28, Y0 + 33.2, "missed", 8),
       rect(LX, Y0 + 38.7, 23, 9.6, FALSE_FILL, rx=2, stroke=FALSE_LINE), t(LX + 11.5, Y0 + 46, "ACR", 7, INK, anchor="middle"),
       t(LX + 28, Y0 + 46.2, "false", 8)]

H = 330
svg = ['<?xml version="1.0" encoding="UTF-8"?>',
       f'<svg id="Figure3" xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="0 0 518.4 {H}" width="518.4pt" height="{H}pt">',
       '<rect x="0" y="0" width="518.4" height="%d" fill="#fff"/>' % H,
       '<g id="A_complex_network">', *A_, '</g>', '<g id="B_nullclines">', *B_, '</g>',
       '<g id="C_acr_table">', *C_, '</g>', '</svg>']
open(sys.argv[1], "w").write("\n".join(svg) + "\n")
