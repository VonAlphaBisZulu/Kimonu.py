"""Panels C-E of Figure 1, drawn into a copy of the Illustrator draft.

C: K = ker(N), with three of its columns drawn on thumbnails of the metabolic panel
D: complex activities M = A K
E: tangent space T = ker(N diag(r*) S^T Y^T) and its complex signatures Y^T T
Units are points on the 7.2 x 9 in artboard; text 8 pt Arial, panel letters 10 pt bold with their tops as in the draft.

    python fig1CDE_matrices.py <Figure1_draft.svg> <met_network_reversible.svg> <out.svg>
"""
import re, sys
import sympy as sy

fig_in, net_in, out = sys.argv[1:4]
CX = ['A', 'B', 'C', 'D', 'D+E', 'F', 'B+G', 'C+G', 'H', 'C+E', 'A+G', 'I', 'A+E']
RX = [('A', 'B'), ('B', 'A'), ('B', 'C'), ('C', 'B'), ('C', 'D'), ('D+E', 'F'), ('F', 'D+E'), ('F', 'B+G'),
      ('C+G', 'H'), ('H', 'C+G'), ('H', 'C+E'), ('A+G', 'I'), ('I', 'A+G'), ('I', 'A+E')]
SP = list("ABCDEFGHI")
Y = sy.zeros(9, 13); A = sy.zeros(13, 14)
for k, c in enumerate(CX):
    for t in c.split('+'): Y[SP.index(t), k] += 1
for j, (s, p) in enumerate(RX): A[CX.index(s), j] = -1; A[CX.index(p), j] = 1
N = Y * A
# K in a basis read off the drawing: two routes through the network, then the five reversible-pair cycles.
# Any basis of ker(N) would do; these columns are just easy to see on the network.
PAIRS = [(0, 1), (2, 3), (5, 6), (8, 9), (11, 12)]
P1, P2 = [2, 4, 5, 7, 8, 10], [2, 4, 5, 7, 11, 13]
cols = [[1 if j in e else 0 for j in range(14)] for e in (P1, P2)] + [[1 if j in p else 0 for j in range(14)] for p in PAIRS]
K = sy.Matrix(cols).T
assert N * K == sy.zeros(9, 7) and K.rank() == 7
M = A * K
r = K * sy.ones(7, 1)                                  # every reaction positive
G = sy.Matrix([list(Y[:, CX.index(s)].T) for s, _ in RX])
T = sy.Matrix.hstack(*(N * sy.diag(*r) * G).nullspace()).T.rref()[0].T
YT = Y.T * T
SUB = {s for s, _ in RX}

INK, MUTE, FONT = "#1b2430", "#9aa3ad", 'font-family="ArialMT, Arial"'
MOD = {**{c: "#fff2cd" for c in ['B', 'D', 'D+E', 'B+G']}, **{c: "#ffe0ff" for c in ['C+G', 'C+E']},
       **{c: "#e5d9ff" for c in ['A+G', 'A+E']}, **{c: "#5a9bd5" for c in ['A', 'C', 'F', 'H', 'I']}}
KIN = "#c9ecf4"
ACRR_FILL, ACR_FILL = "#d9e7c4", "#f2c14e"
# k1-k3 are drawn on the thumbnails: strong hue for the drawing and the column label, tint for the column band
KC = {"k1": ("#e6862b", "#fbe2cc"), "k2": ("#3a9d5d", "#d5eedd"), "k3": ("#56708c", "#dfe6ee")}  # k3 avoids lavender, the colour of concordance module 3

def t(x, y, s, size=8, fill=INK, bold=False, anchor=None):
    fam = 'font-family="Arial-BoldMT, Arial" font-weight="700"' if bold else FONT
    if anchor == "middle": x -= size * 0.28 * len(s)
    return f'<text transform="translate({x:.2f} {y:.2f})" {fam} font-size="{size}" fill="{fill}"><tspan x="0" y="0">{s}</tspan></text>'
def rect(x, y, w, h, fill, rx=2, stroke=None):
    st = f' stroke="{stroke}" stroke-width="0.6"' if stroke else ""
    return f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" fill="{fill}"{st}/>'
def letter(x, y, L):
    """Panel letter, 10 pt bold (Science). y is the baseline of the 18 pt letter of the draft; the top of the cap height
    stays there (Arial Bold cap height 0.716 em)."""
    return f'<text transform="translate({x} {y - (18 - 10) * 0.716:.3f})" fill="{INK}" font-family="Helvetica-Bold, Helvetica" font-size="10" font-weight="700"><tspan x="0" y="0">{L}</tspan></text>'

def matrix(x0, y0, rows, colnames, data, rh=10, cw=13, lw=26, chip=None, band=None, zero=MUTE, fmt=str):
    """row labels (optionally on coloured chips), column header, numbers; zeros muted.
    band: {column name: (hue, tint)}, a tinted column under a label in the hue"""
    o = []
    band = band or {}
    for j, c in enumerate(colnames):
        if c in band: o.append(rect(x0 + lw + j * cw, y0 - rh + 1, cw, rh * (len(rows) + 1), band[c][1], rx=0))
    for j, c in enumerate(colnames):
        o.append(t(x0 + lw + j * cw + cw / 2, y0, c, 7, band[c][0] if c in band else MUTE, bold=c in band, anchor="middle"))
    for i, r_ in enumerate(rows):
        y = y0 + (i + 1) * rh
        if chip and chip.get(r_):
            o.append(rect(x0, y - rh + 2.2, lw - 3, rh - 1.4, chip[r_]))
        o.append(t(x0 + 2, y, r_, 7.5))
        for j in range(len(colnames)):
            v = data[i][j]; s = fmt(v)
            o.append(t(x0 + lw + j * cw + cw / 2, y, s, 7.5, zero if v == 0 else INK, anchor="middle"))
    return "\n  ".join(o)

# ---- thumbnails of the metabolic panel, one column of K each ------------------------------------
net = open(net_in).read()
full = net[net.index('>', net.index('<svg')) + 1: net.rindex('</svg>')]
full = re.sub(r'<!--.*?-->', '', full)
body = re.sub(r'<g id="R\d+_reverse">.*?</g>', '', full, flags=re.S)                   # reverse directions off
body = re.sub(r'<g isolation="isolate">\s*<g isolation="isolate">\s*<text[^>]*font-size="8".*?</g>\s*</g>', '', body, flags=re.S)  # reaction labels off
# black elements of each reaction, identified by a coordinate they start with
OWN = {"M3.834,65.717": 14, "143.175 162.258": 14, "M41.111,133.368": 14, "3.834 62.166": 14,
       "M14.222,57.228": 12, "M79.748,120.782": 12, "x1=\"47.347\"": 12, "44.494 115.108": 12,
       "M95.877,73.196c-20.24": 8, "80.633 115.182": 8, "M79.818,36.572": 8, "77.475 37.268": 8, "x1=\"95.877\"": 8,
       "M68.936,16.709": 1, "71.348 16.847": 1,
       "M128.256,73.196c27.81": 6, "M128.256,73.196c39.052": 6, "x1=\"122.12\"": 6, "119.753 73.8": 6,
       "M147.687,86.152": 9, "M188.84,95.72": 9, "M188.36,95.18": 9, "x1=\"188.84\"": 9, "194.799 102.105": 9,
       "M211.383,125.74c22.162": 11, "264.154 93.465": 11, "M211.383,125.74c23.753": 11, "164.488 158.581": 11, "x1=\"205.959\"": 11,
       "M176.063,32.461": 5, "176.666 34.901": 5, "M87.838,19.071": 3, "262.085 71.797": 3}
def grow(poly, k):
    """scale an arrowhead polygon about its centroid, so heads stay legible on a thumbnail"""
    m = re.search(r'points="([^"]+)"', poly)
    p = [float(v) for v in m.group(1).split()]; xs, ys = p[0::2], p[1::2]
    cx_, cy_ = sum(xs) / len(xs), sum(ys) / len(ys)
    q = " ".join(f"{cx_ + k * (x - cx_):.3f} {cy_ + k * (y - cy_):.3f}" for x, y in zip(xs, ys))
    return poly[:m.start(1)] + q + poly[m.end(1):]

def thumb(em, colour):
    """the network with the reactions of one K column in colour and all others pale"""
    b = body
    for j in em:
        # a reverse direction drawn as its own grey harpoon (e.g. R2 in the R1/R2 cycle) comes back in colour
        m = re.search(r'<g id="R%d_reverse">.*?</g>' % (j + 1), full, flags=re.S)
        if m:
            g = m.group(0).replace('stroke="#bbb"', f'stroke="{colour}"').replace('fill="#bbb"', f'fill="{colour}"')
            g = re.sub(r'<polygon[^>]*/>', lambda p: grow(p.group(0), 1.9), g.replace('stroke-width="1.2"', 'stroke-width="2"'))
            b = g + b
    for key, rx in OWN.items():
        on = (rx - 1) in em
        pat = re.compile(r'<(path|polygon|line)[^>]*' + re.escape(key) + r'[^>]*/>')
        m = pat.search(b); assert m, key
        el = m.group(0)
        if on:
            el = el.replace('stroke="#000"', f'stroke="{colour}"').replace('stroke-width="1.2"', 'stroke-width="2"')
            if el.startswith('<polygon'):
                el = grow(el, 1.9)
                el = el.replace('/>', f' fill="{colour}"/>') if 'fill=' not in el else el
        else:
            el = el.replace('stroke="#000"', 'stroke="#d5d9de"')
            el = el.replace('/>', ' fill="#d5d9de"/>') if el.startswith('<polygon') and 'fill=' not in el else el
        b = b[:m.start()] + el + b[m.end():]
    return b

out_parts = []
Y0 = 222
KCOLS = [f"k{i + 1}" for i in range(7)]
# C: K, with its first three columns drawn on the network, stacked to its left ---------------------
out_parts.append(letter(4.302, Y0, "C"))
out_parts.append(t(24, Y0 - 2, "K = ker(N)", 8, INK))
S = 0.47
for i, (em, name) in enumerate(((P1, "k1"), (P2, "k2"), (PAIRS[0], "k3"))):
    x, y = 4, Y0 + 8 + i * 90
    out_parts.append(f'<g transform="translate({x} {y}) scale({S})">' + thumb(em, KC[name][0]) + '</g>')
    out_parts.append(t(x + S * 242, y + 8, name, 8, KC[name][0], bold=True))
chipK = {"R1": "#e5d9ff", "R2": "#e5d9ff", "R5": "#fff2cd", "R8": "#fff2cd"}
out_parts.append(matrix(136, Y0 + 22, [f"R{j + 1}" for j in range(14)], KCOLS,
                        [list(K[j, :]) for j in range(14)], chip=chipK, band=KC))
# D: complex activities --------------------------------------------------------------------------
DX = 266
out_parts.append(letter(DX, Y0, "D"))
out_parts.append(t(DX + 18, Y0 - 2, "M = A·K", 8, INK))
out_parts.append(matrix(DX + 2, Y0 + 22, CX, KCOLS, [list(M[i, :]) for i in range(13)], lw=30, chip=MOD, band=KC))
# E: tangent space and complex signatures --------------------------------------------------------
EX = 400
out_parts.append(letter(EX, Y0, "E"))
out_parts.append(t(EX + 18, Y0 - 2, "T = ker(N·diag(r*)·Sᵀ·Yᵀ)", 8, INK))
chipT = {s: ACRR_FILL for s in "ABCFHI"}; chipT["G"] = ACR_FILL
out_parts.append(t(EX + 20, Y0 + 12, "T", 7.5, MUTE))
out_parts.append(matrix(EX + 2, Y0 + 22, SP, ["t1", "t2"], [list(T[i, :]) for i in range(9)], lw=16, chip=chipT))
chipC = {c: KIN for c in CX if YT[CX.index(c), :] == sy.Matrix([[1, 0]]) and c in SUB}
out_parts.append(t(EX + 60, Y0 + 12, "Yᵀ·T", 7.5, MUTE))
out_parts.append(matrix(EX + 52, Y0 + 22, CX, ["t1", "t2"], [list(YT[i, :]) for i in range(13)], lw=26, chip=chipC))

fig = open(fig_in).read()
# artboard cropped to 7 in (504 pt), as in the draft Phil cropped; three thumbnails end just above it
fig = fig.replace('viewBox="0 0 518.4 648"', 'viewBox="0 0 518.4 504"', 1)
halos = []
nodes = re.findall(r'<text transform="translate\(([\d.]+) ([\d.]+)\)" font-family="Arial-BoldMT, Arial" font-size="14"[^>]*><tspan x="0" y="0">([A-I])</tspan>', fig)[:9]
groups = {}
for i in range(9):
    if all(v == 0 for v in T[i, :]): continue
    groups.setdefault(tuple(T[i, :]), []).append(SP[i])
acrr = [g for g in groups.values() if len(g) > 1]
for x, y, L in nodes:
    w = 14 * {"I": .278, "F": .611}.get(L, .72); cx_, cy_ = float(x) + w / 2, float(y) - 5
    fill = ACR_FILL if all(v == 0 for v in T[SP.index(L), :]) else ACRR_FILL if any(L in g for g in acrr) else None
    if fill: halos.append(f'<circle cx="{cx_:.2f}" cy="{cy_:.2f}" r="9.5" fill="{fill}"/>')
key = [f'<circle cx="16" cy="196" r="4.5" fill="{ACRR_FILL}"/>', t(24, 199, "ACRR group", 8),
       f'<circle cx="86" cy="196" r="4.5" fill="{ACR_FILL}"/>', t(94, 199, "ACR", 8)]
fig = fig.replace('<!-- Generator', '<g id="acr_acrr_halos">' + "".join(halos) + '</g>\n  <!-- Generator', 1)
out_parts = key + out_parts
fig = fig.replace('</svg>', '  <g id="panels_C_to_E_proposal">\n  ' + "\n  ".join(out_parts) + '\n  </g>\n</svg>')
open(out, "w").write(fig)
print("K", K.shape, "M", M.shape, "T", T.shape, "r*", list(r))
