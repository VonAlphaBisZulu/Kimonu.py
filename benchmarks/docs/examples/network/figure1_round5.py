"""Figure 1, round 5: dx/dt = N·r = Y·A·r = 0 above B; coupling brackets in F for the true groups; panel G
redrawn with W spelled out. Applied to Phil's hand-edited SVG."""
import re, sys
src, out = sys.argv[1:3]
s = open(src).read()
INK, MUTE, GRAY, KIN, MINT = "#1b2430", "#9aa3ad", "#bbb", "#c9ecf4", "#cdf3df"
REG, BOLD = 'font-family="ArialMT, Arial"', 'font-family="Arial-BoldMT, Arial" font-weight="700"'

# ---- B: insert "N·r = " after "dx/dt = " (same advance widths as panel A's line), shift the rest
W_INS = 79.53 - 56.274
def shift(m):
    x = float(m.group(1)); return f'translate({x + W_INS:.3f} 12.5)' if x >= 341 else m.group(0)
i = s.index('<text transform="translate(341.274 12.5)"'); j = s.index('translate(373.53 12.5)') + 200
seg = re.sub(r'translate\(([\d.]+) 12\.5\)', shift, s[i:j])
ins = "".join(f'<text transform="translate({x:.3f} 12.5)" fill="{INK}" {f} font-size="9"><tspan x="0" y="0"{p}>{t}</tspan></text>\n      '
              for x, f, t, p in [(341.274, BOLD, "N", ""), (347.773, REG, "·", ""), (350.771, BOLD, "r", ""),
                                 (354.273, REG, " = ", ' xml:space="preserve"')])
s = s[:i] + ins + seg + s[j:]

# ---- F: replace the placeholder coupling brackets (grey lines and link icons left of T and Y^T.T)
def in_f(el):
    m = re.search(r'\bd="M\s*(-?[\d.]+)[ ,](-?[\d.]+)', el) or re.search(r'x1="(-?[\d.]+)"\s+y1="(-?[\d.]+)"', el)
    if not m: return False
    x, y = float(m.group(1)), float(m.group(2))
    return 'stroke="#bbb"' in el and 270 <= y <= 420 and (285 <= x <= 300 or 352 <= x <= 362)
link = re.search(r'<path d="M289\.698,302\.901[^>]*/>\s*<path d="M290\.939,305\.344[^>]*/>', s).group(0)
s, n = re.subn(r'<(path|line)\b[^>]*/>', lambda m: "" if in_f(m.group(0)) else m.group(0), s)
LINK_C = (290.318, 304.12)                       # centre of the link icon on its bracket line
def bracket(x0, ys, icon_y, colour=GRAY, tick=5.608, gap=8.5):
    """Phil's style: vertical line through the member rows, a tick at each member, link icon in a gap."""
    L = lambda x1, y1, x2, y2: (f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" fill="none" '
                                f'stroke="{colour}" stroke-linecap="round" stroke-miterlimit="10" stroke-width=".75"/>')
    top, bot = min(ys), max(ys)
    o = [L(x0, top, x0, icon_y - gap), L(x0, icon_y + gap, x0, bot)]
    o += [L(x0 + tick, y, x0, y) for y in ys]
    icon = link.replace('stroke="#bbb"', f'stroke="{colour}"')
    o.append(f'<g transform="translate({x0 - LINK_C[0]:.3f} {icon_y - LINK_C[1]:.3f})">{icon}</g>')
    return "".join(o)
T_off, YT_off = 3.634, 3.761
T_rows = dict(zip("ABCDEFGHI", [292.383 + 10 * k for k in range(9)]))
YT_names = ["A", "B", "C", "D", "D+E", "F", "B+G", "C+G", "H", "C+E", "A+G", "I", "A+E"]
YT_rows = {c: 290.376 + 10 * k for k, c in enumerate(YT_names)}
acrr = [T_rows[c] - T_off for c in "ABCFHI"]                      # one ACRR group; G is ACR on its own
mod1 = [YT_rows[c] - YT_off for c in ["A", "B", "C", "D+E", "F", "B+G", "C+G", "H", "A+G", "I"]]
mod2 = [YT_rows[c] - YT_off for c in ["C+E", "A+E"]]
new_f = (bracket(290.318, acrr, T_rows["D"] - T_off + 5) +
         bracket(354.839, mod1, YT_rows["D"] - YT_off) +
         bracket(349.2, mod2, YT_rows["A+G"] - YT_off, tick=11.247))
s = s.replace("</svg>", f'  <g id="coupling_brackets_F">{new_f}</g>\n</svg>')

# ---- G: redraw with W written out
s = re.sub(r'<g id="panel_G_iterated_route">.*?</g>\n', "", s, flags=re.S)
def t(x, y, txt, size=8, fill=INK, bold=False, anchor=None):
    a = f' text-anchor="{anchor}"' if anchor else ""
    return f'<text x="{x}" y="{y}" {BOLD if bold else REG} font-size="{size}" fill="{fill}"{a} xml:space="preserve">{txt}</text>'
b = lambda txt: f'<tspan {BOLD}>{txt}</tspan>'
def box(x, y, w, h, fill): return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" ry="4" fill="{fill}" stroke="{INK}" stroke-width="1"/>'
def arrow(x1, y1, x2, y2, curve=None):
    sgn = 1 if x2 > x1 else -1
    d = f'M{x1},{y1} Q{curve[0]},{curve[1]} {x2},{y2}' if curve else f'M{x1},{y1} L{x2},{y2}'
    return (f'<path d="{d}" fill="none" stroke="{INK}" stroke-width="1"/>'
            f'<polygon points="{x2},{y2} {x2 - 4.5 * sgn},{y2 - 2.2} {x2 - 4.5 * sgn},{y2 + 2.2}" fill="{INK}"/>')
g = [f'<text transform="translate(208 438)" fill="{INK}" font-family="Helvetica-Bold, Helvetica" font-size="18" font-weight="700">G</text>',
     t(228, 437, "iterated route (Fig. 2)", 8, INK, True),
     box(212, 448, 86, 40, KIN),
     t(255, 460, f'ker({b("W")})', 8, anchor="middle"),
     t(255, 471, f'{b("W")} = {b("N")}·diag({b("r")}*)·{b("S")}ᵀ', 7, anchor="middle"),
     t(255, 481, "r*: a steady-state flux", 6.5, MUTE, anchor="middle"),
     box(340, 448, 76, 40, "#fff"),
     t(378, 460, "log closure", 8, anchor="middle"),
     t(378, 471, f'ln {b("ψ")} = {b("Y")}ᵀ·ln {b("x")}', 7, anchor="middle"),
     t(378, 481, "merge if difference fixed", 6.5, MUTE, anchor="middle"),
     arrow(298, 459, 339, 459), t(319, 455.5, "coupled sets", 6.5, MUTE, anchor="middle"),
     arrow(340, 479, 299, 479, curve=(319.5, 492)), t(319.5, 499, "feedback u_c = u_d", 6.5, MUTE, anchor="middle"),
     arrow(416, 468, 432, 468), t(436, 464, "modules,", 8), t(436, 474, "ACR, ACRR", 8), t(436, 483, "(proven)", 7, MUTE),
     t(212, 513, f'one-step local readout (F): {b("T")}, {b("Y")}ᵀ·{b("T")} (Pascual-Escudero &amp; Feliu)', 7, MUTE),
     t(212, 522, "proven ⊆ true ⊆ local; equal on 30 of 33 models", 7, MUTE)]
s = s.replace("</svg>", '  <g id="panel_G_iterated_route">' + "".join(g) + "</g>\n</svg>")
open(out, "w").write(s)
print("removed placeholder elements:", n)
