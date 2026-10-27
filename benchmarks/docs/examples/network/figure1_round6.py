"""Figure 1, round 6: rows instead of columns in C and E, panel F as the readout of B = ker(V^T), and panel G
with both routes. Applied to Phil's hand-edited SVG (vault, 2026-10-14).

C, E: the k1-k3 column bands go, their coloured header labels stay. Each row chip is pulled across its row:
light chips are widened, saturated ones (balanced complexes, ACR) keep the chip and get a light band behind it.
F: the same row bands; T, t1, t2 become B, b1, b2; header B = ker(V^T) with the general containment and the
equality that holds in this network. G: the iterated route ends in ker((Y.Delta)^T), which contains B; a second
row adds the tangent route, whose nullspace ker(W.Y^T) lies in B. D: "nullspace basis" becomes a plain label."""
import re
import sys

src, out = sys.argv[1:3]
s = open(src, encoding="utf-8").read()
INK, MUTE = "#1b2430", "#9aa3ad"
REG = 'font-family="ArialMT, Arial"'
BOLD = 'font-family="Arial-BoldMT, Arial" font-weight="700"'


def once(old, new):
    global s
    n = s.count(old)
    assert n == 1, (n, old[:80])
    s = s.replace(old, new)


def runs(parts, size):
    """tspans for (text, bold, colour) runs that flow on one baseline"""
    out = []
    for t, b, c in parts:
        sp = ' xml:space="preserve"' if t != t.strip() else ""
        out.append(f'<tspan {BOLD if b else REG} fill="{c}" font-size="{size}"{sp}>{t}</tspan>')
    return "".join(out)


def text(x, y, parts, size, anchor="start"):
    a = f' text-anchor="{anchor}"' if anchor != "start" else ""
    return f'<text transform="translate({x} {y})"{a}>{runs(parts, size)}</text>'


# ---- C, E: drop the k1-k3 column bands (two sets of three, 150 and 140 high)
for x, h in [("46.56", 150), ("59.56", 150), ("72.56", 150), ("181.653", 140), ("194.653", 140), ("207.653", 140)]:
    s, n = re.subn(r'<rect x="%s" y="247\.931" width="13" height="%d"[^>]*/>\s*' % (re.escape(x), h), "", s)
    assert n == 1, (x, n)

# ---- C, E, F: pull each row chip across its row
TINT = {"#5a9bd5": "#d6e6f5", "#f2c14e": "#fcefcf"}     # saturated chips: light band, chip stays
ROW_END = [(22.786, 137.56), (159.102, 272.653), (298.748, 337.75), (364.346, 409.97)]  # chip x -> table end


def pull(m):
    el = m.group(0)
    x = float(m.group(1))
    end = next((e for cx, e in ROW_END if abs(cx - x) < 1e-3), None)
    if end is None:
        return el
    fill = re.search(r'fill="([^"]+)"', el).group(1)
    w = f'{end - x:.3f}'
    if fill in TINT:
        band = re.sub(r'width="[^"]+"', f'width="{w}"', el).replace(f'fill="{fill}"', f'fill="{TINT[fill]}"')
        return band + el
    return re.sub(r'width="[^"]+"', f'width="{w}"', el)


s, n = re.subn(r'<rect x="([\d.]+)" y="([\d.]+)" width="[\d.]+" height="8\.6"[^>]*/>', pull, s)
assert n == 36, n          # 4 in C, 13 in E, 7 + 12 in F

# ---- F: header, title, sub-labels
s, n = re.subn(r'<text transform="translate\(([\d.]+) 255\.965\)"[^>]*>.*?</text>',
               lambda m: m.group(0) if not 350 <= float(m.group(1)) <= 460 else "", s, flags=re.S)
assert n == 10, n
once('<text transform="translate(336.406 239.354)"', '<text transform="translate(336.406 233)"')
header = (text(406, 247.5, [("B", True, INK), (" = ker(", False, INK), ("V", True, INK), ("ᵀ)", False, INK)], 10, "middle")
          + text(406, 258, [("ker(", False, INK), ("W", True, INK), ("·", False, INK), ("Y", True, INK),
                            ("ᵀ) ⊆ ", False, INK), ("B", True, INK), (" ⊆ ker((", False, INK), ("Y", True, INK),
                            ("·Δ)ᵀ)", False, INK), (",  here all equal", False, MUTE)], 7, "middle"))
once("</svg>", header + "\n</svg>")
for x, y in [("312.165", "272.383"), ("383.106", "270.376")]:
    once(f'translate({x} {y})" fill="{INK}" font-family="Arial-BoldMT, Arial" font-size="8" font-weight="700" '
         f'isolation="isolate"><tspan x="0" y="0">T</tspan>',
         f'translate({x} {y})" fill="{INK}" font-family="Arial-BoldMT, Arial" font-size="8" font-weight="700" '
         f'isolation="isolate"><tspan x="0" y="0">B</tspan>')
for k in ("1", "2"):
    s, n = re.subn(r'(<text transform="translate\([\d.]+ 28[02]\.\d+\)" fill="#9aa3ad"[^>]*><tspan x="0" y="0">)t%s<' % k,
                   r'\g<1>b%s<' % k, s)
    assert n == 2, (k, n)

# ---- D: label
once('<tspan x="0" y="0">nullspace basis</tspan></text>',
     '<tspan x="0" y="-5">examples of</tspan><tspan x="0" y="5.5">steady-state fluxes</tspan></text>')

# ---- G: title, iterated output, tangent row
once('<tspan x="0" y="0">iterated route (Fig. 2)</tspan>', '<tspan x="0" y="0">two routes to B (Fig. 2)</tspan>')
g = s.index('<g id="panel_G_iterated_route">')
for t in ("modules,", "ACR, ACRR", "(proven)"):
    s, n = re.subn(r'<g isolation="isolate">\s*<text transform="translate\(436 4\d\d\)"[^>]*><tspan x="0" y="0">%s</tspan></text>\s*</g>'
                   % re.escape(t), "", s)
    assert n == 1, t
once('translate(291.04 499)', 'translate(291.04 496)')
# old footer: "one-step local readout (F): T, Y^T.T ..." and "proven ⊆ true ⊆ local ..."
a = s.index('<text transform="translate(212 513)"', g)
a = s.rindex("<g>", g, a)
b = s.index("local; equal on 30 of 33 models</tspan></text>", a)
b = s.index("</g>", b) + len("</g>")
s = s[:a] + s[b:]
arrow = lambda y: (f'<path d="M416,{y}h16" fill="none" stroke="{INK}"/>'
                   f'<polygon points="432 {y} 427.5 {y - 2.2} 427.5 {y + 2.2} 432 {y}" fill="{INK}"/>')
G_NEW = (
    text(436, 466, [("ker((", False, INK), ("Y", True, INK), ("·Δ)ᵀ) ⊇ ", False, INK), ("B", True, INK)], 8)
    + text(436, 475, [("iterated: sound", False, MUTE)], 7)
    + f'<rect x="212" y="503" width="204" height="15" rx="4" ry="4" fill="#fff" stroke="{INK}"/>'
    + text(314, 513.2, [("tangent space at ", False, INK), ("r", True, INK), ("*", False, INK),
                        ("     ", False, INK), ("W", True, INK), ("·", False, INK), ("Y", True, INK),
                        ("ᵀ·", False, INK), ("h", True, INK), (" = 0", False, INK)], 8, "middle")
    + arrow(510.5)
    + text(436, 509, [("ker(", False, INK), ("W", True, INK), ("·", False, INK), ("Y", True, INK),
                      ("ᵀ) ⊆ ", False, INK), ("B", True, INK)], 8)
    + text(436, 518, [("tangent: complete, local", False, MUTE)], 7)
)
end_g = s.index('<text transform="translate(44.077 231.61)"')
end_g = s.rindex("</g>", g, end_g)
s = s[:end_g] + G_NEW + "\n  " + s[end_g:]

open(out, "w", encoding="utf-8").write(s)
print("ok")
