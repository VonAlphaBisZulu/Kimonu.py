"""Figure 1: the defining equations above panels A and B."""
import re, sys
src, out = sys.argv[1:3]
s = open(src).read()
SHIFT = 24                                           # room for two equation lines under the letters' baselines
s = s.replace('viewBox="0 0 518.4 504"', f'viewBox="0 0 518.4 {504 + SHIFT}"', 1)
open_tag = re.search(r'<svg[^>]*>', s); close = s.rindex('</svg>')
body = s[open_tag.end():close]
# panel letters A and B stay where they are; everything else moves down
for x in ("4.302", "289.912"):
    body = body.replace(f'translate({x} 16.991)', f'translate({x} {16.991 - SHIFT})', 1)
INK = "#1b2430"
def line(x, y, parts, size=9):
    """parts: (text, bold) runs, flowing inline"""
    BOLD, REG = 'font-family="Arial-BoldMT, Arial" font-weight="700"', 'font-family="ArialMT, Arial"'
    sp = "".join(f'<tspan {BOLD if b else REG}>{t}</tspan>' for t, b in parts)
    return f'<text transform="translate({x} {y})" fill="{INK}" font-size="{size}" xml:space="preserve">{sp}</text>'
R, B = False, True
eqs = [
    # A: stoichiometric model, mass balance with irreversible (split) reactions
    line(26, 12.5, [("d", R), ("x", B), ("/d", R), ("t", R), (" = ", R), ("N", B), ("·", R), ("r", B), (" = ", R), ("0", B)]),
    line(26, 24, [("r", B), (" ≥ ", R), ("0", B)]),
    # B: the same network under mass action, rates from complex concentrations
    line(311, 12.5, [("d", R), ("x", B), ("/d", R), ("t", R), (" = ", R), ("Y", B), ("·", R), ("A", B), ("·", R), ("r", B), (" = ", R), ("0", B)]),
    line(311, 24, [("r", B), (" = diag(", R), ("k", B), (")·", R), ("S", B), ("ᵀ·", R), ("ψ", B), (",   ln ", R), ("ψ", B),
                   (" = ", R), ("Y", B), ("ᵀ·ln ", R), ("x", B)]),
]
s = s[:open_tag.end()] + f'\n  <g id="equations">{"".join(eqs)}</g>\n  <g transform="translate(0 {SHIFT})">' + body + '</g>\n' + s[close:]
open(out, "w").write(s)
