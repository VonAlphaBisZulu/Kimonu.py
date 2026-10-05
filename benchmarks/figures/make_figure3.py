"""Figure 3: kinetic modules and ACR, bounded from both sides.

A, B  the 33 benchmark models: iterated route (proven) against tangent route (local ACR), make_hierarchy.py
C-E   the toy network where local and global differ: network, nullclines, ACR per metabolite (figure3.py)
F, G  the tangent route at scale: symbolic local-ACR criterion against one exact nullspace, and float64
      against exact arithmetic (make_tangent_benchmark.py)

The three parts are drawn by their own scripts at the size they get here (text kept as text) and nested
into one SVG with transforms; the panel letters are added here, in Figure 1's style.

    python figures/make_figure3.py      # from benchmarks/
"""
import os, re, subprocess, sys, tempfile

W, GAP = 518.4, 4                         # Figure 1's artboard width; vertical gap between rows (pt)
TOY_SCALE = 0.9                           # the toy row is drawn at 518.4 x 330 pt and shrunk to fit
LETTER = '<text transform="translate({x:.2f} {y:.2f})" fill="#1b2430" font-family="Helvetica-Bold, Helvetica" font-size="18" font-weight="700">{L}</text>'

tmp = tempfile.mkdtemp()
env = dict(os.environ, FIG_LETTERS="0", FIG_OUT=tmp)
subprocess.run([sys.executable, "figures/make_hierarchy.py", tmp], check=True,
               env=dict(env, FIG_HEIGHT_IN="2.3", FIG_NAME="hierarchy", FIG_LEGEND_BOTTOM="1"))
subprocess.run([sys.executable, "figures/make_tangent_benchmark.py"], check=True,
               env=dict(env, FIG_HEIGHT_CM="5.3", FIG_NAME="tangent"))
subprocess.run([sys.executable, "docs/examples/network/figure3.py", f"{tmp}/toy.svg"], check=True)

def part(path):
    """(width, height, body) of an SVG file, root element stripped"""
    s = open(path).read()
    root = re.search(r"<svg\b[^>]*>", s)
    w, h = map(float, re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', root.group(0)).groups())
    return w, h, s[root.end():s.rindex("</svg>")]

rows, y, letters = [], 0.0, []
# A, B: hierarchy
w, h, body = part(f"{tmp}/hierarchy.svg"); s = W / w
rows.append(f'<g id="AB_hierarchy" transform="translate(0 {y:.2f}) scale({s:.4f})">{body}</g>')
letters += [LETTER.format(x=2, y=y + 16, L="A"), LETTER.format(x=W / 2 - 12, y=y + 16, L="B")]
y += h * s + GAP
# C, D, E: toy network; its own 18 pt letters are replaced by ours
w, h, body = part(f"{tmp}/toy.svg")
body = re.sub(r'<text[^>]*font-size="18"[^>]*>.*?</text>', "", body)
body = re.sub(r'<rect x="0" y="0" width="518.4" height="\d+" fill="#fff"/>', "", body)
rows.append(f'<g id="CDE_toy" transform="translate(0 {y:.2f}) scale({TOY_SCALE})">{body}</g>')
letters += [LETTER.format(x=2, y=y + 17, L="C"), LETTER.format(x=262 * TOY_SCALE, y=y + 17, L="D"),
            LETTER.format(x=2, y=y + 231 * TOY_SCALE + 2, L="E")]
y += h * TOY_SCALE + GAP
# F, G: tangent-route benchmark
w, h, body = part(f"{tmp}/tangent.svg"); s = W / w
rows.append(f'<g id="FG_tangent" transform="translate(0 {y:.2f}) scale({s:.4f})">{body}</g>')
letters += [LETTER.format(x=2, y=y + 14, L="F"), LETTER.format(x=0.51 * W, y=y + 14, L="G")]
y += h * s

H = round(y + 2, 1)
svg = ['<?xml version="1.0" encoding="UTF-8"?>',
       f'<svg id="Figure3" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" '
       f'viewBox="0 0 {W} {H}" width="{W}pt" height="{H}pt">',
       f'<rect x="0" y="0" width="{W}" height="{H}" fill="#fff"/>', *rows, '<g id="letters">', *letters, '</g>', '</svg>']
os.makedirs("figures/main", exist_ok=True)
open("figures/main/Figure3_draft.svg", "w").write("\n".join(svg) + "\n")
subprocess.run(["rsvg-convert", "-d", "300", "-p", "300", "-o", "figures/main/Figure3_draft.png", "figures/main/Figure3_draft.svg"], check=True)
print(f"Figure 3: {W} x {H} pt")
