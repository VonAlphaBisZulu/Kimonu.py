"""Figure 2: kinetic modules and ACR, bounded from both sides.

A-E   the toy network where local and global differ: metabolite view, complex view, nullclines, the two
      steady states on the network, ACR per metabolite (figure2.py, which draws its own letters A-E)
F, G  the 33 benchmark models: iterated route (proven) against tangent route (local ACR), make_hierarchy.py
H, I  the tangent route at scale: symbolic local-ACR criterion against one exact nullspace, and how many of
      the exact ACR/ACRR relations float64 still finds (make_tangent_benchmark.py)

The three parts are drawn by their own scripts at the size they get here (text kept as text) and nested
into one SVG with transforms; letters F-I are added here, in Figure 1's style.

    python figures/make_figure2.py      # from benchmarks/
"""
import os, re, subprocess, sys, tempfile

W, GAP = 518.4, 4                         # Figure 1's artboard width; vertical gap between rows (pt)
LETTER = '<text transform="translate({x:.2f} {y:.2f})" fill="#1b2430" font-family="Helvetica-Bold, Helvetica" font-size="18" font-weight="700">{L}</text>'

tmp = tempfile.mkdtemp()
env = dict(os.environ, FIG_LETTERS="0", FIG_OUT=tmp)
subprocess.run([sys.executable, "figures/make_hierarchy.py", tmp], check=True,
               env=dict(env, FIG_HEIGHT_IN="2.3", FIG_NAME="hierarchy", FIG_LEGEND_BOTTOM="1"))
subprocess.run([sys.executable, "figures/make_tangent_benchmark.py"], check=True,
               env=dict(env, FIG_HEIGHT_CM="5.3", FIG_NAME="tangent"))
subprocess.run([sys.executable, "docs/examples/network/figure2.py", f"{tmp}/toy.svg"], check=True)

def part(path):
    """(width, height, body) of an SVG file, root element stripped"""
    s = open(path).read()
    root = re.search(r"<svg\b[^>]*>", s)
    w, h = map(float, re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', root.group(0)).groups())
    return w, h, s[root.end():s.rindex("</svg>")]

rows, y, letters = [], 0.0, []
# A-E: toy network, drawn at final width with its own letters
w, h, body = part(f"{tmp}/toy.svg")
body = re.sub(r'<rect x="0" y="0" width="518.4" height="\d+" fill="#fff"/>', "", body)
rows.append(f'<g id="AE_toy" transform="translate(0 {y:.2f})">{body}</g>')
y += h + GAP
# F, G: hierarchy
w, h, body = part(f"{tmp}/hierarchy.svg"); s = W / w
rows.append(f'<g id="FG_hierarchy" transform="translate(0 {y:.2f}) scale({s:.4f})">{body}</g>')
letters += [LETTER.format(x=2, y=y + 16, L="F"), LETTER.format(x=W / 2 - 12, y=y + 16, L="G")]
y += h * s + GAP
# H, I: tangent-route benchmark
w, h, body = part(f"{tmp}/tangent.svg"); s = W / w
rows.append(f'<g id="HI_tangent" transform="translate(0 {y:.2f}) scale({s:.4f})">{body}</g>')
letters += [LETTER.format(x=2, y=y + 14, L="H"), LETTER.format(x=0.51 * W, y=y + 14, L="I")]
y += h * s

H = round(y + 2, 1)
svg = ['<?xml version="1.0" encoding="UTF-8"?>',
       f'<svg id="Figure3" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" '
       f'viewBox="0 0 {W} {H}" width="{W}pt" height="{H}pt">',
       f'<rect x="0" y="0" width="{W}" height="{H}" fill="#fff"/>', *rows, '<g id="letters">', *letters, '</g>', '</svg>']
os.makedirs("figures/main", exist_ok=True)
open("figures/main/Figure2_draft.svg", "w").write("\n".join(svg) + "\n")
subprocess.run(["rsvg-convert", "-d", "300", "-p", "300", "-o", "figures/main/Figure2_draft.png", "figures/main/Figure2_draft.svg"], check=True)
print(f"Figure 2: {W} x {H} pt")
