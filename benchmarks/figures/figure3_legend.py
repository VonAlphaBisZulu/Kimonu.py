"""The shared legend of Figure 3A, B without colons, edited in place on the vault Figure 3.

"decided: ..." and "undecided: ..." take a comma, which has the width of the colon in Arial, so nothing moves. "+n: ..."
loses its colon, so the end of that line moves left by the width of one colon. Only the legend band may change, which
is checked by rendering both files.

    python figures/figure3_legend.py <vault Figure3.svg> <out.svg>     # from benchmarks/
"""
import os, sys, tempfile
import numpy as np
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_finish

src, dest = sys.argv[1], sys.argv[2]
text = open(src, encoding="utf-8").read()
i = text.index('<g id="panel_legend_AB"'); j = text.index("</g>", i)
band = text[i:j]
for a, b in (("decided: iterated = tangent", "decided, iterated = tangent"),
             ("undecided: from iterated route", "undecided, from iterated route"),
             ("+n: candidates that require", "+n candidates that require")):
    assert band.count(a) == 1, a
    band = band.replace(a, b)
new = text[:i] + band + text[j:]
open(dest, "w", encoding="utf-8").write(new)

tmp = tempfile.mkdtemp()
a, b = os.path.join(tmp, "a.png"), os.path.join(tmp, "b.png")
svg_finish._render(text, a); svg_finish._render(new, b)
A = np.asarray(Image.open(a).convert("L"), dtype=int); B = np.asarray(Image.open(b).convert("L"), dtype=int)
k = A.shape[1] / 518.4
ch = np.abs(A - B) > 40
lo, hi = int(151.5 * k), int(163.5 * k)                   # the legend band between the rows A, B and C, D
outside = int(ch[:lo].sum() + ch[hi:].sum())
print(f"pixels changed: legend band {int(ch[lo:hi].sum())}, elsewhere {outside}")
assert outside == 0
