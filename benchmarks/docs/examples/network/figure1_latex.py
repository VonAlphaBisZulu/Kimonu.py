"""Figure 1: the formulas typeset by LaTeX (Computer Modern) and drawn as outlines.

Replaces the Arial formula runs of Phil's Illustrator export by paths, so the figure needs no math font: Arial
has no glyph for the containment signs, and Illustrator substituted Myriad Pro. Each formula is set by
pdflatex in its own box at the size of the run it replaces, converted to SVG by pdftocairo, and its glyphs
are inlined as plain paths. Words next to a formula ("here all equal", "tangent space at") stay Arial text.
A formula keeps the left edge of the runs it replaces, or their centre where the label is centred.

    python figure1_latex.py <in.svg> <out.svg>       # needs pdflatex (amsmath, bm) and pdftocairo
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latex_text import replace

INK, MUTE = "#1b2430", "#9aa3ad"

# (baseline y, x range of the runs replaced, size, anchor, parts); a part is ("tex", source) or ("text", words)
F = [
    (12.5, (25, 86), 9, "left", [("tex", r"d\mathbf{x}/dt = \mathbf{N}\mathbf{r} = \mathbf{0}")]),
    (12.5, (310, 405), 9, "left", [("tex", r"d\mathbf{x}/dt = \mathbf{N}\mathbf{r} = \mathbf{Y}\mathbf{A}\mathbf{r} = \mathbf{0}")]),
    (24.0, (25, 46), 9, "left", [("tex", r"\mathbf{r} \geq \mathbf{0}")]),
    (24.0, (310, 376), 9, "left", [("tex", r"\mathbf{r} = \mathrm{diag}(\mathbf{k})\,\mathbf{S}^\top\bm{\psi}")]),
    (37.26, (310, 366), 9, "left", [("tex", r"\ln\bm{\psi} = \mathbf{Y}^\top\ln\mathbf{x}")]),
    (244.4, (63, 112), 10, "center", [("tex", r"\mathbf{K} = \ker(\mathbf{N})")]),
    (244.4, (206, 250), 10, "center", [("tex", r"\mathbf{M} = \mathbf{A}\mathbf{K}")]),
    (247.5, (380, 435), 10, "center", [("tex", r"\mathbf{B} = \ker(\mathbf{V}^\top)")]),
    (258.0, (338, 490), 7, "center", [("tex", r"\ker(\mathbf{W}\mathbf{Y}^\top) \subseteq \mathbf{B} \subseteq \ker\big((\mathbf{Y}\Delta)^\top\big)"),
                                       ("text", ", here all equal", MUTE)]),
    (270.4, (372, 392), 8, "center", [("tex", r"\mathbf{Y}^\top\mathbf{B}")]),
    (460.0, (242, 269), 8, "center", [("tex", r"\ker(\mathbf{W})")]),
    (466.0, (435, 495), 8, "left", [("tex", r"\ker\big((\mathbf{Y}\Delta)^\top\big) \supseteq \mathbf{B}")]),
    (471.0, (226, 285), 7, "center", [("tex", r"\mathbf{W} = \mathbf{N}\,\mathrm{diag}(\mathbf{r}^*)\,\mathbf{S}^\top")]),
    (471.0, (356, 402), 7, "center", [("tex", r"\ln\bm{\psi} = \mathbf{Y}^\top\ln\mathbf{x}")]),
    (496.0, (290, 360), 6.5, "center", [("text", "feedback ", MUTE), ("tex", r"u_c = u_d", MUTE)]),
    (509.0, (435, 495), 8, "left", [("tex", r"\ker(\mathbf{W}\mathbf{Y}^\top) \subseteq \mathbf{B}")]),
    (513.2, (254, 372), 8, "center", [("text", "tangent space at ", INK), ("tex", r"\mathbf{r}^*"), ("gap", 8),
                                       ("tex", r"\mathbf{W}\mathbf{Y}^\top\mathbf{h} = \mathbf{0}")]),
]
src, out = sys.argv[1:3]
s, n = replace(open(src, encoding="utf-8").read(), F)
open(out, "w", encoding="utf-8").write(s)
print(f"{len(F)} formulas, {n} runs replaced")
