"""Supplementary Text S1: the example network of Figure 1 (EnvZ-OmpR) matrix by matrix, as Markdown with LaTeX.

Every matrix is computed by Kimonu's own functions, exactly, in the species and complex names of Figure 1
(A = EnvZ-ADP, B = EnvZ, C = EnvZ-ATP, D = EnvZ-P, E = OmpR, F = EnvZ-P.OmpR, G = OmpR-P, H = EnvZ-ATP.OmpR-P,
I = EnvZ-ADP.OmpR-P). The interior flux is r* = K 1, as in Figure 1F.

    python figures/textS1_walkthrough.py > textS1_walkthrough.md      # from scripts/
"""
import os, sys
from fractions import Fraction
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src")); sys.path.insert(0, os.path.join(ROOT, "tests"))
from networks import matrices, ENVZ, ENVZ_SPECIES
from kimonu.linalg import matrix_product, nullspace
from kimonu.readout import row_readout
from kimonu import crn_from_matrices, free_directions, readout, steady_state_map, tangent_route

NAME = dict(zip(ENVZ_SPECIES, "ABCDEFGHI"))
A, Y, N, cx = matrices(ENVZ, ENVZ_SPECIES)
crn = crn_from_matrices(A, Y)
CX = ["+".join(NAME[s] for s in c.split("+")) for c in cx]
RX = [f"R{j + 1}" for j in range(A.shape[1])]
SP = [NAME[s] for s in ENVZ_SPECIES]


def frac(v):
    v = Fraction(v)
    return str(v.numerator) if v.denominator == 1 else rf"{'-' if v < 0 else ''}\tfrac{{{abs(v.numerator)}}}{{{v.denominator}}}"


def dense(rows, ncols):
    """a row form (or rows as dicts) -> list of lists of Fractions"""
    return [[Fraction(dict(r).get(j, 0)) for j in range(ncols)] for r in rows]


def tex(M, rows, cols, name):
    """A matrix with row and column labels, as a LaTeX array."""
    head = " & ".join(rf"\scriptstyle\text{{{c}}}" for c in cols)
    body = r" \\ ".join(rf"\scriptstyle\text{{{r}}} & " + " & ".join(frac(v) for v in row) for r, row in zip(rows, M))
    return (rf"$$\mathbf{{{name}}} = \begin{{array}}{{r|{'r' * len(cols)}}} & {head} \\ \hline {body} \end{{array}}$$")


def as_rows(Mnp):
    if hasattr(Mnp, "toarray"):
        Mnp = Mnp.toarray()
    return [{j: Fraction(int(v)) for j, v in enumerate(r) if v} for r in np.asarray(Mnp)]


out = []
p = out.append
p("### S1.1 The network and its matrices\n")
p(f"Nine species, {len(CX)} complexes and {len(RX)} reactions. Each column of $\\mathbf{{N}}$, $\\mathbf{{A}}$ is a reaction; "
  "$\\mathbf{N} = \\mathbf{Y}\\mathbf{A}$.\n")
p(tex(dense(as_rows(N), len(RX)), SP, RX, "N"))
p(tex(dense(as_rows(Y), len(CX)), SP, CX, "Y"))
p(tex(dense(as_rows(A), len(RX)), CX, RX, "A"))

# flux coupling: K = ker N
Krows = nullspace(crn.N)
d = max((j for r in Krows for j, _ in r), default=-1) + 1
p("\n### S1.2 Flux coupling: $\\mathbf{K} = \\ker(\\mathbf{N})$\n")
p(f"An exact basis with {d} columns. Zero rows are blocked reactions (none), proportional rows are fully coupled reactions.\n")
p(tex(dense(Krows, d), RX, [f"k{j + 1}" for j in range(d)], "K"))
zr, cl, _ = row_readout(Krows)
p("Fully coupled: " + "; ".join("{" + ", ".join(RX[i] for i in g) + "}" for g in cl if len(g) > 1) + ".\n")

# concordance: M = A K
Mrows = matrix_product(crn.A, Krows)
p("\n### S1.3 Balanced and concordant complexes: $\\mathbf{M} = \\mathbf{A}\\mathbf{K}$\n")
p(tex(dense(Mrows, d), CX, [f"k{j + 1}" for j in range(d)], "M"))
zr, cl, _ = row_readout(Mrows)
p("Balanced (zero rows): " + ", ".join(CX[i] for i in zr) + ". Concordance modules (proportional rows): "
  + "; ".join("{" + ", ".join(CX[i] for i in g) + "}" for g in cl if len(g) > 1) + ".\n")

# iterated route at r* = K 1
flux = [sum((v for _, v in r), Fraction(0)) for r in Krows]
p("\n### S1.4 The iterated route at $\\mathbf{r}^* = \\mathbf{K}\\mathbf{1}$\n")
p("$\\mathbf{r}^* = (" + ", ".join(frac(v) for v in flux) + ")^\\top$.\n")
W, nc = steady_state_map(crn, flux, False)
p(tex(dense(W, nc), SP, CX, r"W = \mathbf{N}\,\mathrm{diag}(\mathbf{r}^*)\,\mathbf{S}^\top"))
equal, rnd = [], 0
while True:
    rnd += 1
    Wrows = nullspace(W + [tuple(sorted(((c, 1), (e, -1)))) for c, e in equal], nc)
    dw = max((j for r in Wrows for j, _ in r), default=-1) + 1
    _, cls, _ = row_readout(Wrows)
    sets = [set(g) for g in cls]
    p(f"\n**Round {rnd}.** $\\ker(\\mathbf{{W}})$" + (" with the feedback rows $u_c - u_d$" if equal else "")
      + f" has {dw} dimensions. Equal rows give the coupled sets "
      + "; ".join("{" + ", ".join(CX[i] for i in sorted(s)) + "}" for s in sets if len(s) > 1) + ".\n")
    if rnd == 1:
        p(tex(dense(Wrows, dw), CX, [f"w{j + 1}" for j in range(dw)], r"\ker(W)"))
    Lrows = free_directions(crn, [s for s in sets if s])
    dl = max((j for r in Lrows for j, _ in r), default=-1) + 1
    p(f"Log closure: $\\ker((\\mathbf{{Y}}\\Delta)^\\top)$ has {dl} dimensions.\n")
    p(tex(dense(Lrows, dl), SP, [f"l{j + 1}" for j in range(dl)], r"\ker((Y\Delta)^\top)"))
    found = readout(crn, Lrows, detach_uptakes=False)
    merged, acr = found.modules, set(found.acr)
    p("Readout: kinetic modules " + "; ".join("{" + ", ".join(CX[i] for i in sorted(s)) + "}" for s in merged if len(s) > 1)
      + (", ACR " + ", ".join(SP[i] for i in sorted(acr)) if acr else ", no ACR yet") + ".\n")
    if {frozenset(s) for s in merged} == {frozenset(s) for s in sets}:
        break
    equal = [(min(s), c) for s in merged for c in s if c != min(s)]

# final readout and tangent route
R = readout(crn, free_directions(crn, merged), detach_uptakes=False)
p("\n### S1.5 Final readout\n")
p("Kinetic modules: " + "; ".join("{" + ", ".join(CX[i] for i in sorted(s)) + "}" for s in R.modules if len(s) > 1)
  + ". ACR: " + ", ".join(SP[i] for i in R.acr) + ". ACRR groups: "
  + "; ".join("{" + ", ".join(SP[i] for i in g) + "}" for g in R.acrr_groups) + ".\n")
T = tangent_route(crn)
p(f"The tangent route, $\\ker(\\mathbf{{W}}\\mathbf{{Y}}^\\top)$ at a random interior point, gives the same modules, ACR and ACRR: "
  f"{'yes' if sorted(map(sorted, (m for m in T.modules if len(m) > 1))) == sorted(map(sorted, (m for m in R.modules if len(m) > 1))) else 'NO'}.\n")
print("\n".join(out))
