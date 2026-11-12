"""Toric check on the original (unsplit) model, the most binomial an intermediate-eliminated core could be.

Eliminating enzyme intermediates turns each mechanism back into one overall reaction with a rational rate;
distinct enzymes keep distinct rate terms, so the core cannot be closer to binomial than the lumped network
with monomial rates. Boundary species become the zero complex; blocked reactions are removed by one LP.

    python kinetic_variants/lumped_toric.py <model>
"""
import os, sys, xml.etree.ElementTree as ET, numpy as np, scipy.sparse as sp
from fractions import Fraction
from scipy.optimize import linprog
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [B, os.path.join(os.path.dirname(B), 'src')]
from kinetic_variants.toric_check import check
from kimonu.linalg import exact
M = sys.argv[1]
U = os.environ.get('UPSTREAM_ALGORITHM', os.path.join(B, 'reference', 'Upstream_Algorithm'))   # bash reference/fetch_upstream.sh --full
root = ET.parse(f'{U}/Models/original/{M}.xml').getroot()
ns = {'s': root.tag.split('}')[0].strip('{')}
rx = []
for r in root.iter(f"{{{ns['s']}}}reaction"):
    rev = r.get('reversible', 'false') == 'true'
    def side(tag):
        out = {}
        for lst in r:
            if lst.tag.endswith(tag):
                for s in lst:
                    out[s.get('species')] = out.get(s.get('species'), 0) + Fraction(s.get('stoichiometry', '1'))
        return out
    sub, pro = side('listOfReactants'), side('listOfProducts')
    rx.append((sub, pro)); 
    if rev: rx.append((pro, sub))
# boundary species (_b / boundaryCondition) are dropped: their complexes become the zero complex
bnd = {s.get('id') for s in root.iter(f"{{{ns['s']}}}species") if s.get('boundaryCondition') == 'true'}
rx = [({k: v for k, v in a.items() if k not in bnd}, {k: v for k, v in b.items() if k not in bnd}) for a, b in rx]
mets = sorted({m for a, b in rx for m in list(a) + list(b)})
mi = {m: i for i, m in enumerate(mets)}; q = len(rx); m_ = len(mets)
Nf = np.zeros((m_, q))
for j, (a, b) in enumerate(rx):
    for k, v in a.items(): Nf[mi[k], j] -= float(v)
    for k, v in b.items(): Nf[mi[k], j] += float(v)
# unblocked reactions: one LP maximising sum t, t <= v, t <= 1 gives the full support of the flux cone
Asp = sp.csr_matrix(Nf)
c = np.concatenate([np.zeros(q), -np.ones(q)])
A_ub = sp.hstack([-sp.identity(q), sp.identity(q)]).tocsr()
lp = linprog(c, A_ub=A_ub, b_ub=np.zeros(q), A_eq=sp.hstack([Asp, sp.csr_matrix((m_, q))]).tocsr(), b_eq=np.zeros(m_),
             bounds=[(0, None)] * q + [(0, 1)] * q, method="highs")
keep = [j for j in range(q) if lp.x[q + j] > 0.5]
rx = [rx[j] for j in keep]
mets = sorted({m for a, b in rx for m in list(a) + list(b)}); mi = {m: i for i, m in enumerate(mets)}; q = len(rx); m_ = len(mets)
cx = []; cidx = {}
def cid(side):
    key = tuple(sorted(side.items()))
    if key not in cidx: cidx[key] = len(cx); cx.append(key)
    return cidx[key]
Ae, Ye, Ne = [], [], []
for j, (a, b) in enumerate(rx):
    Ae += [(cid(a), j, -1), (cid(b), j, 1)]
    for k, v in a.items(): Ne.append((mi[k], j, -v))
    for k, v in b.items(): Ne.append((mi[k], j, v))
for c_, key in enumerate(cx):
    for k, v in key: Ye.append((mi[k], c_, float(v)))
A = sp.csr_matrix(([v for _, _, v in Ae], ([i for i, _, _ in Ae], [j for _, j, _ in Ae])), shape=(len(cx), q), dtype=np.int64)
Y = sp.csr_matrix(([v for _, _, v in Ye], ([i for i, _, _ in Ye], [j for _, j, _ in Ye])), shape=(m_, len(cx)))
Ndict = {}
for i, j, v in Ne: Ndict[(i, j)] = Ndict.get((i, j), 0) + v
N = exact(((i, j, v) for (i, j), v in Ndict.items() if v != 0), (m_, q))
print(f"{M} original: {m_} metabolites, {q} unblocked reactions, {len(cx)} complexes", flush=True)
rounds = check(A, Y, N)
f, l = rounds[0], rounds[-1]
print(f"   initial: dim {f[0]}, classes {f[1]} -> {f[1] - f[0]} relations | after feedback: dim {l[0]}, classes {l[1]} -> {l[1] - l[0]} relations")
