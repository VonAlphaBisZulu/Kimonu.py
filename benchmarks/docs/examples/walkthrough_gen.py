"""Build nullspace_walkthrough.html: the Langary et al. (2025) example network, directions split,
walked from A to ACR in matrices. Run with the straindesign env and PYTHONPATH=<sd>:<repo>/src."""
import numpy as np, scipy.sparse as sp, sympy as sy, kimonu, os
from kimonu.kinetic import _Topology, _sccs, upstream_set, classical_deficiency
HERE = os.path.dirname(os.path.abspath(__file__))
cx = ['A','B','C','D','D+E','F','B+G','C+G','H','C+E','A+G','I','A+E']
rx = [('A','B'),('B','A'),('B','C'),('C','B'),('C','D'),('D+E','F'),('F','D+E'),('F','B+G'),('C+G','H'),('H','C+G'),('H','C+E'),('A+G','I'),('I','A+G'),('I','A+E')]
mets = list('ABCDEFGHI'); R = [f'R{j+1}' for j in range(14)]
A = np.zeros((13,14),int); Y = np.zeros((9,13),int)
for j,(s,p) in enumerate(rx): A[cx.index(s),j] = -1; A[cx.index(p),j] = 1
for j,c in enumerate(cx):
    for x in c.split('+'): Y[mets.index(x),j] = 1
N = Y @ A
# kernel with the five reversible pairs first, then the two elementary modes
pairs = [(0,1),(2,3),(5,6),(8,9),(11,12)]
Kcols = []
for a,b in pairs:
    v = np.zeros(14,int); v[a]=v[b]=1; Kcols.append(v)
for e in ([2,4,5,7,8,10],[2,4,5,7,11,13]):
    v = np.zeros(14,int); v[e]=1; Kcols.append(v)
K = np.array(Kcols).T; knames = ['R1+R2','R3+R4','R6+R7','R9+R10','R12+R13','EFM 1','EFM 2']
assert not (N @ K).any() and np.linalg.matrix_rank(K)==14-np.linalg.matrix_rank(N)
M = A @ K
As, Ys = sp.csr_matrix(A), sp.csr_matrix(Y)
res = kimonu.concordance_from_matrices(As, sp.csr_matrix(N).astype(np.int64)); bal = set(res.balanced)
topo = _Topology(As); succ = topo.succ
nm = lambda s: sorted((cx[i] for i in s), key=cx.index)
steps = []
for mod in res.modules:
    cur = set(mod)|bal; p1 = []
    while True:
        entry = {c for c in cur if any(any(s not in cur for s in topo.substrates[j]) for j in topo.producers[c])}
        if not entry: break
        p1 += nm(entry); cur -= entry
    t = [c for c in _sccs(succ,cur) if not any(succ[v]-c for v in c)]; p2 = []
    for c in t: p2 += nm(c); cur -= c
    steps.append((nm(mod), p1, p2, nm(cur)))
U1 = [cx.index(c) for c in steps[0][3]]; k = sy.symbols('k1:15'); LU = sy.zeros(len(U1))
for j,(s,p) in enumerate(rx):
    if cx.index(s) in U1:
        i = U1.index(cx.index(s)); LU[i,i] -= k[j]
        if cx.index(p) in U1: LU[U1.index(cx.index(p)),i] += k[j]
det = sy.factor(LU.det())
lam = [int(M[c, 5]) for c in U1]   # activity direction of U1's complexes along EFM 1 (same class along EFM 2)
# S3-4 signatures on the three upstream sets
ups = [s[3] for s in steps]
def leftnull(mods):
    YD = sy.Matrix.hstack(*[sy.Matrix(Y[:,cx.index(m[0])]-Y[:,cx.index(o)]) for m in mods for o in m[1:]])
    return YD, sy.Matrix.hstack(*YD.T.nullspace())
YD0, L0 = leftnull(ups)
sig = [tuple(int(v) for v in L0.T*sy.Matrix(Y[:,cx.index(u[0])])) for u in ups]
delta = classical_deficiency(Ys, As, topo.succ, topo.pred)
final = nm(upstream_set(set().union(*[set(m) for m in res.modules])|bal, topo))
YD, L = leftnull([final]); acr = kimonu.acr_acrr(Ys, [set(cx.index(c) for c in final)], mets)
cls = {c:'bal' for c in nm(bal)}
for n,m in enumerate(res.modules): cls.update({c:f'c{n+1}' for c in nm(m)})
labn = {'bal':'balanced','c1':'class 1','c2':'class 2','c3':'class 3'}

def tbl(Mx, rows, cols, rowcls=None, note=None, fmt=lambda v: str(int(v)), zero=lambda v: v==0):
    h = '<div class="mx"><table><tr><th></th>'+''.join(f'<th>{c}</th>' for c in cols)+('<th></th>' if note else '')+'</tr>'
    for i,r in enumerate(rows):
        c = rowcls.get(r,'') if rowcls else ''
        h += f'<tr class="{c}"><th>{r}</th>'+''.join(f'<td class="{"z" if zero(v) else ""}">{fmt(v)}</td>' for v in np.array(Mx,dtype=object)[i])+(f'<td class="nt">{note(r)}</td>' if note else '')+'</tr>'
    return h+'</table></div>'
svg = open(os.path.join(HERE,'net_fork.svg')).read()
lu_html = tbl([[sy.latex(v).replace('k_{','k').replace('}','') if v!=0 else '0' for v in LU.row(i)] for i in range(LU.rows)], steps[0][3], steps[0][3], fmt=str, zero=lambda v: v=='0')
Lrows = {mets[i]: tuple(int(v) for v in L.row(i)) for i in range(9)}
lcls = {m: ('acr' if all(v==0 for v in Lrows[m]) else ('grp' if Lrows[m]==(0,1) else '')) for m in mets}
steps_html = ''.join(f'<tr><td>{", ".join(m)}</td><td>{", ".join(p1) or "–"}</td><td>{", ".join(p2) or "–"}</td><td><b>{", ".join(u)}</b></td></tr>' for m,p1,p2,u in steps)
sig_html = ''.join(f'<tr><td>{", ".join(u)}</td><td>{u[0]}</td><td>({", ".join(map(str,s))})</td></tr>' for u,s in zip(ups,sig))
css = open(os.path.join(HERE,'walkthrough.css')).read()
html = f'''<title>Nullspace Walkthrough</title>
<style>{css}</style>
<div class="wrap">
<header style="display:flex;flex-direction:column;gap:8px">
<div class="eyebrow">Supplementary walkthrough · draft</div>
<h1>Nullspace walkthrough</h1>
<p>The example network of Langary et al. (2025), from the incidence matrix to absolute concentration robustness, one matrix per step. Every reaction direction is its own column (14 reactions, 13 complexes, 9 metabolites). That is the paper's convention and the one the kinetics needs: under mass action R1 runs at k<sub>1</sub>[A] and R2 at k<sub>2</sub>[B], and the fact that their net flux is always zero is not a dead end to discard but the statement k<sub>1</sub>[A] = k<sub>2</sub>[B], a fixed concentration ratio. Merging the pair into one net reaction keeps the flux picture and loses exactly that (section 7).</p>
</header>

<section><h2>1 · The network and its incidence matrix A</h2>
<p>A records which complex each reaction consumes (−1) and produces (+1): the graph on the left, written as a matrix. A reversible pair is drawn once but occupies two columns. A alone cannot decide balance; that needs what the complexes are made of.</p>
<div class="pair"><div class="card"><svg viewBox="60 60 880 610" role="img" aria-label="Example network"><defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,1L9,5L0,9z" fill="var(--edge)"/></marker></defs>{svg}</svg></div>
<div class="card">{tbl(A,cx,R)}</div></div></section>

<section><h2>2 · Composition Y and the stoichiometric matrix N = Y·A</h2>
<p>Y lists the metabolites in each complex. Multiplying it onto A folds complex bookkeeping back onto metabolites and gives the N of the paper.</p>
<div class="pair"><div class="card">{tbl(Y,mets,cx)}</div><div class="card">{tbl(N,mets,R)}</div></div></section>

<section><h2>3 · Steady states: K = ker(N)</h2>
<p>N has rank 7, so ker(N) has 14 − 7 = 7 dimensions and every steady state is v = K·b. Five columns are the reversible pairs: they are ker(A), cycles in the complex graph that leave every complex where it was. The other two are the elementary modes, steady states only because Y folds complexes back into metabolites.</p>
<div class="card">{tbl(K,R,knames)}</div></section>

<section><h2>4 · Complex activity over all steady states: M = A·K</h2>
<p>Row c of M is the net formation rate of complex c along each steady-state direction. A zero row is a <b>balanced</b> complex. Rows proportional to each other, with any sign, are <b>concordant</b>; a maximal group of them is a concordance module, shown by colour. Scaling each row by its first nonzero entry turns every module into one identical vector, so the modules are read off in one pass with no pairwise test. The reversible-pair columns are all zero: they never affect balance or concordance.</p>
<div class="card">{tbl(M,cx,knames,cls,lambda r: labn[cls[r]])}</div></section>

<section><h2>5 · From concordance to kinetic modules</h2>
<p><b>Mass action in matrix form.</b> Give each complex a concentration ψ<sub>c</sub> = ∏<sub>s</sub> x<sub>s</sub><sup>Y<sub>sc</sub></sup>, so ln ψ = Yᵀ ln x. Each reaction runs at its rate constant times its substrate complex's ψ, v = diag(k)·Sᵀψ with S the substrate part of A. Complex activity is then A·v = L·ψ with L = A·diag(k)·Sᵀ, a weighted Laplacian of the complex graph.</p>
<p><b>Upstream sets.</b> Each concordance module is joined with all balanced complexes and pruned on the graph read from A: Phase I removes, to a fixed point, every complex produced from a substrate outside the set (autonomy); Phase II removes the terminal strong components (feeding).</p>
<div class="card mx"><table class="plain"><tr><th>Concordance module (+ balanced A, C, F, H, I)</th><th>Phase I: entry complexes</th><th>Phase II: terminal</th><th>Upstream set U</th></tr>{steps_html}</table></div>
<p><b>Why an upstream set fixes concentration ratios.</b> Autonomy means every reaction producing a complex of U starts in U, so the activities of U depend on ψ<sub>U</sub> alone: A<sub>U</sub>·v = L<sub>U</sub>·ψ<sub>U</sub>. From M, those activities point in one fixed direction λ<sub>U</sub> scaled by a single number s: zero for the balanced members, the class-1 pattern for the others. Feeding means every complex of U drains out of it, which makes L<sub>U</sub> invertible (a Laplacian with leakage). So ψ<sub>U</sub> = s·L<sub>U</sub><sup>−1</sup>λ<sub>U</sub>: all complex concentrations in U move together with s, in ratios fixed by the rate constants. For U = {{{", ".join(steps[0][3])}}}:</p>
<div class="pair"><div class="card">{lu_html}</div><div class="card eq" style="display:flex;flex-direction:column;gap:8px">
<div>λ<sub>U</sub> = ({", ".join(map(str,lam))})  along EFM 1</div><div>det L<sub>U</sub> = {str(det).replace("*","·")} ≠ 0</div>
<div style="color:var(--muted);font-family:system-ui,sans-serif;font-size:14px">Phase II is what makes this work: a terminal class inside U would have no leakage, its block of L<sub>U</sub> would be a closed Laplacian, and the determinant would vanish.</div></div></div>
<p><b>Merging.</b> Upstream sets that share a complex are one module; these three share none. Proposition S3-4 merges two sets whose representatives' composition difference lies in the span of the differences already fixed, equivalently whose signatures Lᵀy agree, with L = ker(YΔᵀ) over the three sets:</p>
<div class="card mx"><table class="plain"><tr><th>Upstream set</th><th>Representative</th><th>Signature Lᵀy</th></tr>{sig_html}</table></div>
<p>The sets through H and through I share a signature and merge, which also merges their two concordance modules. The network has deficiency δ = {delta} and is not weakly reversible; one concordance merge leaves δ − 1 = 1, and Lemma S4-5 then places every unbalanced complex in one module. Its upstream set is the single kinetic module</p>
<div class="card eq">{{{", ".join(final)}}}</div>
<p>with D, B+G, C+E and A+E, the terminal complexes, left as singletons.</p></section>

<section><h2>6 · ACR and ACRR: the left nullspace of YΔ</h2>
<p>Inside a kinetic module every ratio ψ<sub>c</sub>/ψ<sub>d</sub> is fixed, i.e. (y<sub>c</sub> − y<sub>d</sub>)·ln x is the same in every positive steady state. Collect these differences as the columns of YΔ. A metabolite has ACR when its unit vector lies in their span, which is a zero row of L = ker(YΔᵀ); two metabolites have ACRR when their rows of L are equal.</p>
<div class="pair"><div class="card">{tbl(YD,mets,[c.replace(" ","") for c in [f"{final[0]}−{o}" for o in final[1:]]])}</div>
<div class="card">{tbl(L,mets,["ℓ1","ℓ2"],{m:lcls[m] for m in mets},lambda m: {"acr":"ACR","grp":"equal rows"}.get(lcls[m],""))}</div></div>
<p>G has a zero row: its concentration is the same in every positive steady state, the Shinar–Feinberg result for phosphorylated OmpR. A, B, C, F, H and I share one row, so every ratio among them is fixed: {len(acr.acrr)} ACRR pairs.</p></section>

<section><h2>7 · What merging the directions would have cost</h2>
<div class="card" style="display:flex;flex-direction:column;gap:8px;font-size:14px">
<div>Merging each forward/reverse pair and cutting the directions no steady state uses leaves 8 reactions. R1/R2 disappears (zero net flux), and complex A with it. Balance (C, F, H, I) and concordance (the same three classes) come out unchanged: reversible pairs only ever add zero columns to M.</div>
<div>Kinetic modules do not. In the full network R1/R2 puts complex A into the first upstream set together with C, which fixes the difference A − C. That is exactly what makes the H-branch and the I-branch signatures equal (C+G − A+G = C − A), so they merge, the deficiency drops to one, and Lemma S4-5 builds the module of nine. Without R1/R2 complex A is gone, the three signatures differ, nothing merges, and three small modules remain, {{B, C, D+E, F}}, {{C+G, H}} and {{A+G, I}}, with no ACR at all. Net flux decides balance and concordance; the individual reaction directions decide the kinetics.</div>
<div style="color:var(--muted)">The benchmark networks keep every reversible reaction as two irreversible columns (1798 pairs on iCN900), so this loss never occurs there.</div></div></section>
</div>'''
open(os.path.join(HERE,'nullspace_walkthrough.html'),'w').write(html)
print('written', len(html)//1024, 'KB; final module', final, 'ACR', acr.acr)
