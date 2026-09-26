"""Figure 1: the Langary et al. (2025) example as metabolic network and as complex network, with the
nullspace readouts under each: flux coupling and elementary modes under the metabolic view, balance,
concordance and the kinetic module under the complex view."""
import numpy as np, os, sys, sympy as sy
OUT = sys.argv[1]
cx = ['A','B','C','D','D+E','F','B+G','C+G','H','C+E','A+G','I','A+E']
rx = [('A','B'),('B','A'),('B','C'),('C','B'),('C','D'),('D+E','F'),('F','D+E'),('F','B+G'),('C+G','H'),('H','C+G'),('H','C+E'),('A+G','I'),('I','A+G'),('I','A+E')]
mets = list('ABCDEFGHI'); R = [f'R{j+1}' for j in range(14)]
A = np.zeros((13,14),int); Y = np.zeros((9,13),int)
for j,(s,p) in enumerate(rx): A[cx.index(s),j] = -1; A[cx.index(p),j] = 1
for j,c in enumerate(cx):
    for x in c.split('+'): Y[mets.index(x),j] = 1
N = Y @ A
Kc = []
for a,b in [(0,1),(2,3),(5,6),(8,9),(11,12)]:
    v = np.zeros(14,int); v[a]=v[b]=1; Kc.append(v)
for e in ([2,4,5,7,8,10],[2,4,5,7,11,13]):
    v = np.zeros(14,int); v[e]=1; Kc.append(v)
K = np.array(Kc).T; assert not (N@K).any()
M = A @ K
# flux coupling: proportional rows of K
def groups(rows):
    g = {}
    for i,r in enumerate(rows):
        nz = [x for x in r if x]
        if not nz: continue
        g.setdefault(tuple(sy.Rational(int(x), int(nz[0])) for x in r), []).append(i)
    return [v for v in g.values() if len(v) > 1]
fc = groups(K)
cls = {'A':'bal','C':'bal','F':'bal','H':'bal','I':'bal','B':'c1','D':'c1','D+E':'c1','B+G':'c1','C+G':'c2','C+E':'c2','A+G':'c3','A+E':'c3'}
kin = {'A','B','C','D+E','F','C+G','H','A+G','I'}
col = {'bal':'#D9DDE3','c1':'#9EC4E8','c2':'#F6C08A','c3':'#A8D8A0'}
ink = '#1B2430'; mut = '#5B6673'; edge = '#4B5563'
W, H = 1760, 1340
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="Helvetica, Arial, sans-serif">',
       f'<rect width="{W}" height="{H}" fill="white"/>',
       '<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,1L9,5L0,9z" fill="#4B5563"/></marker></defs>',
       '<style>.edges path{fill:none;stroke:#4B5563;stroke-width:1.8}.edges .stem{stroke-width:3}.sp{fill:#E3EEF9;stroke:#377EB8;stroke-width:1.8}.hub{fill:#FDE7D3;stroke:#D9731A;stroke-width:1.8}.sl{font:600 18px Helvetica,Arial,sans-serif;text-anchor:middle;fill:#1B2430}.rl{font:13px Helvetica,Arial,sans-serif;fill:#374151;paint-order:stroke;stroke:white;stroke-width:4px}</style>']
def label(x, y, t, sub):
    svg.append(f'<text x="{x}" y="{y}" font-size="30" font-weight="700" fill="{ink}">{t}</text><text x="{x+34}" y="{y}" font-size="20" fill="{mut}">{sub}</text>')
# A: metabolic network (curved-fork drawing), scaled into the left half
label(20, 44, 'A', 'Metabolic network: metabolites and reactions')
body = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'net_fork.svg')).read()
svg.append(f'<g transform="translate(-10,40) scale(0.95)">{body}</g>')
# B: complex network, 4 linkage classes
label(900, 44, 'B', 'Complex network: complexes and reactions')
rows = [['A','B','C','D'], ['D+E','F','B+G'], ['C+G','H','C+E'], ['A+G','I','A+E']]
rl = [[('R1','R2'),('R3','R4'),('R5',None)], [('R6','R7'),('R8',None)], [('R9','R10'),('R11',None)], [('R12','R13'),('R14',None)]]
x0, y0, dx, dy, bw, bh = 930, 110, 205, 150, 120, 48
pos = {}
for i,row in enumerate(rows):
    for j,c in enumerate(row): pos[c] = (x0 + j*dx, y0 + i*dy)
for i,row in enumerate(rows):
    for j,(f,b) in enumerate(rl[i]):
        (xa,ya),(xb,yb) = pos[row[j]], pos[row[j+1]]
        if b:
            svg.append(f'<path d="M{xa+bw+4},{ya+bh/2-8} H{xb-6}" stroke="{edge}" stroke-width="1.8" marker-end="url(#ah)"/><path d="M{xb-4},{ya+bh/2+8} H{xa+bw+6}" stroke="{edge}" stroke-width="1.8" marker-end="url(#ah)"/>')
            svg.append(f'<text x="{(xa+bw+xb)/2}" y="{ya+bh/2-16}" font-size="14" text-anchor="middle" fill="{mut}">{f}</text><text x="{(xa+bw+xb)/2}" y="{ya+bh/2+30}" font-size="14" text-anchor="middle" fill="{mut}">{b}</text>')
        else:
            svg.append(f'<path d="M{xa+bw+4},{ya+bh/2} H{xb-6}" stroke="{edge}" stroke-width="1.8" marker-end="url(#ah)"/><text x="{(xa+bw+xb)/2}" y="{ya+bh/2-10}" font-size="14" text-anchor="middle" fill="{mut}">{f}</text>')
for c,(x,y) in pos.items():
    k = c in kin
    svg.append(f'<rect x="{x}" y="{y}" width="{bw}" height="{bh}" rx="8" fill="{col[cls[c]]}" stroke="{ink if k else "#9AA3AE"}" stroke-width="{3 if k else 1.4}" {"" if k else "stroke-dasharray=\"5 4\""}/><text x="{x+bw/2}" y="{y+bh/2+7}" font-size="20" font-weight="600" text-anchor="middle" fill="{ink}">{c}</text>')
ly = y0 + 3*dy + bh + 26
items = [(col['bal'],'balanced'),(col['c1'],'concordance module 1'),(col['c2'],'concordance module 2'),(col['c3'],'concordance module 3')]
for n,(c,t) in enumerate(items):
    xx = x0 + (n % 2) * 360; yy = ly + (n // 2) * 30
    svg.append(f'<rect x="{xx}" y="{yy}" width="22" height="18" rx="4" fill="{c}"/><text x="{xx+30}" y="{yy+15}" font-size="16" fill="{ink}">{t}</text>')
svg.append(f'<rect x="{x0}" y="{ly+62}" width="22" height="18" rx="4" fill="white" stroke="{ink}" stroke-width="3"/><text x="{x0+30}" y="{ly+77}" font-size="16" fill="{ink}">kinetic module (bold outline)</text>')
svg.append(f'<rect x="{x0+360}" y="{ly+62}" width="22" height="18" rx="4" fill="white" stroke="#9AA3AE" stroke-dasharray="5 4"/><text x="{x0+390}" y="{ly+77}" font-size="16" fill="{ink}">terminal, not in a module</text>')
# matrix helper
def matrix(x, y, Mx, rlab, clab, cw=34, rh=26, rowfill=None, colhl=None, title=None):
    if title: svg.append(f'<text x="{x}" y="{y-14}" font-size="17" fill="{ink}">{title}</text>')
    for j,c in enumerate(clab):
        fill = '#FFF3C4' if colhl and j in colhl else 'none'
        svg.append(f'<rect x="{x+70+j*cw}" y="{y}" width="{cw}" height="{rh*(len(rlab)+1)}" fill="{fill}"/>')
        svg.append(f'<text x="{x+70+j*cw+cw/2}" y="{y+18}" font-size="13" text-anchor="middle" fill="{mut}">{c}</text>')
    for i,r in enumerate(rlab):
        yy = y + rh*(i+1)
        if rowfill and r in rowfill: svg.append(f'<rect x="{x}" y="{yy}" width="64" height="{rh-2}" rx="4" fill="{rowfill[r]}"/>')
        svg.append(f'<text x="{x+6}" y="{yy+18}" font-size="15" fill="{ink}">{r}</text>')
        for j in range(len(clab)):
            v = int(Mx[i][j])
            svg.append(f'<text x="{x+70+j*cw+cw/2}" y="{yy+18}" font-size="15" text-anchor="middle" fill="{ink if v else "#B4BCC6"}">{v}</text>')
# C: K with flux coupling and elementary modes
label(20, 780, 'C', 'Steady states K = ker(N): flux coupling and elementary modes')
kl = ['R1+R2','R3+R4','R6+R7','R9+R10','R12+R13','EFM1','EFM2']
fcf = {}
palette = ['#C9B6E4','#F4D35E']
for n,g in enumerate(fc):
    for i in g: fcf[R[i]] = palette[n % 2]
matrix(40, 820, K, R, ['k1','k2','k3','k4','k5','e1','e2'], cw=62, rh=24, rowfill=fcf, colhl={5,6}, title=None)
nx = 40 + 70 + 7*62 + 30
svg.append(f'<text x="{nx}" y="850" font-size="16" fill="{ink}">k1–k5: reversible pairs</text><text x="{nx}" y="874" font-size="16" fill="{ink}">e1, e2: elementary modes</text>')
svg.append(f'<text x="{nx}" y="922" font-size="16" fill="{ink}">Fully coupled reactions</text><text x="{nx}" y="944" font-size="16" fill="{ink}">(proportional rows):</text>')
for n,g in enumerate(fc):
    svg.append(f'<rect x="{nx}" y="{960+n*28}" width="22" height="18" rx="4" fill="{palette[n % 2]}"/><text x="{nx+30}" y="{975+n*28}" font-size="16" fill="{ink}">{{{", ".join(R[i] for i in g)}}}</text>')
# D: M = A·K
label(900, 780, 'D', 'Complex activity M = A·K: balance and concordance')
Mr = [M[cx.index(c)] for c in cx]
matrix(930, 820, [r[5:] for r in Mr], cx, ['e1','e2'], cw=60, rh=24, rowfill={c: col[cls[c]] for c in cx})
tx = 930 + 70 + 2*60 + 40
notes = ['Zero rows: balanced (A, C, F, H, I).', 'Proportional rows: concordant.', 'Columns k1–k5 are zero for every', 'complex and are omitted.', '', 'Kinetic module (panel B, bold):', 'fixed ratios of complex concentrations', 'ψ under mass action, from', 'ker(N·diag(v*)·Sᵀ) and ker(YΔᵀ).', 'G has ACR.']
for n,t in enumerate(notes): svg.append(f'<text x="{tx}" y="{850+n*26}" font-size="16" fill="{ink}">{t}</text>')
svg.append('</svg>')
open(os.path.join(OUT,'Figure1_draft.svg'),'w').write('\n'.join(svg)); print('fully coupled', [[R[i] for i in g] for g in fc])
