import sys, json, numpy as np
sys.path.insert(0,'/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks'); sys.path.insert(0,'/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/src')
import kimonu
from compat import upstream
from run_kimonu_crn import triplets_to_sparse, to_exact
sys.path.insert(0, '/home/sadreyev_lab/pschneider/pipelines/2026-08_Kimonu/Kimonu.py/benchmarks/kinetic_variants')
from measure_phases import load, tentative
from compat.upstream import _Topology
for M in sys.argv[1:]:
    crn, nc, A, Y, S = load(M); topo=_Topology(A); succ,pred=topo.succ,topo.pred; nonr={c for c in range(nc) if not succ[c]}
    res=kimonu.concordance_from_matrices(A,S); bal=set(res.balanced)
    seeds=[set(m)-bal for m in res.modules]; seeds=[m for m in seeds if m]
    tsets=[t for t in (tentative(m|bal,succ,pred,nonr)[0] for m in seeds) if t]
    merged=[]
    for g in tsets:
        hit=[x for x in merged if x&g]
        for x in hit: merged.remove(x); g=g|x
        merged.append(g)
    gt=max((len(x) for x in merged),default=0); nt=sum(1 for x in merged if len(x)>1)
    kY=upstream.kinetic_modules(A,res.modules,res.balanced,Y=Y,nullity=res.nullity); k0=upstream.kinetic_modules(A,res.modules,res.balanced,Y=None)
    print(f"{M:26s} 4-phase overlap-only: giant {gt:4d} mods {nt:4d} | 2-phase overlap-only: giant {k0.giant_module:4d} mods {k0.n_modules:4d} | 2-phase +coupling (Kimonu default): giant {kY.giant_module:4d} mods {kY.n_modules:4d} | ref R giant per figdata: {json.load(open('figdata.json'))['kimonu'][M]['giant_4phase']}", flush=True)
