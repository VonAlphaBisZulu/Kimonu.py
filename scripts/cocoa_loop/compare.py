"""Compare COCOA's kinetic loop (on Kimonu's exact concordance) with Kimonu's iterated route, per model.

Partitions are compared on complexes keyed by composition. "Kimonu coarser" means every COCOA module
lies inside one Kimonu module; "restricted to substrates" drops complexes that are no reaction's
substrate (Kimonu keeps those as singletons by construction, COCOA can attach them by Proposition S3-4).
ACR/ACRR are compared on free metabolites (those before E_1, the enzyme species of the elementary-step
split); ACRR pairs between two ACR metabolites are dropped on both sides, as Kimonu reports them.

    python compare.py <dir with kimonu_*.json and cocoa_*.json> <models.txt> <out.json>
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import records                                         # paths of the records (KIMONU_RECORDS)


def relation(P, Q, universe):
    """'identical', 'P coarser', 'Q coarser' or 'incomparable' for two partitions given as module lists."""
    def label(parts):
        lab = {}
        for i, g in enumerate(parts):
            for c in g:
                lab[c] = i
        return {c: lab.get(c, ("s", c)) for c in universe}
    lp, lq = label(P), label(Q)
    q_in_p = all(len({lp[c] for c in g}) == 1 for g in Q)          # every Q-module inside one P-module
    p_in_q = all(len({lq[c] for c in g}) == 1 for g in P)
    return "identical" if q_in_p and p_in_q else "kimonu coarser" if q_in_p else "cocoa coarser" if p_in_q else "incomparable"


def restrict(P, keep):
    return [r for r in ([c for c in g if c in keep] for g in P) if len(r) > 1]


def free_sets(acr, acrr, free):
    acr = {m for m in acr if m in free}
    acrr = {tuple(sorted(p)) for p in acrr if p[0] in free and p[1] in free and not (p[0] in acr and p[1] in acr)}
    return acr, acrr


def main(d, models, out):
    rows = []
    for M in [l.strip() for l in open(models) if l.strip()]:
        kf, cf = f"{d}/kimonu_{M}.json", f"{d}/cocoa_{M}.json"
        if not (os.path.exists(kf) and os.path.exists(cf)):
            rows.append({"model": M, "status": "missing"}); continue
        K, C = json.load(open(kf)), json.load(open(cf))
        crn = json.load(open(records.network(M)))
        mets = crn["metabolites"]; free = set(mets[:mets.index("E_1")] if "E_1" in mets else mets)
        universe = set(K["balanced"]) | {c for g in K["concordance_modules"] for c in g}
        subs = set(K["substrate"])
        row = {"model": M, "cocoa_seconds": C["seconds"], "cocoa_n": len(C["modules"]),
               "cocoa_giant": max(map(len, C["modules"]), default=0)}
        ca, cr = free_sets(C["acr"], C["acrr"], free)
        row.update(cocoa_acr=len(ca), cocoa_acrr=len(cr))
        for mode in ("attached", "detached"):
            P = K[mode]["modules"]
            ka, kr = free_sets(K[mode]["acr"], K[mode]["acrr"], free)
            Ps, Cs = restrict(P, subs), restrict(C["modules"], subs)
            row[mode] = {"n": len(P), "giant": max(map(len, P), default=0),
                         "relation": relation(P, C["modules"], universe),
                         "relation_substrates": relation(Ps, Cs, subs),
                         "acr": len(ka), "acrr": len(kr),
                         "acr_only_kimonu": sorted(ka - ca), "acr_only_cocoa": sorted(ca - ka),
                         "acrr_only_kimonu": len(kr - cr), "acrr_only_cocoa": len(cr - kr),
                         "acrr_only_cocoa_examples": sorted(cr - kr)[:3]}
        rows.append(row)
    json.dump(rows, open(out, "w"), indent=1)
    print(f"{'model':28s} {'COCOA n/giant':>13s} {'Kim att n/giant':>15s} {'rel(att)':>15s} {'rel subs':>15s} {'Kim det':>9s} {'rel(det) subs':>15s}  ACR c/ka/kd  ACRR c/ka/kd")
    for r in rows:
        if r.get("status"):
            print(f"{r['model']:28s} missing"); continue
        a, t = r["attached"], r["detached"]
        print(f"{r['model']:28s} {r['cocoa_n']:>6d}/{r['cocoa_giant']:<6d} {a['n']:>7d}/{a['giant']:<7d} {a['relation']:>15s} {a['relation_substrates']:>15s} "
              f"{t['n']:>4d}/{t['giant']:<4d} {t['relation_substrates']:>15s}  {r['cocoa_acr']}/{a['acr']}/{t['acr']}  {r['cocoa_acrr']}/{a['acrr']}/{t['acrr']}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
