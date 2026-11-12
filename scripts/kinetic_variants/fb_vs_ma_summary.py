"""Totals and per-model table of fb_vs_ma.py over the 33 networks, for both uptake modes.

    python kinetic_variants/fb_vs_ma_summary.py [csv]   -> prints; with "csv" also writes summary.csv to
results/kimonu/flux_balance_vs_mass_action/ of KIMONU_OUT (records.py), from the same folder of KIMONU_RECORDS
"""
import csv, json, os, statistics, sys
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, B)
import records
D = records.read(records.FB_VS_MA)
models = [l.strip() for l in open(f"{B}/kinetic_variants/models33.txt") if l.strip()]


def load(M, mode):
    f = f"{D}/{M}{'' if mode == 'detached' else '.attached'}.json"
    return json.load(open(f)) if os.path.exists(f) else None


def row(d):
    """One flat record per model and mode. Tangent values fall back to the iterated ones where the routes
    were not run separately (attached uptakes without a gap)."""
    it, tan = d["rates_iterated"], d["rates_tangent"] or d["rates_iterated"]
    ci, ct = d["complexes_iterated"], d["complexes_tangent"] or d["complexes_iterated"]
    c, cb = ci["concordant"], ci["concordant_or_balanced"]
    return dict(
        model=d["model"], uptakes=d["uptakes"], reactions=d["n_reactions"], complexes=d["n_complexes"],
        sources=d["n_sources"], gap=d["gap"],
        # reactions
        fc_pairs=d["flux_coupling"]["pairs"], fc_members=d["flux_coupling"]["members"],
        same_substrate_pairs=d["shared_substrate"]["pairs"],
        fc_or_same_substrate_pairs=d["flux_coupling_and_shared_substrate"]["pairs"],
        ma_pairs=it["pairs"], ma_members=it["members"], ma_largest=it["largest"],
        ma_beyond_fc=it["beyond_flux_coupling"], ma_beyond_fc_and_substrate=it["beyond_flux_coupling_and_shared_substrate"],
        ma_fixed_rate=it["fixed_rate"],
        ma_pairs_tangent=tan["pairs"], ma_fixed_rate_tangent=tan["fixed_rate"], rates_routes_agree=d["rates_routes_agree"],
        # complexes
        conc_pairs=c["flux_balance"], conc_members=ci["members_concordance"], balanced=d["concordance"]["balanced"],
        km_pairs=c["mass_action"], km_members=ci["members_modules"], km_pairs_tangent=ct["concordant"]["mass_action"],
        both_pairs=c["both"], only_conc_pairs=c["only_flux_balance"], only_km_pairs=c["only_mass_action"],
        both_pairs_with_balanced=cb["both"], km_members_concordant=ci["members_both"],
        km_members_balanced=ci["balanced_in_modules"],
        conc_refines_km=c["flux_balance_refines_mass_action"])


def med(v):
    return statistics.median(v) if v else float("nan")


if __name__ == "__main__":
    rows = []
    for mode in ("detached", "attached"):
        R = [row(d) for d in (load(M, mode) for M in models) if d]
        rows += R
        if not R:
            continue
        S = lambda k: sum(r[k] for r in R)
        print(f"\n=== {mode} uptakes, {len(R)} of {len(models)} networks ===")
        print("reactions", S("reactions"), "| complexes", S("complexes"), "| sources", S("sources"))
        print("-- reaction pairs with a fixed rate ratio")
        print(f"flux balance (full coupling)            {S('fc_pairs'):>12,}")
        print(f"  + shared substrate complex            {S('fc_or_same_substrate_pairs'):>12,}")
        print(f"mass action, iterated route             {S('ma_pairs'):>12,}   tangent {S('ma_pairs_tangent'):,}")
        print(f"  beyond flux coupling                  {S('ma_beyond_fc'):>12,}")
        print(f"  beyond flux coupling and substrate    {S('ma_beyond_fc_and_substrate'):>12,}")
        print(f"reactions in a class: flux balance {S('fc_members'):,}, mass action {S('ma_members'):,} "
              f"({100 * S('ma_members') / S('reactions'):.1f} % of reactions; flux balance "
              f"{100 * S('fc_members') / S('reactions'):.1f} %)")
        print(f"fixed rates: iterated {S('ma_fixed_rate'):,} in {sum(r['ma_fixed_rate'] > 0 for r in R)} networks, "
              f"tangent {S('ma_fixed_rate_tangent'):,}")
        print("routes agree on rate classes in", sum(bool(r["rates_routes_agree"]) for r in R), "of", len(R),
              "| not:", [r["model"] for r in R if not r["rates_routes_agree"]])
        ratio = [r["ma_pairs"] / r["fc_pairs"] for r in R if r["fc_pairs"]]
        print(f"pairs mass action / flux balance per network: median {med(ratio):.1f}, "
              f"min {min(ratio):.1f}, max {max(ratio):.1f}")
        ratio = [r["ma_pairs"] / r["fc_or_same_substrate_pairs"] for r in R if r["fc_or_same_substrate_pairs"]]
        print(f"pairs mass action / (flux coupling + substrate) per network: median {med(ratio):.2f}, "
              f"min {min(ratio):.2f}, max {max(ratio):.2f}")
        print("-- complex pairs")
        print(f"concordant (flux balance)               {S('conc_pairs'):>12,}")
        print(f"kinetic modules, iterated route         {S('km_pairs'):>12,}   tangent {S('km_pairs_tangent'):,}")
        print(f"  in both                               {S('both_pairs'):>12,}")
        print(f"  only concordant                       {S('only_conc_pairs'):>12,}")
        print(f"  only kinetic                          {S('only_km_pairs'):>12,}")
        print(f"  in both, balanced taken as one group  {S('both_pairs_with_balanced'):>12,}")
        print(f"balanced complexes {S('balanced'):,} | complexes in kinetic modules {S('km_members'):,}, of these "
              f"{S('km_members_balanced'):,} balanced and {S('km_members_concordant'):,} in a concordance module")
        share = [100 * r["both_pairs"] / r["km_pairs"] for r in R if r["km_pairs"]]
        print(f"share of kinetic pairs that are concordant, per network: median {med(share):.1f} %, "
              f"min {min(share):.1f} %, max {max(share):.1f} %")
        share = [100 * r["both_pairs"] / r["conc_pairs"] for r in R if r["conc_pairs"]]
        print(f"share of concordant pairs that are kinetic, per network: median {med(share):.2f} %, "
              f"min {min(share):.2f} %, max {max(share):.2f} %")
        print("concordance refines the kinetic modules in", sum(bool(r["conc_refines_km"]) for r in R), "networks")
    if "csv" in sys.argv[1:] and rows:
        with open(f"{D}/summary.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        print(f"\nwrote {D}/summary.csv ({len(rows)} rows)")
