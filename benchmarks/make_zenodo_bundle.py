#!/usr/bin/env python3
"""Assemble the Zenodo archive: the records behind the paper, in the layout the figure scripts expect.

    python make_zenodo_bundle.py <list of recorded inputs> <out dir>      # from benchmarks/

The list comes from figures/record_inputs.py, run over every figure script, so the archive holds what the figures
read. EXTRA adds records that the text quotes but no figure reads. Paths under benchmarks/ keep their place, paths
under scratch go to results/. Code is not copied, it is on GitHub.
"""
import glob, hashlib, os, shutil, sys

listing, out = sys.argv[1:3]
SCRATCH = "/scratch/sadreyev_lab/pschneider/"
HERE = os.getcwd() + os.sep
EXTRA = {
    "timing, attached tangent route": SCRATCH + "kimonu_timing_seeds/tanatt/*",
    "timing, preprocessing": SCRATCH + "kimonu_timing_seeds/fva/*.log",
    "genericity certificate": SCRATCH + "kimonu_cert_modp/*.json",
    "uniform random flux": SCRATCH + "kimonu_sz_flux/*",
    "second-order test": SCRATCH + "kimonu_second_order/*",
    "SVD basis check": SCRATCH + "kimonu_ortho_rref/*",
    "sign condition on the gap networks": "kinetic_variants/injectivity_results/*",
    "two-point criterion": SCRATCH + "kimonu_twopoint/*",
    "both routes, both uptake modes": "kinetic_variants/uptake_modes/*",
    "both routes, input of Figure 3A, B": "kinetic_variants/hierarchy_detached/*",
    "coverage shares": "kinetic_variants/member_counts/*",
    "flux coupling under mass action": "kinetic_variants/fb_vs_ma/*",
    "ACR and ACRR lists per uptake mode": SCRATCH + "kimonu_uptake_modes/lists/*",
    "open questions on the gap networks": SCRATCH + "kimonu_gaproots/*.json",
    "three seeds of the reference flux": SCRATCH + "kimonu_seeds/*.json",
    "modules of the reference-equivalent mode": SCRATCH + "kimonu_trip_2026-11-02/refmode_modules/*.json",
    "uniform random flux, 33 networks": SCRATCH + "kimonu_sz_flux_33/*.json",
    "FVA time": SCRATCH + "kimonu_timing_seeds/fva/fva_*.log",
    "questions the gap comes down to": "kinetic_variants/reduce*",
    "exhaustive reference, all records": "reference_final/**/*",
    "sampled reference, all records": "reference_final_sampled/**/*",
    "networks before blocked-reaction removal": "sweep/*_fixed.crn.json",
    "one LP against FVA": "kinetic_variants/one_lp_results.log",
}
paths = {p.strip() for p in open(listing) if p.strip()}
missing_extra = []
for name, pat in EXTRA.items():
    hit = [p for p in glob.glob(pat, recursive=True) if os.path.isfile(p)]
    if not hit:
        missing_extra.append(name)
    paths |= {os.path.abspath(p) for p in hit}
paths = {p for p in paths if not p.endswith((".py", ".sbatch", ".sh")) and "figures/main" not in p and "/out/" not in p}

manifest, size = [], 0
for p in sorted(paths):
    if not os.path.isfile(p):
        continue
    rel = "results/" + p[len(SCRATCH):] if p.startswith(SCRATCH) else p[len(HERE):] if p.startswith(HERE) else None
    if rel is None:
        continue
    dst = os.path.join(out, rel)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(p, dst)
    size += os.path.getsize(p)
    manifest.append(f"{hashlib.sha256(open(p, 'rb').read()).hexdigest()}  {rel}")
open(os.path.join(out, "SHA256SUMS"), "w").write("\n".join(manifest) + "\n")
print(f"{len(manifest)} files, {size / 1e6:.0f} MB; no files found for: {missing_extra}")
