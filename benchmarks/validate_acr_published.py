#!/usr/bin/env python
"""Validate Kimonu's ACR/ACRR against the reference implementation's published results.

The reference ships, per model, the elementary-step-split network it analysed and the metabolites
and metabolite pairs it reported:

    Results/MetSingle/MetSingle_<model>_concordant_<variant>.csv     ACR, all species
    Results/MetSingle/M_MetSingle_<model>_...csv                     ACR, free metabolites only
    Results/MetDouble/MetDouble_<model>_...csv                       ACRR, all species
    Results/MetDouble/M_MetDouble_<model>_...csv                     ACRR, free metabolites only

"Free" excludes enzymes and enzyme-substrate species, which elementary-step splitting appends after
the original metabolites; the reference takes everything before the species named `E_1`. The same
cut is applied here, so the two sides are counting the same thing.

Sets are compared, not counts: agreeing on how many ACR metabolites there are while disagreeing on
which ones would not be agreement.

Usage:
    python benchmarks/validate_acr_published.py <upstream_repo> [model ...] [--variant fixed|random]
"""
import argparse
import csv
import os
import sys
import time
from fractions import Fraction
from math import gcd

import h5py
import numpy as np
from scipy import sparse

import straindesign as sd

import kimonu
from compat import upstream
from kimonu.io import exact_matrix


def read_sparse(group):
    nrows = int(group.attrs["MATLAB_sparse"])
    data = np.asarray(group["data"]).ravel()
    ir = np.asarray(group["ir"]).ravel().astype(np.int64)
    jc = np.asarray(group["jc"]).ravel().astype(np.int64)
    return sparse.csc_matrix((data, ir, jc), shape=(nrows, len(jc) - 1))


def read_strings(fh, dataset):
    """MATLAB cell array of char, as stored in a v7.3 file."""
    out = []
    for ref in np.asarray(dataset).ravel():
        chars = np.asarray(fh[ref]).ravel()
        out.append("".join(chr(c) for c in chars))
    return out


def to_exact(M):
    """Exact rational form of ``M`` for the kernel computation."""
    M = M.tocoo()
    return exact_matrix(zip(M.row, M.col, M.data), M.shape)



def read_published(path):
    """Read one of the reference's published result files.

    Two formats are in play and they are easy to confuse:

    * ``MetSingle_*`` / ``MetDouble_*`` come from R's ``write.table``: whitespace separated, every
      field quoted, a leading empty field for the row name, and a **header line** naming the columns
      ``"V1"`` (and ``"V2"`` for pairs). Reading that header as data invents a species called V1.
    * ``M_MetSingle_*`` / ``M_MetDouble_*`` come from MATLAB's ``writecell``: comma separated, no
      quotes, no header.
    """
    if not os.path.exists(path):
        return None

    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        return []

    quoted_header = lines[0].replace('"', "").split() in (["V1"], ["V1", "V2"])
    entries = []
    if quoted_header:
        for line in lines[1:]:
            fields = [f.strip('"') for f in line.split() if f.strip('"')]
            if not fields:
                continue
            entries.append(fields[0] if len(fields) == 1 else tuple(sorted(fields)))
    else:
        for row in csv.reader(lines):
            fields = [c.strip().strip('"').strip() for c in row if c.strip()]
            if not fields:
                continue
            entries.append(fields[0] if len(fields) == 1 else tuple(sorted(fields)))
    return entries


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("upstream")
    ap.add_argument("models", nargs="*")
    ap.add_argument("--variant", default="fixed", choices=["fixed", "random"])
    args = ap.parse_args()

    model_dir = os.path.join(args.upstream, "Models", "models_with_elementary_steps")
    single_dir = os.path.join(args.upstream, "Results", "MetSingle")
    double_dir = os.path.join(args.upstream, "Results", "MetDouble")

    models = args.models or sorted(
        f[: -len(f"_pre_balanced_{args.variant}.mat")] for f in os.listdir(model_dir)
        if f.endswith(f"_pre_balanced_{args.variant}.mat"))

    header = (f"{'model':<28} {'complexes':>9} {'giant':>6} {'ACR K/P':>12} {'ACRR K/P':>14} "
              f"{'sec':>7}  agreement")
    print(header)
    print("-" * len(header))

    for name in models:
        mpath = os.path.join(model_dir, f"{name}_pre_balanced_{args.variant}.mat")
        if not os.path.exists(mpath):
            continue
        try:
            with h5py.File(mpath, "r") as fh:
                key = next(k for k in fh if k.startswith("model_elementary"))
                m = fh[key]
                A = read_sparse(m["A"]).tocsr().astype(np.int64)
                Y = read_sparse(m["Y"])
                S = read_sparse(m["S"])
                mets = read_strings(fh, m["mets"])
                met_names = read_strings(fh, m["metNames"])

            # The unfiltered files are written from metNames, the free ones from mets ids. Map the
            # former back, and note that the mapping is many-to-one: duplicate display names are
            # common, so an unfiltered comparison is inherently coarser than a free one.
            name_to_id = {}
            for i, display in enumerate(met_names):
                name_to_id.setdefault(display, mets[i])

            # Free metabolites are those before the first enzyme species, as the reference cuts it.
            try:
                free_upto = mets.index("E_1")
            except ValueError:
                free_upto = len(mets)
            free = set(mets[:free_upto])

            N = to_exact(S)
            t0 = time.perf_counter()
            conc = kimonu.concordance_from_matrices(A, N)
            nullity = sd.sparse_nullspace(N).shape[1]
            kin = upstream.kinetic_modules(A, conc.modules, conc.balanced, Y=Y, nullity=nullity)
            found = kimonu.acr_acrr(Y, kin.modules, mets)
            seconds = time.perf_counter() - t0

            k_acr = {a for a in found.acr if a in free}
            k_acrr = {tuple(sorted(p)) for p in found.acrr if p[0] in free and p[1] in free}

            p_acr = read_published(os.path.join(
                single_dir, f"M_MetSingle_{name}_concordant_{args.variant}.csv"))
            p_acrr = read_published(os.path.join(
                double_dir, f"M_MetDouble_{name}_concordant_{args.variant}.csv"))
            if p_acr is None or p_acrr is None:
                print(f"{name:<28} (no published result)")
                continue
            p_acr, p_acrr = set(p_acr), {tuple(sorted(p)) for p in p_acrr if isinstance(p, tuple)}

            bits = []
            if k_acr != p_acr:
                bits.append(f"ACR +{len(k_acr - p_acr)}/-{len(p_acr - k_acr)}")
            if k_acrr != p_acrr:
                bits.append(f"ACRR +{len(k_acrr - p_acrr)}/-{len(p_acrr - k_acrr)}")
            verdict = "exact match" if not bits else "differs: " + ", ".join(bits)

            print(f"{name:<28} {A.shape[0]:>9} {kin.giant_module:>6} "
                  f"{f'{len(k_acr)}/{len(p_acr)}':>12} {f'{len(k_acrr)}/{len(p_acrr)}':>14} "
                  f"{seconds:>7.2f}  {verdict}")
        except Exception as exc:
            print(f"{name:<28} FAILED: {type(exc).__name__}: {exc}"[:150])
        sys.stdout.flush()


if __name__ == "__main__":
    main()
