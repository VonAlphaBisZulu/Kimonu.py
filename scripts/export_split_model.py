#!/usr/bin/env python
"""Convert a reference elementary-split model (.mat, v7.3) into the shared crn.json format.

Carries the real metabolite ids across, so results can be compared against the reference's
published ACR/ACRR files without a name mapping in between.

Usage:
    python scripts/export_split_model.py <upstream_repo> <model> [--variant fixed|random]
"""
import argparse
import json
import os

import h5py
import numpy as np

from validate_acr_published import read_sparse, read_strings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("upstream")
    ap.add_argument("model")
    ap.add_argument("--variant", default="fixed")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    path = os.path.join(args.upstream, "Models", "models_with_elementary_steps",
                        f"{args.model}_pre_balanced_{args.variant}.mat")
    with h5py.File(path, "r") as fh:
        m = fh[next(k for k in fh if k.startswith("model_elementary"))]
        A = read_sparse(m["A"]).tocoo()
        Y = read_sparse(m["Y"]).tocoo()
        S = read_sparse(m["S"]).tocoo()
        mets = read_strings(fh, m["mets"])
        rxns = read_strings(fh, m["rxns"])
        lb = np.asarray(m["lb"]).ravel().tolist()
        ub = np.asarray(m["ub"]).ravel().tolist()

    def triplets(M):
        return [[int(i) + 1, int(j) + 1, float(v)] for i, j, v in zip(M.row, M.col, M.data)]

    out = {
        "n_complexes": A.shape[0], "n_reactions": A.shape[1], "n_metabolites": Y.shape[0],
        "reactions": rxns, "metabolites": mets, "lb": lb, "ub": ub,
        "A": triplets(A), "Y": triplets(Y), "S": triplets(S),
    }
    dest = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "published",
                                    f"{args.model}_{args.variant}.crn.json")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "w") as fh:
        json.dump(out, fh)
    print(f"{args.model}: {A.shape[0]} complexes, {A.shape[1]} reactions, {Y.shape[0]} metabolites "
          f"-> {dest}")


if __name__ == "__main__":
    main()
