#!/usr/bin/env python
"""Compare Kimonu against the reference implementation's *published* results, model by model.

Reads the elementary-step-split models the reference ships (HDF5 .mat, via h5py) and the committed
concordance results beside them (v7 .mat, via scipy), runs Kimonu on the same decomposition, and
compares the partitions.

Complexes are identified by index: both sides come from the same model file, so there is nothing to
match up and no naming scheme to trust.

Pair counts are computed combinatorially rather than by enumerating pairs — the giant module of a
genome-scale network has thousands of members. The refinement check is what establishes the
direction of any disagreement: if every published module sits inside a single Kimonu module, the
published partition is strictly coarser-grained nowhere and its concordant pairs are a subset.

Usage:
    python benchmarks/compare_to_published.py <upstream_repo> [model ...] [--variant fixed|random]
"""
import argparse
import os
import sys
import time
from fractions import Fraction
from math import gcd

import h5py
import numpy as np
import scipy.io as sio
from scipy import sparse

import kimonu


def read_sparse(group):
    """Rebuild a MATLAB sparse matrix stored in HDF5 (CSC: data, ir, jc)."""
    nrows = int(group.attrs["MATLAB_sparse"])
    data = np.asarray(group["data"]).ravel()
    ir = np.asarray(group["ir"]).ravel().astype(np.int64)
    jc = np.asarray(group["jc"]).ravel().astype(np.int64)
    ncols = len(jc) - 1
    return sparse.csc_matrix((data, ir, jc), shape=(nrows, ncols))



def components(adjacency, n, exclude):
    """Connected components of a symmetric boolean sparse matrix, skipping excluded vertices."""
    label = np.zeros(n, dtype=np.int64)
    adjacency = adjacency.tocsr()
    current = 0
    for start in range(n):
        if label[start] != 0 or start in exclude:
            continue
        current += 1
        stack = [start]
        label[start] = current
        while stack:
            v = stack.pop()
            for w in adjacency.indices[adjacency.indptr[v]:adjacency.indptr[v + 1]]:
                if label[w] == 0 and w not in exclude:
                    label[w] = current
                    stack.append(int(w))
    return label


def pair_count(labels, exclude):
    sizes = {}
    for i, l in enumerate(labels):
        if l == 0 or i in exclude:
            continue
        sizes[l] = sizes.get(l, 0) + 1
    return sum(s * (s - 1) // 2 for s in sizes.values()), sum(1 for s in sizes.values() if s > 1), \
        max(sizes.values()) if sizes else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("upstream")
    ap.add_argument("models", nargs="*")
    ap.add_argument("--variant", default="fixed", choices=["fixed", "random"])
    args = ap.parse_args()

    model_dir = os.path.join(args.upstream, "Models", "models_with_elementary_steps")
    result_dir = os.path.join(args.upstream, "Results", "concordant")
    models = args.models
    if not models:
        models = sorted(f[: -len(f"_pre_balanced_{args.variant}.mat")]
                        for f in os.listdir(model_dir)
                        if f.endswith(f"_pre_balanced_{args.variant}.mat"))

    print(f"{'model':<32} {'complexes':>9} {'balanced K/P':>14} {'modules K/P':>13} "
          f"{'pairs K/P':>17} {'giant K/P':>13} {'K sec':>7}  verdict")
    for name in models:
        mpath = os.path.join(model_dir, f"{name}_pre_balanced_{args.variant}.mat")
        rpath = os.path.join(result_dir, f"{name}_concordant_{args.variant}.mat")
        if not (os.path.exists(mpath) and os.path.exists(rpath)):
            continue

        try:
            with h5py.File(mpath, "r") as fh:
                key = next(k for k in fh if k.startswith("model_elementary"))
                m = fh[key]
                A = read_sparse(m["A"]).tocsr().astype(np.int64)
                Y = read_sparse(m["Y"])
            net = kimonu.crn_from_matrices(A, Y)          # N = Y·A, the S of the model
            nc = A.shape[0]

            t0 = time.perf_counter()
            res = kimonu.concordance(net)
            secs = time.perf_counter() - t0

            k_label = np.zeros(nc, dtype=np.int64)
            for idx, group in enumerate(res.modules, start=1):
                for i in group:
                    k_label[i] = idx
            k_balanced = set(res.balanced)
            k_pairs, k_modules, k_giant = pair_count(k_label, k_balanced)

            pub = sio.loadmat(rpath)
            p_balanced = set(int(x) - 1 for x in pub["B"].ravel())
            CC = sparse.csr_matrix(pub["CC"])
            E = (CC == 2).astype(bool) + (CC == -1).astype(bool)
            E = (E + E.T).tocsr()
            p_label = components(E, nc, p_balanced)
            p_pairs, p_modules, p_giant = pair_count(p_label, p_balanced)

            if k_balanced == p_balanced:
                bal = "same"
            else:
                bal = f"differ +{len(k_balanced - p_balanced)}/-{len(p_balanced - k_balanced)}"

            refines = True
            for l in set(p_label):
                if l == 0:
                    continue
                members = [i for i in range(nc) if p_label[i] == l and i not in p_balanced]
                if len(members) < 2:
                    continue
                if len({k_label[i] for i in members}) != 1 or any(i in k_balanced for i in members):
                    refines = False
                    break
            verdict = "published refines Kimonu (strict subset of pairs)" if refines \
                else "PARTITIONS CROSS"
            if k_pairs == p_pairs and refines and k_modules == p_modules:
                verdict = "identical"

            print(f"{name:<32} {nc:>9} {f'{len(k_balanced)}/{len(p_balanced)}':>14} "
                  f"{f'{k_modules}/{p_modules}':>13} {f'{k_pairs}/{p_pairs}':>17} "
                  f"{f'{k_giant}/{p_giant}':>13} {secs:>7.2f}  {bal}; {verdict}")
        except Exception as exc:              # keep going; report which model and why
            print(f"{name:<32} {'-':>9} {'-':>14} {'-':>13} {'-':>17} {'-':>13} {'-':>7}  "
                  f"FAILED: {type(exc).__name__}: {exc}"[:200])
        sys.stdout.flush()


if __name__ == "__main__":
    main()
