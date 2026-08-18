#!/usr/bin/env python
"""Run Kimonu on a pre-decomposed network (crn.json), rather than on a cobra model.

Used for the reference implementation's own elementary-step-split models, which already carry A
and Y. Analysing them directly means the comparison against the published numbers is made on the
paper's own input, with no preprocessing of ours in between.

Usage:
    python benchmarks/run_kimonu_crn.py <crn.json> <out.json> [--repeats N]
"""
import argparse
import json
import time
from fractions import Fraction
from math import gcd

import numpy as np
from scipy import sparse

import kimonu


def triplets_to_sparse(t, nrows, ncols, dtype=float):
    if not t:
        return sparse.csr_matrix((nrows, ncols), dtype=dtype)
    arr = np.asarray(t, dtype=float)
    return sparse.csr_matrix((arr[:, 2].astype(dtype), (arr[:, 0].astype(int) - 1,
                                                        arr[:, 1].astype(int) - 1)),
                             shape=(nrows, ncols), dtype=dtype)


def to_integer_rows(M):
    """Scale each row of M to integers, exactly.

    Row scaling by a positive factor leaves ker(M) unchanged, so this keeps the kernel computation
    in integer arithmetic without altering the answer.
    """
    M = M.tocsr()
    rows, cols, vals = [], [], []
    for i in range(M.shape[0]):
        start, end = M.indptr[i], M.indptr[i + 1]
        if start == end:
            continue
        fracs = [Fraction(str(v)) for v in M.data[start:end]]
        scale = 1
        for f in fracs:
            scale = scale * f.denominator // gcd(scale, f.denominator)
        for j, f in zip(M.indices[start:end], fracs):
            rows.append(i)
            cols.append(int(j))
            vals.append(int(f * scale))
    return sparse.csr_matrix((vals, (rows, cols)), shape=M.shape, dtype=np.int64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("crn")
    ap.add_argument("out")
    ap.add_argument("--repeats", type=int, default=2)
    args = ap.parse_args()

    with open(args.crn) as fh:
        crn = json.load(fh)

    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]
    A = triplets_to_sparse(crn["A"], nc, nr, dtype=np.int64)
    Y = triplets_to_sparse(crn["Y"], nm, nc)
    S = to_integer_rows(triplets_to_sparse(crn["S"], nm, nr))

    times = []
    for _ in range(args.repeats):
        t0 = time.perf_counter()
        result = kimonu.concordance_from_matrices(A, S)
        times.append(time.perf_counter() - t0)

    module_of = {}
    for module_id, group in enumerate(result.modules, start=1):
        for i in group:
            module_of[i] = module_id if len(group) > 1 else -1
    for i in result.balanced:
        module_of[i] = 0

    Ycsc = Y.tocsc()
    complexes = []
    for j in range(nc):
        start, end = Ycsc.indptr[j], Ycsc.indptr[j + 1]
        comp = {crn["metabolites"][int(r)]: float(v)
                for r, v in zip(Ycsc.indices[start:end], Ycsc.data[start:end])}
        complexes.append({"id": f"complex_{j + 1}", "module": module_of[j], "composition": comp})

    out = {"tool": "Kimonu.py", "seconds": min(times), "seconds_cold": times[0],
           "repeats": args.repeats, "n_reactions": nr, "n_metabolites": nm,
           "complexes": complexes}
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)

    n_nonsingleton = sum(1 for g in result.modules if len(g) > 1)
    print(f"{nc} complexes, {result.n_balanced} balanced, {n_nonsingleton} modules, "
          f"{result.n_concordant_pairs} pairs, {min(times):.2f}s")


if __name__ == "__main__":
    main()
