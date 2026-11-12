"""Equivalence of two Kimonu commits on the 33 networks: every result of every route, compared exactly.

    PYTHONPATH=<src of one commit> python kinetic_variants/refactor_equivalence.py dump <model> <out.pkl>
    python kinetic_variants/refactor_equivalence.py compare <dir A> <dir B> [models.txt]

dump builds the unblocked network as uptake_modes.py does (kimonu.crn_from_matrices on unblocked/<model>_fixed.crn.json,
metabolites by index) and runs concordance, then with detached and with attached uptakes (seed 0) the iterated route,
the tangent route, certify, and kinetic flux coupling by both routes. Every field of every result is pickled: modules,
ACR, ACRR, the counts, K and M, B, the dimension, the flux r*, the rounds, the certificate's open relations and the
coupling groups with their fixed rates. The restructuring of the package after commit ed9240e renamed some fields
(giant_module to largest_module, classes to groups, class_of to group_of, n_classes to n_groups, largest_class to
largest_group), and dump stores them under the earlier names, so that two dumps compare field by field. Fields that
only the restructured package has (RouteResult.acrr_groups, ConcordanceResult.activity_ratio) are checked against
the fields they derive from instead: the pairs within the ACRR groups must be the ACRR pairs, and every row of M must
be its activity ratio times the first row of its module.

compare loads <dir A>/<model>.pkl and <dir B>/<model>.pkl for every model (all 33 of models33.txt by default), prints
one line per model and every field that differs, and exits with status 1 if any does. The run behind the archive
README's statement is refactor_equivalence.sbatch.
"""
import json
import os
import pickle
import sys
import time

import numpy as np
from scipy import sparse

B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RENAMED = {"largest_module": "giant_module", "groups": "classes", "group_of": "class_of", "n_groups": "n_classes",
           "largest_group": "largest_class"}
NEW_ONLY = {"acrr_groups", "activity_ratio"}


def network(model):
    """The unblocked network of ``model`` as uptake_modes.py builds it."""
    import kimonu
    crn = json.load(open(f"{B}/unblocked/{model}_fixed.crn.json"))
    nc, nr, nm = crn["n_complexes"], crn["n_reactions"], crn["n_metabolites"]

    def matrix(t, rows, cols, dtype):
        arr = np.asarray(t, dtype=float)
        return sparse.csr_matrix((arr[:, 2].astype(dtype), (arr[:, 0].astype(int) - 1, arr[:, 1].astype(int) - 1)),
                                 shape=(rows, cols), dtype=dtype)
    return kimonu.crn_from_matrices(matrix(crn["A"], nc, nr, np.int64), matrix(crn["Y"], nm, nc, float))


def fields(result, checks):
    """{field: value} of a result namedtuple under the earlier field names, numpy arrays as lists, nested results as
    dicts. Fields only the restructured package has are checked (``checks``) and left out."""
    out = {}
    for name, value in result._asdict().items():
        if name in NEW_ONLY:
            checks.append(_check(name, value, result))
            continue
        if isinstance(value, np.ndarray):
            value = value.tolist()
        elif hasattr(value, "_asdict"):
            value = fields(value, checks)
        out[RENAMED.get(name, name)] = value
    return out


def _check(name, value, result):
    """(field, passed) for a field that only the restructured package has."""
    if name == "acrr_groups":
        pairs = sorted((a, b) for g in value for k, a in enumerate(g) for b in g[k + 1:])
        return name, pairs == sorted(result.acrr) and all(g == sorted(g) for g in value)
    first = {c: m[0] for m in result.modules for c in m}
    return name, set(value) == set(first) and all(
        list(result.M[c]) == [(j, value[c] * v) for j, v in result.M[first[c]]] for c in value)


def dump(model, path):
    import kimonu
    net = network(model)
    checks, out, seconds = [], {"model": model, "kimonu": kimonu.__file__}, {}
    t = time.time()
    out["concordance"] = fields(kimonu.concordance(net), checks)
    seconds["concordance"] = time.time() - t
    for mode, detach in (("detached", True), ("attached", False)):
        for key, call in (("iterated", lambda: kimonu.iterated_route(net, detach_uptakes=detach)),
                          ("tangent", lambda: kimonu.tangent_route(net, detach_uptakes=detach)),
                          ("certify", lambda: kimonu.certify(net, detach_uptakes=detach)),
                          ("coupling_iterated", lambda: kimonu.kinetic_flux_coupling(net, "iterated", detach)),
                          ("coupling_tangent", lambda: kimonu.kinetic_flux_coupling(net, "tangent", detach))):
            t = time.time()
            out[f"{mode}/{key}"] = fields(call(), checks)
            seconds[f"{mode}/{key}"] = time.time() - t
    out["checks"] = checks
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(out, f)
    print(model, kimonu.__file__, f"{sum(seconds.values()):.0f} s", {k: round(v) for k, v in seconds.items()},
          "checks", "pass" if all(ok for _, ok in checks) else [c for c in checks if not c[1]], flush=True)


def differences(a, b, where=""):
    """Paths of the fields in which a and b differ."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = [f"{where}/{k}: only in one" for k in sorted(set(a) ^ set(b), key=str)]
        for k in a.keys() & b.keys():
            out += differences(a[k], b[k], f"{where}/{k}")
        return out
    return [] if a == b and type(a) is type(b) or a == b and {type(a), type(b)} <= {list, tuple} else [where]


def compare(dir_a, dir_b, models_file):
    models = open(models_file).read().split()
    failed = 0
    for model in models:
        a, b = (pickle.load(open(f"{d}/{model}.pkl", "rb")) for d in (dir_a, dir_b))
        assert a["kimonu"] != b["kimonu"], f"{model}: both dumps come from {a['kimonu']}"
        bad_checks = [c for c in a["checks"] + b["checks"] if not c[1]]
        diff = differences({k: v for k, v in a.items() if k not in ("kimonu", "checks")},
                           {k: v for k, v in b.items() if k not in ("kimonu", "checks")})
        n_fields = sum(len(v) for k, v in a.items() if isinstance(v, dict))
        status = "identical" if not diff and not bad_checks else "DIFFERENT"
        failed += status != "identical"
        print(f"{model:<30} {status}  {n_fields} fields in {sum(isinstance(v, dict) for v in a.values())} results, "
              f"{len(a['checks']) + len(b['checks'])} derived-field checks", flush=True)
        for d in diff + [f"check {c[0]} failed" for c in bad_checks]:
            print("   ", d)
    print(f"{len(models) - failed} of {len(models)} networks identical")
    return failed


if __name__ == "__main__":
    if sys.argv[1] == "dump":
        dump(sys.argv[2], sys.argv[3])
    else:
        sys.exit(1 if compare(sys.argv[2], sys.argv[3],
                              sys.argv[4] if len(sys.argv) > 4 else f"{B}/kinetic_variants/models33.txt") else 0)
