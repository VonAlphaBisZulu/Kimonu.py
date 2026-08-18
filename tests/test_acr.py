"""ACR/ACRR tests against the analytically known EnvZ-OmpR answer.

EnvZ-OmpR is the canonical Shinar-Feinberg example: the phosphorylated response regulator Yp has
absolute concentration robustness and nothing else does. The test supplies the kinetic modules
directly, so it pins the ACR layer independently of how the modules are derived.
"""
import numpy as np
import pytest
from scipy import sparse

import kimonu

# EnvZ-OmpR as COCOA constructs it: 9 species, 13 complexes.
METABOLITES = ["X", "XD", "XDYp", "XT", "XTYp", "Xp", "XpY", "Y", "Yp"]

# complex index -> composition, in the model's complex order
COMPLEXES = [
    {"XD": 1},                 # 0
    {"X": 1},                  # 1
    {"XTYp": 1},               # 2
    {"XT": 1, "Yp": 1},        # 3
    {"XT": 1, "Y": 1},         # 4
    {"XD": 1, "Yp": 1},        # 5
    {"XDYp": 1},               # 6
    {"XD": 1, "Y": 1},         # 7
    {"XT": 1},                 # 8
    {"Xp": 1},                 # 9
    {"Xp": 1, "Y": 1},         # 10
    {"XpY": 1},                # 11
    {"X": 1, "Yp": 1},         # 12
]

# The kinetic modules of this network; the large one is the giant module.
KINETIC_MODULES = [{0, 1, 2, 3, 5, 6, 8, 10, 11}, {9}, {4}, {7}, {12}]


def _Y():
    index = {m: i for i, m in enumerate(METABOLITES)}
    rows, cols, vals = [], [], []
    for j, comp in enumerate(COMPLEXES):
        for met, coeff in comp.items():
            rows.append(index[met])
            cols.append(j)
            vals.append(coeff)
    return sparse.csr_matrix((vals, (rows, cols)),
                            shape=(len(METABOLITES), len(COMPLEXES)), dtype=np.int64)


def test_acr_is_exactly_yp():
    """Yp has ACR and no other species does."""
    result = kimonu.acr_acrr(_Y(), KINETIC_MODULES, METABOLITES)
    assert result.acr == ["Yp"]


def test_acrr_pairs():
    """The ratio-robust pairs are the ones the exhaustive path reports."""
    result = kimonu.acr_acrr(_Y(), KINETIC_MODULES, METABOLITES)
    assert len(result.acrr) == 15
    # ACR species are excluded from the pair list: a constant ratio between two constants is not
    # an independent finding.
    assert all("Yp" not in pair for pair in result.acrr)


def test_no_modules_means_no_acr():
    """With no coupled complexes there is nothing to span, so nothing is robust."""
    result = kimonu.acr_acrr(_Y(), [{i} for i in range(len(COMPLEXES))], METABOLITES)
    assert result.acr == []
    assert result.acrr == []
    assert result.n_columns == 0


# --- end to end: concordance -> kinetic modules -> ACR, from the network alone -----------------

# EnvZ-OmpR incidence: (complex, reaction, +-1), in the model's own order.
REACTIONS = [
    (0, 1), (1, 0), (1, 8), (8, 1), (8, 9), (10, 11), (11, 10), (11, 12),
    (3, 2), (2, 3), (2, 4), (5, 6), (6, 5), (6, 7),
]


def _A():
    rows, cols, vals = [], [], []
    for j, (substrate, product) in enumerate(REACTIONS):
        rows += [substrate, product]
        cols += [j, j]
        vals += [-1, 1]
    return sparse.csr_matrix((vals, (rows, cols)),
                             shape=(len(COMPLEXES), len(REACTIONS)), dtype=np.int64)


def test_kinetic_modules_reproduce_the_giant_module():
    """The nine non-terminal complexes form one kinetic module.

    Reaching it needs the deficiency rule: coupling alone cannot merge the set holding XD with the
    one holding XD + Yp, because that merge is equivalent to knowing Yp has ACR.
    """
    A, Y = _A(), _Y()
    concordance = [{0, 2, 6, 8, 11}, {1, 9, 10, 12}, {3, 4}, {5, 7}]
    balanced = {0, 2, 6, 8, 11}
    result = kimonu.kinetic_modules(A, concordance, balanced, Y=Y)

    assert result.deficiency == 2
    assert result.n_modules == 1
    assert result.giant_module == 9
    giant = next(m for m in result.modules if len(m) > 1)
    assert giant == {0, 1, 2, 3, 5, 6, 8, 10, 11}


def test_end_to_end_acr_is_yp():
    """From the incidence matrix alone, the pipeline finds Yp and only Yp."""
    A, Y = _A(), _Y()
    concordance = [{0, 2, 6, 8, 11}, {1, 9, 10, 12}, {3, 4}, {5, 7}]
    modules = kimonu.kinetic_modules(A, concordance, {0, 2, 6, 8, 11}, Y=Y)
    found = kimonu.acr_acrr(Y, modules.modules, METABOLITES)
    assert found.acr == ["Yp"]
    assert len(found.acrr) == 15


def test_without_Y_the_coarser_merges_are_skipped():
    """Omitting Y reproduces overlap-only merging, which is what the reference does."""
    A, Y = _A(), _Y()
    concordance = [{0, 2, 6, 8, 11}, {1, 9, 10, 12}, {3, 4}, {5, 7}]
    result = kimonu.kinetic_modules(A, concordance, {0, 2, 6, 8, 11})
    assert result.giant_module < 9
