"""Small networks with known answers, shared by the tests, and builders of their matrices A, Y and N and their CRN.

Reactions are written as (substrate, product) complexes, each a multiset of species joined by "+", with "0" for the
zero complex. Complexes are numbered in the order in which the reactions first use them.
"""
from fractions import Fraction

import numpy as np
from scipy import sparse

import kimonu


def matrices(reactions, species):
    """(A, Y, N, complexes) of a network written as (substrate, product) pairs."""
    complexes = []
    for s, p in reactions:
        for c in (s, p):
            if c not in complexes:
                complexes.append(c)
    A = np.zeros((len(complexes), len(reactions)), dtype=np.int64)
    for j, (s, p) in enumerate(reactions):
        A[complexes.index(s), j] = -1
        A[complexes.index(p), j] = 1
    Y = np.zeros((len(species), len(complexes)), dtype=np.int64)
    for k, c in enumerate(complexes):
        if c == "0":
            continue
        for tok in c.split("+"):
            Y[species.index(tok), k] += 1
    return sparse.csr_matrix(A), sparse.csr_matrix(Y), sparse.csr_matrix(Y @ A), complexes


def network(reactions, species, names=True):
    """(CRN, complexes) of a network written as (substrate, product) pairs. Results name the metabolites by
    ``species``, or by their index when names is False."""
    A, Y, _, complexes = matrices(reactions, species)
    return kimonu.crn_from_matrices(A, Y, metabolite_ids=species if names else None), complexes


# EnvZ-OmpR with ADP-bound EnvZ (XD), as used by Langary et al. (2025) and in Figure 1. Yp has ACR for every choice of
# rate constants, and reaching it needs the feedback of the iterated route.
ENVZ = [("XD", "X"), ("X", "XD"), ("X", "XT"), ("XT", "X"), ("XT", "Xp"),
        ("Xp+Y", "XpY"), ("XpY", "Xp+Y"), ("XpY", "X+Yp"),
        ("XT+Yp", "XTYp"), ("XTYp", "XT+Yp"), ("XTYp", "XT+Y"),
        ("XD+Yp", "XDYp"), ("XDYp", "XD+Yp"), ("XDYp", "XD+Y")]
ENVZ_SPECIES = ["XD", "X", "XT", "Xp", "Y", "XpY", "Yp", "XTYp", "XDYp"]

# uptake of a, a <-> b, b -> 0. Attached, a is fed from the zero complex at constant rate.
OPEN = [("0", "a"), ("a", "b"), ("b", "a"), ("b", "0")]

# Figure 2: X is made and degraded, dimerizes into an inactive kinase C, and active P activates C cooperatively.
# Reactions R1 to R6 in the order of the Figure 2 legend and of the rate constants K. At k = K the network has the
# steady states X = 10, (C, P) = (1, 4) and (4, 1).
KINASE = [("0", "X"), ("X", "0"), ("X+X", "C"), ("C+P+P", "P+P+P"), ("C", "0"), ("P", "0")]
KINASE_SPECIES = ["X", "C", "P"]
K = [20, 1, Fraction(1, 20), Fraction(1, 4), 1, 1]

# Figure 1, the EnvZ-OmpR network with its metabolites labeled A to I
FIG1 = [("A", "B"), ("B", "A"), ("B", "C"), ("C", "B"), ("C", "D"), ("D+E", "F"), ("F", "D+E"), ("F", "B+G"),
        ("C+G", "H"), ("H", "C+G"), ("H", "C+E"), ("A+G", "I"), ("I", "A+G"), ("I", "A+E")]

# 0 <-> A, 2A <-> 3A: three steady states x_A = 1, 2, 3 for suitable rate constants
CUBIC = [("0", "A"), ("A", "0"), ("A+A", "A+A+A"), ("A+A+A", "A+A")]

TOYS = {"envz": (ENVZ, ENVZ_SPECIES), "open": (OPEN, ["a", "b"]), "kinase": (KINASE, KINASE_SPECIES),
        "fig1": (FIG1, list("ABCDEFGHI")), "cubic": (CUBIC, ["A"])}
