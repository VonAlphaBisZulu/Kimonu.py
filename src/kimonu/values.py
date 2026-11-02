"""What the fixed log-combinations are fixed at: ACR values, ACRR ratios and module ψ ratios at given k.

The routes decide which combinations vᵀ ln x take the same value at every positive steady state; this
module computes that value for given rate constants k, exactly. It reruns the iterated route of
kinetic.py at k itself instead of at a generic flux:

1. Linear: every positive steady state x for k has its complex monomials ψ = ψ(x) in the kernel of

       W_k = N diag(k) Sᵀ,

   the map of kinetic.py with k in place of r* (steady_state_map). When rows c and d of an exact basis of
   ker(W_k) are proportional, row c = λ row d, every kernel vector has ψ_c = λ ψ_d, so the ratio of the
   two monomials is λ, a rational number for rational k. No steady state has to be found.
2. Log: within a coupled set, ln ψ_c - ln ψ_d = (y_c - y_d)ᵀ ln x = ln λ. A fixed combination v found
   by the readout lies in the span of these composition differences, v = Σ_t a_t (y_c_t - y_d_t), so

       vᵀ ln x = Σ_t a_t ln λ_t,   e^(vᵀ ln x) = Π_t λ_t^a_t.

   The coefficients a_t are rational, so the value is a rational number times a root of one (an ACR
   metabolite whose value is fixed through 2X is a square root: 0 -> X, 2X -> 0 gives
   x = (k_1 / 2 k_2)^(1/2)). :class:`Value` holds it in that form, exactly.
3. Feedback: a merge found in log space is the linear fact ψ_c = μ ψ_d with μ the value of
   y_c - y_d; it goes into W_k as a row, as u_c = u_d does in kinetic.py, and the steps repeat until
   nothing merges.

What this proves holds at the given k itself, not only for generic k, and for generic k it is what
kinetic_modules finds. The attached zero complex has ψ_0 = 1, which is how absolute values (ACR) arise:
a complex coupled to the zero complex has the value of its monomial fixed, not only its ratio. Detached,
only ratios that hold whatever the uptake rates are reported.

It also tells when k has no positive steady state at all: a complex whose row of ker(W_k) is zero, a
negative ratio, or two derivations of one combination that disagree. A ValueError says which.

The reverse direction needs no computation. For an interior steady-state flux r and any x > 0, the rate
constants k_j = r_j / ψ_s(j)(x) have x as a positive steady state with flux r (rate_constants). Each
route picks one such r*, the ``flux`` of its result, at x* = 1, so its relations are found at k = r*
with every fixed combination equal to zero there; choosing x instead moves every ACR value to x_i,
every ACRR ratio to x_i / x_j and every module ratio to ψ_c(x) / ψ_d(x). These are the rate constants
to take to a kinetic model when given flux and concentrations should be a steady state, with the values
of the relations prescribed.
"""
from collections import namedtuple
from fractions import Fraction
from math import floor, lcm, log, exp

from kimonu.acr import free_directions, readout
from kimonu.io import substrate_complexes, zero_complexes
from kimonu.kinetic import steady_state_map
from kimonu.linalg import exact, null_space_rows

from scipy import sparse

Values = namedtuple("Values", ["acr", "acrr", "modules", "rounds"])


class Value:
    """An exact positive number r · b^(1/q): a rational r times the q-th root of a rational b.

    ``rational`` is the number as a Fraction when it is rational (q = 1), else None; float() evaluates
    it. The form is canonical: b is not a perfect p-th power for any p dividing q, so two Values are equal
    exactly when their numbers are.
    """
    __slots__ = ("r", "b", "q")

    def __init__(self, logs):
        """From {base: exponent}, the number Π base^exponent with rational bases and exponents."""
        r, fractional = Fraction(1), {}
        for base, e in logs.items():
            if e == 0 or base == 1:
                continue
            if e < 0:                                 # b^(-1/2) as (1/b)^(1/2)
                base, e = 1 / base, -e
            n = floor(e)
            r *= base ** n
            if e != n:
                fractional[base] = e - n
        q = lcm(*(e.denominator for e in fractional.values())) if fractional else 1
        b = Fraction(1)
        for base, e in fractional.items():
            b *= base ** int(e * q)
        for p in range(2, q + 1):                     # take every exact root out of b
            while q % p == 0:
                root = _root(b, p)
                if root is None:
                    break
                b, q = root, q // p
        if q == 1:
            r, b = r * b, Fraction(1)
        self.r, self.b, self.q = r, b, q

    @property
    def rational(self):
        return self.r if self.q == 1 else None

    def __float__(self):
        return exp(_ln(self.r) + _ln(self.b) / self.q)

    def __eq__(self, other):
        if isinstance(other, Value):
            return (self.r, self.b, self.q) == (other.r, other.b, other.q)
        return self.q == 1 and self.r == other

    def __hash__(self):
        return hash((self.r, self.b, self.q))

    def __repr__(self):
        if self.q == 1:
            return str(self.r)
        root = f"({self.b})^(1/{self.q})"
        return root if self.r == 1 else f"{self.r}·{root}"


def _ln(x):
    """ln of a positive Fraction, also when numerator or denominator are too large for a float"""
    return log(x.numerator) - log(x.denominator)


def _root(x, p):
    """The exact p-th root of a positive Fraction, or None when it is irrational"""
    out = []
    for n in (x.numerator, x.denominator):
        r = n
        if n > 1:                                     # integer Newton from above
            r = 1 << -(-n.bit_length() // p)
            while True:
                s = ((p - 1) * r + n // r ** (p - 1)) // p
                if s >= r:
                    break
                r = s
        if r ** p != n:
            return None
        out.append(r)
    return Fraction(*out)


def _add(a, b, c=1):
    """a + c·b for log values {base: exponent}"""
    out = dict(a)
    for base, e in b.items():
        out[base] = out.get(base, 0) + c * e
        if out[base] == 0:
            del out[base]
    return out


class _Span:
    """The span of composition differences with known log values, in echelon form.

    Each column is a vector v (a {metabolite: Fraction} dict) with the value of vᵀ ln x as a log value
    {base: exponent}. Columns are kept reduced against the earlier pivots, so a vector is reduced by
    subtracting the pivots in the order they were added. A new column that reduces to zero repeats a
    combination already known; its two values must agree, or no positive steady state exists.
    """

    def __init__(self):
        self.pivots = []                              # [(metabolite, vector, log value)], vector[metabolite] = 1

    def _reduce(self, v, value):
        v, value = dict(v), dict(value)
        for met, p, pv in self.pivots:
            c = v.get(met)
            if c:
                for i, x in p.items():
                    v[i] = v.get(i, 0) - c * x
                    if v[i] == 0:
                        del v[i]
                value = _add(value, pv, -c)
        return v, value

    def add(self, v, value):
        v, value = self._reduce(v, value)
        if not v:
            self._repeated(value)
            return
        met = min(v)
        c = v[met]
        self.pivots.append((met, {i: x / c for i, x in v.items()}, {b: e / c for b, e in value.items()}))

    def _repeated(self, value):
        """value: the difference of two derivations of one combination, which must be zero"""
        if Value(value) != 1:
            raise ValueError("no positive steady state for these rate constants: two derivations of one "
                             "fixed combination give different values")

    def logs(self, v):
        """vᵀ ln x as a log value {base: exponent}, or None when v is not in the span"""
        rest, value = self._reduce(v, {})
        return None if rest else {b: -e for b, e in value.items()}

    def value(self, v):
        """e^(vᵀ ln x) as a Value, for v in the span"""
        logs = self.logs(v)
        assert logs is not None, "not a fixed combination"
        return Value(logs)


def _ratio_classes(rows, n):
    """Classes of proportional rows among the first n, with the ratio of every member to the first one.

    A zero row is a complex that every kernel vector, and so every steady state, sets to zero, and a
    negative ratio one that no positive vector satisfies: either way there is no positive steady state.
    """
    classes = {}
    for c, row in enumerate(rows[:n]):
        if not row:
            raise ValueError(f"no positive steady state for these rate constants: complex {c} is forced to zero")
        nz = sorted(row.items())
        classes.setdefault(tuple((col, v / nz[0][1]) for col, v in nz), []).append((c, nz[0][1]))
    out = []
    for members in classes.values():
        (first, f0), *_ = members
        ratios = {c: f / f0 for c, f in members}      # ψ_c / ψ_first
        if any(x < 0 for x in ratios.values()):
            raise ValueError(f"no positive steady state for these rate constants: complexes {sorted(ratios)} "
                             "would need monomials of opposite sign")
        out.append(ratios)
    return out


def _column(Ycsc, c):
    return {int(i): Fraction(str(v)) for i, v in zip(Ycsc.indices[Ycsc.indptr[c]:Ycsc.indptr[c + 1]],
                                                    Ycsc.data[Ycsc.indptr[c]:Ycsc.indptr[c + 1]])}


def _difference(Ycsc, c, d):
    """y_c - y_d"""
    out = _column(Ycsc, c)
    for i, v in _column(Ycsc, d).items():
        out[i] = out.get(i, 0) - v
        if out[i] == 0:
            del out[i]
    return out


def relation_values(A, Y, N, k, detach_uptakes=False, metabolite_ids=None):
    """Exact values of the ACR, ACRR and kinetic-module relations at the rate constants k.

    Args:
        A, Y, N: the network, as for :func:`kimonu.kinetic_modules`.
        k: one positive rate constant per reaction, in reaction order. Ints and Fractions are exact;
            a float is read as the decimal it prints as.
        detach_uptakes: as for :func:`kimonu.kinetic_modules`.
        metabolite_ids: optional names for ``acr`` and ``acrr``; indices are used when omitted.

    Returns a :class:`Values` with every value a :class:`Value`:
        ``acr``: {metabolite: x_i} for the metabolites with ACR;
        ``acrr``: one (reference, {metabolite: x_j / x_reference}) per ACRR group, the reference being
        the group's first metabolite;
        ``modules``: one (reference, {complex: ψ_c / ψ_reference}) per kinetic module with two or more
        complexes, the reference being its first complex;
        ``rounds``: as in :class:`kimonu.KineticResult`.

    Raises ValueError when k admits no positive steady state, as far as the relations show it.
    """
    n, q = A.shape[0], A.shape[1]
    k = [x if isinstance(x, (int, Fraction)) else Fraction(str(x)) for x in k]
    if len(k) != q or any(x <= 0 for x in k):
        raise ValueError(f"k needs one positive rate constant for each of the {q} reactions")
    skip = zero_complexes(Y) if detach_uptakes else set()
    Ycsc = sparse.csc_matrix(Y)
    substrates = set(substrate_complexes(A).values()) - skip
    equal, rounds = [], 0                             # feedback rows ψ_c = μ ψ_d as (c, d, μ)
    while True:
        rounds += 1
        # 1) Linear: proportional rows of ker(W_k), with the extra rows of the feedback
        W, (m, n_columns) = steady_state_map(A, Y, N, k, detach_uptakes)
        for t, (c, d, mu) in enumerate(equal):
            W[(m + t, c)], W[(m + t, d)] = Fraction(1), -mu
        rows = null_space_rows(exact(((i, j, v) for (i, j), v in W.items()), (m + len(equal), n_columns)))
        classes = [{c: x for c, x in g.items() if c not in skip} for g in _ratio_classes(rows, n)]
        classes = [g for g in classes if g]
        # 2) Log: the composition differences within each class, with ln of their ratios, span the fixed
        #    combinations; the readout of kinetic.py merges the classes they connect
        span = _Span()
        for g in classes:
            first = min(g)
            for c in sorted(g):
                if c != first:                        # y_first - y_c, at ln ψ_first / ψ_c
                    span.add(_difference(Ycsc, first, c), {g[first] / g[c]: Fraction(1)})
        found = readout(Y, free_directions(Y, [set(g) for g in classes]), skip)
        # 3) Feedback: the merges that no class already holds, with their values. A complex that is no
        #    reaction's substrate has a zero column in W_k, so a row for it would fix only its own ψ.
        class_of = {c: t for t, g in enumerate(classes) for c in g}
        new = []
        for module in found.modules:
            inner = sorted(module & substrates)
            first = inner[0] if inner else None
            for c in inner[1:]:
                if class_of[c] != class_of[first]:
                    mu = span.value(_difference(Ycsc, first, c))                # ψ_first / ψ_c
                    if mu.rational is None:
                        raise NotImplementedError(f"the feedback needs the irrational ratio {mu} "
                                                  f"of complexes {first} and {c}")
                    new.append((first, c, mu.rational))
        if not new:
            break
        equal += new

    names = metabolite_ids if metabolite_ids is not None else list(range(Y.shape[0]))
    acr = {names[i]: span.value({i: Fraction(1)}) for i in found.acr}
    acrr = [(names[g[0]], {names[j]: span.value({j: Fraction(1), g[0]: Fraction(-1)}) for j in g[1:]})
            for g in found.acrr]
    modules = []
    for module in sorted(found.modules, key=min):
        if len(module) > 1:
            first = min(module)
            modules.append((first, {c: span.value(_difference(Ycsc, c, first)) for c in sorted(module - {first})}))
    return Values(acr=acr, acrr=acrr, modules=modules, rounds=rounds)


def rate_constants(A, Y, flux, x=None):
    """Rate constants k_j = r_j / ψ_s(j)(x) for which x is a positive steady state with flux r.

    Args:
        A, Y: the network.
        flux: an interior steady-state flux r (every entry positive, N r = 0), e.g. the ``flux`` of a
            :class:`kimonu.KineticResult` or :class:`kimonu.TangentResult`.
        x: positive concentrations, one per metabolite; all ones when omitted, which gives k = r.

    At these k every relation the routes find takes its value at x: x_i for an ACR metabolite, x_j / x_i
    for an ACRR pair, ψ_d(x) / ψ_c(x) within a module. Exact when x is rational and the substrate
    compositions are integers (or x_i = 1); a fractional power of x_i ≠ 1 is irrational and given as
    a float.
    """
    if x is None:
        return list(flux)
    Ycsc = sparse.csc_matrix(Y)
    x = [v if isinstance(v, (int, Fraction)) else Fraction(str(v)) for v in x]
    k = []
    for j, c in sorted(substrate_complexes(A).items()):
        psi, inexact = Fraction(1), 1.0                # ψ_c(x) = psi · inexact
        for i, y in _column(Ycsc, c).items():
            if x[i] == 1:
                continue
            if y.denominator == 1:
                psi *= Fraction(x[i]) ** int(y)
            else:
                inexact *= float(x[i]) ** float(y)
        k.append(flux[j] / psi if inexact == 1.0 else float(flux[j] / psi) / inexact)
    return k
