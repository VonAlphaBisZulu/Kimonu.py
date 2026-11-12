"""Values of the relations, beyond the main text: exact values at given k, k for given values, formulas in k.

The routes decide which combinations vᵀ ln x take the same value at every positive steady state. This module
computes that value for given rate constants k, exactly (:func:`relation_values`). It reruns the iterated route
(:mod:`kimonu.iterated`) at k itself instead of at a generic flux.

1. Mass balance. Every positive steady state x for k has its complex concentrations ψ = ψ(x) in the nullspace of
   W_k = N diag(k) Sᵀ, the W of Eq. 3 with k in place of r* (:func:`kimonu.steady_state_map`). When rows c and d of
   an exact basis of ker(W_k) are proportional, row c = λ row d (:func:`kimonu.readout.row_readout`), every vector of
   the nullspace has ψ_c = λ ψ_d, so the ratio of the two complex concentrations is λ, a rational number for rational
   k. No steady state has to be found.
2. Composition. Within a coupled set, ln ψ_c - ln ψ_d = (y_c - y_d)ᵀ ln x = ln λ. A fixed combination v found by
   the readout lies in the span of these composition differences, v = Σ_t a_t (y_c_t - y_d_t), so
   vᵀ ln x = Σ_t a_t ln λ_t and e^(vᵀ ln x) = Π_t λ_t^a_t. The coefficients a_t are rational, so the value is a
   rational number times a root of one. An ACR metabolite whose value is fixed through 2X is a square root, as
   0 -> X, 2X -> 0 gives x = (k_1 / 2 k_2)^(1/2). :class:`Value` holds it in that form, exactly.
3. Feedback. A merge found on the log scale is the linear fact ψ_c = μ ψ_d with μ the value of y_c - y_d. It goes
   into W_k as a row, as u_c = u_d does in the iterated route, and the steps repeat until nothing merges.

What this proves holds at the given k itself, not only for generic k, and for generic k it is what the iterated
route finds. The attached zero complex has ψ_∅ = 1, which is how absolute values (ACR) arise, since a complex
coupled to the zero complex has its concentration fixed, not only its ratio. Detached, only ratios that hold
whatever the uptake rates are reported. The computation also tells when k has no positive steady state at all, a
complex whose row of ker(W_k) is zero, a negative ratio, or two derivations of one combination that disagree. A
ValueError says which.

The reverse direction needs no computation. For an interior steady-state flux r and any x > 0, the rate constants
k_j = r_j / ψ_σ(j)(x) have x as a positive steady state with flux r (:func:`rate_constants`). Each route picks one
such r*, the ``flux`` of its result, at x* = 1, so its relations are found at k = r* with every fixed combination
equal to zero there. Choosing x instead moves every ACR value to x_i, every ACRR ratio to x_j / x_i and every module
ratio to ψ_d(x) / ψ_c(x). These are the rate constants to take to a kinetic model when given flux and
concentrations should be a steady state, with the values of the relations prescribed.

Formulas in k (:func:`binomial_formulas`). When two reactions j and l are fully coupled, v_j = ρ v_l for every
steady-state flux (proportional rows of ker N), so mass action gives k_j ψ_σ(j) = ρ k_l ψ_σ(l) at every positive
steady state, for every k, and ln ψ_σ(j) - ln ψ_σ(l) = ln ρ + ln k_l - ln k_j. These facts span composition
differences as in step 2, and every relation whose vector lies in that span has a formula c Π_j k_j^e_j, a
:class:`Monomial`. The facts need only the exact nullspace of N, as concordance does, so this runs at genome scale.
A detached uptake runs at k_j s_j with a free source s_j and gives no fact. A fact that repeats a combination
already spanned is a constraint instead, a monomial in k that equals one at every k with a positive steady state.
Two uptakes from the attached zero complex whose fluxes are fully coupled, for instance, need k_j = ρ k_l. Other
relations proven by the iterated route have values that are products of powers of rational functions of k, which
the monomials do not cover. A relation that only the tangent route reports can have a value that is no such product
at all. On the kinase of Fig. 2, X has ACR with x = (√(b² + 8ac) - b) / 4c at both steady states, with a, b and c
the rate constants of 0 -> X, X -> 0 and 2X -> C, the root of a quadratic.
"""
from collections import namedtuple
from fractions import Fraction
from math import floor, lcm, log, exp

from kimonu.linalg import fraction, nullspace
from kimonu.network import as_crn, substrate_complexes, zero_complexes
from kimonu.readout import free_directions, readout, row_readout
from kimonu.steady_state import interior_flux, steady_state_map

Values = namedtuple("Values", ["acr", "acrr", "modules", "rounds"])
Values.__doc__ = """The exact values of the relations at given rate constants, each a :class:`Value`.

acr      {metabolite: x_i} for the metabolites with ACR
acrr     one (reference, {metabolite: x_j / x_reference}) per ACRR group, the reference its first metabolite
modules  one (reference, {complex: ψ_c / ψ_reference}) per kinetic module of two or more complexes, the reference
         its first complex
rounds   rounds of mass balance and composition, as in :class:`kimonu.RouteResult`
"""

Formulas = namedtuple("Formulas", ["acr", "acrr", "modules", "constraints"])
Formulas.__doc__ = """Monomial formulas in k for the relations, in the layout of :class:`Values`.

acr, acrr, modules  a :class:`Monomial` for every relation that flux coupling fixes, None for the others
constraints         monomials that equal one at every k with a positive steady state
"""


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


class Monomial:
    """c · Π_j k_j^e_j: an exact constant c (a :class:`Value`) times a monomial in the rate constants.

    ``constant`` is c, ``exponents`` is {reaction: e_j} with rational e_j. Called with one rate constant
    per reaction, it returns its exact Value there.
    """
    __slots__ = ("constant", "exponents", "_logs")

    def __init__(self, logs):
        """From a log value whose bases are Fractions (constants) and ("k", j) for rate constant j"""
        self._logs = {b: e for b, e in logs.items() if not isinstance(b, tuple)}
        self.exponents = dict(sorted((b[1], e) for b, e in logs.items() if isinstance(b, tuple) and e))
        self.constant = Value(self._logs)

    def __call__(self, k):
        logs = dict(self._logs)
        for j, e in self.exponents.items():
            base = fraction(k[j])
            logs[base] = logs.get(base, 0) + e
        return Value(logs)

    def format(self, names=None):
        """as text, with rate constant j written k_<names[j]> (k<j> when names are omitted)"""
        parts = [] if self.constant == 1 else [repr(self.constant)]
        for j, e in self.exponents.items():
            k = f"k_{names[j]}" if names is not None else f"k{j}"
            parts.append(k if e == 1 else f"{k}^{e}" if e.denominator == 1 else f"{k}^({e})")
        return "·".join(parts) or "1"

    def __repr__(self):
        return self.format()


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


class _Facts(_Span):
    """The span of the flux-coupling facts; a repeated combination is a constraint on k"""

    def __init__(self):
        super().__init__()
        self.constraints = []

    def _repeated(self, value):
        if any(isinstance(b, tuple) for b in value):
            self.constraints.append(Monomial(value))
        else:                                         # the constants agree, since k = r* satisfies every fact
            assert Value(value) == 1


def _ratio_groups(rows, n):
    """Groups of proportional rows among the first n, with the ratio of every member to the first one.

    A zero row is a complex that every vector of the nullspace, and so every steady state, sets to zero, and a
    negative ratio one that no positive vector satisfies: either way there is no positive steady state.
    """
    zero, groups, ratio = row_readout(rows[:n])
    if zero:
        raise ValueError(f"no positive steady state for these rate constants: complex {zero[0]} is forced to zero")
    out = []
    for group in groups:
        ratios = {c: ratio[c] for c in group}         # ψ_c / ψ_first
        if any(x < 0 for x in ratios.values()):
            raise ValueError(f"no positive steady state for these rate constants: complexes {sorted(ratios)} "
                             "would need concentrations of opposite sign")
        out.append(ratios)
    return out


def _difference(crn, c, d):
    """y_c - y_d as {metabolite: value}"""
    out = dict(crn.Yt[c])
    for i, v in crn.Yt[d]:
        out[i] = out.get(i, 0) - v
        if out[i] == 0:
            del out[i]
    return out


def relation_values(network, k, detach_uptakes=True):
    """Exact values of the ACR, ACRR and kinetic-module relations at the rate constants k.

    The iterated route runs at k itself (see the module docstring), so the modules, ACR and ACRR are those that hold
    at k, and each comes with its exact value.

    Args:
        network: a :class:`kimonu.CRN` without blocked reactions, or a cobra model, which is preprocessed first
            (:func:`kimonu.as_crn`). Pass the CRN when calling several functions, so that preprocessing runs once.
        k: one positive rate constant per reaction, in reaction order. Ints and Fractions are exact;
            a float is read as the decimal it prints as.
        detach_uptakes: as for :func:`kimonu.iterated_route`.

    Returns a :class:`Values`, with every value a :class:`Value` and metabolites named by crn.metabolite_ids.
    Raises ValueError when k admits no positive steady state, as far as the relations show it.
    """
    crn = as_crn(network)
    values = _values(crn, k, detach_uptakes)
    name = crn.metabolite_ids
    return values._replace(acr={name[i]: x for i, x in values.acr.items()},
                           acrr=[(name[r], {name[j]: x for j, x in ratios.items()}) for r, ratios in values.acrr])


def _values(crn, k, detach_uptakes):
    """:func:`relation_values` with metabolites as indices"""
    n, q = crn.A.shape
    k = [fraction(x) for x in k]
    if len(k) != q or any(x <= 0 for x in k):
        raise ValueError(f"k needs one positive rate constant for each of the {q} reactions")
    skip = zero_complexes(crn) if detach_uptakes else set()
    substrates = set(substrate_complexes(crn).values()) - skip
    equal, rounds = [], 0                             # feedback rows ψ_c = μ ψ_d as (c, d, μ)
    while True:
        rounds += 1
        # 1) Linear: proportional rows of ker(W_k), with the extra rows of the feedback
        W, n_columns = steady_state_map(crn, k, detach_uptakes)
        W += [tuple(sorted(((c, 1), (d, -mu)))) for c, d, mu in equal]
        rows = nullspace(W, n_columns)
        groups = [{c: x for c, x in g.items() if c not in skip} for g in _ratio_groups(rows, n)]
        groups = [g for g in groups if g]
        # 2) Log: the composition differences within each group, with ln of their ratios, span the fixed
        #    combinations; the readout merges the groups they connect
        span = _Span()
        for g in groups:
            first = min(g)
            for c in sorted(g):
                if c != first:                        # y_first - y_c, at ln ψ_first / ψ_c
                    span.add(_difference(crn, first, c), {g[first] / g[c]: Fraction(1)})
        found = readout(crn, free_directions(crn, [set(g) for g in groups]), detach_uptakes)
        # 3) Feedback: the merges that no group already holds, with their values. A complex that is no
        #    reaction's substrate has a zero column in W_k, so a row for it would fix only its own ψ.
        group_of = {c: t for t, g in enumerate(groups) for c in g}
        new = []
        for module in found.modules:
            inner = sorted(module & substrates)
            first = inner[0] if inner else None
            for c in inner[1:]:
                if group_of[c] != group_of[first]:
                    mu = span.value(_difference(crn, first, c))                 # ψ_first / ψ_c
                    if mu.rational is None:
                        raise NotImplementedError(f"the feedback needs the irrational ratio {mu} "
                                                  f"of complexes {first} and {c}")
                    new.append((first, c, mu.rational))
        if not new:
            break
        equal += new

    acr = {i: span.value({i: Fraction(1)}) for i in found.acr}
    acrr = [(g[0], {j: span.value({j: Fraction(1), g[0]: Fraction(-1)}) for j in g[1:]}) for g in found.acrr_groups]
    modules = []
    for module in sorted(found.modules, key=min):
        if len(module) > 1:
            first = min(module)
            modules.append((first, {c: span.value(_difference(crn, c, first)) for c in sorted(module - {first})}))
    return Values(acr=acr, acrr=acrr, modules=modules, rounds=rounds)


def rate_constants(network, flux, x=None):
    """Rate constants k_j = r_j / ψ_s(j)(x) for which x is a positive steady state with flux r.

    Args:
        network: a :class:`kimonu.CRN`, or a cobra model, which is preprocessed first (:func:`kimonu.as_crn`).
        flux: an interior steady-state flux r (every entry positive, N r = 0), e.g. the ``flux`` of a
            :class:`kimonu.RouteResult`.
        x: positive concentrations, one per metabolite; all ones when omitted, which gives k = r.

    At these k every relation the routes find takes its value at x: x_i for an ACR metabolite, x_j / x_i
    for an ACRR pair, ψ_d(x) / ψ_c(x) within a module. Exact when x is rational and the substrate
    compositions are integers (or x_i = 1); a fractional power of x_i ≠ 1 is irrational and given as
    a float.
    """
    if x is None:
        return list(flux)
    crn = as_crn(network)
    x = [fraction(v) for v in x]
    k = []
    for j, c in sorted(substrate_complexes(crn).items()):
        psi, inexact = Fraction(1), 1.0                # ψ_c(x) = psi · inexact
        for i, y in crn.Yt[c]:
            if x[i] == 1:
                continue
            if y.denominator == 1:
                psi *= Fraction(x[i]) ** int(y)
            else:
                inexact *= float(x[i]) ** float(y)
        k.append(flux[j] / psi if inexact == 1.0 else float(flux[j] / psi) / inexact)
    return k


def binomial_formulas(network, detach_uptakes=True, seed=0, draw="lp", bits=64):
    """Monomial formulas in k for the relations that flux coupling fixes, and the constraints it puts on k.

    Args:
        network: a :class:`kimonu.CRN` without blocked reactions, or a cobra model, which is preprocessed first
            (:func:`kimonu.as_crn`). Pass the CRN when calling several functions, so that preprocessing runs once.
        detach_uptakes: as for :func:`kimonu.iterated_route`.
        seed, draw, bits: of the interior flux r* at which the relations are found (:func:`relation_values` at
            k = r*, :func:`kimonu.interior_flux`).

    Returns a :class:`Formulas` in the layout of :class:`Values`, with a :class:`Monomial` for every relation that
    flux coupling fixes and None for the others, and the constraints flux coupling puts on k (see the module
    docstring).
    """
    crn = as_crn(network)
    found = _values(crn, interior_flux(crn, seed, draw, bits), detach_uptakes)
    facts = _flux_facts(crn, detach_uptakes)

    def formula(vector):
        logs = facts.logs(vector)
        return None if logs is None else Monomial(logs)

    acr, acrr, modules = _assemble(_targets(crn, found), formula, crn.metabolite_ids)
    return Formulas(acr=dict(acr), acrr=acrr, modules=modules, constraints=facts.constraints)


def _flux_facts(crn, detach_uptakes):
    """The span of the flux-coupling facts, with the constraints they put on k"""
    substrate = substrate_complexes(crn)
    skip = zero_complexes(crn) if detach_uptakes else set()
    # fully coupled reactions: proportional rows of ker N, each with its ratio to the first of its group
    _, groups, ratio = row_readout(nullspace(crn.N))
    facts = _Facts()
    for group in groups:
        members = [j for j in group if substrate[j] not in skip]
        if not members:
            continue
        j, *rest = members
        for l in rest:                                # v_l = ρ v_j: k_l ψ_s(l) = ρ k_j ψ_s(j), ρ = ratio[l] / ratio[j]
            c, d = substrate[l], substrate[j]
            facts.add(_difference(crn, c, d),
                      {ratio[l] / ratio[j]: Fraction(1), ("k", j): Fraction(1), ("k", l): Fraction(-1)})
    return facts


def _targets(crn, found):
    """(kind, key, vector) of every relation of a Values result in indices"""
    out = [("acr", i, {i: Fraction(1)}) for i in found.acr]
    out += [("acrr", (r, j), {j: Fraction(1), r: Fraction(-1)}) for r, ratios in found.acrr for j in ratios]
    out += [("module", (c, d), _difference(crn, d, c)) for c, ratios in found.modules for d in ratios]
    return out


def _assemble(found, formula, names):
    """acr, acrr and modules from the relations (kind, key, vector) and formula(vector) for each"""
    acr, acrr, modules = {}, {}, {}
    for kind, key, vector in found:
        f = formula(vector)
        if kind == "acr":
            acr[names[key]] = f
        elif kind == "acrr":
            acrr.setdefault(names[key[0]], {})[names[key[1]]] = f
        else:
            modules.setdefault(key[0], {})[key[1]] = f
    return list(acr.items()), list(acrr.items()), list(modules.items())
