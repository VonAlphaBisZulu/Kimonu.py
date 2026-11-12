"""Formulas in k for the relations: binomial ones from flux coupling at any size, all of them on small networks.

relation_values gives the value of every relation the iterated route proves at one k. As a function of
k, each value is Π_t λ_t(k)^a_t, with λ_t the ratio of two rows of ker(W_k), a rational function of k.
So every relation the iterated route proves has a formula. What limits it is size: λ_t comes from a
kernel over the field Q(k) of rational functions in all rate constants, and its numerator and
denominator can have as many terms as there are spanning trees in the part of the network that fixes
it. Two ways are offered.

Binomial values (binomial_formulas). When two reactions j and l are fully coupled, v_j = ρ v_l for every
steady-state flux (proportional rows of ker N), mass action gives k_j ψ_s(j) = ρ k_l ψ_s(l) at every
positive steady state, for every k. So

    ln ψ_s(j) - ln ψ_s(l) = ln ρ + ln k_l - ln k_j,

a composition difference whose value is the logarithm of a monomial in k. These facts span composition
differences as in values.py, and every relation whose vector lies in that span has the formula
c Π_j k_j^e_j that the reduction returns. The facts hold at every k, not only for generic k. They need
only the exact null space of N, as concordance does, so this runs at genome scale. A detached uptake
runs at k_j s_j with a free source s_j and gives no fact.

A fact that repeats a combination already spanned gives a constraint instead: a monomial in k that
equals one at every k with a positive steady state. Two uptakes from the attached zero complex whose
fluxes are fully coupled, for instance, need k_j = ρ k_l.

All values (symbolic_formulas). The derivation of relation_values carried out over Q(k) with sympy:
every relation the iterated route proves, as an expression in k. Exact, and for small networks only.

What neither gives is a formula for a relation the iterated route does not prove. On the kinase of
Fig. 2 (:mod:`kimonu.tangent`), X has ACR: with a, b and c the rate constants of 0 -> X, X -> 0 and
2X -> C, the balance of X alone is a - b x - 2 c x² = 0, so x = (√(b² + 8ac) - b) / 4c at both steady
states. That is the root of a quadratic, not a product of powers of rational functions of k, and no
linear step can produce it.
"""
from collections import namedtuple
from fractions import Fraction

from kimonu.io import substrate_complexes, zero_complexes
from kimonu.linalg import entries, null_space_rows
from kimonu.readout import free_directions, readout
from kimonu.steady_state import interior_flux, rate_columns
from kimonu.values import Value, _Span, _difference, relation_values

Formulas = namedtuple("Formulas", ["acr", "acrr", "modules", "constraints"])


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
            base = k[j] if isinstance(k[j], (int, Fraction)) else Fraction(str(k[j]))
            logs[Fraction(base)] = logs.get(Fraction(base), 0) + e
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


def binomial_formulas(crn, detach_uptakes=True, seed=0, draw="lp", bits=64):
    """Monomial formulas in k for the relations that flux coupling fixes, and the constraints it puts on k.

    Args:
        crn: the network, as for :func:`kimonu.iterated_route`.
        detach_uptakes: as for :func:`kimonu.iterated_route`.
        seed, draw, bits: of the interior flux r* at which the relations are found (relation_values at k = r*,
            :func:`kimonu.interior_flux`).

    Returns a :class:`Formulas` in the layout of :class:`kimonu.Values`, with a :class:`Monomial` for
    every relation that flux coupling fixes and None for the others: ``acr`` {metabolite: formula},
    ``acrr`` one (reference, {metabolite: formula of x_j / x_reference}) per group, ``modules`` one
    (reference, {complex: formula of ψ_c / ψ_reference}) per module. ``constraints`` are monomials
    that equal one at every k with a positive steady state.
    """
    found = relation_values(crn._replace(metabolite_ids=list(range(len(crn.metabolite_ids)))),
                            interior_flux(crn, seed, draw, bits), detach_uptakes)
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
    # fully coupled reactions: proportional rows of ker N, each with its ratio to the first of its class
    classes = {}
    for j, row in enumerate(null_space_rows(crn.N)):
        first = row[0][1]
        classes.setdefault(tuple((col, v / first) for col, v in row), []).append((j, first))
    facts = _Facts()
    for members in classes.values():
        members = [(j, f) for j, f in members if substrate[j] not in skip]
        if not members:
            continue
        (j, f), *rest = members
        for l, g in rest:                             # v_l = (g / f) v_j: k_l ψ_s(l) = (g / f) k_j ψ_s(j)
            c, d = substrate[l], substrate[j]
            facts.add(_difference(crn, c, d), {g / f: Fraction(1), ("k", j): Fraction(1), ("k", l): Fraction(-1)})
    return facts


class _SymbolicSpan(_Span):
    """The span with log values over Q(k); a repeated combination agrees for generic k"""

    def _repeated(self, value):
        pass


def symbolic_formulas(crn, detach_uptakes=True):
    """Every relation the iterated route proves, with its value as a sympy expression in k.

    The derivation of :func:`kimonu.relation_values` over the field Q(k) of rational functions in the
    rate constants, which are the positive symbols k_<reaction id>. Needs sympy; the kernels over Q(k) grow fast, so this is for small networks.

    Where flux coupling constrains k (see :func:`binomial_formulas`), generic k has no positive steady
    state, so each constraint is solved for one rate constant first, and the formulas hold on the k
    that satisfy them. Returns a :class:`Formulas` in the layout of :class:`kimonu.Values` with sympy
    expressions; ``constraints`` are these substitutions as equations k_j = expression.
    """
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix

    n = len(crn.Yt)
    k = [sp.Symbol(f"k_{r}", positive=True) for r in crn.reaction_ids]
    # the constraints of flux coupling, each solved for its last rate constant with exponent ±1
    solved = {}
    for constraint in _flux_facts(crn, detach_uptakes).constraints:
        if constraint.constant.rational is None or any(e.denominator != 1 for e in constraint.exponents.values()):
            continue
        expr = sp.Rational(*constraint.constant.rational.as_integer_ratio()) * sp.Mul(
            *(k[j] ** int(e) for j, e in constraint.exponents.items()))
        powers = sp.powsimp(expr.subs(solved)).as_powers_dict()
        unit = [s_ for s_ in k if powers.get(s_) in (1, -1)]
        if unit:
            s_, e = unit[-1], powers[unit[-1]]
            value = (expr.subs(solved) / s_ ** e) ** -e
            solved = {key: v.subs(s_, value) for key, v in solved.items()}
            solved[s_] = value
    F = sp.QQ.frac_field(*k)
    rate = [F.from_sympy(sp.sympify(solved.get(s_, s_))) for s_ in k]
    skip = zero_complexes(crn) if detach_uptakes else set()
    substrate = substrate_complexes(crn)
    substrates = set(substrate.values()) - skip
    by_reaction = {}
    triples, (m, _) = entries(crn.N)
    for i, j, v in triples:
        by_reaction.setdefault(j, []).append((i, F.convert(sp.Rational(v.numerator, v.denominator))))
    columns, n_columns = rate_columns(crn, detach_uptakes)       # as steady_state_map

    def power(logs):
        out = F.one
        for b, e in logs.items():
            out *= b ** int(e)
        return out

    equal, rounds = [], 0
    while True:
        rounds += 1
        W = [[F.zero] * n_columns for _ in range(m + len(equal))]
        for j, c in substrate.items():
            for i, v in by_reaction.get(j, ()):
                W[i][columns[j]] += rate[j] * v
        for t, (c, d, mu) in enumerate(equal):
            W[m + t][c], W[m + t][d] = F.one, -mu
        basis = DomainMatrix(W, (len(W), n_columns), F).nullspace().to_list()
        groups = {}
        for c in range(n):
            if c in skip:
                continue
            row = [b[c] for b in basis]
            first = next((v for v in row if v), None)
            if first is None:
                raise ValueError(f"complex {c} is zero at every steady state for generic k: the rate constants "
                                 "are constrained beyond what flux coupling shows")
            groups.setdefault(tuple(v / first for v in row), {})[c] = first
        classes = [{c: f / g[min(g)] for c, f in g.items()} for g in groups.values()]
        span = _SymbolicSpan()
        for g in classes:
            first = min(g)
            for c in sorted(g):
                if c != first:
                    span.add(_difference(crn, first, c), {g[first] / g[c]: Fraction(1)})
        found = readout(crn, free_directions(crn, [set(g) for g in classes]), detach_uptakes)
        class_of = {c: t for t, g in enumerate(classes) for c in g}
        new = []
        for module in found.modules:
            inner = sorted(module & substrates)
            for c in inner[1:]:
                if class_of[c] != class_of[inner[0]]:
                    logs = span.logs(_difference(crn, inner[0], c))
                    if any(e.denominator != 1 for e in logs.values()):
                        raise NotImplementedError(f"the feedback needs a root of a ratio of complexes {inner[0]} and {c}")
                    new.append((inner[0], c, power(logs)))
        if not new:
            break
        equal += new

    def formula(vector):
        logs = span.logs(vector)
        return sp.factor(sp.Mul(*(F.to_sympy(b) ** sp.Rational(e.numerator, e.denominator) for b, e in logs.items())))

    acr_found = [("acr", i, {i: Fraction(1)}) for i in found.acr]
    acrr_found = [("acrr", (g[0], j), {j: Fraction(1), g[0]: Fraction(-1)}) for g in found.acrr_groups for j in g[1:]]
    module_found = [("module", (min(mod), d), _difference(crn, d, min(mod)))
                    for mod in sorted(found.modules, key=min) if len(mod) > 1 for d in sorted(mod - {min(mod)})]
    acr, acrr, modules = _assemble(acr_found + acrr_found + module_found, formula, crn.metabolite_ids)
    return Formulas(acr=dict(acr), acrr=acrr, modules=modules,
                    constraints=[sp.Eq(s_, v) for s_, v in solved.items()])
