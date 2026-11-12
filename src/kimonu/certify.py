"""Both routes at one reference flux: where they agree, the list of relations is correct and complete (Eq. 4).

The iterated route reports only true relations and can miss some, and the tangent route reports all true relations
and can add false couplings, so their nullspaces bound the true B from both sides,

    ker((YΔ)ᵀ)  ⊇  B  ⊇  ker(W Yᵀ).                                                              (Eq. 4)

Where the two outer nullspaces give the same kinetic modules, ACR and ACRR, both equal the true relations for
generic rate constants. Neither route can show this alone, since the tangent route cannot tell a false coupling
from a true one and the iterated route cannot tell whether it missed a relation. What the tangent route reports
beyond the iterated route stays open. Each open relation asks whether a second steady state exists that no path of
steady states connects to the reference, as in the self-activating kinase of Fig. 2.

:func:`certify` draws one r* and runs both routes there, so they are compared at the same point. It is not the
default analysis, since the tangent route's nullspace is denser and takes longer (up to about 8 minutes on the
largest genome-scale networks at the default LP point, and about ten times that with the uniform draw).
"""
import warnings
from collections import namedtuple

from kimonu.iterated import iterated_route
from kimonu.steady_state import interior_flux
from kimonu.tangent import tangent_route

Certificate = namedtuple("Certificate", ["certified", "modules", "acr", "acrr", "open_pairs", "open_acr", "open_acrr",
                                         "iterated", "tangent", "flux"])
Certificate.__doc__ = """The comparison of the two routes at one reference flux.

certified   True when both routes give the same kinetic modules, ACR and ACRR, so the lists are correct and complete
modules, acr, acrr  the proven relations, those of the iterated route
open_pairs  sorted complex pairs (c, d), c < d, that the tangent route couples and the iterated route does not
open_acr    metabolites with ACR by the tangent route only
open_acrr   metabolite pairs with ACRR by the tangent route only
iterated, tangent   the two :class:`kimonu.RouteResult`
flux        the reference flux r* both routes ran at
"""


def _coupled_pairs(modules):
    return {(a, b) for g in modules for a in g for b in g if a < b}


def certify(A, Y, N, detach_uptakes=True, seed=0, draw="lp", bits=64, metabolite_ids=None):
    """Run the iterated and the tangent route at one reference flux and tell whether their lists agree.

    Args:
        A, Y, N, detach_uptakes, seed, draw, bits, metabolite_ids: as for :func:`kimonu.iterated_route`. One r* is drawn
            with :func:`kimonu.interior_flux` and passed to both routes.

    Returns a :class:`Certificate`. A gap between the routes is not an error. It gives a warning that says how many
    relations stay open.
    """
    flux = interior_flux(N, seed, draw, bits)
    it = iterated_route(A, Y, N, detach_uptakes, metabolite_ids=metabolite_ids, flux=flux)
    tan = tangent_route(A, Y, N, detach_uptakes, metabolite_ids=metabolite_ids, flux=flux)
    pairs_it, pairs_tan = _coupled_pairs(it.modules), _coupled_pairs(tan.modules)
    acr_it, acr_tan, acrr_it, acrr_tan = set(it.acr), set(tan.acr), set(it.acrr), set(tan.acrr)
    if not (pairs_it <= pairs_tan and acr_it <= acr_tan and acrr_it <= acrr_tan):
        warnings.warn("the iterated route reports a relation the tangent route does not, which happens only at a "
                      "finely tuned reference flux. Draw another one with a different seed.", stacklevel=2)
    open_pairs, open_acr, open_acrr = sorted(pairs_tan - pairs_it), sorted(acr_tan - acr_it), sorted(acrr_tan - acrr_it)
    certified = pairs_it == pairs_tan and acr_it == acr_tan and acrr_it == acrr_tan
    if not certified:
        warnings.warn(f"the routes differ, so {len(open_pairs)} complex pairs, {len(open_acr)} ACR and "
                      f"{len(open_acrr)} ACRR relations stay open. The iterated route's lists are proven but may be "
                      "incomplete.", stacklevel=2)
    return Certificate(certified=certified, modules=it.modules, acr=it.acr, acrr=it.acrr, open_pairs=open_pairs,
                       open_acr=open_acr, open_acrr=open_acrr, iterated=it, tangent=tan, flux=flux)
