"""One call from a cobra model to the relations: the default mode and its alternatives.

:func:`analyze` runs the preprocessing that :func:`kimonu.concordance` also runs (:func:`kimonu.prepare`: split
reversible reactions, remove blocked reactions, find pinned rates, build N = Y·A) and then one route. The default
mode is the iterated route with detached uptakes and seed 0 for the reference flux.
"""
from kimonu.certify import certify
from kimonu.io import CRN, prepare
from kimonu.iterated import iterated_route
from kimonu.reference import reference_route
from kimonu.tangent import tangent_route

ROUTES = ("iterated", "tangent", "reference", "both")


def analyze(model, route="iterated", detach_uptakes=True, seed=0, draw="lp", bits=64, solver=None):
    """Kinetic modules, ACR and ACRR of a cobra model, or of a :class:`kimonu.CRN` from :func:`kimonu.prepare`.

    Args:
        model: a cobra model, which is prepared first, or a CRN that :func:`kimonu.prepare` returned. Pass the CRN
            to read the complexes and reactions behind the indices of the result (``crn.complexes``,
            ``crn.reaction_ids``).
        route: "iterated" (default, :func:`kimonu.iterated_route`), proven relations; "tangent"
            (:func:`kimonu.tangent_route`), every true relation and possibly false couplings; "reference"
            (:func:`kimonu.reference_route`), the definitions of the reference implementation of Langary et al.
            (2025); or "both" (:func:`kimonu.certify`), both mass-action routes at one reference flux and whether
            they agree.
        detach_uptakes, seed, draw, bits: as for :func:`kimonu.iterated_route`. The reference route uses none of
            them.
        solver: solver for FVA, which :func:`kimonu.prepare` needs only when a positive lower bound forces flux.

    Returns a :class:`kimonu.RouteResult`, or a :class:`kimonu.Certificate` for route="both". Metabolites are named
    by their ids, complexes by their index in ``crn.complexes``.
    """
    if route not in ROUTES:
        raise ValueError(f"route must be one of {', '.join(ROUTES)}, not {route!r}")
    crn = model if isinstance(model, CRN) else prepare(model, solver)
    A, Y, N, ids = crn.A, crn.Y, crn.N, crn.metabolite_ids
    if route == "iterated":
        return iterated_route(A, Y, N, detach_uptakes, seed, draw, bits, metabolite_ids=ids)
    if route == "tangent":
        return tangent_route(A, Y, N, detach_uptakes, seed, draw, bits, metabolite_ids=ids)
    if route == "reference":
        return reference_route(A, Y, N, metabolite_ids=ids)
    return certify(A, Y, N, detach_uptakes, seed, draw, bits, metabolite_ids=ids)
