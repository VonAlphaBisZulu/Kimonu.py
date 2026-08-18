"""Kimonu — exact, deterministic concordance and kinetic-module analysis of metabolic networks.

Pure Python (no Julia). Depends on SENUS (exact sparse nullspace) and cobra.

Core idea: concordance is flux coupling lifted one matrix product up. With N = Y·A,
K = ker(N) (exact, sparse, via SENUS), and M = A·K, concordant complexes are the
proportional rows of M and balanced complexes are its zero rows — exact linear algebra,
no LP, no sampling. See ROADMAP.md.
"""
from kimonu.io import build_crn
from kimonu.concordance import concordance, concordance_from_matrices, ConcordanceResult
from kimonu.kinetic import kinetic_modules, KineticResult
from kimonu.acr import acr_acrr, AcrResult

__version__ = "0.0.1"
__all__ = ["build_crn", "concordance", "concordance_from_matrices", "ConcordanceResult",
           "kinetic_modules", "KineticResult", "acr_acrr", "AcrResult", "__version__"]
