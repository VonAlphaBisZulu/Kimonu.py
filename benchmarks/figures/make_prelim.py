"""Preliminary figures F1-F5 from figdata.json. See FINDINGS.md for what each arm is and how far
it can be trusted. Run from benchmarks/: python figures/make_prelim.py"""
# figdata.json schema:
#   models, n_complexes[m]
#   kimonu[m]: balanced, conc, giant, giant_4phase, acr, acrr, acr_free, acrr_free, seconds, sec_conc, sec_kin, sec_acr
#   cocoa_false[m] / cocoa_true[m]: list of {seed, balanced, conc, giant, acr, acrr, acr_free, acrr_free, seconds}
#   ref_exh[m]: balanced, conc, seconds, sec_conc, and where the R stage ran: giant, acr, acrr, acr_free, acrr_free
#   ref_sampled[m]: list of {seed, balanced, conc, seconds}
# Free counts: metabolite index < index of "E_1" in the model's own metabolite list.
