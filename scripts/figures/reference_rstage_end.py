"""End of each R stage of the exhaustive reference, for fig. S1: the modification time of its ACR log.

    python figures/reference_rstage_end.py      # from benchmarks/; writes figures/reference_rstage_end.json

The R stage of the reference (reference_final/acr/run_chain.sh) writes no end stamp, and its ACR log was last written
when the stage ended. fig4_figS1_S2_benchmark.py takes the stage as the time from the start stamp in the chain log to
that moment. A zip keeps file times only to 2 s, so the Zenodo archive cannot carry them, and they are recorded here
(seconds since the epoch, as os.path.getmtime gives them) from the frozen records in the checkout, once. A file that
exists already must agree with the logs, so the record is never rewritten from copies whose times were rounded.
"""
import glob, json, os

OUT = "figures/reference_rstage_end.json"
new = {os.path.basename(f)[:-len(".acr.log")]: os.path.getmtime(f) for f in sorted(glob.glob("reference_final/acr/*.acr.log"))}
assert new, "no reference_final/acr/*.acr.log here (run from benchmarks/ of the checkout with the frozen records)"
if os.path.exists(OUT):
    old = json.load(open(OUT))
    assert old == new, f"{OUT} differs from the modification times of the ACR logs (copies with rounded times?)"
    print(f"{OUT} agrees with the {len(new)} ACR logs")
else:
    json.dump(new, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}, {len(new)} ACR logs")
