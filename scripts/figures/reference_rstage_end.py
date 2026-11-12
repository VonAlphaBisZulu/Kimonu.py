"""End of each R stage of the exhaustive reference, for fig. S1: the modification time of its ACR log.

    python figures/reference_rstage_end.py      # from scripts/; writes reference_rstage_end.json to figures/source_data/

The R stage of the reference (reference/run_chain.sh) writes no end stamp, and its ACR log was last written when the
stage ended. fig4_figS1_S2_benchmark.py takes the stage as the time from the start stamp in the chain log to that
moment. A zip keeps file times only to 2 s, so the Zenodo archive cannot carry them, and they are recorded here
(seconds since the epoch, as os.path.getmtime gives them) from the frozen records (results/reference/exhaustive/acr/
of KIMONU_RECORDS, copied with their times), once. A file that exists already must agree with the logs, so the record
is never rewritten from copies whose times were rounded.
"""
import glob, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import records

OLD, OUT = records.read(records.SOURCE_DATA, "reference_rstage_end.json"), records.write(records.SOURCE_DATA, "reference_rstage_end.json")
logs = records.read(records.REF_EXHAUSTIVE, "acr", "*.acr.log")
new = {os.path.basename(f)[:-len(".acr.log")]: os.path.getmtime(f) for f in sorted(glob.glob(logs))}
assert new, f"no {logs}"
if os.path.exists(OLD):
    old = json.load(open(OLD))
    assert old == new, f"{OLD} differs from the modification times of the ACR logs (copies with rounded times?)"
    print(f"{OLD} agrees with the {len(new)} ACR logs")
    if OUT != OLD:
        json.dump(old, open(OUT, "w"), indent=1)
else:
    json.dump(new, open(OUT, "w"), indent=1)
    print(f"wrote {OUT}, {len(new)} ACR logs")
