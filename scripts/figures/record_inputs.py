"""Run a figure script and record every data file it reads, which gives the list for the archive.

    python figures/record_inputs.py <list file> <script> [args...]      # from scripts/
Appends to <list file> the data files the script opened for reading, without code, fonts and caches, one per line:
"records/<path>" for a file of the records (KIMONU_RECORDS, records.py), "repo/<path>" for a file of the repository
(the source snapshot of the archive) and the absolute path of any other file under the repository, the records or the
cluster folders below, which make_zenodo_bundle.py and zenodo_bundle.sh count as read outside the archive. Files the
script reads back from KIMONU_OUT (its own output) or from the temporary folder are not inputs.
"""
import os, runpy, sys, tempfile

listing, script = sys.argv[1:3]
sys.argv = sys.argv[2:]
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import records

real = lambda p: os.path.realpath(p) + os.sep
REC, REPO, OUT = real(records.RECORDS), real(records.REPO), real(records.OUT)
CLUSTER = ("/scratch/sadreyev_lab/pschneider/", "/home/sadreyev_lab/pschneider/pipelines/")   # earlier results, the vault
TMP = real(tempfile.gettempdir())
WRITE = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_APPEND
SKIP = ("__pycache__", "font_cache", "/fontlist-v", "/.git/")   # caches of Python, the SVG fonts and matplotlib
seen = set()


def hook(event, args):
    # io.open reports the mode, FileIO and os.open report None and the flags; a file opened for writing is no input
    if event == "open" and isinstance(args[0], (str, bytes)) and (args[1] in ("r", "rb", "rt") or args[1] is None and not args[2] & WRITE):
        p = os.path.realpath(args[0] if isinstance(args[0], str) else args[0].decode())
        if p.endswith((".py", ".pyc", ".ttf", ".otf", ".woff2", ".afm")) or any(s in p for s in SKIP):
            return
        if p.startswith(REC):
            seen.add("records/" + p[len(REC):])
        elif OUT != REC and p.startswith(OUT) or p.startswith(TMP):
            return
        elif p.startswith(REPO):
            seen.add("repo/" + p[len(REPO):])
        elif p.startswith(CLUSTER):
            seen.add(p)


sys.addaudithook(hook)
sys.path.insert(0, os.path.dirname(os.path.abspath(script)))
try:
    runpy.run_path(script, run_name="__main__")
finally:
    with open(listing, "a") as f:
        for p in sorted(seen):
            f.write(p + "\n")
    print(f"{script}: {len(seen)} inputs", file=sys.stderr)
