"""Run a figure script and record every data file it reads, which gives the list for the archive.

    python figures/record_inputs.py <list file> <script> [args...]      # from benchmarks/
Appends to <list file> the data files the script opened for reading, without code and fonts. Files under benchmarks/
are written relative to it, files under the scratch results as absolute paths. make_zenodo_bundle.py checks that the
archive holds every file of the list.
"""
import os, runpy, sys

listing, script = sys.argv[1:3]
sys.argv = sys.argv[2:]
HERE = os.getcwd() + os.sep
ROOTS = (HERE, "/scratch/sadreyev_lab/pschneider/")
WRITE = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_APPEND
SKIP = ("__pycache__", "font_cache", "/fontlist-v")            # caches of Python, the SVG fonts and matplotlib
seen = set()


def hook(event, args):
    # io.open reports the mode, FileIO and os.open report None and the flags; a file opened for writing is no input
    if event == "open" and isinstance(args[0], (str, bytes)) and (args[1] in ("r", "rb", "rt") or args[1] is None and not args[2] & WRITE):
        p = os.path.abspath(args[0] if isinstance(args[0], str) else args[0].decode())
        if p.startswith(ROOTS) and not p.endswith((".py", ".pyc", ".ttf", ".otf", ".woff2", ".afm")) and not any(s in p for s in SKIP):
            seen.add(p[len(HERE):] if p.startswith(HERE) else p)


sys.addaudithook(hook)
sys.path.insert(0, os.path.dirname(os.path.abspath(script)))
try:
    runpy.run_path(script, run_name="__main__")
finally:
    with open(listing, "a") as f:
        for p in sorted(seen):
            f.write(p + "\n")
    print(f"{script}: {len(seen)} inputs", file=sys.stderr)
