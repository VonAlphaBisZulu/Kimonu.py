"""Run a figure script and record every data file it reads, which gives the list for the archive.

    python figures/record_inputs.py <list file> <script> [args...]
Appends to <list file> the absolute paths opened for reading under benchmarks/ or the scratch results, without code.
"""
import os, runpy, sys

listing, script = sys.argv[1:3]
sys.argv = sys.argv[2:]
ROOTS = (os.getcwd() + os.sep, "/scratch/sadreyev_lab/pschneider/")
seen = set()


def hook(event, args):
    if event == "open" and isinstance(args[0], (str, bytes)) and (args[1] in (None, "r", "rb", "rt") or args[1] is None):
        p = os.path.abspath(args[0] if isinstance(args[0], str) else args[0].decode())
        if p.startswith(ROOTS) and not p.endswith((".py", ".pyc", ".ttf", ".otf", ".woff2", ".afm")) and "__pycache__" not in p and "font_cache" not in p:
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
