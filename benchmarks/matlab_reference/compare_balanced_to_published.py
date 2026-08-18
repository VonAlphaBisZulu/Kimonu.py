#!/usr/bin/env python
"""Compare a rerun balanced set against the reference's committed result.

Octave's `load` fails on some of the committed .mat files with a zlib error where scipy reads them
without complaint, so the published side is read here rather than in Octave. Sets are compared, not
just counts.

Usage:
    python compare_balanced_to_published.py <published.mat> <rerun_indices.txt>
"""
import sys

import scipy.io as sio


def main():
    published_path, rerun_path = sys.argv[1], sys.argv[2]
    d = sio.loadmat(published_path)
    published = set(int(x) for x in d["B"].ravel())
    with open(rerun_path) as fh:
        rerun = set(int(x) for x in fh.read().split())
    print(f"published {len(published)}, rerun {len(rerun)}")
    print(f"identical: {published == rerun}")
    if published != rerun:
        print(f"  only in rerun: {sorted(rerun - published)[:20]}")
        print(f"  only in published: {sorted(published - rerun)[:20]}")


if __name__ == "__main__":
    main()
