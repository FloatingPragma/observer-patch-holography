"""Compare a freshly executed full run with the committed measured evidence.

Exact arrays and integer fields must be identical. Derived floating readbacks,
kernels and statistics allow 5e-10 absolute/relative rounding tolerance. ZIP
compression bytes and the hash chains binding those bytes are not portable
numeric invariants; each trace is independently authenticated by verify.py.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def same(a, b, path="root"):
    if type(a) is not type(b):
        raise ValueError("reproduction type: "+path)
    if isinstance(a, dict):
        if set(a) != set(b):
            raise ValueError("reproduction keys: "+path)
        for key in a:
            same(a[key], b[key], path+"."+key)
    elif isinstance(a, list):
        if len(a) != len(b):
            raise ValueError("reproduction list length: "+path)
        for index, (x, y) in enumerate(zip(a, b)):
            same(x, y, path+f"[{index}]")
    elif isinstance(a, float):
        if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a, b, rel_tol=5e-10, abs_tol=5e-10):
            raise ValueError("reproduction float: "+path)
    elif a != b:
        raise ValueError("reproduction exact field: "+path)


def arrays(a, b):
    with np.load(a) as old, np.load(b) as new:
        if set(old.files) != set(new.files):
            raise ValueError("reproduction array census")
        for key in old.files:
            x, y = old[key], new[key]
            if x.shape != y.shape or x.dtype != y.dtype:
                raise ValueError("reproduction array shape or dtype: "+key)
            ok = (np.all(np.isfinite(x)) and np.all(np.isfinite(y))
                  and np.allclose(x, y, atol=5e-10, rtol=5e-10)) if x.dtype.kind == "f" else np.array_equal(x, y)
            if not ok:
                raise ValueError("reproduction array values: "+key)


def compare(fresh, original):
    read = lambda base, name: json.loads((base/name).read_text())
    old, new = read(original, "trace/trace.json"), read(fresh, "trace/trace.json")
    same(old["binding"], new["binding"])
    same(old["events"], new["events"])
    same(len(old["chunks"]), len(new["chunks"]))
    for a, b in zip(old["chunks"], new["chunks"]):
        same({k: v for k, v in a.items() if k not in ("sha256", "previous", "chain")},
             {k: v for k, v in b.items() if k not in ("sha256", "previous", "chain")})
        arrays(original/"trace"/a["file"], fresh/"trace"/b["file"])
    for name in ("controls.json", "provenance.json", "q13.json", "q13_controls.json"):
        a, b = read(original, name), read(fresh, name)
        if name == "q13.json":
            a.pop("trace_sha256")
            b.pop("trace_sha256")
        if name == "q13_controls.json":
            for family in a["families"] + b["families"]:
                family.pop("trace_sha256")
        same(a, b, name)
    for name in ("kernels.npz", "interval_members.npz", "q13_reads.npz",
                 "q13_control_d1_reads.npz", "q13_control_d2_reads.npz"):
        arrays(original/name, fresh/name)
    print("SUPPORT_WIRING_FULL_RUN_REPRODUCED", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fresh", type=Path)
    p.add_argument("--original", type=Path, default=ROOT/"evidence/support_wiring_776")
    a = p.parse_args()
    compare(a.fresh, a.original)
