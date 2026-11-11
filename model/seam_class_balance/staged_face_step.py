"""driver for the second reading (A at the midpoint of the face step), engine engine_face_step.

Usage: python -m model.seam_class_balance.staged_face_step <tasks.json> <process index> <number of processes> <output directory>
Task fields as in staged.py, plus mode: "fm" (time-like plaquette above the seam with weight 1, the faces of the seam
in that time slice with weight 1/2) or "el" (the time-like plaquette alone; equals staged.py). Reads no measured value.
"""
import json
import math
import os
import sys
import time

import numpy as np

from . import complexes as C
from . import engine as E
from . import engine_face_step as F
from .staged import NBLOCK, estimate


def twist_set(cx, mode):
    P0 = cx["P0"]
    tw = np.zeros(cx["nplaq"])
    tw[P0] = 1.0
    twp = [P0]
    if mode == "fm":
        l = int(cx["pl"][P0, 3]) if cx["name"].startswith("cub") else int(cx["pl"][P0, 0])
        for a in range(cx["lp_ptr"][l], cx["lp_ptr"][l + 1]):
            p = int(cx["lp_idx"][a])
            if p != P0 and cx["ptype"][p] == 0:
                tw[p] = 0.5
                twp.append(p)
    return tw, np.array(twp, np.int64)


def run_task(t):
    t0 = time.time()
    cx = C.build(t["cx"])
    kind, par = int(t["kind"]), float(t["par"])
    N, ncyc, nloc, ntherm, nhit = int(t["N"]), int(t["ncyc"]), int(t["nloc"]), int(t["ntherm"]), int(t.get("nhit", 4))
    E.seed_rng(int(t["seed"]))
    U, nf = E.cold(cx)
    g2 = 4.0 / par if kind == 0 else par
    a = (cx["pl"], cx["ps"], cx["pn"], cx["lp_ptr"], cx["lp_idx"], cx["lp_pos"])
    links = np.arange(cx["nlink"])
    tw, twp = twist_set(cx, t.get("mode", "fm"))
    tw0 = np.zeros(cx["nplaq"])
    step = 0.5 * math.sqrt(g2) / 2.0
    # step size from the acceptance of single hits on the untwisted lattice (through the main engine)
    sites = np.arange(cx["nsite"])
    full = E.args_of(cx)
    for _ in range(24):
        acc, _n = E.sweep(U, nf, links, sites, *full, kind, par, -1, 0.0, 0.0, 0.0, nhit, step, 0.3)
        r = acc / (links.size * nhit)
        step = min(1.5, step * min(1.5, max(0.6, (r + 0.1) / 0.55)))
    near = set()
    for p in twp:
        near |= set(C.shells(cx, int(p), int(t.get("nshell", 1))).tolist())
    near = np.array(sorted(near), np.int64)
    out = {"task": t, "complex": {k: cx[k] for k in ("name", "nsite", "nlink", "nplaq")}, "step": step,
           "near_links": int(near.size), "twisted_plaquettes": twp.tolist(), "twist_weights": tw[twp].tolist(),
           "g2": g2, "x": 4.0 * math.pi / g2, "S": 2.0 * math.pi ** 2 / g2}
    for _ in range(int(t.get("pretherm", 100))):
        F.sweep(U, links, *a, kind, par, tw0, 0.0, nhit, step)
    nodes = np.zeros((N + 1, NBLOCK, 7))
    ahalf = 0.5 / N
    for k in range(N + 1):
        nodes[k] = F.run_node(U, links, near, *a, kind, par, tw, twp, k / N, nhit, step, ntherm, ncyc, nloc, ahalf, NBLOCK)
    out.update(estimate(nodes, N))
    out["nodes"] = nodes.tolist()
    out["seconds"] = time.time() - t0
    return out


if __name__ == "__main__":
    if not E.AVAILABLE:
        raise SystemExit("the staged Monte Carlo needs numba")
    tasks = json.load(open(sys.argv[1]))
    i, n, od = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    os.makedirs(od, exist_ok=True)
    order = list(enumerate(tasks))
    if len(sys.argv) > 5 and sys.argv[5] == "rev":      # a helper process takes the same share from the end of the list
        order = order[::-1]
    for j, t in order:
        if j % n != i:
            continue
        fn = os.path.join(od, t["name"] + ".json")
        if os.path.exists(fn):
            continue
        r = run_task(t)
        json.dump(r, open(fn + ".tmp", "w"))
        os.replace(fn + ".tmp", fn)
        print("%s lnA %.4f +- %.4f (fwd %.4f bwd %.4f TI %.4f) %.0fs" % (
            t["name"], r["lnA"], r["lnA_err"], r["lnA_forward"], r["lnA_backward"], r["lnA_TI"], r["seconds"]), flush=True)
