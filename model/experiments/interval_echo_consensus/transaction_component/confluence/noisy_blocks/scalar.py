"""Single-launch held branches using the original routed scalar meter."""

import argparse
import gzip
from itertools import product
import json
from math import pi
from pathlib import Path

import numpy as np

from ....seam_word_meter import apparatus as seam
from ..apparatus import COMMIT, FAR, MASKS, Experiment, layout
from .instrument import NOMINAL

HERE = Path(__file__).resolve().parent
MODES = (2, 1, 0, 2, 0, 1, 1, 0)


def scalar_quarter(device, branch):
    clock = device.m.rows[14][0]
    source = tuple(device.loads[j] for j in clock)
    if source != (2, 0):
        raise ArithmeticError("Scalar clock was not reset")
    full = seam.pulse(pi / 3)[:, seam.LAW.SRC.index(source)]
    chosen = np.kron(seam.projectors()[0 if branch else 2], np.eye(3)) @ full
    weight = float(np.linalg.norm(chosen) ** 2)
    np.testing.assert_allclose(weight, (1 if branch else 3) / 4, atol=seam.TOL, rtol=0)
    # Each projected local output is sharp; publish those two physical loads.
    idx = int(np.argmax(abs(chosen)))
    residual = np.linalg.norm(np.delete(chosen, idx))
    if residual > seam.TOL:
        raise ArithmeticError("Clock branch not sharp")
    a, b = clock
    device.loads[a], device.loads[b] = seam.LAW.SRC[idx]
    seen = 1 - device.read_bit(clock, COMMIT)
    if seen:
        device.loads[a], device.loads[b] = device.loads[b], device.loads[a]
    if device.read_bit(clock, FAR) != 1:
        raise ArithmeticError("Scalar clock return failed")
    return seen


def run(boundary, flags, modes):
    words = tuple(
        tuple(b ^ flags[3 * side + i] for i in range(3))
        for side, b in enumerate(boundary)
    )
    dev = Experiment(boundary, words)
    dev.initialize(dev.m.rows[14][0], 1)
    steps = [dev.checkpoint("start")]
    selectors = []
    for (side, m), chosen in zip(NOMINAL, modes, strict=True):
        start = len(dev.raw)
        first = scalar_quarter(dev, chosen != 0)
        second = scalar_quarter(dev, chosen == 2) if first else None
        mode = (1 + second) if first else 0
        selectors.append([start, len(dev.raw), mode])
        boundary_bit, old, old_rank = dev.state(side, 6 * side, "prepare")
        if dev.state(side, COMMIT, "revalidate") != (boundary_bit, old, old_rank):
            raise ArithmeticError("Guard revalidation failed")
        mask = MASKS[m]
        guards = dev.guard(side, old_rank, mask, COMMIT)
        if guards != [int(old[i] != boundary_bit) for i in mask]:
            raise ArithmeticError("Independent material guard failed")
        target = list(old)
        for i in mask:
            target[i] = (boundary_bit, old[i], 1 - boundary_bit)[mode]
        new_rank = sum(x != boundary_bit for x in target)
        p = layout(side)
        for pair, x in zip(p["banks"][new_rank], target, strict=True):
            dev.write(pair, x)
        if old_rank != new_rank:
            for a, b in dev.m.move(p["selector"][old_rank], p["selector"][new_rank]):
                dev.loads[a], dev.loads[b] = dev.loads[b], dev.loads[a]
        steps.append(dev.checkpoint("commit"))
    return dict(
        boundary=boundary,
        initial_flags=flags,
        branch=modes,
        selectors=selectors,
        checkpoints=steps,
        observations=dev.observations,
        raw=dev.raw,
    )


def report():
    probes = [run(b, (1, 0, 1, 0, 1, 0), MODES) for b in product((0, 1), repeat=2)]
    probes.extend(run((0, 1), (v,) * 6, MODES) for v in (0, 1))
    return dict(
        provenance="Selected sharp native quantum outcomes. Every pulse probability is checked against "
        "its routed ternary effect and each branch has nonzero weight.",
        runs=probes,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--output", type=Path)
    options = cli.parse_args()
    data = (json.dumps(report(), separators=(",", ":")) + "\n").encode()
    if options.output:
        with options.output.open("xb") as f:
            f.write(gzip.compress(data, mtime=0))
    elif gzip.decompress((HERE / "scalar.json.gz").read_bytes()) != data:
        raise ValueError("Scalar held read failed")
