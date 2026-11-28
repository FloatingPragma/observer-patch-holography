"""Count native exterior oscillator ticks in the ordinary and paired-record source trials."""

import argparse
from hashlib import sha256
import json
from math import asin, atan2, ceil, floor, hypot, log, pi, sqrt

import numpy as np
from scipy import sparse
import sympy as sp

from . import settled_midpoint_force as force

HERE = force.HERE
SEED = 28122827
SHOTS = 2**24
TIMES = (0.25, 1, 4, 16, 64)
ARMS = ("empty", "ordinary", "fresh", "shared")
SETTINGS = (2 * len(ARMS) + 1) * len(TIMES) * 4
EPS = sqrt(log(2 * SETTINGS / 0.01) / (2 * SHOTS)) + 1e-12


def laplacian(graph):
    adjacency, degree = graph
    return (sparse.diags(degree, dtype=int) - adjacency.astype(int)).tocsr()


def clock(lap, cell):
    sites = 12 * cell + np.arange(12)
    inside = lap[sites][:, sites]
    outside = np.setdiff1d(np.arange(lap.shape[0]), sites)
    gateways = np.unique(lap[outside][:, sites].nonzero()[1]).tolist()
    local = sp.Matrix(inside.toarray())
    bank = []
    multiplicities = []
    for energy in (4, 6):
        rows = (local - energy * sp.eye(12)).col_join(sp.eye(12)[gateways, :])
        solutions = rows.nullspace()
        multiplicities.append(len(solutions))
        v = solutions[0] * sp.ilcm(*[k.q for k in solutions[0]])
        full = np.zeros(lap.shape[0], int)
        full[sites] = np.array(v, dtype=int).ravel()
        np.testing.assert_array_equal(lap @ full, energy * full)
        bank.append(full)
    levels = np.column_stack(bank)
    assert int(levels[:, 0] @ levels[:, 1]) == 0
    port = int(np.argmax(abs(levels[:, 0] * levels[:, 1])))
    assert levels[port].prod() != 0
    levels[:, 1] *= np.sign(levels[port].prod())
    return (
        levels,
        port,
        {
            "multiplicities_at_4_and_6": multiplicities,
            "port": port,
            "levels": levels[sites].T.tolist(),
        },
    )


def held_tape(amplitudes, rng):
    # For start (u + i^j v)/sqrt(2), read only occupation and its complement.
    u, v = amplitudes
    chances = [abs(u + c * v) ** 2 / 2 for c in (1, 1j, -1, -1j)]
    hits = [int(rng.binomial(SHOTS, p)) for p in chances]
    return [[n, SHOTS - n] for n in hits]


def estimate(tapes, product):
    band = [-pi / TIMES[0], pi / TIMES[0]]
    phases = []
    for time, tape in zip(TIMES, tapes, strict=True):
        counts = [x[0] / SHOTS for x in tape]
        x = (counts[0] - counts[2]) / (2 * product)
        y = (counts[1] - counts[3]) / (2 * product)
        delta = sqrt(2) * EPS / (product * hypot(x, y))
        assert delta < 1
        width = asin(delta) + 1e-12
        angle = atan2(y, x)
        first = ceil((time * band[0] - angle - width) / (2 * pi))
        last = floor((time * band[1] - angle + width) / (2 * pi))
        assert first == last
        angle += 2 * pi * first
        band = [
            max(band[0], (angle - width) / time),
            min(band[1], (angle + width) / time),
        ]
        assert band[0] <= band[1]
        phases.append({"phase": angle, "error": width, "frequency": band})
    return {"counts": tapes, "phases": phases}


def difference(a, b):
    return [a[0] - b[1], a[1] - b[0]]


def run():
    rng = np.random.default_rng(SEED)
    source = json.loads((HERE / "settled_midpoint_force.json").read_bytes())
    board = force.federation(0)
    cal, _ = force.routes.calibrate(force.routes.codebook()[0])
    active = force.routes.graph(board, cal)
    adj, deg = active(-1, -1, 0, 0)
    base = laplacian((adj, deg))
    near, far = source["potential_edges"][0]["owners"]
    banks = {c: clock(base, c) for c in range(board.cells)}
    gaps = {str(c): b[2] for c, b in banks.items()}
    group = np.column_stack([banks[c][0] for c in (near, far)])
    for i, j in source["pairs"]:
        for k, w in np.ndindex(9, 9):
            lap = laplacian(active(i, j, k + 1, w + 1))
            np.testing.assert_array_equal(lap @ group, group * [4, 6, 4, 6])
    ordinary = banks[0][0][:, 0]
    for q, interaction in force.meter.contact_terms().items():
        h, x, y = force.meter.two_loads(base * 0.5, adj, interaction, q)
        h2 = (2 * h).astype(int)
        assert (2 * h - h2).nnz == 0
        modes = ordinary[x, None] * group[y, :]
        if q == 1:
            modes += ordinary[y, None] * group[x, :]
        np.testing.assert_array_equal(h2 @ modes, modes * [8, 10, 8, 10])
    readings = {}
    for c in (near, far):
        levels, port, _ = banks[c]
        meters = levels / np.linalg.norm(levels, axis=0)
        u, v = meters[port]
        reads = {}
        for arm in ARMS:
            tapes = [held_tape((u, v * np.exp(-1j * t)), rng) for t in TIMES]
            reads[arm] = estimate(tapes, u * v)
        readings[str(c)] = reads
    change = laplacian(active(near, -1, 9, 0))
    eig, basis = np.linalg.eigh(change.toarray() * 0.5)
    levels, port, _ = banks[near]
    coeff = basis.T @ (levels / np.linalg.norm(levels, axis=0))
    control = [
        held_tape((basis[port] * np.exp(-1j * t * eig)) @ coeff, rng) for t in TIMES
    ]
    near_base = readings[str(near)]["empty"]["counts"]
    sensitivity = max(
        (abs(changed[0] - blank[0]) / SHOTS, j, k)
        for j, (new, old) in enumerate(zip(control, near_base, strict=True))
        for k, (changed, blank) in enumerate(zip(new, old, strict=True))
    )
    assert sensitivity[0] > 2 * EPS
    shifts = {}
    for arm in ARMS:
        far_read, near_read = [
            readings[str(c)][arm]["phases"][-1]["frequency"] for c in (far, near)
        ]
        shifts[arm] = difference(far_read, near_read)
    contrast = {
        s: difference(shifts[s], shifts[b])
        for s, b in (("ordinary", "empty"), ("shared", "fresh"))
    }
    return {
        "seed": SEED,
        "shots_per_setting": SHOTS,
        "number_of_settings": SETTINGS,
        "simultaneous_probability_error": EPS,
        "times": TIMES,
        "owners": [near, far],
        "fixed_modes_per_cell": gaps,
        "exterior_sharp_pair_configurations": len(source["pairs"]) * 81,
        "material_clock_reads": readings,
        "far_minus_near": shifts,
        "change_of_clock_difference": contrast,
        "near_local_word": {
            "word": 9,
            "counts": control,
            "largest_absolute_occupation_change": sensitivity[0],
            "difference_halfwidth": 2 * EPS,
            "time": TIMES[sensitivity[1]],
            "phase_index": sensitivity[2],
        },
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("MIDPOINT_CLOCKS.md", "settled_midpoint_force.json")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opts = cli.parse_args()
    z = run()
    target = HERE / "midpoint_clocks.json"
    z["times"] = list(z["times"])
    if opts.write:
        target.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(target.read_bytes()), z)
    print(
        z["change_of_clock_difference"],
        z["near_local_word"]["largest_absolute_occupation_change"],
    )
