"""Change the flight law to test sensitivity to a joint-head potential on formed stock."""

import argparse
from hashlib import sha256
import json

import numpy as np
from scipy import sparse

from . import formed_tides, tidal

FORMATION = formed_tides.FORMATION
RECEIPT = tidal.HERE / "joint_control.json"
MISMATCH = np.diag([0.0, 1.0, 1.0, 0.0])


def measured_mismatch(arm, shots):
    read = next(row for row in arm["modular"]["held_reads"] if row["exposure"] == 4)
    counts = np.asarray(read["joint_and_product_counts"])
    if not np.all(counts.sum(axis=(1, 2)) == shots):
        raise ValueError("Incomplete formed-stock count")
    # W=1 on sharp outcomes 01 and 10. Four means across both formed stocks.
    mismatch = counts[:, 0, 1] + counts[:, 1, 0]
    radius = 2 * np.sqrt(np.log(2 * 4 / 0.01) / (2 * shots))
    return {
        "mismatch_counts_formed_product": mismatch.tolist(),
        "product_minus_formed": float((mismatch[1] - mismatch[0]) / shots),
        "simultaneous_99_radius": float(radius),
    }


def run():
    data = json.loads(formed_tides.SOURCE_COUNTS.read_text())
    bank = next(row for row in data["experiments"] if row["cells"] == 80)
    board = FORMATION.federation(1)
    route = FORMATION.region(board)
    radii, paths, launch = tidal.chart(board, route)
    if list(route) != bank["route"]:
        raise ValueError("Flight and sources differ")
    potential = sparse.kron(
        sparse.diags(radii.astype(float) ** 2), MISMATCH, format="csr"
    )
    formed = {}
    for arm in bank["arms"]:
        s = arm["charge"]
        a, b = FORMATION.starting_coins(route, s)
        initial = np.zeros(board.n * 4, complex)
        initial[4 * route[0] + 2 * a + b] = 1
        native_source = FORMATION.hamiltonian(board, route, s, 8)
        wave = FORMATION.flow(native_source, initial, [4.0])[0]
        source = FORMATION.reduce(wave)
        np.testing.assert_allclose(np.trace(source), 1, atol=3e-12, rtol=0)
        source /= np.trace(source)
        product = np.kron(*FORMATION.marginals(source))
        stocks = np.stack([source, product])
        expectation_change = float(np.trace((product - source) @ MISMATCH).real)
        observed = measured_mismatch(arm, data["shots"])
        if (
            abs(observed["product_minus_formed"] - expectation_change)
            > observed["simultaneous_99_radius"]
        ):
            raise ValueError("Sharp measurement of the introduced potential differs")
        probes = {}
        for q in (-1, 1):
            native = tidal.front.native.generator(board, route, q, False) / 2
            for path in paths:
                for i, j in zip(path[1:-1], path[2:], strict=True):
                    block = native[4 * i : 4 * i + 4, 4 * j : 4 * j + 4].toarray()
                    np.testing.assert_array_equal(block, -np.eye(4) / 2)
            readings = {}
            for g in (0, 1, -1):
                full, independent = formed_tides.release(
                    native + g * potential, launch, stocks, radii
                )
                # Rail r=2→3 minus rail r=1→2, separately at each end of the bridge.
                response = (independent - full)[[1, 3]]
                np.testing.assert_allclose(
                    response, -g * expectation_change, atol=2e-12, rtol=0
                )
                readings[str(g)] = response.tolist()
            probes[str(q)] = readings
        formed[str(s)] = {
            "mismatch_change_population": expectation_change,
            "mismatch_change_observed": observed,
            "relative_acceleration_change_by_probe_and_coupling": probes,
        }
    return {
        "model_change": "imposed nonlocal potential g times squared seam distance times the sharp-head mismatch projector",
        "formation_shots_per_endpoint": data["shots"],
        "source_counts_sha256": sha256(
            formed_tides.SOURCE_COUNTS.read_bytes()
        ).hexdigest(),
        "paths": paths,
        "formed": formed,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with RECEIPT.open("x") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        tidal.verify(result, json.loads(RECEIPT.read_text()))
    print(json.dumps(result["formed"], indent=2))
