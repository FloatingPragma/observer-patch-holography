"""Join the free collar formation, its independent source counts and the release read."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from . import tidal

HERE = tidal.HERE
RECEIPT = HERE / "formed_tides.json"
FORMATION = tidal.front.optics.ranging.formation
SOURCE_COUNTS = Path(FORMATION.__file__).with_name("results.json")


def source_read(arm, shots):
    meter = arm["modular"]
    n0 = np.asarray(meter["calibration_counts"])
    held = next(row for row in meter["held_reads"] if row["exposure"] == 4)
    draws = np.asarray(held["joint_and_product_counts"])
    if n0.sum() != shots or not np.all(draws.sum(axis=(1, 2)) == shots):
        raise ValueError("Incomplete source counts")
    p = n0 / shots
    weights = np.log(np.outer(p.sum(axis=1), p.sum(axis=0)) / p)
    contrast = np.sum((draws[1] - draws[0]) * weights) / shots
    np.testing.assert_allclose(contrast, -held["measured_nats"], atol=1e-13, rtol=0)
    return float(contrast), held["radius_99_nats"]


def release(h, launch, states, radii):
    # Assemble the probability-acceleration observable before tracing the head
    # density. Use an arbitrary scalar shift to lower cancellation error.
    h = h - h.diagonal().mean() * sparse.eye(h.shape[0], format="csr")
    records = []
    for ray in launch.T:
        psi = np.kron(ray[:, None], np.eye(4))
        v = h @ psi
        rv = np.repeat(radii, 4)[:, None]
        second = psi.conj().T @ (rv * (h @ v))
        observable = 2 * v.conj().T @ (rv * v) - second - second.conj().T
        records.append(np.einsum("ij,sji->s", observable, states).real)
    return tidal.two_rays(np.stack(records, axis=1))


def run():
    raw = json.loads(SOURCE_COUNTS.read_text())
    bank = next(e for e in raw["experiments"] if e["cells"] == 80)
    board = FORMATION.federation(1)
    route = FORMATION.region(board)
    if list(route) != bank["route"]:
        raise ValueError("Source and position experiment use different marked contacts")
    radii, paths, launch = tidal.chart(board, route)
    result = {}
    for arm in bank["arms"]:
        s = arm["charge"]
        coin = FORMATION.starting_coins(route, s)
        h = tidal.front.native.generator(board, route, s, False) / 2
        initial = np.zeros(h.shape[0])
        initial[4 * route[0] + 2 * coin[0] + coin[1]] = 1
        wave = expm_multiply(-4j * h, initial, traceA=-4j * h.diagonal().sum())
        f = FORMATION.reduce(wave)
        np.testing.assert_allclose(np.trace(f), 1, atol=3e-12, rtol=0)
        f = f / np.trace(f)
        p = np.kron(*FORMATION.marginals(f))
        states = np.stack([f, p])
        eigenvalues = np.linalg.eigvalsh(states)
        if eigenvalues.min() <= 0:
            raise ValueError("A source logarithm needs a common faithful support")
        score = np.diff(tidal.transfer.meter.head.station.log_blocks(states), axis=0)[0]
        jeffreys = float(np.trace((p - f) @ score).real)
        observed, bound = source_read(arm, raw["shots"])
        if jeffreys < observed - bound:
            raise ValueError(
                "The computed full score violates the measured coarse lower bound"
            )
        reads = {}
        for q in (-1, 1):
            native = tidal.front.native.generator(board, route, q, False)
            tidal.check_contacts(native, route)
            full, independent = release(native / 2, launch, states, radii)
            reads[str(q)] = (independent - full).tolist()
        result[str(s)] = dict(
            least_formed_eigenvalue=float(eigenvalues[0].min()),
            quantum_score_change=jeffreys,
            sharp_score_change=observed,
            sharp_radius=bound,
            full_score_lower_99=observed - bound,
            radial_acceleration_changes=reads,
        )
    return dict(
        source_sha256=sha256(SOURCE_COUNTS.read_bytes()).hexdigest(),
        paths=paths,
        formed=result,
    )


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
