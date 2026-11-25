"""Fit a collar record at one appointment, then test its predictions on separate appointments."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space

from . import measure, observer_records

HERE = Path(__file__).resolve().parent


def instrument_rows(ref, ordinary, basis, board, route, rings, time):
    variations = dict(ordinary=ordinary, **{str(i): d for i, d in enumerate(basis)})
    records = {}
    for q in (-1, 1):
        stresses = measure.stress_effects(board, route, rings, q, time)
        forward, reverse = [
            measure.focusing(board, route, aperture, q, ref, variations, time)
            for aperture in (rings, rings[::-1])
        ]
        for name, sign in [("flux", -1), ("work", 1)]:
            ordinary_read = measure.station.pair(stresses[name], ordinary)
            ordinary_f = forward["ordinary"] + sign * reverse["ordinary"]
            gain = ordinary_f / ordinary_read
            observed = np.array(
                [forward[str(i)] + sign * reverse[str(i)] for i in range(len(basis))]
            )
            predicted = gain * measure.station.pair(basis, stresses[name])
            records[f"{q}_{name}"] = (observed, predicted)
    return records


def run(write):
    ref, prior, *_ = measure.jacobi.beam.aperture.source_and_changes()
    _, _, basis, coeff, _ = measure.select_source(ref)
    m = measure.station.radar_menus()[6]
    dv = measure.station.range_gradient(measure.station.pair(m, ref))[1]
    gradient = np.array([dv @ measure.station.pair(m, b) for b in basis])
    api = measure.jacobi.beam.front.optics.pairing.formation
    board = api.federation(1)
    route = api.region(board)
    rings = measure.jacobi.beam.screens(board, route)

    def panel(time):
        return instrument_rows(ref, prior["ordinary"], basis, board, route, rings, time)

    training = panel(4)
    constraints = np.array(
        [gradient] + [observed - forecast for observed, forecast in training.values()]
    )
    constraints /= np.linalg.norm(constraints, axis=1)[:, None]
    free = null_space(constraints, rcond=1e-10)
    projected = free @ (free.T @ coeff)
    length = np.linalg.norm(projected)
    if length < 1e-10:
        raise ArithmeticError("Training leaves no detectable modular variation")
    vector = projected / length
    np.testing.assert_allclose(constraints @ vector, 0, atol=2e-12)
    direction = np.einsum("i,ibjk->bjk", vector, basis)
    for h in (2.0**-i for i in range(7, 28)):
        states = [ref + s * h * direction for s in (-1, 1)]
        if min(np.linalg.eigvalsh(s).min() for s in states) <= 0:
            continue
        preps = [measure.station.preparation(s) for s in states]
        if min(min(p["min_kept"], p["min_dephased"]) for p in preps) > 0:
            break
    else:
        raise ArithmeticError("Training fit not preparable")
    for x, p in zip(states, preps, strict=True):
        assert p["max_error"] < 1e-11
        np.testing.assert_allclose(
            measure.marginals(x), measure.marginals(ref), atol=1e-12
        )

    def apply(rows):
        out = {}
        for k, (o, p) in rows.items():
            observed, predicted = float(o @ vector), float(p @ vector)
            assert abs(predicted) > 1e-8, "Unresolved predicted response"
            out[k] = dict(
                measured=observed,
                predicted=predicted,
                ratio=observed / predicted,
                within_one_percent=bool(abs(observed / predicted - 1) <= 0.01),
            )
        return out

    # The source, marginals and both endpoint recipes are fixed before held acquisition.
    result = dict(
        training=apply(training),
        joint_reading_coefficients=vector.tolist(),
        step=h,
        endpoint_preparations=preps,
        singular_values=np.linalg.svd(constraints, compute_uv=False).tolist(),
        fitted_constraints_rank=18 - free.shape[1],
        feasible_dimension=free.shape[1],
        modular_slope=length,
        pauli_count=measure.meter_read(states, basis, coeff, length, h, 2**36, 190816),
        held={},
    )
    for t, seed in [(3.5, 190814), (4.5, 190815)]:
        expected = apply(panel(t))
        counts = observer_records.collect(direction, h, t, seed)
        decoded = observer_records.decode(counts)
        counts["decoded"] = decoded
        path = HERE / f"held_{t}.json"
        if write:
            path.write_text(json.dumps(counts, indent=2, allow_nan=False) + "\n")
        else:
            assert decoded == json.loads(path.read_text())["decoded"]
        result["held"][str(t)] = dict(
            population=expected, recorded_vacancy_flux=decoded
        )
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    p.add_argument("--counts-only", action="store_true")
    args = p.parse_args()
    if args.counts_only:
        for t in (3.5, 4.5):
            file = json.loads((HERE / f"held_{t}.json").read_text())
            assert observer_records.decode(file) == file["decoded"]
            print(
                t,
                file["decoded"]["predicted_odd_focusing"],
                file["decoded"]["measured_odd_focusing"],
            )
    else:
        result = run(args.write)
        path = HERE / "search.json"
        if args.write:
            path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        else:
            saved = json.loads(path.read_text())
            np.testing.assert_allclose(
                result["joint_reading_coefficients"],
                saved["joint_reading_coefficients"],
                rtol=1e-9,
                atol=1e-11,
            )
            for t, held in result["held"].items():
                for k, r in held["population"].items():
                    np.testing.assert_allclose(
                        r["ratio"],
                        saved["held"][t]["population"][k]["ratio"],
                        rtol=1e-8,
                    )
        print(
            json.dumps(
                {k: v["recorded_vacancy_flux"] for k, v in result["held"].items()},
                indent=2,
            )
        )
