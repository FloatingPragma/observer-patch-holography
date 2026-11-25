"""Prepare the strongest volume-matched joint modular source and test its directional response."""

import argparse
import json
from pathlib import Path

import numpy as np
from mpmath import mp
from scipy.sparse.linalg import expm_multiply

from matched_einstein_limit.source_implication.directional_stress import bridge
from matched_einstein_limit.source_implication.full_algebra import station
from matched_einstein_limit.source_implication.jacobi_screen import measure as jacobi

HERE = Path(__file__).resolve().parent
SHOTS = 2**26


def marginals(state):
    blocks = state.reshape(2, 2, 2, 2, 2)
    return np.trace(blocks, axis1=2, axis2=4), np.trace(blocks, axis1=1, axis2=3)


def joint_basis():
    axes = station.phase.AXES[1:]
    basis = np.zeros((18, 2, 4, 4), complex)
    for b in range(2):
        for a in range(3):
            for j in range(3):
                basis[9 * b + 3 * a + j, b] = np.kron(axes[a], axes[j]) / 2
    np.testing.assert_allclose(np.einsum("ibjk,lbkj->il", basis, basis), np.eye(18))
    return basis


def select_source(ref):
    basis = joint_basis()
    left, right = marginals(ref)
    product = np.array(
        [np.kron(a, z) / np.trace(a) for a, z in zip(left, right, strict=True)]
    )
    modular = station.log_blocks(product) - station.log_blocks(ref)
    coefficients = station.pair(basis, modular)
    menu = station.radar_menus()[6]
    ranges = station.pair(menu, ref)
    area, volume = station.range_gradient(ranges)
    range_basis = np.array([station.pair(menu, d) for d in basis])
    normal = range_basis @ volume
    assert np.linalg.norm(normal) > 1e-10
    projection = coefficients - normal * (normal @ coefficients) / (normal @ normal)
    length = np.linalg.norm(projection)
    unit = projection / length
    direction = np.einsum("i,ibjk->bjk", unit, basis)
    np.testing.assert_allclose(marginals(direction), 0, atol=1e-13)
    np.testing.assert_allclose(station.pair(modular, direction), length, atol=1e-12)
    np.testing.assert_allclose(
        station.pair(basis, station.log_blocks(product)), 0, atol=1e-12
    )
    dv = normal @ unit
    np.testing.assert_allclose(dv, 0, atol=1e-12)
    for h in (2.0**-i for i in range(7, 28)):
        endpoints = [ref + s * h * direction for s in (-1, 1)]
        if min(np.min(np.linalg.eigvalsh(x)) for x in endpoints) <= 0:
            continue
        preps = [station.preparation(x) for x in endpoints]
        if min(min(p["min_kept"], p["min_dephased"]) for p in preps) <= 0:
            continue
        assert max(p["max_error"] for p in preps) < 1e-11
        break
    else:
        raise ArithmeticError("No native two-sided preparation realized")
    for state in endpoints:
        np.testing.assert_allclose(marginals(state), marginals(ref), atol=1e-12, rtol=0)
    info = dict(
        basis_order="parity even, odd; left and right X,Y,Z; product divided by two",
        basis_coefficients=unit.tolist(),
        modular_slope=length,
        product_modular_slope=station.pair(-station.log_blocks(product), direction),
        matched_survey_volume_slope=dv,
        survey_area_slope=area @ range_basis.T @ unit,
        step=h,
        endpoint_preparations=preps,
        endpoint_parked_odd=[float(np.trace(x[1]).real) for x in endpoints],
        endpoint_survey_volume=[
            station.geometry(station.pair(menu, x))[1] for x in endpoints
        ],
    )
    return direction, endpoints, basis, coefficients, info


def meter_read(endpoints, basis, coefficients, slope, h, shots, seed):
    rng = np.random.default_rng(seed)
    banks = []
    estimates = []
    for state in endpoints:
        readings, counts = [], []
        for i, effect in enumerate(basis):
            b = i // 9
            mass = np.trace(state[b]).real
            mean = station.pair(effect, state)
            p = [mass / 2 + mean, mass / 2 - mean, 1 - mass]
            n = rng.multinomial(shots, p)
            counts.append(n.tolist())
            readings.append((n[0] - n[1]) / (2 * shots))
        banks.append(counts)
        estimates.append(coefficients @ readings)
    observed = estimates[1] - estimates[0]
    radius = (
        2 * np.sum(abs(coefficients)) * np.sqrt(np.log(2 * 36 / 0.01) / (2 * shots))
    )
    expected = 2 * h * slope
    assert abs(observed - expected) <= radius
    return dict(
        seed=seed,
        shots_per_setting=shots,
        outcomes=["plus", "minus", "other_parity"],
        minus_plus_banks=banks,
        expected_difference=expected,
        sampled_difference=observed,
        simultaneous_99pct_half_width=radius,
        positive_resolved=bool(observed > radius),
    )


def stress_effects(board, route, rings, q, time=4):
    target = jacobi.beam.indices(rings[0] + rings[1])
    effects = {k: [] for k in ("flux", "work")}
    for odd in (False, True):
        h = bridge.native.generator(board, route, q, odd)
        cells = bridge.operators(board, route, q, odd)
        edge = route[1:3]
        operators = dict(
            flux=1j * bridge.current(board, cells, edge) / 8,
            work=bridge.bridges(board, [edge]) / 2,
        )
        w = np.zeros((h.shape[0], len(target)))
        w[target] = np.eye(len(target))
        ht = -(0.5j * time) * h
        w = expm_multiply(ht, w, traceA=ht.diagonal().sum())
        for name, op in operators.items():
            obs = [v.conj().T @ (op @ v) for v in np.hsplit(w, 10)]
            effects[name].append(np.mean(obs, axis=0))
    return {k: np.array(v) for k, v in effects.items()}


def focusing(board, route, rings, q, ref, changes, time=4):
    effects = jacobi.signal_jets(board, route, rings, q, time)
    bearings = np.column_stack((np.ones(5), jacobi.beam.BEARINGS))

    def moments(x):
        return (
            bearings.T
            @ np.array([jacobi.beam.probabilities(e, x) for e in effects])
            @ bearings
        )

    baseline = moments(ref)
    jets = jacobi.matrix_jets(baseline)

    def observable(m):
        return -jacobi.optical_matrices(jacobi.matrix_jets(m))[2].trace()

    slopes = {}
    for name, d in changes.items():
        change = moments(d)
        delta = jacobi.response_jets(baseline, change)
        slope = -jacobi.deformation_response(jets, delta)[1].trace()
        # Near a caustic, subtracting two large double-precision focus readings
        # loses a small source slope. Evaluate the two measured-moment mixtures
        # in extended precision before taking their centered difference.
        with mp.workdps(48):
            to_mp = np.vectorize(lambda x: mp.mpf(float(x)), otypes=[object])
            center, displacement = to_mp(baseline), to_mp(change)
            step = mp.mpf("1e-12")
            high = observable(center + displacement * step)
            low = observable(center - displacement * step)
            centered = float((high - low) / (2 * step))
        np.testing.assert_allclose(centered, slope, rtol=1e-8, atol=1e-8)
        slopes[name] = slope
    return slopes


def run():
    ref, prior, _, _ = jacobi.beam.aperture.source_and_changes()
    d, endpoints, basis, coefficients, source = select_source(ref)
    source["independent_pauli_read"] = meter_read(
        endpoints,
        basis,
        coefficients,
        source["modular_slope"],
        source["step"],
        SHOTS,
        190811,
    )
    source["precision_pauli_read"] = meter_read(
        endpoints,
        basis,
        coefficients,
        source["modular_slope"],
        source["step"],
        2**36,
        190812,
    )
    changes = dict(ordinary=prior["ordinary"], published=prior["unpairing"], selected=d)
    formation = jacobi.beam.front.optics.pairing.formation
    board = formation.federation(1)
    route = formation.region(board)
    rings = jacobi.beam.screens(board, route)
    by_probe = {}
    published = json.loads(
        (Path(jacobi.__file__).parent / "source_response.json").read_text()
    )["species"]
    for q in (-1, 1):
        stress = stress_effects(board, route, rings, q)
        forward, reverse = [
            focusing(board, route, screen, q, ref, changes)
            for screen in (rings, rings[::-1])
        ]
        rows = {}
        for channel, sign, panel, key in (
            (
                "flux",
                -1,
                "forward_minus_reverse",
                "held_focusing_over_current_prediction",
            ),
            (
                "work",
                1,
                "forward_plus_reverse",
                "held_focusing_over_whole_seam_work_prediction",
            ),
        ):
            f = {name: forward[name] + sign * reverse[name] for name in changes}
            s = {
                name: station.pair(stress[channel], value)
                for name, value in changes.items()
            }
            gain = f["ordinary"] / s["ordinary"]
            lo, hi = map(float, published[str(q)][panel][key][1:-1].split(","))
            assert lo <= f["published"] / (gain * s["published"]) <= hi
            forecast = gain * s["selected"]
            rows[channel] = dict(
                ordinary_focusing=f["ordinary"],
                ordinary_stress=s["ordinary"],
                ordinary_gain=gain,
                selected_stress=s["selected"],
                predicted_focusing=forecast,
                observed_focusing=f["selected"],
                observed_over_predicted=f["selected"] / forecast,
                within_one_percent=bool(abs(f["selected"] / forecast - 1) <= 0.01),
            )
        by_probe[str(q)] = rows
    return dict(source=source, probes=by_probe)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = run()
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        saved = json.loads(path.read_text())
        for q, rows in out["probes"].items():
            for key, reading in rows.items():
                for k, v in reading.items():
                    np.testing.assert_allclose(
                        v, saved["probes"][q][key][k], atol=1e-8, rtol=1e-8
                    )
        np.testing.assert_allclose(
            out["source"]["basis_coefficients"], saved["source"]["basis_coefficients"]
        )
        for bank in ("independent_pauli_read", "precision_pauli_read"):
            assert (
                out["source"][bank]["minus_plus_banks"]
                == saved["source"][bank]["minus_plus_banks"]
            )
    print(
        json.dumps(
            dict(
                slope=out["source"]["modular_slope"],
                step=out["source"]["step"],
                meter=out["source"]["independent_pauli_read"]["sampled_difference"],
                radius=out["source"]["independent_pauli_read"][
                    "simultaneous_99pct_half_width"
                ],
                probes=out["probes"],
            ),
            indent=2,
        )
    )
