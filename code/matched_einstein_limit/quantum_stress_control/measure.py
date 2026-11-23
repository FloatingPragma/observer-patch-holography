"""Measure a single Maxwell quantum's integrated stress and retarded propagation on the same box."""
import argparse
from hashlib import sha256
from itertools import combinations, product
import json
import os
from pathlib import Path
import sys

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from model.candidate_laws.interval_echo.openwave.photon_identity.quantum_transfer import measure as maxwell

PROTOCOL = json.loads((HERE / "protocol.json").read_text())


def face(box, point, axes):
    def advance(x, a):
        return tuple(v + (i == a) for i, v in enumerate(x))
    a, b = axes
    route = [point, advance(point, a), advance(advance(point, a), b), advance(point, b), point]
    return [box.lookup[x, y] for x, y in zip(route[:-1], route[1:], strict=True)]


def loop_vector(box, path):
    return sum(sign * box.q[e] for e, sign in path)


def physical_rails(box):
    electric = box.frame * np.sqrt(box.freq / 2)
    edge_axes = np.argmax(np.array([np.subtract(head, tail) for tail, head in box.edges]), axis=1)
    magnetic, face_axes = [], []
    for point in product(range(box.n + 1), repeat=3):
        for a, b in combinations(range(3), 2):
            if point[a] < box.n and point[b] < box.n:
                magnetic.append(loop_vector(box, face(box, point, (a, b))))
                face_axes.append(3 - a - b)
    return ((electric, edge_axes), (np.array(magnetic), np.array(face_axes)))


def squared_outcomes(max_occ, initial):
    creation = np.diag(np.sqrt(np.arange(1, max_occ + 1)), k=-1)
    values, vectors = eigh(creation + creation.T)
    return values**2, vectors[initial]**2


def field_stress(n, rails, state, omega, rng):
    # Each local quadrature can be aligned with one canonical mode. The component
    # orthogonal to it carries one spectator quantum, leaving cap-1 slots.
    cap = maxwell.CAP
    vac_values, vac_probs = squared_outcomes(cap, 0)
    one_values, one_probs = squared_outcomes(cap, 1)
    spectator_values, spectator_probs = squared_outcomes(cap - 1, 0)
    shots = PROTOCOL["enrollments_per_quadrature"]
    energy = float(abs(state)**2 @ omega)
    result = []
    for units in PROTOCOL["source_stress_times_in_side_units"]:
        exact, observed, variance, records = np.zeros(3), np.zeros(3), np.zeros(3), []
        phase = state * np.exp(-1j * omega * (units * n))
        for vectors, axes in rails:
            norms = np.einsum("ij,ij->i", vectors, vectors)
            overlaps = abs(vectors @ phase)**2 / norms
            assert np.min(overlaps) >= 0 and np.max(overlaps) <= 1 + 1e-12
            rail_groups = []
            for axis in range(3):
                # Weighted sufficient statistics, since individual rail normalizations differ.
                sq_sum, sq2_sum = np.zeros(2), np.zeros(2)
                for norm, u in zip(norms[axes == axis], overlaps[axes == axis], strict=True):
                    n_aligned = rng.binomial(shots, np.clip(u, 0, 1))
                    vac_counts = rng.multinomial(shots, vac_probs)
                    one_counts = rng.multinomial(n_aligned, one_probs)
                    spec_counts = rng.multinomial(shots - n_aligned, spectator_probs)
                    vac_sq = norm * vac_values
                    one_sq = norm * one_values
                    spec_sq = norm * spectator_values
                    sq_sum += [vac_counts @ vac_sq, one_counts @ one_sq + spec_counts @ spec_sq]
                    sq2_sum += [vac_counts @ vac_sq**2,
                                one_counts @ one_sq**2 + spec_counts @ spec_sq**2]
                    means = np.array([vac_counts @ vac_sq, one_counts @ one_sq + spec_counts @ spec_sq]) / shots
                    moments = np.array([vac_counts @ vac_sq**2, one_counts @ one_sq**2 + spec_counts @ spec_sq**2]) / shots
                    variance[axis] += (moments - means**2).sum() / (shots - 1)
                exact[axis] += float(2 * (norms[axes == axis] @ overlaps[axes == axis]))
                observed[axis] += (sq_sum[1] - sq_sum[0]) / shots
                rail_groups.append({"rail_count": int(np.count_nonzero(axes == axis)),
                                    "sum_of_squared_outcomes_vacuum_photon": sq_sum.tolist(),
                                    "sum_of_fourth_power_outcomes_vacuum_photon": sq2_sum.tolist()})
            records.append(rail_groups)
        trace, trace_se = float(observed.sum() / 2), float(np.sqrt(variance.sum()) / 2)
        np.testing.assert_allclose(exact.sum() / 2, energy, atol=1e-10, rtol=0)
        assert trace > 5 * trace_se
        result.append({"time_in_side_units": units, "field_square_groups_electric_magnetic": records,
                       "normal_ordered_diagonal_stress": (energy - exact).tolist(),
                       "record_diagonal_stress": (trace - observed).tolist(),
                       "record_diagonal_standard_errors": [trace_se] * 3,
                       "record_trace": trace, "trace_standard_error": trace_se})
    return result


def compressed_states(directions):
    frame, levels, operators = maxwell.coordinate_bank(directions)
    vac = np.zeros(len(levels), complex)
    vac[0] = 1
    one = np.zeros_like(vac)
    one[levels == 1] = frame.conj().T @ directions["state"]
    np.testing.assert_allclose(np.linalg.norm(one), 1, atol=1e-12)
    return (vac, one), operators


def susceptibilities(source, receiver, state):
    (vac, one), bank = compressed_states({"source": source, "receiver": receiver, "state": state})
    s, r = bank["source"], bank["receiver"]
    retarded = 1j * (r @ s - s @ r)
    return np.array([np.vdot(x, retarded @ x).real for x in (vac, one)])


def finite_kicks(box, state, path, receiver, rng):
    # Use the parent's ordered addressed Wilson phases for finite amplitude.
    directions = {e: box.q[e] for e, _ in path}
    directions.update({"state": state, "receiver": receiver})
    states, bank = compressed_states(directions)
    values, read = eigh(bank["receiver"].toarray())
    shots = PROTOCOL["enrollments_per_kick_sign"]
    results = []
    for eps in PROTOCOL["finite_kicks"]:
        pairs, errors, slopes, records = [], [], [], []
        for initial in states:
            means, sample, variances, draw = [], [], [], []
            for sign in (1, -1):
                evolved = initial.copy()
                for edge, orient in path:
                    evolved = expm_multiply((1j * sign * eps * orient) * bank[edge], evolved, traceA=0)
                probs = abs(read.conj().T @ evolved)**2
                counts = rng.multinomial(shots, probs / probs.sum())
                obs = counts @ values / shots
                means.append(probs @ values)
                sample.append(obs)
                variances.append((counts @ values**2 / shots - obs**2) / (shots - 1))
                draw.append(counts.tolist())
            slopes.append((means[0] - means[1]) / (2 * eps))
            pairs.append((sample[0] - sample[1]) / (2 * eps))
            errors.append(np.sqrt(sum(variances)) / (2 * eps))
            records.append(draw)
        results.append({"kick": eps, "slopes_vacuum_photon": slopes, "record_slopes_vacuum_photon": pairs,
                        "record_slope_standard_errors": errors, "receiver_outcome_values": values.tolist(),
                        "counts_state_by_kick_sign_by_outcome": records})
    return results


def measure(n, rng):
    box = maxwell.Box(n)
    _, levels, rays, waves = maxwell.circuit_bank(box, 0)
    emitted, state = maxwell.source_record(1, levels, rays, waves)
    spectrum = abs(state)**2
    counts = rng.multinomial(PROTOCOL["enrollments_per_spectrometer"], spectrum / spectrum.sum())
    energy = float(spectrum @ box.freq)
    rails = physical_rails(box)
    stress = field_stress(n, rails, state, box.freq, rng)
    channels, finite = [], []
    for axes in combinations(range(3), 2):
        p = (n // 2, n // 2, 0)
        path = face(box, p, axes)
        src = loop_vector(box, path)
        for z in (1, n - 1):
            dst = loop_vector(box, face(box, (*p[:2], z), axes))
            for units in PROTOCOL["times_in_side_units"]:
                t = n * units
                recv = dst * np.exp(1j * box.freq * t)
                pair = susceptibilities(src, recv, state)
                formula = -2 * np.vdot(recv, src).imag
                ctrl = -2 * np.vdot(dst * np.exp(1j * 1.2 * box.freq * t), src).imag
                np.testing.assert_allclose(pair, formula, atol=1e-11, rtol=0)
                channels.append({"face_axes": list(axes), "receiver_z": z, "time": t,
                                 "vacuum_photon": pair.tolist(), "changed_clock": float(ctrl)})
                if axes == (0, 2) and z == 1 and units == 0.5:
                    finite = finite_kicks(box, state, path, recv, rng)
    max_diff = max(abs(a["vacuum_photon"][0] - a["vacuum_photon"][1]) for a in channels)
    control = max(abs(a["vacuum_photon"][0] - a["changed_clock"]) for a in channels)
    assert control > 1e-5 and finite
    return {"side": n, "transverse_modes": len(box.freq), "cochain_hash": box.hash,
            "pointer_zero_one_quantum_probability": float(emitted[0, 1]),
            "one_quantum_probability_given_pointer_zero": float(emitted[0, 1] / emitted[0].sum()),
            "spectrometer_frequencies": box.freq.tolist(), "spectrometer_photon_counts": counts.tolist(),
            "spectrometer_vacuum_zero_energy_counts": PROTOCOL["enrollments_per_spectrometer"],
            "spectrometer_record_mean_energy": float(counts @ box.freq / counts.sum()),
            "normal_ordered_energy": energy, "stress_records": stress, "weak_optical_channels": channels,
            "weak_response_maximum_difference": float(max_diff), "maximum_clock_control_difference": control,
            "adjacent_xz_finite_kicks_at_half_side_time": finite}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    args = parser.parse_args()
    rng = np.random.default_rng(PROTOCOL["seed"])
    runs = []
    for n in PROTOCOL["sides"]:
        runs.append(measure(n, rng))
    packet = {"protocol_sha256": sha256((HERE / "protocol.json").read_bytes()).hexdigest(),
              "parent_code_sha256": sha256(Path(maxwell.__file__).read_bytes()).hexdigest(),
              "runs": runs}
    args.output.write_text(json.dumps(packet, indent=2) + "\n")
    for run in runs:
        print(run["side"], run["normal_ordered_energy"], run["weak_response_maximum_difference"],
              run["maximum_clock_control_difference"])


if __name__ == "__main__":
    main()
