"""Read phase-sensitive matter currents invisible to a static slow-sector modular frame."""

import argparse
import hashlib
from itertools import combinations
import json
from types import SimpleNamespace

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from model.experiments.interval_echo_null_balance.directional import native_analyzer
from model.experiments.interval_echo_null_balance.regional import verify
from model.experiments.interval_echo_space.shell_source import prepare
from . import modular


def traceless_basis():
    """Frobenius-orthonormal real coordinates for normalized Hermitian variations."""
    identity = np.eye(4)
    out = []
    for j in range(1, 4):
        diag = identity[:j].sum(axis=0) - j * identity[j]
        out.append(np.diag(diag) / np.linalg.norm(diag))
    for i, j in combinations(range(4), 2):
        entry = np.outer(identity[i], identity[j])
        out.extend(
            ((entry + entry.T) / np.sqrt(2), 1j * (entry - entry.T) / np.sqrt(2))
        )
    return np.array(out)


def signed_effect(mapped, analyzer, edge):
    """Pull a native signed endpoint read back through the four launch columns."""
    rows = analyzer @ mapped[list(edge)]
    return rows.conj().T @ (np.array([1.0, -1.0])[:, None] * rows)


def sample_endpoints(state, analyzer, edge, rng, shots=modular.SHOTS):
    p = abs(analyzer @ state[list(edge)]) ** 2
    chances = np.r_[p, max(0.0, 1 - p.sum())]
    counts = rng.multinomial(shots, chances / chances.sum())
    estimate = (counts[0] - counts[1]) / shots
    variance = (p.sum() - (p[0] - p[1]) ** 2) / shots
    return float(estimate), float(variance), counts.tolist()


def neighbor_energy(sphere, owner):
    """Allocate each signed-contrast bond equally to its two endpoints."""
    ports = set(range(12 * owner, 12 * (owner + 1)))
    q = np.zeros((sphere.n, sphere.n))
    seams = []
    for a, b in sphere.seams:
        share = (int(a in ports) + int(b in ports)) / 2
        if share:
            q[a, a] += share * (candidate.KAPPA + 0.5)
            q[b, b] += share * (candidate.KAPPA + 0.5)
            q[a, b] = q[b, a] = -0.5 * share
            seams.append((a, b))
    for port in ports:
        q[port, port] += candidate.KAPPA * (6 - sphere.degree[port])
    return q, seams


def measure_energy(state, instrument, rng):
    q, seams = instrument
    diagonal = np.diag(q)
    active = np.flatnonzero(diagonal)
    weights = diagonal[active]
    p = abs(state[active]) ** 2
    n = 2**22
    counts = rng.multinomial(n, np.r_[p, max(0.0, 1 - p.sum())])
    estimate = weights @ counts[:-1] / n
    variance = weights**2 @ p - (weights @ p) ** 2
    histograms = [counts.tolist()]
    for a, b in seams:
        plus, minus = (abs(state[a] + sign * state[b]) ** 2 / 2 for sign in (1, -1))
        hist = rng.multinomial(n, [plus, minus, max(0.0, 1 - plus - minus)])
        histograms.append(hist.tolist())
        estimate += q[a, b] * (hist[0] - hist[1]) / n
        variance += q[a, b] ** 2 * (plus + minus - (plus - minus) ** 2)
    return float(estimate), float(variance / n), histograms


def source_frame(directions):
    # Fix the triplet coordinates using existing physical port labels.
    axes = []
    for direction in directions:
        rem = direction.copy()
        for prior in axes:
            rem -= (prior @ rem) * prior
        if np.linalg.norm(rem) > 1e-8:
            axes.append(rem / np.linalg.norm(rem))
        if len(axes) == 3:
            break
    return np.column_stack(
        (np.ones(12) / np.sqrt(12), directions @ np.array(axes).T / 2)
    )


def measure():
    car, sphere = cell(), federation(0)
    directions = modular.cap.chord_coordinates()
    basis = source_frame(directions)
    generator = candidate.response(car)
    energies = np.diag(basis.T @ generator @ basis)
    np.testing.assert_allclose(generator @ basis, basis * energies, atol=3e-13)
    gap = energies[1] - energies[0]
    np.testing.assert_allclose(energies[1:], energies[1], atol=3e-13)
    phases = np.array([1.0, -1j, -1j, -1j])
    parent, _, _ = modular.common_reference()
    stars, charges, _ = modular.regional_operators(parent)
    initial = basis.T @ charges @ basis
    delayed = phases.conj()[None, :, None] * initial * phases[None, None, :]
    analyzer = native_analyzer()
    edges = car.seams
    wedges = np.array([np.cross(basis[a, 1:], basis[b, 1:]) for a, b in edges])
    triples = list(combinations(range(len(edges)), 3))
    minimum = np.linalg.svd(wedges[np.array(triples)], compute_uv=False)[:, -1]
    picks = list(triples[int(np.argmax(np.round(minimum, 12)))])
    chosen = [edges[i] for i in picks]
    calibrators = np.array(
        [*initial, *delayed, *(signed_effect(basis, analyzer, e) for e in chosen)]
    )
    axes = traceless_basis()

    def features(effects):
        return np.einsum("pij,kji->pk", effects, axes).real

    design = features(calibrators)
    spectrum = np.linalg.svd(design, compute_uv=False)
    inverse = np.linalg.pinv(design)
    counts = [
        int(np.linalg.matrix_rank(design[:i], tol=1e-9)) for i in (12, 24, 25, 26, 27)
    ]
    if counts != [9, 12, 13, 14, 15]:
        raise ArithmeticError(f"Unexpected ranks {counts}")
    carrier = SimpleNamespace(car=sphere, active=sphere.seams)
    launch_map = np.zeros((sphere.n, 4))
    launch_map[:12] = basis
    # Aperture opens for the held experiment. Propagate the same four possible
    # input columns across the full free board, including every external seam.
    full_h = candidate.generator(sphere)
    held_exposure = 2.0
    flight = expm_multiply(-1j * held_exposure * full_h, launch_map)
    held_edges = [e for e in edges if e not in chosen] + [
        e for e in sphere.external if e[0] < 12
    ]
    neighbors = [b // 12 for a, b in sphere.external if a < 12]
    energy_meters = [neighbor_energy(sphere, owner) for owner in neighbors]
    all_region_q = sum(neighbor_energy(sphere, k)[0] for k in range(sphere.cells))
    np.testing.assert_allclose(all_region_q, full_h.toarray(), atol=1e-12)
    exterior = flight[12:].conj().T @ flight[12:]
    held_ops = np.array(
        [
            *(signed_effect(flight, analyzer, e) for e in held_edges),
            exterior,
            *(flight.conj().T @ q @ flight for q, _ in energy_meters),
        ]
    )
    kernel = features(held_ops) @ inverse
    baseline = np.trace(calibrators, axis1=1, axis2=2).real / 4
    held_baseline = np.trace(held_ops, axis1=1, axis2=2).real / 4
    sources = (
        np.array([0, 1, 1j, 0]) / np.sqrt(2),
        np.array([0, 1, -1j, 0]) / np.sqrt(2),
        np.array([1 + 1j, 2, -1j, 1]) / np.sqrt(8),
    )

    def calibration(prepared, cohort, key, shots):
        obs, var, digests, seam_counts = [], [], [], []
        for delay, phase in enumerate((np.ones(4), phases)):
            local = basis @ (phase * (basis.T @ prepared[:12]))
            for pole in range(12):
                val, noise, digest = modular.read_star(
                    local,
                    stars[pole],
                    charges[pole],
                    np.random.default_rng([key, cohort, delay, pole]),
                    shots,
                )
                obs.append(val)
                var.append(noise)
                digests.append(digest)
        for i, edge in enumerate(chosen):
            val, noise, records = sample_endpoints(
                prepared,
                analyzer,
                edge,
                np.random.default_rng([key + 1, cohort, i]),
                shots,
            )
            obs.append(val)
            var.append(noise)
            seam_counts.append(records)
        return obs, var, digests, seam_counts

    results = []
    for cohort, source in enumerate(sources):
        prepared, gates = prepare(carrier, launch_map @ source, root=0)
        obs, var, digests, seam_counts = calibration(
            prepared, cohort, 76511, modular.SHOTS
        )
        prediction = held_baseline + kernel @ (obs - baseline)
        true_initial = np.einsum("i,pij,j->p", source.conj(), calibrators, source).real
        ideal_prediction = held_baseline + kernel @ (true_initial - baseline)
        exact = np.einsum("i,pij,j->p", source.conj(), held_ops, source).real
        np.testing.assert_allclose(ideal_prediction, exact, rtol=0, atol=1e-11)
        state = expm_multiply(-1j * held_exposure * full_h, prepared)
        reads = [
            sample_endpoints(
                state, analyzer, e, np.random.default_rng([76513, cohort, i])
            )
            for i, e in enumerate(held_edges)
        ]
        p = float(np.sum(abs(state[12:]) ** 2))
        hits = int(np.random.default_rng([76514, cohort]).binomial(modular.SHOTS, p))
        reads.append(
            (
                hits / modular.SHOTS,
                p * (1 - p) / modular.SHOTS,
                [hits, modular.SHOTS - hits],
            )
        )
        for j, meter in enumerate(energy_meters):
            reads.append(
                measure_energy(state, meter, np.random.default_rng([76515, cohort, j]))
            )
        observed = np.array([s[0] for s in reads])
        se = np.sqrt(kernel**2 @ var + [s[1] for s in reads])
        precise_obs, precise_var, precise_digests, precise_seams = calibration(
            prepared, cohort, 76521, 2**26
        )
        energy_kernel = kernel[-len(neighbors) :]
        energy_prediction = held_baseline[-len(neighbors) :] + energy_kernel @ (
            precise_obs - baseline
        )
        precise_reads = [
            measure_energy(state, meter, np.random.default_rng([76523, cohort, j]))
            for j, meter in enumerate(energy_meters)
        ]
        energy_observed = np.array([read[0] for read in precise_reads])
        energy_var = np.array([read[1] for read in precise_reads])
        energy_se = np.sqrt(energy_var + energy_kernel**2 @ precise_var)
        precision = dict(
            predicted=energy_prediction.tolist(),
            measured=energy_observed.tolist(),
            measured_se=np.sqrt(energy_var).tolist(),
            difference_se=energy_se.tolist(),
            max_standardized_error=float(
                max(abs(energy_observed - energy_prediction) / energy_se)
            ),
            count_digest=hashlib.sha256(
                json.dumps(
                    [precise_digests, precise_seams, *[s[2] for s in precise_reads]]
                ).encode()
            ).hexdigest(),
        )
        results.append(
            dict(
                preparation=gates,
                precision_neighbor_energy=precision,
                charge_count_digests=digests,
                signed_calibration_counts=seam_counts,
                signed_calibration_exact=true_initial[24:].tolist(),
                neighboring_energy_exact=exact[-len(neighbors) :].tolist(),
                held_prediction=prediction.tolist(),
                held_measurement=observed.tolist(),
                held_standard_error=se.tolist(),
                held_counts=[s[2] for s in reads[: len(held_edges) + 1]],
                neighbor_energy_count_digests=[
                    hashlib.sha256(json.dumps(s[2]).encode()).hexdigest()
                    for s in reads[-len(neighbors) :]
                ],
                largest_held_standard_error_units=float(
                    max(abs(observed - prediction) / se)
                ),
            )
        )
    # All current-free analysis settings have real effects; the stationary
    # pure triplet pair is conjugate at every service exposure.
    pair = [basis @ s for s in sources[:2]]
    difference = np.outer(pair[0], pair[0].conj()) - np.outer(pair[1], pair[1].conj())
    np.testing.assert_allclose(difference.real, 0, rtol=0, atol=1e-13)
    signed = np.array([signed_effect(basis, analyzer, e) for e in edges])
    contrasts = np.array(
        [np.einsum("i,pij,j->p", s.conj(), signed, s).real for s in sources[:2]]
    )
    strongest = int(np.argmax(np.round(abs(contrasts[0] - contrasts[1]), 12)))
    return dict(
        cells=sphere.cells,
        shots_per_setting=modular.SHOTS,
        selected_current_seams=chosen,
        held_seams=held_edges,
        subsequent_held_reads=[
            "outside birth cell",
            "neighbor energies in the listed order",
        ],
        neighboring_cells=neighbors,
        shots_per_energy_setting=2**22,
        precision_calibration_shots=2**26,
        beat_gap=float(gap),
        quarter_beat_exposure=float(np.pi / (2 * gap)),
        held_exposure=held_exposure,
        trace_zero_ranks=counts,
        inverse_smallest_singular_value=float(spectrum[-1]),
        stationary_pair=dict(
            current_seam=edges[strongest],
            signed_responses=contrasts[:, strongest].tolist(),
            real_density_difference_max=float(np.max(abs(difference.real))),
        ),
        sources=[
            "circular triplet +",
            "circular triplet -",
            "complex scalar-triplet superposition",
        ],
        rows=results,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = modular.HERE / "matter_current.json"
    # JSON roundtrip normalizes native port tuple/list representations.
    data = json.loads(json.dumps(measure(), allow_nan=False))
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        verify(json.loads(path.read_text()), data)
    print(
        data["trace_zero_ranks"],
        data["stationary_pair"],
        [row["largest_held_standard_error_units"] for row in data["rows"]],
    )
