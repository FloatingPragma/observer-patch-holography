"""Read ring energy, transport currents and charge-probe coarse Ricci curvature on one field."""

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

from ..field_geometry_equilibrium.measure import charge_sector, ring_laplacian, search


HERE = Path(__file__).resolve().parent
THETA = np.pi / 8
SHOTS = 2**32
SEED = 517812
SIGNS = (-1, 1)


def components(patch, words):
    index = {int(w): i for i, w in enumerate(words)}
    pairs = [[] for _ in patch.faces]
    for i, w in enumerate(words):
        for f, end in patch.ring_moves(int(w)):
            j = index[end]
            if i < j:
                pairs[f].append((i, j))
    operators = []
    for edges in pairs:
        a, b = np.array(edges).T
        data = np.repeat([0.5, 0.5, -0.5, -0.5], len(a))
        matrix = coo_matrix(
            (data, (np.r_[a, b, a, b], np.r_[a, b, b, a])), shape=(len(words),) * 2
        ).tocsr()
        assert (matrix @ matrix - matrix).nnz == 0
        operators.append(matrix)
    return operators


def stress(ops, h, states):
    action = np.array([p @ states for p in ops])
    np.testing.assert_allclose(action.sum(axis=0), h @ states, atol=1e-13, rtol=0)
    energies = np.einsum("wi,fwi->if", states.conj(), action).real
    currents = -2 * np.einsum("awi,bwi->iab", action.conj(), action).imag
    time_derivative = (
        2 * np.einsum("wi,fwi->if", (-1j * (h @ states)).conj(), action).real
    )
    np.testing.assert_allclose(
        time_derivative, currents.sum(axis=1), atol=2e-15, rtol=0
    )
    np.testing.assert_allclose(
        energies, np.einsum("fwi,fwi->if", action.conj(), action).real, atol=2e-14
    )
    return energies, currents


def kernels(patch, words, vacuum, excitations):
    calibration = np.zeros((2, 12, 12))
    vacuum_flight = np.zeros_like(calibration)
    sources = np.zeros((2, 30, 3, 12, 12))
    c, s = np.cos(THETA), np.sin(THETA)
    for sign_id, sign in enumerate(SIGNS):
        for port in range(12):
            patch.compensator(port, sign)
            ports, h = charge_sector(patch, words, port, sign)
            start = np.zeros((len(ports), 31), complex)
            start[: len(words)] = np.column_stack([vacuum, excitations])
            series = expm_multiply(
                -1j * h,
                start,
                start=0,
                stop=3 / 16,
                num=4,
                traceA=-1j * h.diagonal().sum(),
            )
            calibration[sign_id, port] = np.bincount(
                ports, abs(series[2][:, 0]) ** 2, minlength=12
            )
            quiet, excited = np.split(series[3], [1], axis=1)
            for dest in range(12):
                mask = ports == dest
                a, b = quiet[mask], excited[mask]
                vacuum_flight[sign_id, port, dest] = np.sum(abs(a) ** 2)
                low = np.sum(abs(c * a - 1j * s * b) ** 2, axis=0)
                high = np.sum(abs(c * a + 1j * s * b) ** 2, axis=0)
                sources[sign_id, :, 0, port, dest] = low
                sources[sign_id, :, 1, port, dest] = high
                sources[sign_id, :, 2, port, dest] = (low + high) / 2
    for array in (calibration, vacuum_flight, sources):
        np.testing.assert_allclose(array.sum(axis=-1), 1, atol=1e-12, rtol=0)
    return calibration, vacuum_flight, sources


def metric(vacuum_cal):
    profile = vacuum_cal.mean(axis=0)
    return np.sum(abs(profile[:, None, :] - profile[None, :, :]), axis=2) / 2


class Transport:
    def __init__(self, d):
        self.d = d
        self.pairs = list(combinations(range(len(d)), 2))
        delta = np.array([np.eye(len(d))[i] - np.eye(len(d))[j] for i, j in self.pairs])
        self.constraints = np.r_[delta, -delta][:, 1:]
        weights = np.array([d[i, j] for i, j in self.pairs])
        if weights.min() <= 0:
            raise ValueError("Vacuum distances do not separate enrolled ports")
        self.limits = np.r_[weights, weights]

    def curvature(self, kernel, edges):
        result = []
        for i, j in edges:
            difference = kernel[i] - kernel[j]
            solution = linprog(
                -difference[1:],
                A_ub=self.constraints,
                b_ub=self.limits,
                bounds=(None, None),
                method="highs-ds",
            )
            if not solution.success:
                raise ArithmeticError(solution.message)
            result.append(1 + solution.fun / self.d[i, j])
        return np.array(result)


def signs_and_pooled(table):
    return np.concatenate([table, table.mean(axis=0, keepdims=True)])


def multinomial_rows(table, rng):
    rows = np.array(
        [rng.multinomial(SHOTS, row / row.sum()) for row in table.reshape(-1, 12)]
    )
    return rows.reshape(table.shape)


def scored_law(matrix, waves, rng):
    values, vec = eigh(matrix.toarray())
    levels = np.round(values, 12)
    spectrum, idx = np.unique(levels, return_inverse=True)
    prob = np.array([np.bincount(idx, abs(vec.conj().T @ w) ** 2) for w in waves])
    trials = [rng.multinomial(SHOTS, p / p.sum()) for p in prob]
    return dict(
        eigenvalues=spectrum.tolist(),
        probabilities=prob.tolist(),
        counts=[c.tolist() for c in trials],
        expected=(prob @ spectrum).tolist(),
        measured=(np.array(trials) @ spectrum / SHOTS).tolist(),
    )


def run():
    patch = search.Patch()
    words = np.array(patch.field_component())
    g = np.ones(len(words)) / np.sqrt(len(words))
    h = ring_laplacian(patch, words)
    ops = components(patch, words)
    np.testing.assert_allclose([p @ g for p in ops], 0, atol=1e-14)
    z = np.column_stack([np.where(words & patch.bit[e], -1, 1) for e in patch.edges])
    v = expm_multiply(-0.5j * h, g[:, None] * z)
    waves = [
        np.cos(THETA) * g[:, None] + 1j * sign * np.sin(THETA) * v for sign in SIGNS
    ]
    low, high = [stress(ops, h, state) for state in waves]
    np.testing.assert_allclose(low[0], high[0], atol=1e-14, rtol=0)
    np.testing.assert_allclose(low[1], high[1], atol=1e-14, rtol=0)
    energy, currents = low
    calibration, quiet, kernels_by_sign = kernels(patch, words, g, v)
    distance = metric(calibration)
    transport = Transport(distance)
    edges = list(patch.edges)
    face_edges = np.array(
        [
            [edges.index(tuple(sorted((a, b)))) for a, b in search.around(face)]
            for face in patch.faces
        ]
    )
    baseline = np.array(
        [transport.curvature(k, edges) for k in signs_and_pooled(quiet)]
    )
    kernels_with_pooled = signs_and_pooled(kernels_by_sign)
    curvature = np.array(
        [
            [
                [transport.curvature(arm, edges) for arm in branches]
                for branches in by_head
            ]
            for by_head in kernels_with_pooled
        ]
    )
    delta = curvature - baseline[:, None, None, :]
    faces = delta[:, :, 2, face_edges].mean(axis=-1)
    hold = np.arange(30) % 3 == 2
    train_e = energy[~hold]
    coupling = np.sum(train_e * faces[2, ~hold]) / np.sum(train_e**2)
    forecast = coupling * energy[hold]
    residuals = [
        np.linalg.norm(response[hold] - forecast) / np.linalg.norm(response[hold])
        for response in faces
    ]
    correlations = [
        np.corrcoef(energy[hold].ravel(), response[hold].ravel())[0, 1]
        for response in faces
    ]
    effects = curvature[:, :, 1] - curvature[:, :, 0]
    contender = np.where(hold[:, None], abs(effects[2]), -1.0)
    head, ray = (int(x) for x in np.unravel_index(contender.argmax(), contender.shape))
    best_edge = [edges[ray]]
    rng = np.random.default_rng(SEED)
    cal_counts = multinomial_rows(calibration, rng)
    d_read = metric(cal_counts / SHOTS)
    reader = Transport(d_read)
    quiet_counts = multinomial_rows(quiet, rng)
    held_counts = multinomial_rows(kernels_by_sign[:, head], rng)
    held_k = np.array(
        [
            [reader.curvature(arm, best_edge)[0] for arm in arms]
            for arms in signs_and_pooled(held_counts) / SHOTS
        ]
    )
    # A union over 120 empirical twelve-category rows, accounting for both metric and kernel reads.
    tv_bound = np.sqrt(np.log(120 * (2**12 - 2) / 0.01) / (2 * SHOTS))
    edge_price = d_read[best_edge[0]]
    # Each transported separation changes by <= 2 eps for the two marginal prices and
    # <= 2 eps from the calibration. The shared denominator has price error <= 2 eps.
    err = (8 * tv_bound + 2 * tv_bound * np.max(abs(held_k[:, 1] - held_k[:, 0]))) / (
        edge_price - 2 * tv_bound
    )
    er = [rng.binomial(SHOTS, np.clip(e[head], 0, 1)) for e in (low[0], high[0])]
    j_pair = np.unravel_index(abs(currents[head]).argmax(), (20, 20))
    a, b = map(int, j_pair)
    j_read = scored_law(
        1j * (ops[a] @ ops[b] - ops[b] @ ops[a]), [w[:, head] for w in waves], rng
    )
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        field_words=len(words),
        head_order=edges,
        face_order=patch.faces,
        angle=THETA,
        energy_unit="ring-projector eigenvalue",
        energy_profiles=energy.tolist(),
        energy_currents=currents.tolist(),
        vacuum_calibration_kernels=calibration.tolist(),
        vacuum_deadline_kernels=quiet.tolist(),
        vacuum_distances=distance.tolist(),
        vacuum_ray_curvatures=baseline.tolist(),
        ray_curvatures=curvature.tolist(),
        source_arms=["negative phase", "positive phase", "incoherent"],
        probe_order=["vacancy", "excess", "charge-even"],
        phase_flights=kernels_by_sign.tolist(),
        local_fit=dict(
            trained_head_indices=np.flatnonzero(~hold).tolist(),
            held_head_indices=np.flatnonzero(hold).tolist(),
            coefficient=float(coupling),
            relative_held_errors=np.array(residuals).tolist(),
            correlations=correlations,
            held_curvature_forecast=forecast.tolist(),
            face_response=faces.tolist(),
        ),
        selected=dict(
            head_index=head,
            head=edges[head],
            ray_index=ray,
            ray=edges[ray],
            phase_differences=effects[:, head, ray].tolist(),
            measured_phase_differences=(held_k[:, 1] - held_k[:, 0]).tolist(),
            conservative_simultaneous_99_percent_radius=float(err),
            all_ring_click_counts=np.array(er).tolist(),
            current_faces=[a, b],
            current=j_read,
            curvature_estimates=held_k.tolist(),
            held_vacuum_calibration_counts=cal_counts.tolist(),
            held_vacuum_deadline_counts=quiet_counts.tolist(),
            held_source_counts=held_counts.tolist(),
        ),
        field_dwell_exposure=0.5,
        calibration_exposure=1 / 8,
        ray_exposure=3 / 16,
        shots=SHOTS,
        seed=SEED,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    encoded = json.dumps(
        result, default=lambda v: v.tolist(), allow_nan=False, separators=(",", ":")
    )
    with args.output.open("x") as f:
        f.write(encoded + "\n")
