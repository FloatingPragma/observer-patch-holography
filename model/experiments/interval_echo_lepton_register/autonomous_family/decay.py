"""Weak-current spectral records of the same effective charged field."""

from math import log, pi, sqrt

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.sparse import block_diag, csr_matrix, diags
from scipy.sparse.linalg import expm_multiply


def phase_space(x):
    if x >= 1.0:
        return 0.0
    if x == 0.0:
        return 1.0
    if x > 0.9:
        # Use the positive integral when the polynomial terms nearly cancel.
        integral = quad(lambda u: (1 - u) ** 2 * u**2 / (x + (1 - x) * u), 0, 1)[0]
        return 12 * (1 - x) ** 5 * integral
    return 1 - 8 * x + 8 * x**3 - x**4 - 12 * x**2 * log(x)


def current_spin_trace(parent, daughter, pauli):
    """Compare explicit left-spin sums and the invariant form at a recoiling point."""
    q = (parent**2 - daughter**2) / (2 * parent + sqrt(2 * (parent**2 + daughter**2)))
    p_nu = np.array([q, q, 0, 0])
    p_anti = np.array([q, 0, 0, q])
    p_daughter = np.array([np.hypot(daughter, sqrt(2) * q), -q, 0, -q])
    p_parent = np.array([parent, 0, 0, 0])
    assert_allclose(p_daughter + p_nu + p_anti, p_parent, atol=1e-15)
    identity = np.eye(2)

    def left_density(momentum):
        return momentum[0] * identity - np.einsum("i,ijk->jk", momentum[1:], pauli)

    neutral = left_density(p_nu)
    antineutral = left_density(p_anti)
    charged = left_density(p_daughter)
    rest = left_density(p_parent) / 2
    bar = np.concatenate((identity[None, :, :], -pauli))
    metric = np.array([1, -1, -1, -1])
    trace = sum(
        metric[i]
        * metric[j]
        * np.trace(charged @ bar[i] @ antineutral @ bar[j])
        * np.trace(neutral @ bar[i] @ rest @ bar[j])
        for i in range(4)
        for j in range(4)
    )
    expected = 8 * np.dot(p_parent * metric, p_anti) * np.dot(p_daughter * metric, p_nu)
    assert expected > 0
    assert_allclose(trace, expected, atol=0, rtol=2e-13)
    return float(trace.real)


def first_pair_hamiltonian(energies, weights, bins, coupling):
    """One Hermitian operator for all parent tags and all one-pair daughter bins."""
    high = energies[-1]
    nodes, quadrature = np.polynomial.legendre.leggauss(8)
    reservoirs = []
    for daughter in energies:
        edges = np.linspace(daughter, 2 * high, bins + 1)
        midpoints = (edges[:-1] + edges[1:]) / 2
        half_width = (edges[1:] - edges[:-1]) / 2
        points = midpoints[:, None] + half_width[:, None] * nodes
        # A common ultraviolet form factor equals one at all three parent energies.
        u = np.clip(points / high - 1, 0, 1)
        taper = 1 - 3 * u**2 + 2 * u**3
        f = np.array([phase_space((daughter / energy) ** 2) for energy in points.flat])
        density = points**5 * f.reshape(points.shape) * taper**2 / (3072 * pi**4)
        reservoirs.append((midpoints, np.sqrt(half_width * (density @ quadrature))))

    blocks = []
    dim = 1 + len(energies) * bins
    for index, parent in enumerate(energies):
        diagonal = np.concatenate(([parent], *(grid for grid, _ in reservoirs)))
        amplitude = coupling * np.concatenate(
            [
                sqrt(weights[index] * weights[daughter]) * coefficient
                for daughter, (_, coefficient) in enumerate(reservoirs)
            ]
        )
        vertex = csr_matrix(
            (amplitude, (np.zeros(len(amplitude), int), np.arange(1, dim))),
            shape=(dim, dim),
        )
        blocks.append(diags(diagonal, format="csr") + vertex + vertex.T)
    return block_diag(blocks, format="csr")


def time_records(energies, weights, bins, coupling):
    dim = 1 + len(energies) * bins
    h = first_pair_hamiltonian(energies, weights, bins, coupling)
    parents = np.zeros((h.shape[0], len(energies)))
    parents[np.arange(len(energies)) * dim, np.arange(len(energies))] = 1
    times = np.linspace(0, 128, 5)
    traces = expm_multiply(-1j * h, parents, start=0, stop=128, num=len(times))
    banks = []
    for t, columns in zip(times, traces, strict=True):
        assert_allclose(np.sum(abs(columns) ** 2, axis=0), 1, atol=2e-12)
        energy_read = np.einsum("ij,ij->j", columns.conj(), h @ columns).real
        assert_allclose(energy_read, energies, atol=2e-12)
        for tag, vector in enumerate(columns.T):
            outside = np.r_[vector[: tag * dim], vector[(tag + 1) * dim :]]
            assert_allclose(outside, 0, atol=0)
            probability = abs(vector[tag * dim : (tag + 1) * dim]) ** 2
            daughters = probability[1:].reshape(len(energies), bins)
            banks.append(
                dict(
                    time=float(t),
                    parent=tag,
                    survival=float(probability[0]),
                    daughter_probabilities=daughters.sum(axis=1).tolist(),
                    energy_expectation=float(energy_read[tag]),
                )
            )
    # Retain the resolved bins for every channel at the last time.
    final = [
        abs(traces[-1, f * dim + 1 : (f + 1) * dim, f]) ** 2
        for f in range(len(energies))
    ]
    return dict(
        bins=bins,
        weak_coefficient=coupling,
        banks=banks,
        last_daughter_bins=[p.reshape(len(energies), bins).tolist() for p in final],
    )


def current_words(current):
    """Lift the two flavor currents to CAR words, with an idle positive charge."""
    flavor = np.eye(3, dtype=int)
    charges = np.array([-1] * 3 + [0] * 6 + [1])
    lepton = np.concatenate((flavor, flavor, -flavor, np.zeros((1, 3), int)))
    idle = 1 << 9

    def apply(state, mode, create):
        bit = 1 << mode
        if bool(state & bit) == create:
            return state, 0
        return state ^ bit, (-1) ** ((state & (bit - 1)).bit_count())

    records = []
    for parent in range(3):
        for daughter in range(3):
            state = idle | (1 << parent)
            amplitude = current[parent, parent] * current[daughter, daughter].conj()
            for mode, create in (
                (parent, False),
                (3 + parent, True),
                (6 + daughter, True),
                (daughter, True),
            ):
                state, sign = apply(state, mode, create)
                amplitude *= sign
            occupied = [i for i in range(10) if state >> i & 1]
            assert occupied == [daughter, 3 + parent, 6 + daughter, 9]
            assert charges[occupied].sum() == 0
            np.testing.assert_array_equal(lepton[occupied].sum(axis=0), flavor[parent])
            assert_allclose(abs(amplitude), 1, atol=2e-14)
            records.append(
                dict(
                    parent=parent,
                    charged=daughter,
                    neutral=parent,
                    antineutral=daughter,
                    car_sign=float(amplitude.real),
                )
            )
    return records


def weak_records(energies, current, pauli):
    """Use the measured internal current with the declared local left-current vertex."""
    weights = abs(current.diagonal()) ** 2
    assert_allclose(current, np.diag(current.diagonal()), atol=2e-14)
    widths = np.zeros((len(energies), len(energies)))
    spin_traces = []
    for parent, M in enumerate(energies):
        for daughter, m in enumerate(energies):
            if m < M:
                x = (m / M) ** 2
                invariant = quad(
                    lambda s: (1 - s) ** 2 * (s - x) ** 2 / s,
                    x,
                    1,
                    epsabs=1e-14,
                    epsrel=1e-12,
                )[0]
                assert_allclose(12 * invariant, phase_space(x), atol=2e-11)
                widths[parent, daughter] = (
                    weights[parent]
                    * weights[daughter]
                    * M**5
                    * invariant
                    / (128 * pi**3)
                )
                spin_traces.append(
                    dict(
                        parent=parent,
                        charged=daughter,
                        neutral=parent,
                        antineutral=daughter,
                        amplitude_squared=current_spin_trace(M, m, pauli),
                    )
                )

    finite = time_records(energies, weights, bins=96, coupling=1)
    finer = time_records(energies, weights, bins=192, coupling=1)
    cut = time_records(energies, weights, bins=96, coupling=0)
    cut_maximum = max(sum(row["daughter_probabilities"]) for row in cut["banks"])
    assert cut_maximum == 0.0
    maximum_drift = max(
        abs(
            np.array(left["daughter_probabilities"]) - right["daughter_probabilities"]
        ).max()
        for left, right in zip(finite["banks"], finer["banks"], strict=True)
    )
    return dict(
        branch_order="light, middle, heavy",
        widths_divided_by_coupling_squared=widths.tolist(),
        higher_to_middle_over_higher_to_light=float(widths[2, 1] / widths[2, 0]),
        signed_current_words=current_words(current),
        spin_trace_checks=spin_traces,
        original_start_probabilities=finite,
        finer_time_records=finer,
        maximum_bin_refinement_difference=float(maximum_drift),
        cut_daughter_probability=cut_maximum,
    )
