"""Reconstruct the finite angular map independently from the published detector banks."""

from itertools import combinations_with_replacement, product
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "records.json").read_text())


def test_actual_apertures_and_endpoint_diamonds():
    e = DATA["endpoints"]
    delta = np.array(e["deltas"])
    np.testing.assert_array_equal(delta, -delta[e["opposites"]])
    assert len(e["pair_disagreements"]) == 15 and all(
        n == 0 for _, _, n in e["pair_disagreements"]
    )
    assert len(e["broken_owners"]) == 8
    assert len(e["address_witnesses"]) == 125 and all(
        i != v for i, v in enumerate(e["address_witnesses"])
    )
    for word, (start, finish) in zip(e["aperture_words"], e["arrivals"], strict=True):
        v = delta[word].sum(axis=0)
        assert start == (-v % 5) @ [25, 5, 1] and finish == (v % 5) @ [25, 5, 1]
    assert len({pair[1] for pair in e["arrivals"]}) == 32


def multiply(a, b):
    # Integer coefficients in Z[phi], phi²=phi+1.
    return (a[0] * b[0] + a[1] * b[1], a[0] * b[1] + a[1] * b[0] + a[1] * b[1])


def test_exact_orbit_designs_through_degree_five():
    e = DATA["endpoints"]
    raw = np.array(
        [
            [(int(k), 0) if abs(k) <= 1 else (0, int(k // 3)) for k in v]
            for v in e["deltas"]
        ]
    )
    all_rays = np.stack([raw[word].sum(axis=0) for word in e["aperture_words"]])
    for orbit in (all_rays[:12], all_rays[12:]):
        keys = {tuple(row.ravel()) for row in orbit}
        assert all(tuple(-row.ravel()) in keys for row in orbit)
        radius = sum(
            (np.array(multiply(a, a)) for a in orbit[0]), start=np.zeros(2, dtype=int)
        )
        for row in orbit:
            np.testing.assert_array_equal(
                sum(np.array(multiply(a, a)) for a in row), radius
            )
        for degree in (2, 4):
            radial = radius if degree == 2 else np.array(multiply(radius, radius))
            for indices in combinations_with_replacement(range(3), degree):
                terms = []
                for row in orbit:
                    value = (1, 0)
                    for k in indices:
                        value = multiply(value, row[k])
                    terms.append(value)
                # Uniform-sphere second/fourth tensor using integer identities.
                total = np.sum(terms, axis=0)
                if any(indices.count(k) % 2 for k in indices):
                    np.testing.assert_array_equal(total, [0, 0])
                else:
                    divisor = (
                        3 if degree == 2 else (5 if len(set(indices)) == 1 else 15)
                    )
                    np.testing.assert_array_equal(divisor * total, len(orbit) * radial)


def test_polarimeter_all_preparations_and_basis_settings():
    d = np.array(DATA["directions"])
    z = np.array([0.0, 0.0, 1.0])

    def axes(n):
        phi = np.arctan2(n[1], n[0])
        t = np.array([-np.sin(phi), np.cos(phi), 0])
        return np.column_stack((np.cross(t, n), t, n))

    frames = [axes(n) for n in d]
    h, v = np.eye(3)[:, :2].T
    diagonal = (h + v) / np.sqrt(2)
    right = (h + 1j * v) / np.sqrt(2)
    inputs = np.stack((h, v, diagonal, right), axis=1)
    read = [
        np.eye(3),
        np.column_stack((diagonal, (-h + v) / np.sqrt(2), z)),
        np.column_stack((right, (-h + 1j * v) / np.sqrt(2), z)),
    ]
    with np.load(HERE / "detectors.npz") as bank:
        recorded = bank["tomography"] / DATA["polarization_probability_denominator"]
        for i, j in product(range(32), repeat=2):
            state = frames[j].T @ frames[i] @ inputs
            expected = np.stack(
                [abs(basis.conj().T @ state).T ** 2 for basis in read], axis=1
            )
            np.testing.assert_allclose(recorded[i, j], expected, rtol=0, atol=5.3e-10)


def test_raw_gains_counts_and_conditional_normalization():
    d = np.array(DATA["directions"])
    a = DATA["angular"]
    dots = d @ d.T
    h6 = (
        (231 * dots[:12] ** 6 - 315 * dots[:12] ** 4 + 105 * dots[:12] ** 2 - 5) / 16
    ).mean(axis=0)
    c = np.r_[np.full(12, 1 / 12), np.full(20, -1 / 20)]
    q = (3 * dots[0] ** 2 - 1) / 2
    with np.load(HERE / "detectors.npz") as bank:
        z = (
            bank["tomography"][:, :, :2, 0]
            / DATA["polarization_probability_denominator"]
        )
        collision = z[:, :, :, :2].sum(axis=-1).mean(axis=2) / 32
        count = bank["counts"]
    assert np.all(count.sum(axis=-1) == a["originals_per_stratum"])
    rates = count[:, :, :, :2].sum(axis=-1) / a["originals_per_stratum"]
    assert (
        np.max(abs(rates - z[:, :, :, :2].sum(axis=-1))) < a["simultaneous_max_error"]
    )
    sampled = c @ (rates.mean(axis=2) / 32) @ h6 / (c @ h6)
    assert abs(sampled - a["counted_raw_six_gain"]) < 1e-11
    assert abs(sampled) < a["simultaneous_gain_halfwidth"]
    for name, bypass in (("ballistic", 1), ("half_bypass", 0.5), ("collision", 0)):
        k = np.eye(32) * bypass + (1 - bypass) * collision
        survival = k.sum(axis=1).mean()
        stored = a["transfer"][name]
        assert abs(stored["survival"] - survival) < 1e-9
        for key, input_, score in (("six", c, h6), ("two", q / 32, q)):
            g = input_ @ k @ score / (input_ @ score)
            np.testing.assert_allclose(
                stored[key], [g, g / survival], rtol=0, atol=1e-9
            )
    stokes = np.array(a["quadrupole_stokes"])
    np.testing.assert_allclose(stokes[:, 0], q / 15, atol=2e-12)
    np.testing.assert_allclose(
        np.hypot(stokes[:, 1], stokes[:, 2]), (1 - dots[0] ** 2) / 10, atol=2e-12
    )
    np.testing.assert_allclose(stokes[:, 3], 0, atol=2e-12)


def test_second_encounter_from_first_encounter_tomography():
    directions = np.array(DATA["directions"])
    e = DATA["endpoints"]
    q = (3 * (directions @ directions[0]) ** 2 - 1) / 2
    with np.load(HERE / "detectors.npz") as bank:
        density_reads = (
            bank["tomography"][:, :, :2].mean(axis=2)
            / DATA["polarization_probability_denominator"]
        )
        saved = bank["second_kernels"] / DATA["polarization_probability_denominator"]
        arrivals = bank["second_stations"]
    result = np.zeros_like(saved)
    for j, n in enumerate(directions):
        phi = np.arctan2(n[1], n[0])
        v = np.array([-np.sin(phi), np.cos(phi), 0])
        transverse = np.stack((np.cross(v, n), v))
        components = directions @ transverse.T
        # Reconstruct the surviving two-dimensional density from first-hit
        # tomography. The second real normal tests the real off-diagonal part.
        h, vv = density_reads[:, j, 0, :2].T
        u = density_reads[:, j, 1, 0] - density_reads[:, j, 1, 1]
        trace = h + vv
        x, y = components.T
        result[0] += (
            trace[:, None]
            - np.outer(h, x * x)
            - np.outer(vv, y * y)
            - np.outer(u, x * y)
        ) / 32**2
        result[1] += np.outer(trace, 1 - (x * x + y * y) / 2) / 32**2
        middle = np.array(np.unravel_index(e["arrivals"][j][1], (5, 5, 5)))
        delta = np.array(e["deltas"])
        expected = [
            ((middle + delta[word].sum(axis=0)) % 5) @ [25, 5, 1]
            for word in e["aperture_words"]
        ]
        np.testing.assert_array_equal(arrivals[j], expected)
    np.testing.assert_allclose(result, saved, rtol=0, atol=5.5e-10)
    for matrix, reference in zip(result, (0.07, 0.01), strict=True):
        np.testing.assert_allclose(matrix.sum(axis=1), 4 / 9, rtol=0, atol=1e-9)
        assert abs((q @ matrix @ q) / (q @ q) / (4 / 9) - reference) < 1e-9
