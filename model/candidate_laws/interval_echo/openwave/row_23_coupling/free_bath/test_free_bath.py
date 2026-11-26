"""Independent microscopic corner Hamiltonians and unconditional signed source reads."""

import gzip
import json
from itertools import combinations, product
import numpy as np
import pytest
from scipy.linalg import expm
from model.executable_model.carrier import cell
from model.executable_model import federation
from . import measure, native


def microscopic(board, pins, holes, k):
    free = [p for p in range(board.n) if p not in pins]
    basis = list(combinations(free, holes))
    index = {w: i for i, w in enumerate(basis)}
    h0, h1 = native.LAW.local_term_parts()
    raw = h0 + k * h1
    h = np.zeros((len(basis), len(basis)))
    for idx, word in enumerate(basis):
        levels = np.ones(board.n, dtype=int)
        levels[list(word)] = 0
        for p, v in pins.items():
            levels[p] = v
        h[idx, idx] = k * ((6 - board.degree) @ (levels**2))
        for a, b in board.seams:
            col = native.LAW.SRC.index((levels[a], levels[b]))
            h[idx, idx] += raw[col, col]
            if a not in pins and b not in pins and levels[a] != levels[b]:
                alt = tuple(
                    sorted((set(word) - {a, b}) | ({a} if levels[b] == 0 else {b}))
                )
                swap = native.LAW.SRC.index((levels[b], levels[a]))
                h[index[alt], idx] += raw[swap, col]
    np.testing.assert_allclose(h, h.T, atol=2e-13)
    return h, basis


@pytest.mark.parametrize("k", [0, 8])
def test_complete_guarded_hamiltonians_and_spectral_moments(k):
    b = cell()
    home = 0
    d = native.distances(b, home)
    for shell in (1, 2, 3):
        probe = int(np.flatnonzero(d == shell)[0])
        mats = {}
        for a, c in product((-1, 1), repeat=2):
            h, _ = microscopic(b, {home: 1 + a, probe: 1 + c}, 5, k)
            mats[a, c] = h
        signed = sum(a * c * h for (a, c), h in mats.items()) / 4
        gain = 2 * k + 0.25
        wanted = -gain if shell == 1 else 0
        np.testing.assert_allclose(signed, wanted * np.eye(len(h)), atol=1e-11)
        if shell > 1:

            def variance(matrix):
                # Full spectral variance of the uniform density on the 252-state block.
                return (
                    np.sum(matrix * matrix) / len(matrix)
                    - (np.trace(matrix) / len(matrix)) ** 2
                )

            direct = sum(a * c * variance(z) for (a, c), z in mats.items()) / 4
            a = set(b.neighbours[home])
            c = set(b.neighbours[probe])
            m = b.n - 2
            predicted = (
                2
                * gain**2
                * (
                    (0.25 + 0.25 / (m - 1)) * len(a & c)
                    - 0.25 * len(a) * len(c) / (m - 1)
                )
            )
            assert abs(direct - predicted) < 2e-10


def test_dilute_susceptibility_against_full_free_evolution():
    b = cell()
    for k in (0, 8):
        for load in (0, 2):
            dense, basis = microscopic(b, {0: load}, 1, k)
            vac, _ = microscopic(b, {0: load}, 0, k)
            got, free = native.vacancy_operator(b, k, {0: load})
            np.testing.assert_array_equal(free, [w[0] for w in basis])
            np.testing.assert_allclose(
                got.toarray(), dense - vac[0, 0] * np.eye(len(dense)), atol=2e-12
            )
            init = native.uniform_control(b, 0)
            u = expm(-1j * 0.913 * dense)
            np.testing.assert_allclose(
                abs(native.propagation(got, init, 0.913)) ** 2,
                abs(u @ init) ** 2,
                atol=3e-12,
            )
            np.testing.assert_allclose(
                u @ u.conj().T / len(free), np.eye(len(free)) / len(free), atol=2e-12
            )


def test_all_original_starts_and_radius_window():
    d = json.loads(gzip.decompress(measure.OUTPUT.read_bytes()))
    for size, b in d["boards"].items():
        for age, row in b["source"].items():
            assert sum(row["profile"]["counts"]) == measure.SHOTS
            assert sum(row["susceptibility_address_counts"]) == measure.SHOTS
            low, high = row["sampler_rms_interval"]
            assert low <= row["native_rms"] <= high
            for key, z in row["employed"].items():
                assert sum(z["clock"]["counts"]) == measure.SHOTS
                assert abs(z["signed_residual"]) <= z["residual_sampler_halfwidth"]
                if z["in_radius_window"]:
                    lo, hi = z["ratio_to_calibrated_contact"]
                    assert 0.9 < lo <= 1 <= hi < 1.1
            joint = row["joint_susceptibility"]
            for key, read in joint["cohorts"].items():
                assert sum(read["counts"]) == measure.SUS_SHOTS
                start = 0
                for location, width in zip(
                    joint["locations"], joint["shell_counts_per_location"], strict=True
                ):
                    mass = row["profile"]["probability"][location]
                    probs = b["susceptibility"][str(location)]["cohorts"][key][
                        "probability"
                    ]
                    np.testing.assert_allclose(
                        read["probability"][start : start + width],
                        mass * np.array(probs),
                        atol=3e-13,
                    )
                    start += width
                assert abs(sum(read["probability"]) - 1) < 3e-12
        for state in b["susceptibility"].values():
            for key, z in state["readings"].items():
                assert (
                    abs(z["first_shell_half_difference"] - z["native_half_difference"])
                    <= z["sampler_halfwidth"]
                )
                if key.endswith("_mixed"):
                    assert z["native_half_difference"] == 0
            stat = state["readings"]["8_1_coherent"]
            assert (
                abs(stat["first_shell_half_difference"])
                > 10 * stat["sampler_halfwidth"]
            )
        for r in b["fluctuation"].values():
            for z in r["readings"].values():
                assert (
                    abs(z["signed_energy_variance"] - z["exact"])
                    <= z["sampler_halfwidth"]
                )
    for level in range(3):
        b = federation(level)
        assert max(b.degree) == 6 and min(b.degree) == 5
        e = native.contacts(8).diagonal().reshape(3, 3)
        np.testing.assert_allclose(
            (e[2] - e[0]) / 2, 16.25 * np.array([1, 0, -1]), atol=1e-13
        )
