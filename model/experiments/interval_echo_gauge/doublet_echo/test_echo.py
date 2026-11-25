"""Checks independent of the forward-only amplitude reduction and phase inference."""

import hashlib
import json
import numpy as np
from scipy.sparse.linalg import expm_multiply
from scipy.stats import poisson
from .certificate import certificate
from .measure import (
    HERE,
    SHOTS,
    J,
    apparatus,
    endpoint,
    source_words,
    stable_binomial,
    interval,
)


def test_native_round_trip_by_separate_exponentials():
    board, sector, addresses, ports, home, frame = apparatus()
    generator = sector.hamiltonian(8)
    initial = np.zeros((sector.P, 2), complex)
    initial[home] = frame
    first = expm_multiply(-8j * generator, initial)
    local = first[ports]
    middle = np.zeros_like(initial)
    middle[ports[:2]] = J @ local[:2]
    middle[ports[2:]] = 1j * J @ local[2:]
    final = expm_multiply(-8j * generator, middle)
    calculated = frame.T @ final[home]
    expected = (np.linalg.det(local[:2]) + 1j * np.linalg.det(local[2:])) * J
    np.testing.assert_allclose(calculated, expected, atol=2e-20, rtol=1e-8)
    _, transform, alpha = source_words(board, addresses[0], frame)
    tag_minus = (endpoint(addresses[0])["tag"] - np.eye(len(frame))) / 2
    np.testing.assert_allclose(tag_minus, alpha * frame @ frame.T, atol=2e-13)
    assert np.linalg.det(transform) != 0


def test_exact_material_path_traces():
    for k, imag in ((8, 87), (4, 75)):
        assert certificate(k) == {
            "left_nonzero": (17, 85),
            "right_nonzero": (17, 85),
            "relative_imaginary_nonzero": (37, imag),
        }


def test_rare_event_bank_sub_ulp_probability():
    n = np.full(25000, SHOTS, dtype=np.int64)
    draws = stable_binomial(np.random.default_rng(24817), n, 1e-17)
    assert abs(draws.mean() - SHOTS * 1e-17) < 0.12


def test_every_bank_and_held_phase():
    data = json.loads((HERE / "records.json").read_text())
    assert (
        data["protocol_sha256"]
        == hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    )
    for arm in data["arms"]:
        for run in arm["runs"]:
            n = run["batches_per_setting"] * data["batch_size"]
            for kind, rows in run["complete_counts"].items():
                for row in rows:
                    assert sum(row["counts"]) == n
                    assert sum(row["first_batch"]) == SHOTS
                    if kind == "which":
                        assert (
                            list(
                                map(sum, zip(*row["labelled_path_counts"], strict=True))
                            )
                            == row["counts"][:4]
                        )
                    if kind != "ordinary_swap":
                        assert row["counts"][1] == 0
            assert all(v["compatible"] for v in run["forecast"]["held"])
            if run["exposure"] == 32:
                assert run["forecast"]["dark_vs_bright_resolved"]
                assert run["forecast"]["quadratures_resolved"]
                assert (
                    sum(
                        row["counts"][1]
                        for row in run["complete_counts"]["ordinary_swap"]
                    )
                    > 0
                )


def test_every_calibration_and_held_tail_independently():
    # Le Cam bounds the variation distance to Poisson(np) by n*p².
    # At k=0 the exact binomial tail is (1-p)^n.
    data = json.loads((HERE / "records.json").read_text())
    tail = 0.01 / 4096
    for arm in data["arms"]:
        for run in arm["runs"]:
            n = 6 * SHOTS * run["batches_per_setting"]
            for i in range(6):
                k = run["forecast"]["pooled_target_counts"][f"phase{i}"]
                lo, hi = interval(k, n)
                assert lo <= k / n <= hi
                if k:
                    assert poisson.sf(k - 1, n * lo) + n * lo**2 < tail
                    assert poisson.cdf(k, n * hi) + n * hi**2 < tail
                else:
                    assert np.exp(n * np.log1p(-hi)) <= tail * (1 + 1e-13)
