"""Independent physical transfer, counting, collective bound and radial controls."""

import json

import numpy as np
import pytest
from scipy import linalg

from . import measure


def interval(value):
    return tuple(json.loads(value))


def test_rail_contact_and_rotation_invariance():
    mid, targets, *_ = measure.front.inputs()
    mid = np.array(mid, dtype=complex)
    delta = np.array(targets["unpair"], dtype=complex)
    phase = measure.signed_transfers.phase
    axes = phase.AXES
    reads = phase.rail_analyzers(1)
    saved = json.loads((measure.HERE / "results.json").read_text())
    angle_axis = (axes[1] + 2 * axes[2] + 3 * axes[3]) / np.sqrt(14)
    r = measure.signed_transfers.rotation(angle_axis, 1.23)

    def coincidence(source, all_parities, strength):
        m = np.zeros(3)
        for b, x in enumerate(source):
            if b == 1 and not all_parities:
                continue
            a = phase.parity_blocks(1, strength, [1, -1][b])
            h = a[1:, 1:].transpose(0, 2, 1, 3).reshape(4, 4)
            left = linalg.expm(-1j * np.pi * h / 8).reshape(2, 2, 2, 2)[:, :, 0, :]
            right = phase.marked_gate(1, strength)[:, :, 0, :]
            embed = np.einsum("iax,jby->ijabxy", left, right).reshape(4, 4, 4)
            out = np.einsum("phx,xy,qhy->pq", embed, x, embed.conj())
            m += np.einsum("aij,ji->a", reads[13:16], out).real
        return 0.5 + m @ m / 6

    for all_parities in (False, True):
        key = "unconditional_native_rail"
        a, b = (
            (
                saved[key]["collar_first_derivative"],
                saved[key]["collar_quadratic_coefficient"],
            )
            if all_parities
            else (
                saved["coin_Z_right_random_axis_product_derivative"],
                saved["coin_product_quadratic"],
            )
        )
        base = coincidence(mid, all_parities, 8)
        for s in (1 / 8, 1.0):
            moved = mid + s * delta
            p = coincidence(r @ moved @ r.conj().T, all_parities, 8)
            np.testing.assert_allclose(
                p, coincidence(moved, all_parities, 8), atol=1e-15
            )
            lo, hi = np.array(interval(a)) * s + np.array(interval(b)) * s**2
            assert lo < p - base < hi < 0
        # Fixed detectors with zero contrast lose the right-head information.
        np.testing.assert_allclose(
            coincidence(mid, all_parities, 0),
            coincidence(mid + delta, all_parities, 0),
            atol=1e-15,
        )


def test_complete_counts_and_resource_certificates():
    saved = json.loads((measure.HERE / "results.json").read_text())
    assert measure.run() == saved
    for size, arms in saved["positive_source_resource_history"].items():
        for name, row in arms.items():
            n = row["samples_per_cohort"]
            assert all(sum(nums) == n for nums in row["counts"].values())
            d = 2 * row["decision_threshold"]
            assert np.exp(-n * d**2 / 4) < 0.01
            bounds = row["arbitrary_joint_read_raw_lottery"]
            if size == "1/5000000":
                bounds = row["arbitrary_joint_read_installed_angle"]
                assert bounds["tested_copies_required_at_least"] > (
                    2_250_000 if name == "ordinary" else 26_000_000
                )
            assert (bounds["tested_copies_required_at_least"] - 1) * np.log1p(
                interval(bounds["single_copy_chi_square_upper"])[1]
            ) <= np.log1p(0.98**2 / (0.01 * 0.99)) + 1e-12


def test_joint_readout_divergence_and_intensive_speed():
    chi = np.array([[0.7, 0.03j], [-0.03j, 0.3]])
    tau = np.array([[0.53, 0.19], [0.19, 0.47]])
    current = np.array([[0, -0.5j], [0.5j, 0]])

    def divergence(a, b):
        d = b - a
        return np.trace(d @ np.linalg.inv(a) @ d).real

    def times(x, n):
        result = np.ones((1, 1))
        for _ in range(n):
            result = np.kron(result, x)
        return result

    for n in (2, 3, 5):
        a, b = times(chi, n), times(tau, n)
        assert 1 + divergence(a, b) == pytest.approx((1 + divergence(chi, tau)) ** n)
        rng = np.random.default_rng(1953)
        u, _ = np.linalg.qr(rng.normal(size=a.shape) + 1j * rng.normal(size=a.shape))
        p, q = [(u.conj().T @ x @ u).diagonal().real for x in (a, b)]
        assert sum((p - q) ** 2 / p) < divergence(a, b)
        intensive = (
            sum(
                np.kron(np.eye(2**i), np.kron(current, np.eye(2 ** (n - i - 1))))
                for i in range(n)
            )
            / n
        )
        np.testing.assert_allclose(np.linalg.eigvalsh(intensive)[[0, -1]], [-0.5, 0.5])


def test_radial_against_independent_shell_distribution():
    radial = json.loads((measure.HERE / "radial.json").read_text())
    previous = json.loads((measure.HERE.parent / "results.json").read_text())
    ratios = {}
    for q, row in radial["charges"].items():
        charges = previous["tested_station"]["charges"]
        shells = charges[q]["shells"]["4.0"]
        pop = np.array([shells[str(k)]["population"] for k in (1, 2, 3, 4)])
        lo, hi = interval(row["expected_capped_velocity_distance"])
        assert lo < 2 * sum(pop) < hi
        slopes = {}
        for name in ("ordinary", "pairing"):
            score = [shells[str(k)]["matched_log_slopes"][name] for k in (1, 2, 3, 4)]
            slopes[name] = 2 * pop @ score
            lo, hi = interval(row["expected_distance_slopes"][name])
            assert lo < slopes[name] < hi < 0
        ratio = slopes["pairing"] / slopes["ordinary"]
        lo, hi = ratios[q] = interval(row["collar_over_ordinary"])
        assert 0 < lo < ratio < hi
    a, b = ratios["-1"][1], ratios["1"][0]
    assert a < b and (b - a) / (b + a) > 0.14136
    assert interval(radial["common_response_relative_error_at_least"])[
        0
    ] == pytest.approx((b - a) / (b + a))
