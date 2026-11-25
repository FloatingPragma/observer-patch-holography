from fractions import Fraction
import json
from math import comb, isqrt
import random

import numpy as np

from . import counts, measure, receiver as rx


def test_binomial_envelope_and_native_lottery():
    n = 140
    for p in (Fraction(1, 2), Fraction(9, 19)):
        m = ((n + 1) * p).__floor__()
        w = 1 + isqrt(n)
        for k in range(n + 1):
            j = ((k - m) if k >= m else m - 1 - k) // w
            relative = Fraction(comb(n, k), comb(n, m)) * (p / (1 - p)) ** (k - m)
            assert relative <= Fraction(1, 2 ** max(0, j - 1))
        rng = random.Random(4914)
        sample = [counts.random_binomial(n, p, rng) for _ in range(16)]
        assert len(set(sample)) > 1 and all(0 <= x <= n for x in sample)
    a = rx.emit_atoms()
    assert np.array_equal(a[0] + a[6], a[24] + a[30])
    pa, r = (
        np.concatenate(([np.eye(4) / 2], rx.station.iso.basis(False))),
        a.reshape(36, 4, 4),
    )
    frame = np.einsum("aij,bji->ab", pa, r)
    assert np.linalg.svd(frame, compute_uv=False)[-1] > 0.99999


def test_blind_counts_reconstruct_and_round_without_source_labels():
    rx.mp.dps = 80
    rx.iv.dps = 62
    saved = json.loads((rx.HERE / "results.json").read_text())
    rows, labels = measure.regional.read()
    for name, systems in saved["material_banks"].items():
        for i, bank in enumerate(systems):
            keys = saved["endpoint_material_rows"][i]
            assert len(set(map(tuple, keys))) == 31
            assert all(
                labels[a][j][1] in ("shift", "domain_cut", "image_cut") for a, j in keys
            )
            n = bank["trials"]
            observed = [
                Fraction(2 * Fraction(c) * (2 * k - n), n)
                for c, k in zip(bank["score_ranges"], bank["success"])
            ]
            numbers, stock, ver = rx.calibrate(
                keys, rows, observed, (-1 if i == 0 else 1) * rx.SHIFT
            )
            assert sum(numbers) == rx.TOTAL and min(numbers) > int(0.003 * rx.TOTAL)
            assert np.max(abs(stock.sum(axis=0).trace() - 1)) < 1e-14
            assert float(ver["implementation_and_inverse_error_upper"]) < 1e-10
            assert min(np.linalg.eigvalsh(stock).flat) > 0.017
    assert saved["positive_weight_lower_millionths"] > 1000


def test_propagated_hidden_rotation_and_integer_kernel_crosscheck():
    # The exact native replay and marginal counts can use different numerical
    # proposals across platforms; the scientific bounds below are independent.
    s = measure.run()
    assert max(map(float, s["state_error_upper"])) < 1.1e-10
    for q, exp in s["cameras"].items():
        assert float(exp["simultaneous_mismatch_upper_from_material_counts"]) < 3e-12
        assert exp["archived_camera_audit"]["probabilities_checked"] == 200
        assert (
            exp["archived_camera_audit"]["minimum_margin_inside_count_interval"] > 4e-11
        )
        r = exp["sources"]["quarter_turn"]
        frozen, adapt = r["frozen"], r["adaptive"]
        assert frozen["counts"]["simultaneous99_max_pair_mismatch_lower"] > 3e-7
        assert frozen["max_pair_mismatch_lower"] > 4.9e-7
        assert adapt["max_pair_mismatch_upper"] < 3e-14
        assert all(a > 0 for a, b in adapt["native_directional_area_gain"])
