"""Independent finite-word and instrument checks for the winding-reservoir decision."""

import json
from hashlib import sha256
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def independent_small_ring():
    n = 4
    words = np.arange(1 << (3 * n // 2))
    H = np.zeros((len(words), len(words)))
    for x in words:
        for a in range(n):
            b = (a + 1) % n
            left, right = (x >> a) & 1, (x >> b) & 1
            if left == right:
                continue
            y = x ^ (1 << a) ^ (1 << b)
            if a % 2 == 0:
                H[x, x] += 8.5
                H[y, x] -= 0.5
            else:
                h = (x >> (n + a // 2)) & 1
                if h == right:
                    H[x, x] -= 4
                    H[y ^ (1 << (n + a // 2)), x] -= 4
    return words, H


def test_grand_histogram_and_full_local_hamiltonian():
    data = json.loads((HERE / "reservoir.json").read_text())
    p = HERE / "reservoir_counts.npz"
    assert sha256(p.read_bytes()).hexdigest() == data["counts_and_spectra_sha256"]
    words, H = independent_small_ring()
    sector = np.array([int(x % 16).bit_count() for x in words])
    vals = [np.linalg.eigvalsh(H[np.ix_(sector == i, sector == i)]) for i in range(5)]
    np.testing.assert_allclose(
        [arr.min() for arr in vals], data["rings"][0]["sector_minima"]
    )
    Z = np.array([np.exp(-arr / 16).sum() for arr in vals])
    np.testing.assert_allclose(Z, Z[::-1])
    assert np.allclose(H, H[::-1, ::-1])
    with np.load(p) as raw:
        for ring in data["rings"]:
            n = ring["n"]
            for rec in ring["equilibrium"]:
                mu, beta, shift = rec["mu"], rec["beta"], rec["uniform_onsite"]
                y = raw[f"n{n}_b{beta}_mu{mu}_shift{shift}"]
                assert sum(y) == data["trials"]
                freq = y @ np.arange(n + 1) / (sum(y) * n)
                assert (
                    abs(freq - rec["exact_filling"])
                    < rec["simultaneous_filling_radius"]
                )
                assert rec["half_within_bound"] == (mu == shift)
                if n == 4 and beta == 1 / 16:
                    w = Z * np.exp(beta * (mu - shift) * np.arange(5))
                    np.testing.assert_allclose(
                        rec["exact_filling"], w @ np.arange(5) / (4 * sum(w))
                    )


def test_momentum_response_and_each_raw_impulse():
    data = json.loads((HERE / "momentum.json").read_text())
    for ring in data["finite_momentum"]:
        for rec in ring["sectors"]:
            if rec["number"] in {1, ring["n"] // 2 - 1, ring["n"] // 2 + 1}:
                assert (
                    rec["all_same_size"]
                    and rec["max_difference_at_equal_dimension"] < 2e-12
                )
        # Folding the uniform blank hopping band into two-port cells gives
        # the branches 17 ± |cos(k/2)| at cell momentum k.
        cells = ring["n"] // 2
        k = 2 * np.pi * np.arange(cells) / cells
        np.testing.assert_allclose(
            ring["all_blank_one_vacancy"]["ground_by_momentum"],
            17 - np.abs(np.cos(k / 2)),
            atol=2e-12,
            rtol=0,
        )
    winding, blank = data["charge_position_and_restored_heads"]
    for current in [winding, blank]:
        raw = np.array(current["full_word_counts"])
        cat = np.array(current["counts_by_category"])
        assert np.array_equal(raw.sum(axis=1), cat.sum(axis=1))
        assert (raw.sum(axis=1) == 2**22).all()
    assert np.count_nonzero(np.array(winding["counts_by_category"])[:, 2]) == 0
    assert winding["counts_by_category"][3][-1] > 4_180_000
    assert blank["counts_by_category"][3][2] > 4_000_000
    src = json.loads((HERE / "reservoir.json").read_text())["rings"][-1]
    with np.load(HERE / "reservoir_counts.npz") as raw:
        for row in src["local_transfer_spectra"]:
            j = round(row["k"] * 5 / (2 * np.pi))
            E, W = raw[f"add_remove_{row['sign']}_{j}_{row['parity']}"].T
            assert np.all(W >= 0)
            test = W[abs(E - row["kg_frequency"]) < 0.02].sum() / sum(W)
            assert abs(test - row["fraction_within_KG"]) < 1e-12
            assert test < 1e-15 and row["visible_lines"] == 90


def test_neutral_grating_weights_at_registered_sizes():
    study = json.loads((HERE / "neutral_response.json").read_text())
    for ring in study["rings"]:
        cells = ring["n"] // 2
        assert ring["uniform_words"] == 3**cells + 1
        for row in ring["response"]:
            assert row["matrix_residual"] < 1e-8
            a, b = row["channels"]
            np.testing.assert_allclose(
                a["first_visible"], b["first_visible"], atol=1e-8
            )
            for channel in row["channels"]:
                x = sum(line[1] for line in channel["energy_and_fraction"])
                assert abs(x + channel["residual_weight"] - 1) < 1e-8
    short, ten, eighteen, long = study["rings"]
    assert 3.17 < short["response"][-1]["channels"][0]["first_visible"] < 3.18
    assert 0.43 < long["response"][-1]["channels"][0]["first_visible"] < 0.44
    assert 1.90 < eighteen["response"][-1]["channels"][0]["first_visible"] < 1.92
    assert 3.54 < ten["response"][-1]["channels"][0]["first_visible"] < 3.56
