import hashlib
from fractions import Fraction as F
from itertools import combinations, product
import json

import numpy as np

from energy_release import HERE, exact, pair_after, unbounded, stroke_records
from energy_routes import apparatus, public_capacity, rotate


def test_stochastic_energy_contact_from_the_four_state_matrix():
    swap = np.eye(4, dtype=int)[[0, 2, 1, 3]]
    kernel = (3 * np.eye(4, dtype=int) + swap) / 4
    initial = np.array([8, 4, 2, 1]) / 15
    for n in range(8):
        np.testing.assert_allclose(
            initial @ np.linalg.matrix_power(kernel, n),
            list(map(float, pair_after(F(1, 5), F(1, 3), n))),
        )
    native = unbounded(F(1, 5), F(1, 3), (1, 2, 3))
    assert native["cross_covariance"] == [F(1, 48750), F(1, 97500), F(1, 195000)]
    assert native["connected_third"] == -F(17, 848250000)


def test_original_denominator_and_sampled_energy_words():
    result = json.loads((HERE / "energy.json").read_text())
    assert (
        hashlib.sha256((HERE / "energy_counts.npz").read_bytes()).hexdigest()
        == result["counts_sha256"]
    )
    assert (
        hashlib.sha256((HERE / "protocol.json").read_bytes()).hexdigest()
        == result["protocol_sha256"]
    )
    assert F(result["filter_acceptance"]) == F(15, 32) ** 3
    n = 16
    bins = np.arange(64)
    p = exact(F(1, 5), F(1, 3), n, (1, 2, 3))
    probs = np.array(list(map(float, p["words"]))) / float(p["hit"])
    archive = np.load(HERE / "energy_counts.npz")
    for i, board in enumerate(result["boards"]):
        assert json.loads(json.dumps(apparatus(i))) == {
            k: v for k, v in board.items() if k != "sample"
        }
        count = archive[f"L{i}"].sum(axis=0)
        freq = np.zeros(64)
        for parity in range(8):
            final = bins.copy()
            for j in range(3):
                if (parity >> j) & 1:
                    final ^= (((bins >> (2 * j)) ^ (bins >> (2 * j + 1))) & 1) * (
                        3 << (2 * j)
                    )
            # Archive bits are (probe,bath), least significant probe first.
            ids = sum(
                (2 * ((final >> (2 * j)) & 1) + ((final >> (2 * j + 1)) & 1))
                * 4 ** (2 - j)
                for j in range(3)
            )
            np.add.at(freq, ids, count[:, parity])
        s = board["sample"]
        assert freq.sum() == s["delivered"]
        assert (
            sum(map(sum, s["first_filter_leaks"])) + s["timeout"] + s["delivered"]
            == s["launches"]
        )
        assert (
            np.max(
                abs(freq - freq.sum() * probs)
                / np.sqrt(freq.sum() * probs * (1 - probs))
            )
            < 6
        )


def test_independent_count_propagation_and_covariance_bounds():
    bins = json.loads((HERE / "precision.json").read_text())
    settings = json.loads((HERE / "precision_protocol.json").read_text())
    result = json.loads((HERE / "energy.json").read_text())
    for name, v in bins.items():
        theory = result["native"] if name == "native" else result["controls"][name]
        cov = np.array([float(F(v)) for v in theory["delivered"]["cross_covariance"]])
        observed = np.array(v["covariance"])
        half = np.array(v["simultaneous_covariance_halfwidth"])
        assert np.all(abs(observed - cov) < half)
        assert sum(map(sum, v["age_by_terminal_word"])) == v["delivered"]
        assert (
            v["delivered"] + v["timeout"] + sum(map(sum, v["first_filter_leaks"]))
            == settings["launches_per_bank"]
        )
        # The law of total covariance separates age variation from residual
        # within-age correlation directly on the independent endpoint counts.
        counts = np.array(v["age_by_terminal_word"])
        n = counts.sum(axis=1)
        bits = (np.arange(64)[:, None] >> np.array([0, 2, 4])) & 1
        age_means = counts @ bits / n[:, None]
        centered = age_means - np.array(v["probe_upper"])
        between = np.array(
            [
                np.average(centered[:, i] * centered[:, j], weights=n)
                for i, j in combinations(range(3), 2)
            ]
        )
        assert np.max(abs(observed - between)) < max(half)
        if name == "native":
            assert np.all(observed > half)
        else:
            assert np.all(abs(observed) < half)
    assert all(
        sum(map(F, v)) == 1 for v in result["native"]["independent_age_marginals"]
    )


def test_tomography_and_evolution_of_every_region():
    eye = np.eye(2)
    x = eye[::-1]
    z = np.diag([1, -1])
    local = [eye, x, z, z @ x]
    words = list(product(range(2), repeat=3))
    shift = np.zeros((8, 8))
    for j, w in enumerate(words):
        shift[words.index(rotate(w)), j] = 1
    swap01 = np.eye(8)[[words.index((b, a, c)) for a, b, c in words]]
    swap02 = np.eye(8)[[words.index((c, b, a)) for a, b, c in words]]
    np.testing.assert_array_equal(swap02 @ swap01, shift)
    for size in range(4):
        for region in combinations(range(3), size):
            factors = list(
                product(*(range(4) if i in region else [0] for i in range(3)))
            )
            basis = []
            for labels in factors:

                def op(values):
                    a, b, c = [local[i] for i in values]
                    return np.kron(np.kron(a, b), c)

                mat = op(labels)
                for n in range(7):
                    u = np.linalg.matrix_power(shift, n)
                    np.testing.assert_array_equal(u @ mat @ u.T, op(rotate(labels, n)))
                basis.append(mat.ravel())
            assert np.linalg.matrix_rank(basis) == 4**size
    # Equal energy probabilities on the two sharp rails leave their different
    # positions visible in the rotated region.
    np.testing.assert_equal(x.diagonal(), [0, 0])
    assert z[0, 0] != z[1, 1]
    # Dropping the stroke count makes rotations of a word indistinguishable.
    policy = json.loads((HERE / "protocol.json").read_text())
    caps = public_capacity(policy)
    assert [
        caps[k]["M0"] for k in ("collective", "one", "any_two", "count_hidden")
    ] == [8, 2, 4, 4]


def test_clock_counts_actual_read_changes_and_moves_the_three_native_axes():
    t = json.loads((HERE / "energy.json").read_text())["ring_tomography"]
    assert t == stroke_records()
    for name, bank in t["banks"].items():
        for i, item in enumerate(bank, 1):
            assert item["clock_count"] == (i if name == "native" else 0)
            assert item["clock_rail"] == (i % 2 if name == "native" else 0)
            n = t["shots_per_setting"]
            initial = [[1, 0, 0], [0, -1, 0], [0, 0, 1]]
            expected = np.array(rotate(initial, item["clock_count"]))
            observed = 2 * np.array(item["plus_counts"]) / n - 1
            assert np.max(abs(observed - expected)) < 0.08
