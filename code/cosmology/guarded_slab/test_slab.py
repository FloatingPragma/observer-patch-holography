import hashlib
import json
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
from guarded_transfer import HERE, authenticate, refinement, source, worst_word_distance
from guarded_release import compute, energy_control, route


def test_trajectory_frequencies_follow_independent_sharp_word_law():
    k, _ = source()
    receipt = json.loads((HERE / "results.json").read_text())
    with np.load(HERE / "counts.npz") as saved:
        for board in receipt["boards"]:
            counts = saved[f"L{board['level']}"]
            assert counts.shape == (2, board["owners"], 8)
            assert np.all(counts.sum(axis=-1) == board["shots_per_start_per_owner"])
            for start in (0, 1):
                path_probs = []
                for word in product((0, 1), repeat=3):
                    w = (start,) + word
                    path_probs.append(np.prod([k[b, a] for a, b in zip(w, w[1:])]))
                observed = counts[start].sum(axis=0)
                p = np.array(path_probs)
                n = observed.sum()
                assert np.max(abs(observed - n * p) / np.sqrt(n * p * (1 - p))) < 6
    assert refinement()["max_flip_probability_error"] == "0"


def test_authentication_and_rank_do_not_depend_on_free_interleaving():
    log = json.loads((HERE / "results.json").read_text())["log"]
    events, init = log["events"], log["initial"]
    assert (
        authenticate(sorted(events, key=lambda e: -e["register"]), init)
        == log["heights"]
    )
    i = next(i for i, v in enumerate(events) if v["writer"] is not None)
    foreign = next(e for e in events[:i] if e["register"] != events[i]["register"])
    for change in [
        dict(writer=foreign["id"]),
        dict(writer=events[i]["id"]),
        dict(read=1 - events[i]["read"]),
        dict(writer=None),
    ]:
        damaged = [dict(e) for e in events]
        damaged[i].update(change)
        with pytest.raises(ValueError):
            authenticate(damaged, init)


def test_capped_release_by_forward_state_enumeration():
    # Evolve all joint receiver states as a short explicit Markov tree.
    cap = 3
    receipt = compute(cap)
    words = list(product((0, 1), repeat=3))
    dist = {(0, 0, 0): F(1)}

    def step(d, active):
        nxt = {}
        for w, mass in d.items():
            for flips in product((0, 1), repeat=len(active)):
                after = list(w)
                p = mass
                for i, b in zip(active, flips):
                    after[i] ^= b
                    p *= F(1, 4) if b else F(3, 4)
                key = tuple(after)
                nxt[key] = nxt.get(key, 0) + p
        return nxt

    for t in range(1, cap + 1):
        dist = step(dist, (0, 1, 2))
        delivered = dist
        for enrolled in [(0, 1, 2), (1, 2), (2,)]:
            delivered = step(delivered, enrolled)
        observed = [F(3, 4) ** (t - 1) * F(1, 4) * delivered[w] for w in words]
        assert observed == list(map(F, receipt["common_sharp_by_first_change"][t - 1]))
    assert sum(map(F, receipt["independent_sharp"].values())) == 1


def test_marked_ground_transform_is_positive_and_compresses():
    k = np.array([[0.75, 0.25], [0.25, 0.75]])
    t = k
    for _ in range(3):
        t = np.kron(t, k)
    rng = np.random.default_rng(118)
    for keep in range(4):
        read = np.array([bits[keep] for bits in product(range(2), repeat=4)])
        f = rng.normal(size=2)
        np.testing.assert_allclose(t @ f[read], (k @ f)[read])
    assert np.linalg.eigvalsh(t)[0] > 0


def test_exact_worst_word_distance_against_full_rows():
    k = np.array([[0.75, 0.25], [0.25, 0.75]])
    for t in range(4):
        one = np.linalg.matrix_power(k, t)
        full = np.kron(np.kron(one, one), one)
        distances = abs(full - 1 / 8).sum(axis=1) / 2
        np.testing.assert_allclose(distances, float(worst_word_distance(3, t)))


def test_cached_release_and_native_source_bindings():
    release = json.loads((HERE / "release.json").read_text())
    assert release.pop("energy_reader_on_either_rail") == energy_control()
    for age, row in enumerate(release["common_energy_by_first_change"], 1):
        assert list(map(F, row)) == [F(3, 4) ** (age - 1) / 32] * 8
    assert release.pop("routes") == [route(level) for level in range(3)]
    assert release == compute()
    assert release["unbounded_off_diagonal_cov"] == ["3/650", "3/1300", "3/2600"]
    assert release["unbounded_connected_third"] == "51/754000"
    assert sum(map(F, release["common_sharp"])) + F(release["timeout"]) == 1
    receipt = json.loads((HERE / "results.json").read_text())
    assert (
        hashlib.sha256((HERE / "counts.npz").read_bytes()).hexdigest()
        == receipt["counts_sha256"]
    )
    for name, expected in receipt["inputs"].items():
        assert (
            hashlib.sha256((HERE.parents[2] / name).read_bytes()).hexdigest()
            == expected
        )
