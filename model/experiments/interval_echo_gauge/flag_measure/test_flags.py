from fractions import Fraction
import json
from math import log, sqrt

import numpy as np

from .measure import (
    CLASSIFIER,
    EDGES,
    FACES,
    HERE,
    INDEX,
    N,
    SHOTS,
    WORDS,
    calibration_inputs,
    swap_indices,
)


def receipt():
    return json.loads((HERE / "records.json").read_text())


def test_exact_equivariance_and_single_flag_orbit():
    mod, manifest, _, mapping = calibration_inputs()
    _, _, proper, _, _ = mod.load_inputs()
    pinned = {tuple(sorted(f)) for f in mod.canonical_faces(manifest)}
    assert {tuple(sorted(mapping[i] for i in f)) for f in FACES} == pinned
    groups = proper["proper_port_action"]["permutation_rows"]
    flags = {(f, p) for f in pinned for p in f}
    f, p = min(flags)
    assert {(tuple(sorted(g[i] for i in f)), g[p]) for g in groups} == flags
    assert len(flags) == 60 and all(sum(i in f for f in pinned) == 5 for i in range(N))
    for i in range(N):
        marginal = sum(Fraction(1, 60) for _, corner in flags if i == corner)
        assert marginal == Fraction(1, N)


def test_independent_last_insert_position_and_markov_irreducibility():
    for size in range(2, 13):
        pos = [Fraction(0)] * size
        pos[-1] = Fraction(1)
        for i in reversed(range(size - 1)):
            p = Fraction(i + 1, i + 2)
            pos[i], pos[i + 1] = (1 - p) * pos[i] + p * pos[i + 1], p * pos[i] + (
                1 - p
            ) * pos[i + 1]
        assert pos == [Fraction(1, size)] * size
    # Every one-seam interchange has positive mass in a half-exchange sweep.
    visited = {0}
    stack = [0]
    while stack:
        old = stack.pop()
        for u, v in EDGES:
            i = int(swap_indices(u, v)[old])
            if i not in visited:
                visited.add(i)
                stack.append(i)
    assert len(visited) == 660


def test_conditional_face_and_corner_records_with_all_refusals():
    d = receipt()

    def matches(data):
        n = sum(data)
        assert n == SHOTS and len(data) == 660
        table = {f: [data[INDEX[(p, *sorted(set(f) - {p}))]] for p in f] for f in FACES}
        totals = {f: sum(counts) for f, counts in table.items()}
        admitted = sum(totals.values())
        assert 0 < admitted <= n

        def radius(trials):
            return sqrt(log(2 * 1024 / 0.01) / (2 * trials))

        return (
            abs(admitted / n - 1 / 11) <= radius(n)
            and all(
                abs(q / admitted - 1 / 20) <= radius(admitted) for q in totals.values()
            )
            and all(
                q and all(abs(c / q - 1 / 3) <= radius(q) for c in table[f])
                for f, q in totals.items()
            )
        )

    for bank in d["native_and_controls"]:
        assert matches(bank["counts"]) == (bank["speed"] == 1)
        if bank["speed"] == 1:
            assert (
                bank["max_uniform_error"] < 2e-14
                and int(bank["numerator_per_word"]) > 0
            )
    assert all(matches(bank) for bank in d["host_counts"])
    assert matches(d["half_contacts"][-1]["counts"])
    assert d["half_contacts"][-1]["TV"] < 2e-7


def test_registered_metric_and_uniform_full_family_separation():
    d = receipt()["classification"]
    c = json.loads(
        (CLASSIFIER / "port_dual_metric_selector.certificate.json").read_text()
    )
    for family, (a, b) in c["squared_distances_at_measure_point"].items():
        exact = float(Fraction(a)) + sqrt(5) * float(Fraction(b))
        np.testing.assert_allclose(
            d["reference"][family]["squared_distance"], exact, atol=1e-10
        )
        assert (
            d["score_envelopes"][family][0]
            < d["fitted"][family]["squared_distance"]
            < d["score_envelopes"][family][1]
        )
    assert d["score_envelopes"]["G"][1] < min(
        d["score_envelopes"][x][0] for x in ("F", "P")
    )
    for f, s in [("F", -1), ("G", 1)]:
        a, b, _ = d["fitted"][f]["coefficients"]
        assert s * a * b > 0  # nearest point is inside the compact stratum
    assert d["f_control"]["F"]["squared_distance"] < 1e-23
    assert min(d["f_control"][x]["squared_distance"] for x in ("G", "P")) > 90


def test_sharp_port_product_gram_and_independent_dense_probes():
    d = receipt()
    gram = np.zeros((N, N), int)
    for (p, a, b), c in zip(WORDS, d["native_and_controls"][0]["counts"], strict=True):
        if tuple(sorted((p, a, b))) in FACES:
            v = np.eye(N, dtype=int)[d["material_to_pinned"][p]]
            gram += c * np.outer(v, v)
    assert np.array_equal(gram, np.diag(d["classification"]["port_counts"]))
    # A non-coordinate product recovers the corresponding port function norm.
    u = np.arange(N) - 5
    v = u[::-1]
    obs = sum(n * u[k] * v[k] for k, n in enumerate(gram.diagonal()))
    assert u @ gram @ v == obs and gram.trace() == d["classification"]["admissions"]
