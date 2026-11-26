import json
from itertools import permutations, product

import numpy as np
import pytest
from scipy.stats import binom

from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin import (
    apparatus,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release import (
    detectors,
    flight,
    geometry,
)

from .measure import HERE, OWNERS, SIZE, WORDS, analyzers, enrol_matrix, receiver
from .mixture import SUMMED, TAIL, decode, interval, saddle_tail, summed_tail
from .readback import report


def test_ternary_load_read_has_an_intact_local_reference_for_every_possible_word():
    graph = federation(0)
    ring, _, leaves, _, _ = geometry.layout()
    source = apparatus.fixture()
    protected = set(ring) | set(leaves) | set(source["ports"]) | set(source["donor"])
    assert (
        len(detectors.verify()) == 2
    )  # calibrated unit shutter, blank transposition, signed read
    for owner in OWNERS:
        sites = receiver(owner)
        stations = {
            u: min(
                v
                for v in graph.neighbours[u]
                if v // 12 == owner and v not in protected and graph.degree[v] == 5
            )
            for u in sites
        }
        assert all(v in graph.neighbours[u] for u, v in stations.items())
        assert all(v in graph.neighbours[u] for u, v in zip(sites, sites[1:]))
        for w in product(range(3), repeat=3):
            if w.count(0) > 1 or w.count(2) > 1:
                continue
            # At least one unit opens its (1,1) shutter. All three sites
            # communicate by internal blank transpositions. Read at the open
            # unit's private neighbor, then undo; intermediate sites return.
            base = w.index(1)
            for j in range(3):
                stock = [*w, 1]
                route = list(range(j, base, 1 if j < base else -1)) + [base, 3]
                contacts = list(zip(route, route[1:]))
                for i, k in contacts:
                    stock[i], stock[k] = stock[k], stock[i]
                assert stock[3] == w[j]
                for i, k in reversed(contacts):
                    stock[i], stock[k] = stock[k], stock[i]
                assert stock == [*w, 1]


def test_contacts_match_previous_99_state_sector_and_closed_graph():
    states, previous = flight.chart((160, 167, 0), True)
    np.testing.assert_array_equal(enrol_matrix(states, 8, False), previous)
    closed = enrol_matrix(states, 8, True)
    inside = np.array([v // 12 == 13 and e // 12 == 13 for v, e, _ in states])
    assert not closed[np.ix_(inside, ~inside)].any()
    assert not np.array_equal(previous, enrol_matrix(states, 4, False))


def test_local_half_swap_preserves_pure_sector_and_interferes():
    words = tuple(permutations(range(3)))
    assert words == WORDS
    unitary = analyzers()["left"]
    np.testing.assert_allclose(unitary.conj().T @ unitary, np.eye(6), atol=2e-13)
    i, j = words.index((0, 1, 2)), words.index((1, 0, 2))
    sample = (np.eye(6)[:, i] + 1j * np.eye(6)[:, j]) / np.sqrt(2)
    a = abs((unitary @ sample)[[i, j]]) ** 2
    np.testing.assert_allclose(sorted(a), [0, 1], atol=3e-13)


def test_count_limits_hold_the_registered_tail():
    n = 5000
    for k in (1, SUMMED - 1, SUMMED, n // 3, n - SUMMED, n - 1):
        lo, hi = interval(k, n)
        below, above = summed_tail(k, n, lo, True), summed_tail(k, n, hi, False)
        np.testing.assert_allclose(below, binom.sf(k - 1, n, lo), rtol=1e-9)
        np.testing.assert_allclose(above, binom.cdf(k, n, hi), rtol=1e-9)
        np.testing.assert_allclose([below, above], TAIL, rtol=1e-5)
    # At 2^52 starts both tail evaluations agree where the decoder switches.
    lo, hi = interval(SUMMED, SIZE)
    np.testing.assert_allclose(
        [saddle_tail(SUMMED, SIZE, lo, True), saddle_tail(SUMMED, SIZE, hi, False)],
        [summed_tail(SUMMED, SIZE, lo, True), summed_tail(SUMMED, SIZE, hi, False)],
        rtol=1e-5,
    )


def test_both_deposited_count_decoders_and_missing_bank():
    raw = json.loads((HERE / "records.json").read_text())
    s = report(raw)
    assert s["total_banks"] == 324 and s["blocked_clicks"] == 0
    assert s["incompatible_conjugate_bin_pairs"] == [0, 324]
    assert s["resolved_free_contrast_pair_rates"] == [108, 108]
    assert s["resolved_against_incoherent_half_swap"] == [616, 864]
    assert (
        max(
            r["pair_given_both_classes_arrive"]
            for r in raw["ideal_receiver_bounds"]
            if r["key"][4:] == ["character", 1, 2] and not r["key"][1]
        )
        < 2e-11
    )
    raw["banks"].pop()
    with pytest.raises(ValueError):
        report(raw)
    mixture = json.loads((HERE / "mixture.json").read_text())
    s = decode(mixture)
    assert s == dict(
        compatible=1296,
        held_bins=1296,
        pooled_source_exclusions=359,
        pooled_source_bins=432,
    )
    mixture["source_heads"][3]["counts"].pop()
    with pytest.raises(ValueError):
        decode(mixture)


@pytest.mark.parametrize(
    "change", ["duplicate", "replace", "fractional", "exposure", "seed", "source"]
)
def test_mixture_rejects_invalid_enrollment(change):
    raw = json.loads((HERE / "mixture.json").read_text())
    if change == "duplicate":
        raw["stock_clicks"].append(raw["stock_clicks"][0])
    elif change == "replace":
        raw["held"][-1] = raw["held"][0]
    elif change == "fractional":
        raw["held"][0]["counts"][-1] = float(raw["held"][0]["counts"][-1])
    elif change == "exposure":
        raw["starts"] += 1
    elif change == "seed":
        raw["seed_key_prefix"] += 1
    else:
        raw["source_order"].reverse()
    with pytest.raises(ValueError):
        decode(raw)
