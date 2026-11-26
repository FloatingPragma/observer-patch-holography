"""Check the graph classification, material odds and a complete raw-board flight."""

from fractions import Fraction
from itertools import combinations
from math import comb

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from tools.strict_json import load
from . import masks
from .__main__ import ARCHIVE, compare
from .. import open_pair_check as filed


def test_exact_source_ceiling_and_counterfactual_panel():
    from . import admission
    limit = admission.verify()
    assert limit == admission.panel_ceiling([9, 3, 3, 1]) < Fraction(94, 100)
    assert admission.panel_ceiling([0, 1, 0, 0]) == 1 - Fraction(63, 64)**8
    assert admission.panel_ceiling([0, 0, 1, 0]) == admission.panel_ceiling([0, 0, 0, 1]) == 1
    assert admission.panel_ceiling([1, 0, 0, 0]) > Fraction(99, 100)


def test_all_masks_against_cancellation_at_each_destination():
    # Read the raw board. Classify by matched exits instead of multiplying operator columns.
    with np.load(filed.BOARD, allow_pickle=False) as f:
        edges = sorted(tuple(sorted(map(int, e))) for e in f["edges"]
                       if set(e).intersection((156, 157, 166)))
    def bit(p, q):
        return 1 << edges.index(tuple(sorted((p, q))))
    a, b, c, d = 157, 156, 166, 163
    shared = [bit(a, b), bit(c, b), bit(a, d), bit(c, d)]
    tip_only = sum(1 << i for i, e in enumerate(edges) if {a, c}.intersection(e)) - sum(shared)
    center_only = (2**15-1) ^ (tip_only | sum(shared))
    classes = load(ARCHIVE)["census"]["filed_address"]["classes"]
    matches_by_ray = {name: np.zeros(2**15, dtype=bool) for name in classes}
    for mask in range(2**15):
        tips = mask & tip_only == tip_only
        matches = all(bool(mask & x) == bool(mask & y) for x, y in zip(shared[::2], shared[1::2]))
        shut_center = mask & center_only == center_only
        plus = tips and mask & sum(shared) == sum(shared)
        for name, ok in (("one_load_minus", tips and matches), ("one_load_plus", plus),
                         ("pair_minus", tips and matches and shut_center), ("pair_plus", plus and shut_center)):
            matches_by_ray[name][mask] = ok
    mask_ids = np.arange(2**15)
    live = (1 - ((mask_ids[:, None] >> np.arange(15)) & 1)).astype(np.int16)
    board = masks.federation(0)
    for name, expected in matches_by_ray.items():
        words = ((b, a), (b, c)) if name.startswith("pair") else ((a,), (c,))
        sign = -1 if name.endswith("minus") else 1
        np.testing.assert_array_equal(masks.eigenmasks(board, edges, live, words, sign), expected)
        hits = np.flatnonzero(expected)
        assert np.bincount([int(i).bit_count() for i in hits], minlength=16).tolist() == classes[name]["by_write_count_0_to_15"]
        assert hits[(hits & sum(shared)) == 0].tolist() == classes[name]["four_routes_open_masks"]
    # The protected pair has the same integer columns if *every* remote collar changes.
    from ..topology_census import action
    written = set(board.seams) - {e for e in edges if bit(*e) in shared}
    x, z = (b, a), (b, c)
    assert action(board, written, {x: 1, z: -1}) == ({x: 2, z: -2}, {x: 3, z: -3})


def test_fixed_write_budget_and_overlapping_event_odds():
    # Enumerate all finite materials at each write count on a small carrier.
    slots = range(8)
    probs = []
    for k in range(9):
        hits = sum({0, 1, 2} <= set(mask) and set(mask).isdisjoint({3, 4})
                   for mask in combinations(slots, k))
        probs.append(Fraction(hits, comb(8, k)))
    bound = masks.exchangeable_bound(8, 3, 2)
    assert Fraction(bound["sharp_ceiling"]) == max(probs)
    assert bound["best_total_writes"] == [i for i, p in enumerate(probs) if p == max(probs)]
    # Two compatible events share ten writes and four blanks; a conflicting pair
    # shares a contact which must be both open and shut. With two events the
    # second-order bound is equality, including when they share no contacts.
    e = frozenset(range(11)), frozenset(range(11, 15))
    p = Fraction(81, 4**15)
    for other, joint in (
        ((frozenset([*range(10), 15]), e[1]), p/4),
        ((frozenset(range(1, 12)), frozenset(range(12, 16))), 0),
        ((frozenset(range(15, 26)), frozenset(range(26, 30))), p*p),
    ):
        bounds = masks.occurrence_bounds([e, other])
        assert Fraction(bounds["probability_at_least_one_lower"]) == 2*p-joint


def test_inclusive_historical_and_drawn_materials_from_raw_edges():
    archive = load(ARCHIVE)
    unit = archive["units"]
    with np.load(filed.BOARD, allow_pickle=False) as f:
        edges = sorted(tuple(sorted(map(int, e))) for e in f["edges"])
    words = np.array(list(combinations(range(240), 2)))
    a, b = words.T
    ix = np.zeros((240, 240), dtype=int)
    ix[a, b] = ix[b, a] = np.arange(len(words))
    home = [e for e in edges if e[0]//12 == e[1]//12 == 13]
    assert archive["pair_outcomes"] == [*map(list, home), "quiet"]
    # First-gateway and tidal donors predate the pair. A fresh clocked coin
    # history also checks the original formation bit convention end to end.
    bank = archive["ensemble"]["banks"][0]
    bits = np.random.default_rng([bank["seed"], 0]).binomial(1, .5, (630, 2))
    assert np.packbits(bits.ravel(), bitorder="little").tobytes().hex() == bank["coins_hex"]
    tidal = archive["historical"]["tidal_sham"]
    tidal_mask = set(map(tuple, tidal["written"]))
    assert len(tidal_mask) == 160 and all(p//12 == q//12 for p, q in tidal_mask)
    for closed, row, k, t in (
        ({next(e for e in edges if e[0]//12 != e[1]//12)}, archive["historical"]["first_gateway"], .5, 37),
        (tidal_mask, tidal, .5, 37),
        ({e for e, word in zip(edges, bits, strict=True) if all(word)}, bank, 0, 7),
    ):
        sources, targets = [], []
        for p, q in edges:
            if (p, q) not in closed:
                spectator = np.delete(np.arange(240), [p, q])
                sources.append(ix[p, spectator])
                targets.append(ix[q, spectator])
        s, d = np.concatenate(sources), np.concatenate(targets)
        adjacency = sparse.csr_matrix((np.ones(2*len(s)), (np.r_[s, d], np.r_[d, s])),
                                      shape=(len(words),)*2)
        deg = np.asarray(adjacency.sum(axis=1)).ravel()
        exchange = (sparse.diags(deg)-adjacency)/2
        refs = 6-np.bincount(np.ravel(edges), minlength=240)
        h = exchange + sparse.diags(k*(deg+refs[a]+refs[b]), dtype=float)
        launch = np.zeros((len(words), 2))
        launch[[ix[156, 157], ix[156, 166]], [0, 1]] = 1
        evolved = expm_multiply(-1j*t*h, launch)
        amp = evolved[ix[tuple(np.array(home).T)]]
        accepted = np.array([-13, 14])*np.sqrt(float(sum((Fraction(364, 729))**i for i in range(8))))/27
        np.testing.assert_allclose(abs(amp @ accepted)**2, np.array(row["filed_input"][:-1])/unit,
                                   atol=2e-11, rtol=0)
        # Reconstruct the 3-word monitored pulse without the study's acquisition code.
        x, z = int((156, 157) not in closed), int((156, 166) not in closed)
        lap = np.array([[z, -z, 0], [-z, x+z, -x], [0, -x, x]])
        diag = np.diag(lap) + [1, 0, 1]
        w, v = np.linalg.eigh(lap/2+k*np.diag(diag))
        keep = np.diag([1., 0., 1.])
        pulse = (v*np.exp(-2j*np.pi*w/3)) @ v.T
        vec = np.linalg.matrix_power(keep @ pulse, 3)[:, 0][[0, 2]]
        fail = 1-np.vdot(vec, vec).real
        weight = sum(fail**j for j in range(8 if z else 1))
        dist = weight*abs(amp @ vec)**2
        np.testing.assert_allclose(np.r_[dist, 1-dist.sum()], np.array(row["acquired"])/unit,
                                   atol=2e-11, rtol=0)
        swap = expm_multiply(-t*exchange, launch @ (weight*abs(vec)**2))
        dist = swap[ix[tuple(np.array(home).T)]]
        np.testing.assert_allclose(np.r_[dist, 1-dist.sum()], np.array(row["random_swaps"])/unit,
                                   atol=2e-11, rtol=0)
        np.testing.assert_allclose(np.linalg.svd(amp, compute_uv=False)[0]**2, row["optimal_XZ"]/unit,
                                   atol=2e-11, rtol=0)


def test_receipt_rejects_a_different_draw_and_a_changed_receiver():
    row = load(ARCHIVE)["ensemble"]["banks"][0]
    with pytest.raises(ValueError, match="exact"):
        compare(row, {**row, "seed": row["seed"]+1})
    with pytest.raises(ValueError, match="terminal"):
        compare(row, {**row, "acquired": row["acquired"][::-1]})
    invalid = [float("nan"), *row["acquired"][1:]]
    with pytest.raises(ValueError, match="terminal"):
        compare({**row, "acquired": invalid}, row)


@pytest.mark.parametrize("mutation", ["missing_panel", "repeat_seed", "different_history", "shift_loss"])
def test_admission_rejects_invalid_source_histories(monkeypatch, mutation):
    from . import admission
    record = load(ARCHIVE)
    banks = record["ensemble"]["banks"]
    if mutation == "missing_panel":
        banks.clear()
    elif mutation == "repeat_seed":
        banks[1] = banks[0]
    elif mutation == "different_history":
        banks[0]["coins_hex"] = banks[1]["coins_hex"]
    else:
        banks[0]["source_exhaustion"] -= 1
        banks[0]["reset_unavailable"] += 1
    monkeypatch.setattr(admission, "load", lambda _: record)
    with pytest.raises((ValueError, ArithmeticError)):
        admission.verify()
