"""Raw reconstruction on all three carriers. Use pytest -s for discrepancy maxima."""

from copy import deepcopy
from decimal import Decimal as D
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations, product

import numpy as np
import pytest
from scipy.linalg import expm

from . import raw_verify as raw


@pytest.fixture(scope="module")
def receipt():
    return raw.read_json(raw.HERE/"results.json")


def digest(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def test_raw_input_bindings(receipt):
    # Donor git blobs at 44f7680d; CRLF working-tree conversion applies only to text.
    pins = {
        "analyzers.json": "bca73380e3f4960b8b4aa1268023573845e9f40463c8545af35ed550e7f17a69",
        "results.json": "06f003e17f7ad04be6e76b4d48e9cbabc4afbd697711449d5cb7e33467d23447",
        "PROTOCOL.md": "6514822767043ce397a2f313bbfe9e2fcb93c4b1a309d0279298f6f47f3bdb9f",
        "../row_30_arrows/PROTOCOL.md": "5097dd73eb3da0f0aae98113162d14020c607af6b89a498e67b74adcdbf58ec9",
        "staggered.json": "8ea06d8a234a5a63c091ceb9d40a25f951c8209e09c04e8ff624f37cf7ed88bf",
    }
    assert receipt["donor"] == "44f7680df55b4e68b89c301896df3e10237f046d"
    assert receipt["donor_sha256"] == pins == {p: digest(raw.DONOR/p) for p in pins}
    assert digest(raw.DONOR.parent/"row_30_arrows/results.json") == "38c126d0d8c860d52d14748b222b4c1d6006374aa0c8a829ac87cbad981565ba"
    assert receipt["design_sha256"] == digest(raw.HERE/"DESIGN.md")


def test_native_queries_and_wrong_decisions(receipt):
    echo2, ex8, ag4, _, _, mirror = raw.local()
    eye = np.eye(99, dtype=int)
    # Blank swap and blank reference in both orientations, including unused two-load columns.
    for j, (a, b) in enumerate(raw.LOADS):
        np.testing.assert_array_equal(ex8[:, j], 4*(eye[:, j]-eye[:, raw.ix(None, (b, a))]))
        np.testing.assert_array_equal(ag4[:, j], 4*(a-b)**2*eye[:, j])
    swap = expm(-1j*np.pi*(ex8[:9, :9]/8 + 8*ag4[:9, :9]/4))
    np.testing.assert_allclose(abs(swap), eye[:9, [raw.ix(None, pair[::-1]) for pair in raw.LOADS]],
                               atol=2e-14, rtol=0)
    # Omitting the coin reversal while reflecting ports and contrast changes the write.
    wrong = [raw.ix(None if lab is None else (lab[0], -lab[1]), w[::-1])
             for lab, w in product(raw.LABELS, raw.LOADS)]
    assert not np.array_equal(echo2[np.ix_(wrong, wrong)], echo2[np.ix_(mirror, mirror)])
    for sign in (-1, 1):
        ray, refusals = raw.attempt(sign)
        for mounts in product((0, 1), repeat=2):
            assert raw.joint_attempt(sign, mounts) == (ray, refusals)
        assert ray[0] == 0 and ray[2] == sign*ray[1]
        assert sum(x*x for x in ray) == F(1, 8)
        assert [[str(x) for x in row] for row in refusals] == receipt["first_wrong_then_parent_left_right"][str(sign)]
        bad, _ = raw.attempt(sign, omit_last_center=True)
        assert bad[0]**2 == F(1, 16) and sum(x*x for x in bad) != F(1, 8)


def test_two_fuel_records_both_mounts_and_detuning():
    for k in (*raw.STRENGTHS, -128, 128):
        r = np.hypot(1, k/128)
        bare = (np.sin(np.pi*r/4)/r)**2
        for c in (0, 1):
            np.testing.assert_allclose([raw.fuel(k, False, c), raw.fuel(k, True, c)],
                                       [bare, .5], atol=2e-14, rtol=0)
    assert .09 < raw.fuel(128, True, 1)-raw.fuel(128, False, 1) < .11


def test_every_physical_installation_and_order(receipt):
    archive = raw.read_json(raw.DONOR.parent/"row_30_arrows/results.json")["orientations"]
    choices = raw.spectroscopy()
    assert [b["cells"] for b in receipt["boards"]] == [20, 80, 320]
    for level, b in enumerate(receipt["boards"]):
        edges, degree, inc = raw.board(level)
        triangles = [(m, f) for m, f in permutations(inc[0], 2) if max(m, f) < 12 and f in inc[m]]
        assert len(triangles) == b["parking_triangle_orders"] == 10
        assert len(edges) == b["seams"] and len(degree) == 12*b["cells"]
        assert {raw.parking(level, *tri) for tri in triangles} == {b["builder_actions_on_accept"]}
        with pytest.raises(AssertionError, match="arrow points away"):
            raw.parking(level, *triangles[0], wrong_mount=True)
        assert set(b["alternative_orders"]) == {"reverse_local", "descending"}
        for order, old in (("ascending", "ascending"), ("reverse_local", "reversed_cell_labels"), ("descending", "descending")):
            graph = raw.component(level, order)
            if level == 0:
                report = archive[old]
                assert report["paths_per_port"] == np.bincount(graph.sites, minlength=len(degree)).tolist()
                assert report["longest_path"] == max(len(s[1]) for s in graph.states)
                assert report["distinct_coin_words"] == report["reachable_states"] == len(graph.states)
            if order == "ascending":
                assert sorted(set(graph.sites)) == b["physical_port_order"]
                assert len(graph.states) == b["complete_history_dimension"] and b["all_unlisted_ports"] == 0
            else:
                control = b["alternative_orders"][order]
                assert control["reachable_configurations"] == len(graph.states)
                assert control["donor_words_reachable"] == [all(p in graph.words for p in pair) for _, pair in choices]


def test_complete_donor_and_archived_pulse_waits(receipt):
    picks = raw.spectroscopy()
    assert receipt["selection"] == {
        "energy_at_8": 4, "k_squared": 1, "distinct_energies_tested": 31, "spatial_columns_tested": 32,
        "choices": [{"analyzer": i, "paths": [list(w) for w in pair]} for i, pair in picks]}
    maximum = 0.
    for (i, _), pulse in zip(picks, receipt["archived_waits"], strict=True):
        assert pulse["analyzer"] == i
        assert [r["k"] for r in pulse["settings"]] == [0, .5, 1, 2, 4, 8, 16, 64, 128]
        for r in pulse["settings"]:
            v = raw.historical_pulse(i, r["k"])
            maximum = max(maximum, *abs(np.subtract(v, [r["dark_weight"], r["two_port_return_37"]])))
            if r["k"] == .5:
                assert .3 < v[1] < .4
            if r["k"] == 0:
                assert v == (0., 0.)
    print(f"Donor pulse replay: max absolute discrepancy {maximum:.3g}")
    assert maximum < 8e-11


def test_two_port_detector_sees_only_the_donor_subspace_on_refinement():
    for level, rank in enumerate((2, 6, 6)):
        tree = raw.component(level)
        energy, vec = np.linalg.eigh(tree.h2)
        band = vec[:, abs(energy-1) < 1e-9]
        assert band.shape[1] == rank
        effect = band[np.isin(tree.sites, [9, 11])]
        np.testing.assert_allclose(np.linalg.eigvalsh(effect.T @ effect),
                                   [0]*(rank-2)+[1, 1], rtol=0, atol=2e-12)
        if level:
            # Same exact energy and wavelength, zero response at both enrolled ports.
            pair = [tree.words.index((0, 1, 38, 40, port)) for port in (45, 47)]
            ray = np.zeros(len(tree.sites), dtype=int)
            ray[pair] = [1, -1]
            np.testing.assert_array_equal(tree.h2 @ ray, ray)
            np.testing.assert_array_equal(tree.lap @ ray, ray)
            assert not ray[np.isin(tree.sites, [9, 11])].any()


def floors(edges, depth, n):
    rejected = F(7, 8)**n
    clipped = lambda x: float(max(F(0), x))
    return (clipped(1-rejected-F(edges+depth+n, 2**n)),
            clipped(1-rejected-edges*F(3, 4)**n-F(depth+n, 2**n)))


def test_source_ledger_all_caps_and_thresholds(receipt):
    for level, b in enumerate(receipt["boards"]):
        edges = len(raw.board(level)[0])
        for (i, paths), record in zip(raw.spectroscopy(), b["records"], strict=True):
            depth = len(paths[0])-2
            assert (record["donor_column"], record["prefix_hops"]) == (i, depth)
            assert record["histories"] == [list(p) for p in paths]
            assert all(raw.board(level)[1][w[-1]] == 5 for w in (*paths, paths[0][:-1]))
            for n, row in zip(raw.CAPS, record["floor_by_cap"], strict=True):
                assert row["retry_cap"] == n
                bounds = floors(edges, depth, n)
                assert (row["matched"], row["bare_abs_k_at_most_128"]) == bounds
                for p in (.5, *(raw.fuel(k, False, 0) for k in (*raw.STRENGTHS, 128))):
                    for sign in (-1, 1):
                        accept, fail = raw.ledger(edges, depth, n, sign, p)
                        assert accept > 0 and float(accept)+1e-14 >= bounds[int(p != .5)]
                        assert all(x > 0 for key, x in fail.items() if key != "failed_reset")
                        assert (fail["failed_reset"] > 0) == (n > 1)
                        if n == 1 and p == .5:
                            # With no retries this is one fuel at each edge, one
                            # return hop, depth route hops, and one five-query word.
                            exact = F(1, 2**(edges+depth+4))
                            assert abs(accept*D(exact.denominator)-1) < D("1e-25")
            for k, wanted in ((0, 22), (128, (23, 24, 26)[level])):
                p = .5 if k == 0 else raw.fuel(k, False, 0)
                first = next(n for n in range(1, 97) if raw.ledger(edges, depth, n, -1, p)[0] > D(".94"))
                assert first == wanted


def check_arm(data, expected_ports, tips, refusals):
    observed = np.asarray(data["ports"])
    assert observed.shape == expected_ports.shape and np.isfinite(observed).all()
    assert min(observed) >= 0 and data["refused"].keys() == refusals.keys()
    max_rel = 0.
    for key, truth in refusals.items():
        prob = data["refused"][key]
        assert np.isfinite(prob) and prob >= 0
        if truth:
            max_rel = max(max_rel, float(abs(D(str(prob))/truth-1)))
        else:
            assert prob == 0
    delta = max(np.max(abs(observed-expected_ports)), abs(data["home"]-sum(expected_ports[tips])))
    assert np.isfinite(delta) and delta < 8e-12 and max_rel < 5e-12
    assert abs(sum(observed)+sum(data["refused"].values())-1) < 1.5e-11
    return delta, max_rel


def test_all_original_start_port_and_refusal_probabilities(receipt):
    maximum, count = np.zeros(2), 0
    for level, b in enumerate(receipt["boards"]):
        tree = raw.component(level)
        ports = sorted(set(tree.sites))
        read = np.equal.outer(ports, tree.sites).astype(int)
        np.testing.assert_array_equal(read.sum(axis=0), np.ones(len(tree.sites), dtype=int))
        for (_, paths), record in zip(raw.spectroscopy(), b["records"], strict=True):
            depth, tips = len(paths[0])-2, [ports.index(p[-1]) for p in paths]
            appointments = record["appointments"]
            assert [(r["strength"], r["exposure"]) for r in appointments] == list(product(raw.STRENGTHS, raw.TIMES))
            for r in appointments:
                k, t = r["strength"], r["exposure"]
                assert set(r["arms"]) == set(raw.ARMS)
                for name, history in zip(raw.ARMS, raw.flight(level, paths, k, t), strict=True):
                    p = raw.fuel(k, False, 0) if name == "bare_fuel_minus" else .5
                    admitted, fail = raw.ledger(len(raw.board(level)[0]), depth, 96, 1 if name == "plus" else -1, p)
                    dist = float(admitted)*(read @ history)
                    errors = check_arm(r["arms"][name], dist, tips, fail)
                    maximum = np.maximum(maximum, errors)
                    if (k, t) == (8, 37):
                        assert (sum(dist[tips]) > .94) == (name in ("minus", "bare_fuel_minus"))
                        if level == 0 and name == "plus":
                            # Swapping two off-tip records preserves both the headline
                            # and total probability. The full address check must detect it.
                            altered = deepcopy(r["arms"][name])
                            bins = altered["ports"]
                            bins[0], bins[1] = bins[1], bins[0]
                            with pytest.raises(AssertionError):
                                check_arm(altered, dist, tips, fail)
                    count += 1
    print(f"{count} arm/settings: max absolute address or home discrepancy {maximum[0]:.3g}; "
          f"max relative refusal discrepancy {maximum[1]:.3g}")
    assert count == 750


def test_verification_policy_rejects_altered_receipt(receipt, monkeypatch):
    """Exercise the producer CLI's comparison policy only. No scientific oracle uses this import."""
    import sys
    from . import experiment

    sample = receipt["boards"][0]["records"][0]["appointments"][0]["arms"]["minus"]
    for change in ("drop_reset", "zero_reset", "bad_port", "nonfinite"):
        changed = deepcopy(sample)
        if change == "drop_reset":
            del changed["refused"]["failed_reset"]
        elif change == "zero_reset":
            changed["refused"]["failed_reset"] = 0.
        else:
            changed["ports"][0] = .1 if change == "bad_port" else float("nan")
        monkeypatch.setattr(sys, "argv", ["experiment"])
        monkeypatch.setattr(experiment, "run", lambda: sample)
        monkeypatch.setattr(experiment, "load", lambda _: changed)
        with pytest.raises(ValueError, match="Fields differ|Numerical replay differs|Nonfinite"):
            experiment.main()


def test_verification_policy_nonzero_exit_on_physical_gate():
    """Aliasing the two tips gives a zero-norm difference on the physical probability gate."""
    import os
    import subprocess
    import sys

    module = __package__ + ".experiment"
    script = f"""
from {module} import main, terminal
def planted_alias(tree, choices):
    item = choices[0]
    return terminal(tree, [dict(item, paths=[item["paths"][0]]*2)])
main.__globals__["terminal"] = planted_alias
main()
"""
    path = raw.HERE/"results.json"
    before = digest(path)
    run = subprocess.run([sys.executable, "-B", "-c", script], cwd=raw.HERE.parents[4],
                         env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
                         capture_output=True, text=True, timeout=45)
    assert run.returncode != 0 and "Invalid minus conditional law" in run.stderr
    assert digest(path) == before
