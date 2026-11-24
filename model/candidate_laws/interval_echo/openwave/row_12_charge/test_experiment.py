"""Independent native columns and frozen control outcomes for the three instruments."""

import json
from itertools import combinations

import numpy as np

from . import experiment as probe
from . import exterior


def archived(file, k, n, kind, shift=0, sense=1):
    data = json.loads((probe.HERE / file).read_text())
    return next(
        r
        for r in data["rows"]
        if (r["kappa"], r["loads"], r["kind"], r["shift"], r["sense"])
        == (k, n, kind, shift, sense)
    )


def test_projected_two_port_probe_is_native_up_to_common_scalar():
    """Assemble projected 99-state seam diagonals and all allowed transitions literally."""
    kappa = 8
    native = probe.law.local_term_parts()
    full = native[0] + kappa * native[1]
    geometry = probe.geometry()
    belt_edges = probe.belt.ring(10).seams
    north = set(geometry["near_indices"])
    # 0..9 belt, 10/11 poles, 12 far, 13..17 stationary far-end neighbors.
    edges = [
        *belt_edges,
        *((a, 10 if a in north else 11) for a in range(10)),
        (10, 12),
        *((12, v) for v in range(13, 18)),
    ]
    for n in (1, 3, 4):
        box = probe.sector(n)
        for sign in (-1, 1):
            basis = []
            for location in (10, 12):
                for word in box.states:
                    x = np.ones(18, int)
                    x[list(word)] = 2
                    x[location] = 1 + sign
                    basis.append(tuple(x))
            ix = {word: i for i, word in enumerate(basis)}
            projected = np.zeros((len(basis), len(basis)))
            for j, x in enumerate(basis):
                for a, b in edges:
                    sa, sb = x[a], x[b]
                    col = probe.law.SRC.index((sa, sb))
                    projected[j, j] += full[col, col]
                    if sa == sb:
                        continue
                    swapped = list(x)
                    swapped[b], swapped[a] = sa, sb
                    i = ix.get(tuple(swapped))
                    if i is not None:
                        row = probe.law.SRC.index((sb, sa))
                        projected[i, j] += full[row, col]
            _, gain, swap = probe.contact(kappa, sign)
            occ = np.array([len(north.intersection(s)) for s in box.states])
            h = box.hamiltonian(kappa)
            off = np.eye(len(h)) * swap
            proposed = np.block([[h + np.diag(occ * gain), off], [off, h]])
            residual = projected - proposed
            np.testing.assert_allclose(
                residual, np.eye(len(basis)) * residual[0, 0], atol=1e-12, rtol=0
            )


def test_sign_calibration_and_common_control():
    for k, signed in ((8, (15.5, -17)), (0, (-0.5, -1))):
        for iq, expected in enumerate(signed):
            sign = (-1, 1)[iq]
            assert abs(probe.contact(k, sign)[1] - expected) < 1e-12
            for n in (3, 4):
                for t in (0, 1):
                    wall = archived("results.json", k, n, "flight")["probes"][iq]
                    compact = archived("results.json", k, n, "compact")["probes"][iq]
                    assert [c["probabilities"] for c in wall["common"][t]] == [
                        c["probabilities"] for c in compact["common"][t]
                    ]
                    for c in wall["common"][t]:
                        assert sum(c["counts"]) == probe.SHOTS


def test_pole_translation_and_release():
    for n in (3, 4):
        plus = archived("results.json", 8, n, "flight")
        minus = archived("results.json", 8, n, "flight", sense=-1)
        assert plus["arrival_p"][3] > 0.70 and minus["arrival_p"][0] > 0.70
        held = archived("results.json", 8, n, "flight", shift=1)
        for p, other in zip(plus["probes"], held["probes"]):
            for f, fo in zip(p["differential"], other["differential"]):
                for label in (0, 1):
                    assert (
                        abs(
                            f["visibility"][label]
                            - (-1) ** label * fo["visibility"][label]
                        )
                        < 1e-11
                    )
            assert p["mobile"][2]["counts"][0] == [0] * 11
            assert 0 < sum(p["mobile"][3]["counts"][0])
    empty = archived("results.json", 8, 0, "calibration")
    assert (
        abs(sum(empty["probes"][0]["mobile"][3]["probabilities"][0]) - np.sin(1) ** 2)
        < 1e-11
    )


def test_exterior_full_chamber_and_guards():
    kappa = 8
    for n in (1, 3, 4):
        h, number, labels = exterior.chamber(n, kappa)
        indices = [number[s] for s in combinations(range(10), n)]
        closed = h[indices, :][:, indices].toarray() - probe.sector(n).hamiltonian(
            kappa
        )
        np.testing.assert_allclose(
            closed, np.eye(len(closed)) * closed[0, 0], atol=1e-12
        )
        assert np.count_nonzero(labels == 11) == len(
            list(combinations(range(10), n - 1))
        )
    for n in (3, 4):
        flight = archived("exterior.json", kappa, n, "flight")
        for s in flight["settings"]:
            if s["open"]:
                assert abs(complex(*s["visibility"]) - 1) > 0.02
                assert s["gateway_probabilities"][0] > 0.001
            else:
                assert s["visibility"] == [1, 0] and s["gateway_probabilities"] == [
                    0,
                    0,
                ]
            assert all(np.sum(f["counts"]) == probe.SHOTS for f in s["quadratures"])
    standard = archived("exterior.json", kappa, 1, "gateway_standard")
    for s in standard["settings"]:
        if not s["open"]:
            direct = np.exp(-1j * probe.contact(kappa, s["sign"])[1] * s["time"])
            assert abs(complex(*s["visibility"]) - direct) < 1e-11


def test_registered_reference_difference_only_addresses_belt_and_exterior_phase_depends_on_source():
    for n in (1, 3, 4):
        h0, _, labels = exterior.chamber(n, 8, reference=0)
        h1, _, _ = exterior.chamber(n, 8, reference=1)
        # A belt 2 in place of 1 costs 3*k at reference zero, k at reference one.
        # Gateway 2 has no spare reference. The assembly removes a common scalar.
        relative = -16.0 * (labels == 11)
        np.testing.assert_allclose(
            (h0 - h1).toarray(), np.diag(relative - relative.mean()), atol=1e-12
        )
    for n in (3, 4):
        row = archived("exterior_native.json", 8, n, "flight")
        compact = archived("exterior_native.json", 8, n, "compact")
        for s, c in zip(row["settings"], compact["settings"]):
            if s["open"] and s["time"] == 2:
                assert abs(complex(*s["visibility"]) - complex(*c["visibility"])) > 0.01
            if not s["open"]:
                assert s["visibility"] == [1, 0]
    bulk = json.loads((probe.HERE / "exterior_native.json").read_text())
    assert sum(bulk["bulk"]["counts"]) == probe.SHOTS
    assert all(sum(word) - 3 == 1 for word in bulk["bulk_pointer_words"])
