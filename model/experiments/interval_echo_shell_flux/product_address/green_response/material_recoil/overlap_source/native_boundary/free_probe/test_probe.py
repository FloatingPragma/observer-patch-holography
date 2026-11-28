"""Independent physical columns, actual local occupancy, and finite-board response."""

from fractions import Fraction
from hashlib import sha256
import json

import numpy as np

from . import measure as m
from . import time_mean


def test_each_two_load_generator_column_is_the_literal_native_edge_sum():
    board, seam, _, _ = m.make_board()
    h, words, bits = m.generator(board, seam, 2, 8, True)
    assert h.shape == (1104, 1104)
    # Compute columns directly in the native 99-state tensor basis: all
    # unserviced seams are blank, and both binary loads are retained.
    exchange, cost = m.LAW.local_term_parts()
    full = exchange + 8 * cost
    members = [1 + j for j, (_, d) in enumerate(m.LAW.LABELS) if abs(d) == 1]
    possible = [np.flatnonzero(abs(full[:, j]) > 1e-12) for j in range(99)]
    address = {word: i for i, word in enumerate(words)}
    for idx, word in enumerate(words):
        native_edges = []
        for a, b in board.seams:
            left, right = int(a in word), int(b in word)
            source_word = m.LAW.SRC.index((left, right))
            outside = set(word) - {a, b}
            native_edges.append((a, b, source_word, outside))
        for mem in range(4):
            oracle = np.zeros(h.shape[0])
            column = idx * 4 + mem
            oracle[column] = 8 * sum(6 - board.degree[p] for p in word)
            for a, b, src, outside in native_edges:
                written = (a, b) == seam
                local_col = src + 9 * (members[mem] if written else 0)
                for j in possible[local_col]:
                    load = m.LAW.SRC[j % 9]
                    assert max(load) <= 1
                    dst = tuple(
                        sorted(
                            outside
                            | {p for p, v in zip((a, b), load, strict=True) if v}
                        )
                    )
                    if not written:
                        assert j // 9 == 0
                    new_mem = members.index(j // 9) if written else mem
                    oracle[4 * address[dst] + new_mem] += full[j, local_col]
            np.testing.assert_allclose(
                h[:, [column]].toarray().ravel(), oracle, atol=2e-12, rtol=0
            )


def test_reader_uses_physical_occupations_and_both_independent_subtractions():
    data = json.loads((m.HERE / "results.json").read_text())
    b, seam, dist, origins = m.make_board()
    assert tuple(data["times"]) == m.TIMES
    assert (
        data["protocol_sha256"]
        == sha256((m.HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    )
    assert origins == [2, 13] and tuple(seam) == (1, 14)
    for label, rec in data["written"].items():
        assert rec["max_independent_wave_discrepancy"] < 2e-11
        for probe in rec["probes"]:
            contrast = (
                np.array(probe["with_extra_erased"])
                - np.array(rec["no_extra_erased"])
                - np.array(probe["with_extra_paired"])
                + np.array(rec["no_extra_paired"])
            )
            np.testing.assert_allclose(contrast, probe["port_contrast"], atol=1e-14)
            np.testing.assert_allclose(
                contrast @ dist, probe["radial_moment"], atol=1e-13
            )
            for key in ("with_extra_paired", "with_extra_erased"):
                total = np.array(probe[key])
                np.testing.assert_allclose(total.sum(axis=1), 2, atol=8e-12, rtol=0)
                assert np.min(total) > -1e-13 and np.max(total) < 1 + 1e-12
            assert probe["written_source_energy_increment"] == (2 * int(label) + 1) / 8
    native = data["written"]["8"]
    np.testing.assert_allclose(
        native["probes"][0]["radial_moment"],
        native["probes"][1]["radial_moment"],
        atol=2e-12,
    )
    assert native["probes"][0]["radial_moment"][2] > 0.069
    assert data["ordinary"]["probes"][0]["radial_moment"][2] < -0.159
    ordinary, words, _ = m.generator(b, seam, 1, 8, False)
    for p in seam:
        assert ordinary.diagonal()[words.index((p,))] == 51


def test_exact_expansion_on_an_independent_two_state_system_and_real_contact():
    from scipy import sparse

    # H = X swaps a unit from radius zero to radius one with probability sin²t.
    h = sparse.csr_matrix([[0, 1], [1, 0]])
    result = time_mean.coefficients(h, np.array([0, 1]), [0])
    assert result[0, :4].tolist() == [0, 1, Fraction(-1, 3), Fraction(2, 45)]
    data = json.loads((m.HERE / "time_mean.json").read_text())
    assert (
        data["protocol_sha256"]
        == sha256((m.HERE / "TIME_MEAN.md").read_bytes()).hexdigest()
    )
    expected = {
        "native": Fraction(-1445, 2304),
        "no_contrast": Fraction(-5, 2304),
        "ordinary": Fraction(-85, 576),
    }
    for tag, v in data["panels"].items():
        assert max(g["largest_within_group"] for g in v["spectral_groups"]) < 1e-10
        assert min(g["smallest_separation"] for g in v["spectral_groups"]) > 1e-6
        for row in v["rows"]:
            coeff = [Fraction(x) for x in row["taylor_even_coefficients"]]
            assert coeff[:3] == [0, 0, 0] and coeff[3] == expected[tag]
    assert (expected["native"] / (Fraction(17, 8))) / (expected["ordinary"] / 51) == 102


def test_finite_local_counts_resolve_opposite_primary_signs():
    from . import sharp_read

    data = json.loads((m.HERE / "sharp_read.json").read_text())
    assert (
        data["protocol_sha256"]
        == sha256((m.HERE / "SHARP_READ.md").read_bytes()).hexdigest()
    )
    assert data["enrolled_raw_means"] == 2880
    _, _, dist, _ = m.make_board()
    for label, panels in data["results"].items():
        for row in panels:
            raw = np.array(row["counts"])
            assert raw.min() >= 0 and raw.max() <= sharp_read.SHOTS
            reduced = (raw[0] - raw[1] - raw[2] + raw[3]) @ dist / sharp_read.SHOTS
            np.testing.assert_array_equal(reduced, row["radial"])
            if label == "8":
                assert reduced[2] - row["radius_99"] > 0
            else:
                assert reduced[2] + row["radius_99"] < 0
