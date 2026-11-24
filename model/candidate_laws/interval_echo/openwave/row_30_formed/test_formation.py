"""Independent checks on the deposited local service and observer records."""

import json

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_30_formed import experiment as e
from model.executable_model import federation


def receipt():
    return json.loads((e.HERE / "results.json").read_text())


def test_repair_word_and_band_multiplicity():
    for s in receipt()["studies"]:
        level = (20, 80, 320).index(s["cells"])
        board = federation(level)
        edges = np.array(board.seams)
        words = s["formation"]
        bits = np.unpackbits(
            bytes_array(words["coin_pairs_packed_little_bit_first_hex"]),
            bitorder="little",
            count=2 * len(edges),
        ).reshape(-1, 2)
        closed = bits.all(1)
        assert np.count_nonzero(closed) == words["written"]
        writes = np.bincount(edges[closed].ravel(), minlength=board.n)
        assert (
            np.bincount(writes, minlength=7).tolist()
            == words["full_board_band_multiplicities"]
        )
        assert sum(words["full_board_band_multiplicities"]) == board.n
        assert writes[s["sharp_source"]] == s["selected_writes"] == 1
        if level == 0:
            # An independent full-board dense calculation includes all components.
            a, b = edges[~closed].T
            graph = np.zeros((board.n, board.n))
            graph[a, b] = graph[b, a] = 1
            h = (np.diag(graph.sum(1)) - graph) / 2 + np.diag(e.K * (6 - writes))
            energies = np.linalg.eigvalsh(h)
            for r, multiplicity in enumerate(words["full_board_band_multiplicities"]):
                base = e.K * (6 - r)
                assert (
                    sum((energies >= base - 1e-9) & (energies <= base + 6 + 1e-9))
                    == multiplicity
                )


def bytes_array(value):
    return np.frombuffer(bytes.fromhex(value), dtype=np.uint8)


def test_clock_and_joint_return_on_same_domain():
    for s in receipt()["studies"]:
        window, leak = e.band_window(s["selected_writes"])
        r = s["one_write_domain_rest"]
        assert r["spatial_odd_counts"] == 0
        assert len(s["study"]) == len(s["spatial_rereads"]) == 5
        for line, report in zip(
            [*s["study"], r], [*s["spatial_rereads"], r["spatial_reread"]], strict=True
        ):
            counts = np.asarray(line["clock_counts"])
            assert counts.sum() == e.SHOTS
            kept = counts[window]
            two = max(kept[:-1] + kept[1:])
            margin = (
                0.8 * (kept.sum() / e.SHOTS - e.COUNT_EPS - leak)
                - two / e.SHOTS
                - e.COUNT_EPS
            )
            np.testing.assert_allclose(margin, line["single_line_rejection_margin"])
            joint = np.array(report["returned_then_elsewhere_by_clock_bin"])
            assert joint.sum() == e.SHOTS
            p = joint[:, window].sum() / e.SHOTS - e.COUNT_EPS
            q = joint[0, window].sum() / e.SHOTS + e.COUNT_EPS
            rhs = max(0.0, p - leak - (1.0 - p + leak) * np.sqrt(leak)) ** 2
            np.testing.assert_allclose(
                rhs - q, report["joint_single_line_rejection_margin"]
            )
        assert all(
            r["excludes_one_line_in_band_at_simultaneous_99"] for r in s["study"]
        )
        assert r["spatial_reread"]["joint_excludes_one_line_in_band_at_simultaneous_99"]
        assert max(s["exchange_only_same_modes_eigen_residuals"]) < 1e-10
        for row in s["exchange_only_same_modes_clock_counts"]:
            assert sum(row) == e.SHOTS
            assert max(np.add(row, np.roll(row, 1))) > 0.79 * e.SHOTS


def test_unfiltered_denominator_and_independent_flights():
    boards = receipt()["studies"]
    for s in boards:
        for flight in s["flight"]:
            c = flight["radial_counts_then_tag_failure"]
            assert sum(c) == e.SHOTS
            assert flight["two_or_more_cells_all_launches"] == sum(c[2:-1]) / e.SHOTS
            assert flight["two_or_more_cells_given_tag"] == sum(c[2:-1]) / sum(c[:-1])
    first, second = (x for x in boards if x["cells"] == 80)
    assert first["seed"] != second["seed"]

    def remote(flight):
        return sum(flight["radial_counts_then_tag_failure"][2:-1])

    assert remote(first["flight"][1]) - remote(first["flight"][0]) > 60000
    assert max(remote(x) for x in second["flight"]) <= 10
