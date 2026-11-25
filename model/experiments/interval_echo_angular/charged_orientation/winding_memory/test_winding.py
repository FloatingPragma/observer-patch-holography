import json

import numpy as np

from . import measure


def test_full_board_prefix_subspace():
    ring, geometry = measure.apparatus()
    assert len(ring) == 15
    assert geometry["external_contacts"] == 5
    assert len(geometry["owners"]) == 5
    words, terms = measure.sector(len(ring))
    degree = np.full(len(words), 2)
    degree[[0, -1]] = 1
    # Directly generated native matrix versus the analytic path including its end diagonals.
    expected = np.diag(degree) + np.diag(np.ones(15), 1) + np.diag(np.ones(15), -1)
    np.testing.assert_array_equal(sum(terms), -4 * expected)
    base = None
    n = len(ring)
    for hole, heads in words:
        bits = np.array([(heads >> j) & 1 for j in range(n)])
        load = np.array([int(j != hole) for j in range(n)])
        gauss = load - bits + np.roll(bits, 1)
        if base is None:
            base = gauss
        np.testing.assert_array_equal(gauss, base)


def test_observer_controls_and_receipt():
    result = json.loads(json.dumps(measure.run(), default=int))
    saved = json.loads((measure.HERE / "records.json").read_text())
    # Numeric eigensolvers can round at final printed precision on different platforms.
    assert result["geometry"] == saved["geometry"]
    assert result["source_reader"] == [1.0, 0.0]
    assert result["head_sites_needed_for_stationary_winding_coherence"] == 15
    assert result["remote_rank_certificate"]["determinant_lower"] == "1/10"
    assert result["repeated_directed_winding_probability"] == 0
    assert [bank["successes"] for bank in result["addressed_ramsey"]] == [
        measure.SHOTS,
        0,
    ]
    for row, frozen in zip(result["free_records"], saved["free_records"], strict=True):
        for key in (
            "sharp_at_source",
            "sharp_at_target",
            "sharp_source_and_wound",
            "target_dephased",
        ):
            assert abs(row[key]["probability"] - frozen[key]["probability"]) < 2e-11
            assert row[key]["successes"] == frozen[key]["successes"]
        banks = list(
            zip(row["target_two_inputs"], frozen["target_two_inputs"], strict=True)
        )
        for sample, stored in zip(
            row["local_launches"], frozen["local_launches"], strict=True
        ):
            banks.extend((sample[key], stored[key]) for key in sample)
        for sample, stored in banks:
            assert sample["successes"] == stored["successes"]
            np.testing.assert_allclose(
                sample["probability"], stored["probability"], atol=2e-11
            )
            np.testing.assert_allclose(
                sample["interval99"], stored["interval99"], atol=2e-11
            )
    rows = {bank["time"]: bank for bank in result["free_records"]}
    first, second = rows[1]["target_two_inputs"]
    assert first["interval99"][1] < second["interval99"][0]
    assert rows[2]["sharp_source_and_wound"]["interval99"][0] > 0.13
    assert rows[1.5]["remote_pair_determinant"] > 0.10
    assert rows[0]["remote_pair_determinant"] == 0


def test_mixed_charge_formation_and_bulk_control():
    from fractions import Fraction
    from . import mixed_masses as masses

    saved = json.loads((measure.HERE / "mixed_masses.json").read_text())
    assert saved["sweep_cutoff"] == str(masses.LIMIT)
    assert masses.contacts()["even"] == (-32, -32)
    assert masses.mixed_parities() == saved["mixed_parities_selected_motif"]
    # Independent eigensolver vs exact signed-minor counts on both finite and
    # cyclic geometries, including the strict alternating zero mode.
    for word in ("0101010", "010101", "0100001", "0102010"):
        for cyclic in (False, True):
            d, e = masses.coefficients(word, "even", cyclic=cyclic)
            h = np.diag(d) + np.diag(e[: len(d) - 1], 1) + np.diag(e[: len(d) - 1], -1)
            if cyclic:
                h[0, -1] = h[-1, 0] = e[-1]
            vals = np.linalg.eigvalsh(h / 8)
            for q in (Fraction(0), masses.LIMIT, Fraction(51, 10)):
                assert masses.inertia(d, e, q, cyclic) == (
                    np.count_nonzero(vals < float(q) - 1e-10),
                    np.count_nonzero(abs(vals - float(q)) < 1e-10),
                )
    row = masses.periodic("01" * 7 + "00", 9)
    assert row["finite_line"] == saved["periodic_rows"][2]["finite_line"]
    assert row["phase_zero_edge"][0] > 0.034
    assert row["bandwidth"] > 0
    assert row["source_to_detector_pole_norm"] > 0.004
