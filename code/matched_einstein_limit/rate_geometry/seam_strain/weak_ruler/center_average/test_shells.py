"""Rebuild the shell kernels and cross-check the marked pump by sparse evolution."""

import json

import numpy as np

from . import shells


def test_native_shell_sources_and_radar_separation():
    fresh = shells.run()
    saved = json.loads((shells.HERE / "shells.json").read_text())
    for record in (fresh, saved):
        groups = record["groups"]
        for name, trial in groups.items():
            for key in ("fixed_marginals_certificate", "ten_input_certificate"):
                proof = trial[key]
                assert (
                    proof["status"] == "separated"
                    and proof["inverse_defect_upper"] < 1e-5
                )
                assert proof["individually_fixed_volume_area_slope"][0] > 1.1e-6
                assert len(proof["pivot_columns"]) == len(trial["meters"])
                assert (
                    max(max(map(abs, v)) for v in proof["direction"].values())
                    < 1 + 1e-12
                )
            lotteries = trial["same_marginal_lotteries"]
            assert lotteries["area_slope_difference"][0] > 0
            for side in lotteries["probabilities"].values():
                assert all(0 <= lo <= hi <= 1 for lo, hi in side.values())
                assert (
                    sum(r[0] for r in side.values())
                    <= 1
                    <= sum(r[1] for r in side.values())
                )
            if len(trial["meters"]) == 4:
                fit = trial["four_physical_arms"]
                assert fit["arms"] == list(shells.native.ARMS)
                lo, hi = fit["common_origin_minor"]
                assert lo > 2 or hi < -3.6
        trial = groups["zero_and_aperture_remainder"]
        assert (
            trial["fixed_marginals_certificate"][
                "individually_fixed_volume_area_slope"
            ][0]
            > 0.00144
        )
        assert trial["four_physical_arms"]["common_origin_minor"][0] > 2.27
        rings = groups["distance_from_zero_within_aperture"]["meters"][::2]
        assert [r["ports"] for r in rings] == [
            [0],
            [1, 5, 7],
            [8, 9, 38],
            [39, 40, 42, 46, 47],
        ]
        assert groups["graph_ball_about_zero"]["meters"][2]["ports"] == [
            1,
            5,
            7,
            10,
            11,
        ]


def test_actual_shell_operators_add_and_match_independent_native_probe():
    board = shells.source.formation.federation(1)
    roots = shells.source.formation.region(board)
    r = json.loads((shells.HERE / "shells.json").read_text())
    total = json.loads((shells.HERE / "native.json").read_text())["regional"]
    group = r["groups"]["distance_from_zero_within_aperture"]
    wave = (
        shells.live_source.head_channel()
        .reshape(-1, 2, 4, 8)
        .transpose(1, 0, 2, 3)
        .reshape(-1, 8)
    )
    hw = (shells.source.system(board, roots, -1)[0] @ wave) / 40
    axes = dict(
        zip(
            shells.AXES, shells.live_source.qubit_hidden.correlation_axes(), strict=True
        )
    )
    for rr in group["meters"]:
        sites = shells.live_radar.station_rows(board, rr["ports"]).ravel()
        ph = wave[sites].conj().T @ hw[sites]
        material = (
            (ph + ph.conj().T) / 2
            if rr["meter"] == "energy"
            else wave[sites].conj().T @ wave[sites]
        )
        for tag, head in axes.items():
            source = np.kron(np.diag([2 / 3, 1 / 3]), np.asarray(head, complex))
            x = np.trace(material @ source).real
            lo, hi = rr["inputs"][tag]
            assert lo - 2e-11 < x < hi + 2e-11
    for meter in ("number", "energy"):
        for tag in ("reference", *shells.native.ARMS):
            parts = [
                rr["inputs"][tag] for rr in group["meters"] if rr["meter"] == meter
            ]
            low, high = map(sum, zip(*parts, strict=True))
            oldlo, oldhi = total[meter][tag]
            assert low <= oldhi and oldlo <= high


def test_material_lottery_readback_and_response():
    """Recover the physical source difference from the archived preparation settings."""
    r = json.loads((shells.HERE / "shells.json").read_text())
    names, source_columns, _ = shells.jet_lottery.material()
    source_columns = np.array(source_columns, float)
    area = np.array([np.mean(r["area_slopes"][tag]) for tag in shells.AXES])
    for trial in r["groups"].values():
        lottery = trial["same_marginal_lotteries"]
        arms = lottery["probabilities"]
        diff = np.array(
            [np.mean(arms["plus"][n]) - np.mean(arms["minus"][n]) for n in names]
        )
        head = source_columns @ diff
        kernel = trial["fixed_marginals_certificate"]["direction"]
        expected = (
            np.array([np.mean(kernel[tag]) for tag in shells.AXES])
            / lottery["common_scale"]
        )
        np.testing.assert_allclose(head, expected, rtol=0, atol=2e-15)
        for row in trial["meters"]:
            matter = np.array([np.mean(row["inputs"][tag]) for tag in shells.AXES])
            assert abs(matter @ head) < 2e-15
        lo, hi = lottery["area_slope_difference"]
        assert 0 < lo < area @ head < hi
