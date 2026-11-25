import json

import numpy as np
from scipy import sparse

from . import aperture as a


def test_native_cut_continuity_and_independent_short_window_calibrations():
    board = a.source.formation.federation(1)
    for charge in (-1, 1):
        h = a.source.system(board, a.source.formation.region(board), charge)[0]
        identity = sparse.eye(h.shape[0], dtype=np.int64, format="csr")
        for radius, sizes in ((0, (2, 10, 10)), (1, (12, 10, 20))):
            inside, rim, cut = a.enroll(board, radius)
            assert (len(inside), len(rim), len(cut)) == sizes
            assert set(inside).isdisjoint(rim)
            inner_mask = np.zeros(h.shape[0], dtype=np.int64)
            inner_mask[a.live_radar.station_rows(board, inside).ravel()] = 1
            inner = sparse.diags(inner_mask, dtype=np.int64, format="csr")
            for region, targets in (
                (inner, rim),
                (identity - inner, sorted({p for p, q in cut})),
            ):
                monitor = np.zeros_like(inner_mask)
                monitor[a.live_radar.station_rows(board, targets).ravel()] = 1
                absorber = sparse.diags(monitor, dtype=np.int64, format="csr")
                cross = region @ h @ (identity - region)
                # In -d<region>/dt the loss term is zero and all remaining terms
                # come from carrier hops crossing the registered cut.
                assert (absorber @ region).nnz == 0
                assert (region @ h - h @ region - cross + cross.T).count_nonzero() == 0
                active = np.column_stack(cross.nonzero())
                ports = (active // 4) % board.n
                assert {tuple(sorted(map(int, pair))) for pair in ports} == {
                    tuple(sorted(e)) for e in cut
                }
            record = a.short_window(board, charge, inside, cut)
            assert record["census_count"] == len(inside)
            if radius == 1:
                assert record["universal_blank_cut"]
                mat = 400 * 20 * np.eye(8, dtype=np.int64)
                np.testing.assert_equal(record["outward_coefficient_times_400"], mat)
                np.testing.assert_equal(record["inward_coefficient_times_400"], mat)


def test_continuous_detector_includes_arrivals_on_no_click_branches():
    b = a.source.formation.federation(1)
    inner, rim, _ = a.enroll(b, 0)
    h = a.source.system(b, a.source.formation.region(b), -1)[0]
    keep = {
        "crossing": a.live_radar.station_rows(b, inner).ravel(),
        "click": slice(None),
    }
    origin = inner[:1]
    numerical = a.captures(b, h, rim, origin, keep)
    hardware = a.bounded_capture(b, h, rim, origin, keep)
    tol = float(a.effect_error(h))
    for kind in keep:
        np.testing.assert_allclose(
            numerical[kind][origin[0]], hardware[kind][origin[0]], atol=tol, rtol=0
        )
    completed = hardware["crossing"][origin[0]]
    clicked = hardware["click"][origin[0]]
    assert np.linalg.eigvalsh(completed - clicked)[0] > 1e-4
    assert np.linalg.eigvalsh(clicked)[0] > 0 and np.linalg.eigvalsh(completed)[-1] < 1


def test_two_charge_enclosed_equal_census_transfer_and_all_outcomes():
    data = json.loads((a.HERE / "aperture.json").read_text())
    assert data["duration"] == 6 and data["detection_rate"] == 1
    assert data["cells"] == 80 and data["pump_charge"] == data["forming_charge"] == -1
    for k, pop in data["panels"].items():
        assert int(k) in (-1, 1)
        assert [r["radius"] for r in pop] == [0, 1]
        for row in pop:
            c = row["certificate"]
            k = json.loads((a.HERE.parent / "live_modular.json").read_text())[
                "modular_slopes"
            ]
            common = a.common_level_curve(c["slopes"], k)
            assert common["detector_and_entropy_rank"] == 4
            assert common["control_slopes"]["ordinary"][0] > 0
            native = common["finite_native_lotteries"]
            for arm in ("positive", "negative"):
                weights = native[arm].values()
                assert min(v[0] for v in weights) > 0
                assert sum(v[0] for v in weights) <= 1 <= sum(v[1] for v in weights)
            assert native["modular_gap_nats"][0] > 0.00038
            assert native["finite_entropy_gap_nats"][0] > 2e-11
            assert common["retained_entropy_slope_nats"][0] > 0.04
            saved = row["three_equal_readings"]["retained_entropy_slope_nats"]
            assert max(saved[0], common["retained_entropy_slope_nats"][0]) <= min(
                saved[1], common["retained_entropy_slope_nats"][1]
            )
            f = c["high_lottery_weight"]
            assert (
                0 < f[0] <= row["finite_equal_census"]["high_probability"] <= f[1] < 1
            )
            for key in ("outward", "inward", "K", "retained_entropy"):
                lo, hi = c["finite_equal_census"][key]
                assert lo * hi > 0
            stats = row["simulated_counts"]
            n = stats["trials_per_arm_and_setting"]
            assert n == 2**40
            for arm, datum in stats["draws"].items():
                assert 0 < datum["census_click"] < n
                for kind in ("outward", "inward"):
                    clicks, terminal, source = datum[kind]
                    assert (
                        clicks + terminal + source == n
                        and min(clicks, terminal, source) > 0
                    )
                    count = len(
                        {e[0 if kind == "outward" else 1] for e in row["oriented_cut"]}
                    )
                    matched = row["finite_equal_census"]
                    means = row["means"]
                    center = (
                        means[arm][kind]
                        if arm != "matched_lottery"
                        else (
                            (1 - matched["high_probability"])
                            * means[matched["low"]][kind]
                            + matched["high_probability"] * means[matched["high"]][kind]
                        )
                    )
                    fraction = (
                        center * data["isolated_transmission_probability"] / count
                    )
                    assert abs((clicks + terminal) / n - fraction) < 5e-6
            for side in ("left", "right"):
                assert side in c["fixed_census"]
                for kind in ("outward", "inward"):
                    r = c["fixed_census"][side]["transfer"][kind]
                    assert r["decorrelated"][0] * r["decorrelated"][1] > 0
            assert row["short_window"]["universal_blank_cut"] == (row["radius"] == 1)


def test_finite_native_rotation_contrast_realizes_the_solved_tangents():
    from scipy.linalg import expm

    board = a.source.formation.federation(1)
    head, _ = a.source.form(board, a.source.formation.region(board), -1)
    _, deltas, _ = a.live_survey.retained(True)
    for tag, g in a.live_survey.rotations(head).items():
        # The same pre-pump head is rotated in two physically distinct arms.
        np.testing.assert_allclose(g @ g, np.eye(4) / 4, atol=1e-14, rtol=0)
        states = []
        for sign in (1, -1):
            u = expm(-1j * sign * np.pi * g / 2)
            turned = u @ head @ u.conj().T
            for x, y in zip(
                a.source.formation.marginals(turned), a.source.formation.marginals(head)
            ):
                np.testing.assert_allclose(x, y, atol=5e-14, rtol=0)
            states.append(a.live_source.propagate_heads(turned))
        np.testing.assert_allclose(
            (states[0] - states[1]) / 2, deltas[tag], atol=3e-14, rtol=0
        )
