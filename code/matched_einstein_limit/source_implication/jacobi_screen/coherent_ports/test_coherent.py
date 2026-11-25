"""Check the new instruments, equal ensembles, linear phase response, and rank decisions."""

import json

import numpy as np
from mpmath import iv
from scipy.linalg import expm

from ...action_response.four_ticks import parse_interval
from . import certificate, measure


def test_splitter_is_physical_and_preserves_all_source_blocks():
    position = np.eye(5) / 5
    for charge in (-1, 1):
        for sign in (-1, 1):
            h = measure.blank_clock(charge, 8)
            t = measure.DWELL if sign == 1 else 2 * np.pi - measure.DWELL
            g = measure.splitter(charge, sign)
            gate = expm(-1j * t * h)
            np.testing.assert_allclose(
                gate / gate[0, 0], g[:2, :2] / g[0, 0], atol=1e-13
            )
            np.testing.assert_allclose(g @ position @ g.conj().T, position, atol=1e-15)
            np.testing.assert_allclose(np.sum(abs(g) ** 2, axis=1), np.ones(5))
    # Every phase has the same nonsingular, independently calibrated first
    # moment of the actual launch positions on this two-component label.
    x = measure.optics.beam.BEARINGS
    instrument = 2 * x.T @ abs(g) ** 2 @ x / 5
    assert np.linalg.eigvalsh(instrument).min() > 0.75


def test_moving_screen_control_with_known_accelerations():
    rng = np.random.default_rng(24466)
    a, b = rng.normal(size=(2, 2, 2))
    rows = {}
    for tag in ("plus_plus", "minus_minus", "plus_minus"):
        j, d = rng.normal(size=(2, 2, 2))
        rows[tag] = [j, d, a @ d + b @ j]

    def enclose():
        return {
            k: [
                [
                    [(iv.mpf(float(v)) + iv.mpf(["-1e-13", "1e-13"])) for v in r]
                    for r in m
                ]
                for m in jet
            ]
            for k, jet in rows.items()
        }

    # numpy inverses introduce rounded coefficients in the planted data;
    # the intervals cover the rounding necessary to the true linear equation.
    def test():
        return certificate.moving_screen_test(
            {t: np.array(j, dtype=object) for t, j in enclose().items()}
        )

    control = test()
    assert all(not row["certified_full_rank"] for row in control.values())
    rows["plus_minus"][2][0, 0] += 1
    changed = test()
    assert (
        changed["1"]["certified_full_rank"] and not changed["2"]["certified_full_rank"]
    )


def test_gate_before_derivatives_and_independent_stress_elimination():
    cert = json.loads((measure.HERE / "certificate.json").read_text())
    nominal = json.loads((measure.HERE / "results.json").read_text())
    for charge in ["-1", "1"]:
        row = cert["charges"][charge]
        obs = nominal["strength"]["8"][charge]
        assert 0 not in parse_interval(row["orientation_minor"])
        u, v = map(parse_interval, row["required_nonaffinities"])
        assert u.b < v.a or v.b < u.a
        for direction in ("forward", "reverse"):
            assert any(
                p["certified_full_rank"]
                for p in row["moving_screen_tests"][direction].values()
            )
            settings = obs["settings"][direction]
            assert settings["maximum_corrected_probability_jet_error"] < 1e-13
            for key in ("expansion", "focusing"):
                np.testing.assert_allclose(
                    settings["corrected_phase_average"][key],
                    settings["sharp"][key],
                    rtol=1e-8,
                    atol=1e-8,
                )
            for tag in measure.PHASES:
                r = obs["settings"][direction]["channels"][tag]
                p = parse_interval(row["tags"][tag][direction + "_pixel00"])
                assert r["probability00"] in p
                hits, other = r["click_and_complement"]
                n = hits + other
                # Uniform 99% Hoeffding band for all 16 quoted charge,
                # direction and phase setting probabilities.
                half = np.sqrt(np.log(2 * 16 / 0.01) / (2 * n))
                assert abs(hits / n - r["probability00"]) < half
            plus, minus = [
                obs["settings"][direction]["channels"][t]["click_and_complement"][0] / n
                for t in ["plus_plus", "minus_minus"]
            ]
            assert abs(plus - minus) > 2 * half
        for tag in measure.PHASES:
            for key, reading in [
                ("odd_expansion", "expansion"),
                ("odd_clock_focusing", "focusing"),
            ]:
                assert obs["odd"][tag][reading] in parse_interval(row["tags"][tag][key])
        changed = nominal["strength"]["0"][charge]["odd"]
        assert abs(changed["plus_plus"]["focusing"]) < 1
        assert abs(obs["odd"]["plus_plus"]["focusing"]) > 10


def test_delayed_choice_keeps_every_pixel_ensemble_while_separating_the_images():
    from . import delayed

    data = json.loads((measure.HERE / "delayed.json").read_text())
    for row in data["charges"].values():
        for key, value in delayed.odd_minors(row["settings"]).items():
            observed = parse_interval(value)
            # Same outward calculation with the caller's working precision.
            assert observed.a <= parse_interval(row["odd_minors"][key]).b
            assert observed.b >= parse_interval(row["odd_minors"][key]).a
        assert 0 not in parse_interval(
            row["odd_minors"]["minus_minus/plus_minus/minus_plus"]
        )
        for side, reads in row["settings"].items():
            assert all(
                x["certified_full_rank"]
                for x in row["moving_screen_tests"][side].values()
            )
            reference = reads["plus_plus"]["hit_by_launch"]
            for setting in reads.values():
                np.testing.assert_allclose(
                    setting["hit_by_launch"], reference, rtol=1e-13, atol=1e-15
                )
