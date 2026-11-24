"""Independent clock, first-arrival and control checks on the deposited records."""

import json

import numpy as np
from scipy.linalg import expm

from . import experiment as e


def records():
    return json.loads((e.HERE / "results.json").read_text())


def test_native_belt_and_changed_law():
    g = e.geometry()
    assert len(set(g["perimeter"] + [g["parked_unit"]])) == 11
    blocks, inserts = e.calibration()
    for q in e.TAGS:
        free = e.generator(0, q, "sham", blocks, inserts)
        expected = (
            -(np.roll(np.eye(10), 1, axis=1) + np.roll(np.eye(10), -1, axis=1)) / 2
        )
        np.testing.assert_allclose(free, expected, atol=2e-13)
        np.testing.assert_allclose(
            e.generator(0, q, "coin_even", blocks, inserts), free
        )
        cut = e.generator(0, q, "empty_write", blocks, inserts)
        np.testing.assert_allclose(cut[2, [3, 1]], [0, -0.5], atol=2e-13)
        h = e.generator(8, q, "coin_even", blocks, inserts, 1)
        expected_gap = np.sqrt(1025) if q == 1 else 1.0
        np.testing.assert_allclose(np.diff(np.linalg.eigvalsh(h)), expected_gap)


def test_independent_two_leg_first_arrival():
    # Independently propagate subnormalised density matrices of missed detections.
    blocks, inserts = e.calibration()
    report = records()
    for bank in report["banks"][:4]:
        for q in e.TAGS:
            h = e.generator(8, q, bank["source"], blocks, inserts)
            u = expm(-1j * h)
            hits = []
            for start, target in [(0, 5), (5, 0)]:
                rho = np.zeros_like(u)
                rho[start, start] = 1.0
                p = [0.0]
                for _ in range(16):
                    rho = u @ rho @ u.conj().T
                    p.append(rho[target, target].real)
                    rho[target, :] = 0.0
                    rho[:, target] = 0.0
                hits.append(p)
            by16 = np.convolve(*hits)[:17].sum()
            read = bank["tags"][str(q)]["radar"][0]
            assert (
                abs(by16 - read["completed_by_deadline"][1] / report["shots"])
                < e.RADAR_EPS
            )
            assert (
                read["completed_by_deadline"][-1] + read["unreturned"]
                == report["shots"]
            )


def test_saved_controls_and_shared_metric():
    report = records()
    for c in report["reference_clock"]:
        assert (
            c["simultaneous_interval"][0] < c["spacing"] < c["simultaneous_interval"][1]
        )
    for row in report["timing_contrasts"]:
        d = np.array(row["background_subtracted_tag_difference"])
        bound = row["simultaneous_difference_error"]
        if row["contrast"] and row["source"] == "coin_even":
            assert d[2] < -0.5
        elif row["contrast"] and row["source"] == "empty_write":
            assert d[2] > bound
        else:
            assert np.max(abs(d)) < bound
    for bank in report["banks"]:
        for q in ("1", "2"):
            clocks = bank["tags"][q]["clock_bonds"]
            for i, c in enumerate(clocks):
                if c["frequency"] is None:
                    assert bank["source"] == "empty_write" and i == 2
                else:
                    expected = (
                        np.sqrt(1025)
                        if bank["contrast"]
                        and bank["source"] == "coin_even"
                        and q == "1"
                        and i in (1, 2, 7, 8)
                        else 1
                    )
                    assert (
                        c["frequency_enclosure"][0]
                        <= expected
                        <= c["frequency_enclosure"][1]
                    )
