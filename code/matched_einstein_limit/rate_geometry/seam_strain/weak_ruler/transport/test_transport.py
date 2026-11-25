"""Native transports, independent pointer controls and held comparisons."""

import hashlib
import json

import numpy as np

from . import measure as m


def test_rebuild_source_and_event_bounds():
    path = m.HERE / "results.json"
    contents = path.read_bytes()
    assert (
        hashlib.sha256(contents).hexdigest()
        == "3b5e6a0e7172c20a5ae88badb34bf0f898f410454f2efce6fac6e6daa0e9a9a4"
    )
    archive, fresh = json.loads(contents), m.run()
    assert len(fresh["rows"]) == len(archive["rows"]) == 8
    for row, old in zip(fresh["rows"], archive["rows"], strict=True):
        for key in ("charge", "ports", "mark_time", "arrival_time"):
            assert row[key] == old[key]
        for tag in m.ARMS:
            now, recorded = row["arms"][tag], old["arms"][tag]
            assert now["negative_shutter_half_width_lower"] > 3e-5
            for r in (now, recorded):
                assert r["conditional_clock"][1] < -0.024
                counts = r["counts"]
                assert all(sum(v) == m.SAMPLES and min(v) >= 0 for v in counts.values())
                band = m.count_clock(counts)
                assert band == r["joint_from_counts_99_family"]
                assert (
                    band[0]
                    < now["joint_clock"][0]
                    <= now["joint_clock"][1]
                    < band[1]
                    < 0
                )
            if tag == "reference" and row["mark_time"] == 4:
                assert (
                    now["minimum_relative_arrival_change_under_dephasing"][0] > 0.1209
                )
        for r in (row, old):
            transfer = r["fixed_volume_transfer"]["measured_minus_predicted"]
            assert transfer[0] > 0.0027 or transfer[1] < -0.0027


def test_variable_strength_pointer_and_dephased_mobile():
    data = json.loads((m.HERE / "results.json").read_text())
    ref, variation, _ = m.flight.live_survey.retained(True)
    states = {"reference": ref, **{k: ref + variation[k] for k in m.ARMS[1:]}}
    for charge in (-1, 1):
        # Independent floating Chebyshev implementation on the two arms of
        # the actual pointer. These amplitudes include every intermediate
        # port in the complement to the marked site via U_T - V.
        h, coordinates = m.flight.hardware(charge)
        kernels = m.live_cadence.kernel_grid(h, coordinates[:4], np.arange(13) / 2)
        for row in (r for r in data["rows"] if r["charge"] == charge):
            s, t = (int(2 * row[key]) for key in ("mark_time", "arrival_time"))
            a, middle, b = (m.PORTS.index(site) for site in row["ports"])
            full = kernels[t, b, a]
            through = kernels[t - s, b, middle] @ kernels[s, middle, a]
            avoid = full - through
            for tag, rho in states.items():
                bank = row["arms"][tag]

                def p(wave):
                    return np.trace(wave @ rho @ wave.conj().T).real

                expected, missed, arrival = (p(wave) for wave in (through, avoid, full))
                x, y = (np.mean(bank["sharp"][i]) for i in (0, 2))
                np.testing.assert_allclose(
                    [x, y, np.mean(bank["unmarked"][0])],
                    [expected, missed, arrival],
                    rtol=0,
                    atol=7e-13,
                )
                w = (arrival + expected - missed) / 2
                for eps in (1.0, 0.6, 0.1, 0.01):
                    plus, minus = np.sqrt((1 + eps) / 2), np.sqrt((1 - eps) / 2)
                    pp, pm = (
                        p(plus * through + minus * avoid),
                        p(minus * through + plus * avoid),
                    )
                    assert abs((arrival + (pp - pm) / eps) / 2 - w) < 2e-13
                    coherence = np.sqrt(1 - eps**2)
                    assert (
                        abs(
                            pp
                            + pm
                            - (expected + missed + 2 * coherence * (w - expected))
                        )
                        < 2e-13
                    )
                # Threshold dephasing at which the signed clock is zero.
                g = expected / (expected - w)
                assert 0 < g < 1
                new_arrival = expected + missed + 2 * g * (w - expected)
                assert abs((new_arrival / arrival - 1) + 2 * w / arrival) < 2e-12
