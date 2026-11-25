"""Verify physical phase quadratures, dwell derivatives and the clock-rank witness."""

import hashlib
import json

import numpy as np
import sympy as sp
from mpmath import iv
from scipy.linalg import expm

from ...directional_stress.bridge import tangents
from . import counts, measure
from ...action_response.four_ticks import parse_interval


def test_native_dimer_frequency_and_changed_generator():
    blocks = measure.action.native.iso.ranging.parity_blocks
    for q in (-1, 1):
        for k in (7, 8):
            for parity in (1, -1):
                h = (
                    sp.Matrix(
                        2
                        * blocks(q, k, parity)[1:, 1:]
                        .transpose(0, 2, 1, 3)
                        .reshape(4, 4)
                    ).applyfunc(int)
                    / 2
                )
                gap = k + (1 - parity) // 2
                assert h * h == gap * h and h.rank() == 3
                projector = sp.eye(4) - h / gap
                # The far/input block has one accessible coin for each charge.
                transmission = projector[2:, :2]
                effect = 4 * transmission.H * transmission
                assert effect.rank() == 1 and effect**2 == effect
                t = 0.63
                u = expm(-1j * t * np.array(h, dtype=complex))[2:, :2]
                np.testing.assert_allclose(
                    u.conj().T @ u,
                    np.sin(t * gap / 2) ** 2 * np.array(effect, dtype=complex),
                    atol=1e-15,
                )


def test_independent_phase_and_clock_conversions():
    saved = json.loads((measure.HERE / "results.json").read_text())
    native = measure.action.native
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    clock = np.arange(1, 13) / 2
    states = tangents()
    for q in (-1, 1):
        waves = np.array(
            [counts.wave_columns(board, route, q, o) for o in (False, True)]
        )
        for y, port in enumerate(measure.PORTS[2:]):
            a, b = waves[:, :, y, :, 0, :], waves[:, :, y, :, 1, :]
            signals = {
                k: np.einsum("btyi,bij,btyj->t", a, x, b.conj())
                for k, x in states.items()
            }
            for name, z in measure.PHASORS.items():
                omega = 2 * np.angle(complex(z))

                def delay(signal, lag=0, rate=1):
                    time = lag + rate * clock
                    freq = signal * np.exp(1j * omega * time)
                    return np.real(time @ freq / freq.sum())

                row = saved["charges"][str(q)][str(port)][name]

                def check(key, x):
                    lo, hi = json.loads(row[key])
                    assert lo < x < hi

                ref = signals["reference"]
                check("delay", delay(ref))
                np.testing.assert_allclose(
                    delay(ref, lag=0.25), delay(ref) + 0.25, atol=1e-12
                )
                eps = 1e-6
                disp = (delay(ref, rate=1 + eps) - delay(ref, rate=1 - eps)) / (
                    2 * eps
                ) - delay(ref)
                if omega:
                    check("omega_delay_dispersion", disp)
                else:
                    assert abs(disp) < 1e-7
                for label in ("ordinary", "pairing"):
                    d = signals[label]
                    lo, hi = json.loads(row["delay_slopes"][label])
                    assert (
                        lo
                        < (delay(ref + eps * d) - delay(ref - eps * d)) / (2 * eps)
                        < hi
                    )


def test_trial_totals_and_simultaneous_coverage():
    saved = json.loads((measure.HERE / "results.json").read_text())
    data = json.loads((measure.HERE / "counts.json").read_text())
    archived = {q: arm["counts"] for q, arm in data["charges"].items()}
    frozen = json.dumps(archived, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(frozen).hexdigest() == "3b1e2d518007d776ebdd8d97c08508812ed9c104f022c3d5936d136d264ec2e7"
    assert counts.analyze_counts(archived) == data
    laws = counts.detector_probabilities()
    replay = counts.analyze_counts(counts.acquire_counts(laws))
    # Audit every raw bin in both sets against the native probability itself.
    # A union Bernstein bound covers all bins and both complete acquisitions.
    total_bins = sum(np.size(table) for group in laws.values() for table in group.values())
    log_tail = np.log(4 * total_bins / 0.01)
    for observed in (data, replay):
        shots = observed["shots_per_setting"]
        for q, arm in observed["charges"].items():
            for preparation, bins in arm["counts"].items():
                p = np.asarray(laws[q][preparation], dtype=np.longdouble)
                counts_array = np.asarray(bins, dtype=np.longdouble)
                bound = np.sqrt(2 * shots * p * (1-p) * log_tail) + 2 * log_tail / 3
                assert np.all(abs(counts_array - shots*p) <= bound)
        check_trial_coverage(saved, observed)


def check_trial_coverage(saved, data):
    assert (
        data["nominal_affine_clock_witnesses_99pct"]["vacancy"]["nonsingular_subsets"]
        == 6
    )
    assert (
        data["nominal_affine_clock_witnesses_99pct"]["excess"]["nonsingular_subsets"]
        == 4
    )
    shots = data["shots_per_setting"]
    for q, group in data["charges"].items():
        for arr in group["counts"].values():
            np.testing.assert_array_equal(
                np.array(arr).sum(axis=-1), np.full((12, 4), shots)
            )
        for y, settings in group["nominal_intervals_99pct"].items():
            for name, read in settings.items():
                cert = saved["charges"][q][y][name]
                for k in ("delay", "ordinary", "pairing"):
                    a, b = json.loads(read[k])
                    c, d = json.loads(
                        cert["delay"] if k == "delay" else cert["delay_slopes"][k]
                    )
                    assert max(a, c) < min(b, d)


def test_ordinary_calibration_and_free_clock_rank():
    saved = json.loads((measure.HERE / "results.json").read_text())
    group = saved["charges"]
    for q, low, large in [
        ("-1", "zero", "quarter_turn"),
        ("1", "rational_turn", "quarter_turn"),
    ]:
        row = group[q]["38"]
        assert (
            parse_interval(row[low]["collar_over_ordinary"]).b
            < parse_interval(row[large]["collar_over_ordinary"]).a
        )
    rows = {}
    for q, ys in group.items():
        for p, channels in ys.items():
            for n, item in channels.items():
                rows[q, p, n] = [
                    iv.mpf(1),
                    parse_interval(item["delay"]),
                    parse_interval(item["omega_delay_dispersion"]),
                    parse_interval(item["delay_slopes"]["ordinary"]),
                    parse_interval(item["pairing_minus_ordinary_over_volume_ratio"]),
                ]
    for label, data in saved["affine_clock_witnesses"].items():
        if not data["nonsingular_subsets"]:
            continue
        keys = list(map(tuple, data["rows"]))
        slopes = [
            [rows[k][j] - rows[keys[0]][j] for j in range(1, 5)] for k in keys[1:]
        ]
        assert 0 not in measure.determinant(slopes)
    assert saved["affine_clock_witnesses"]["vacancy"]["nonsingular_subsets"] == 5
    assert saved["affine_clock_witnesses"]["excess"]["nonsingular_subsets"] == 4
    assert saved["affine_clock_witnesses"]["both"]["nonsingular_subsets"] > 0
