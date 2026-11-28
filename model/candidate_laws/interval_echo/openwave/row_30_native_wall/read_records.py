"""Reconstruct the Row 30 pilot comparisons from the saved observer counts."""

import json
from math import acos, atan2, cos, pi
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def rate(bank):
    return (bank[0] - bank[1]) / sum(bank[:2])


def energy(record, originals):
    previous = None
    for hold in record["count_banks"]:
        x, y = hold["counts"]
        assert sum(x) == sum(y) == originals
        t = hold["t"]
        angle = atan2(rate(y), rate(x))
        if previous is not None:
            angle += 2 * pi * round((previous * t - angle) / (2 * pi))
        previous = angle / t
    lo, hi = record["phase_frequency_interval"]
    assert lo <= previous <= hi
    return previous


def wave(record, originals):
    banks = record["counts_right_left_unpaired_non_single_wall_upstream_refusals"]
    assert all(sum(b) == originals for b in banks)
    value = sum((b[0] - b[1]) / sum(b[:3]) for b in banks) / 2
    k = pi - acos(value)
    lo, hi = record["measured_wave_number_from_band_minimum"]
    assert lo <= k <= hi
    density = record["held_wall_density_counts_then_other_and_upstream_refused"]
    assert sum(density) == originals
    return k


def compare(xs, ys, rows, result):
    misses = 0
    for x, y, row in zip(xs[1:-1], ys[1:-1], rows[1:-1], strict=True):
        fitted = ys[0] + (ys[-1] - ys[0]) * (x - xs[0]) / (xs[-1] - xs[0])
        residual = y - fitted
        lo, hi = row[result]
        assert lo <= residual <= hi
        misses += int(lo * hi > 0)
    return misses


def check():
    receipt = json.loads((HERE / "records.json").read_text())
    originals = receipt["bank_originals"]
    unit = energy(receipt["bench"]["clock"], originals)
    native_modes = 0
    native_misses = [0, 0, 0]
    with np.load(HERE / "scout_counts.npz", allow_pickle=False) as scout:
        for panel in receipt["runs"]:
            rows = panel["lines"]
            assert len(rows) == panel["dimers"] + 1
            energy(panel["independent_host_phase_read"], originals)
            raw = scout[panel["scout_count_array"]]
            assert raw.shape == (8192, 3) and raw.dtype.kind in "iu"
            assert np.all(raw >= 0) and np.all(
                raw.sum(axis=1) == receipt["scout_originals"]
            )
            measured = [energy(r["line_clock"], originals) / unit for r in rows]
            wave_numbers = [wave(r["spatial_read"], originals) for r in rows]
            recoil = [q**2 for q in measured]
            assert sorted(wave_numbers) == wave_numbers
            for value, row in zip(measured, rows, strict=True):
                lo, hi = row["line_clock"]["frequency_clock_units"]
                assert lo <= value <= hi
            misses = [
                compare(
                    [1 - cos(q) for q in wave_numbers],
                    measured,
                    rows,
                    "endpoint_calibrated_cosine_foil_error",
                ),
                compare(
                    [2 - 2 * cos(q) for q in wave_numbers],
                    recoil,
                    rows,
                    "endpoint_calibrated_lattice_squared_shell_error",
                ),
                compare(
                    [q * q for q in wave_numbers],
                    recoil,
                    rows,
                    "endpoint_calibrated_continuum_squared_shell_error",
                ),
            ]
            if panel["reference_load"] == 0:
                native_modes += len(rows)
                native_misses = [
                    a + b for a, b in zip(native_misses, misses, strict=True)
                ]
            if "direct_phase_reversal" in panel:
                final = {}
                for key, flights in panel["direct_phase_reversal"]["flights"].items():
                    for frame in flights:
                        a, b = frame["sign_arms"]
                        for branch in (a, b):
                            assert (
                                sum(branch["wall_address_counts_then_other"])
                                == originals
                            )
                            assert sum(branch["remote_equal_bond_counts"]) == originals
                        observed = (
                            a["remote_equal_bond_counts"][0]
                            - b["remote_equal_bond_counts"][0]
                        ) / originals
                        lo, hi = frame["remote_sign_difference_interval"]
                        assert lo <= observed <= hi
                    final[key] = observed
                assert final["native"] > 0.3 and abs(final["cut"]) < 0.000003
    assert native_modes == 21 and native_misses == [15, 15, 15]
    print("Verified 21 native phase and wavelength records from saved counts.")
    print("All 15 held interior intervals miss each endpoint-calibrated law.")
    print(
        f"Remote phase-reversal contrast: {final['native']:.7f}; cut {final['cut']:.7f}."
    )


if __name__ == "__main__":
    check()
