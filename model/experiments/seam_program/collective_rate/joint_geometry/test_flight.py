import json
from fractions import Fraction as F
from pathlib import Path

import numpy as np

from model.experiments.seam_program.collective_rate.measure import spectrum
from model.experiments.seam_program.collective_rate.joint_geometry.measure import flight_records
from model.experiments.seam_program.collective_rate.joint_geometry.short_flight import (
    parity_terminals, shell_averages, short_coefficients,
)


def test_exact_signs_and_scalar_rank():
    coeff = short_coefficients(8)
    assert coeff["exchange"] == [0, 0, F(5, 8), F(-623, 32)]
    assert coeff["aligned"] == [0, 0, -1, F(20479, 576)]
    assert coeff["parity"] == [0] * 4
    averages = shell_averages(8)
    _, x, y = averages[1]
    _, u, v = averages[2]
    assert [item[0] for item in averages.values()] == [42, 40, 8]
    assert x*v - y*u == F(-767, 34560)
    changed = short_coefficients(0)
    assert (changed["exchange"][2], changed["aligned"][2]) == (F(-1, 24), 0)


def test_signed_thermal_records_and_material_clock():
    blocks = spectrum(8)
    stored = json.loads(Path(__file__).with_name("results.json").read_text())["by_contrast"]["8"]
    for terminals, name in ((None, "rows"), (parity_terminals(), "parity_route")):
        times = (1/32, 1/8)
        fresh = flight_records(blocks, 8, times, terminals)
        assert fresh["valid_signed_amplitude"] > .001
        for row, archived in zip(fresh["rows"], stored[name][::2], strict=True):
            for field in row:
                np.testing.assert_allclose(row[field], archived[field], rtol=1e-7, atol=9e-14)
            assert np.isclose(row["exposure"], 2*np.arcsin(np.sqrt(row["clock_transfer"])))
            for source in ("exchange", "aligned", "parity"):
                reference = np.array(row["receiver_empty_and_full"])
                change = np.array(row[f"receiver_{source}_derivatives"])
                positive = reference + fresh["valid_signed_amplitude"] * change
                negative = reference - fresh["valid_signed_amplitude"] * change
                assert np.all((positive >= 0) & (positive <= 1) & (negative >= 0) & (negative <= 1))
        row = fresh["rows"][0]
        t = row["exposure"]
        coefficients = short_coefficients(8, half_order=4, terminals=terminals)
        for key in ("exchange", "aligned", "parity"):
            series = sum(float(c)*t**(2*n) for n, c in enumerate(coefficients[key]))
            np.testing.assert_allclose(row[key + "_response"], series, rtol=.004, atol=1e-15)


def test_parity_survives_the_disjoint_sender():
    assert parity_terminals() == (6, 11)
    assert short_coefficients(8, 3, (6, 11))["parity"] == [0, 0, 0, F(-9137, 5760)]
    assert short_coefficients(0, 3, (6, 11))["parity"] == [0, 0, 0, F(-1, 5760)]
    report = json.loads(Path(__file__).with_name("results.json").read_text())["by_contrast"]["8"]
    modular = report["regional_modular"]
    row = report["parity_route"][2]
    forecast = row["aligned_response"] * modular["parity_response"] / (8*modular["aligned_per_energy"])
    assert 290_000 < row["parity_response"] / forecast < 300_000
