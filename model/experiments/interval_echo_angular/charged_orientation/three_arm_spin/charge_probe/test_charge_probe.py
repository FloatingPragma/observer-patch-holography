import json
from itertools import product

import numpy as np

from . import measure, tomography


def test_charge_conjugation_and_the_probe_as_an_independent_interferometer():
    assert measure.conjugation() == 99
    assert measure.physical_probe()["probe_dimer"] == [145, 144]
    for k, q, s in product((4, 8), (-1, 0, 1), (-1, 1)):
        # Retaining the source's initial load projects either branch of a
        # blank seam onto its diagonal return amplitude.
        def amp(n, m):
            if n == m:
                return 1
            return np.exp(-1j * (k * (n - m) ** 2 + 0.5) * measure.DWELL) * np.cos(
                measure.DWELL / 2
            )

        f, n = amp(1 + q, 1), amp(1 + q, 1 + s)
        for angle in (0, -np.pi / 2):
            p = measure.local_assay(q, s, k, angle)
            interfered = n * np.exp(1j * angle)
            np.testing.assert_allclose(
                p[:2],
                [abs(f + interfered) ** 2 / 4, abs(f - interfered) ** 2 / 4],
                atol=1e-14,
            )


def test_four_probe_signs_and_independent_quantum_channels():
    data = json.loads((measure.HERE / "results.json").read_text())
    for row in data["rows"]:
        assert all(sum(raw) == measure.SHOTS for raw in row["counts"])
        if not row["contact"]:
            lo, hi = row["phase99"]
            assert lo <= 0 <= hi
    for r in data["shifts"]:
        prediction = measure.DWELL * (
            (2 * r["kappa"] + 1)
            if r["charge"] == r["probe"]
            else (0.5 - 2 * r["kappa"])
        )
        assert r["shift99"][0] < prediction < r["shift99"][1]
    assert len(data["filtered_probe"]) == 4
    assert all(r["original_bell_witness_99"][1] < 0 for r in data["filtered_probe"])
    for k in (4, 8):
        assert (measure.native_graph(-1, k) != measure.native_graph(1, k)).nnz == 0


def test_prediction_from_ordinary_observer_counts():
    data = json.loads((measure.HERE / "tomography.json").read_text())
    training = tomography.trained_terms()
    assert len(data["observations"]) == 48
    for row in data["observations"]:
        cal = next(
            r
            for r in data["stocked_probe_calibration"]
            if r["charge"] == row["charge"] and r["probe"] == row["probe"]
        )
        factor = tomography.bounds(cal["counts"][0], data["calibration_shots"])
        b, m = tomography.INPUTS[row["input"]], tomography.ANALYZERS[row["analyzer"]]
        ci = tomography.predict(training, b, m) * factor
        held = tomography.bounds(row["counts"][0], data["shots"])
        np.testing.assert_allclose(ci, row["forecast99"], rtol=1e-9)
        assert sum(row["counts"]) == data["shots"]
        assert tomography.overlaps(ci, held)
        assert not tomography.overlaps(
            tomography.predict(training, b, np.zeros(3)) * factor, held
        )
