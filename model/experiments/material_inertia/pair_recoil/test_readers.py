import json

import numpy as np
from scipy import sparse
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import federation
from model.experiments.interval_echo_encounter.conjugation import free_seam
from .. import neutral_defect as nd
from . import clock, force, measure


def test_unfiltered_chebyshev_including_nonorthogonal_sources():
    h = np.array([[34.0, -0.5, 0.0], [-0.5, -14.0, -0.5], [0.0, -0.5, 32.0]])
    sources = np.array([[0.8, 0.1], [0.0, 0.3], [0.6, 0.9]])
    times = (0.01, 1.0, 12.0)
    measured = measure.propagation(sparse.csr_matrix(h), sources, times)
    for t in times:
        np.testing.assert_allclose(
            measured[t],
            expm(-1j * t * (h - 37 * np.eye(3))) @ sources,
            atol=3e-13,
            rtol=0,
        )


def test_independent_enrollment_and_cross_read():
    result = json.loads((measure.HERE / "results.json").read_text())
    pulse = json.loads((measure.HERE / "pulse.json").read_text())
    for survey, calibration in zip(
        result["placements"], pulse["placements"], strict=True
    ):
        J = survey["energy"]["matrix"][1][0]
        for acquisition in survey["energy"]["weak_reference"][:2]:
            time = acquisition["time"]
            a, b = acquisition["counts"]
            shot = a["click"] + a["quiet"]
            contrast = (
                (b["click"] - a["click"])
                / shot
                / (2 * time * measure.EPSILON * np.sqrt(1 - measure.EPSILON**2))
            )
            assert abs(contrast - J) < 4 * acquisition["standard_error"]
        held = next(
            r
            for r in survey["flights"]
            if r["time"] == 12 and r["preparation"] == "filtered"
        )
        assert held["imbalance"]["sampled"] > 3.5 * held["imbalance"]["standard_error"]
        for r in survey["flights"]:
            for group in r["records"]:
                assert np.sum(group["counts"]) == 2**36
        assert calibration["excess"][1]["cell_deflection"]["sampled"] > 0.75
        assert calibration["vacancy"][1]["cell_deflection"]["sampled"] < -0.75
        for row in calibration["sharp"]:
            assert abs(row["intact_imbalance"]["exact"]) < 1e-14
        assert calibration["filtered"][1]["signed_number_change"]["sampled"] < -0.00085


def test_all_cell_connected_energies_from_native_words():
    car = federation(0)
    edges, paths = force.cell_paths(car)
    h = free_seam(nd.KAPPA)

    def energy(occupancy):
        value = nd.KAPPA * ((6 - car.degree) @ (occupancy**2))
        for a, b in car.seams:
            word = law.SRC.index((occupancy[a], occupancy[b]))
            value += h[word, word]
        return value

    for row in json.loads((measure.HERE / "force.json").read_text())["records"]:
        src = []
        for i, gate in enumerate(row["gateways"]):
            occ = np.ones(car.n, dtype=int)
            regular = min(v for v in car.neighbours[gate] if car.degree[v] == 5)
            occ[gate], occ[regular] = (2, 0) if i and row["reverse"] else (0, 2)
            src.append(occ)
        connected = (
            energy(src[0] + src[1] - 1)
            - energy(src[0])
            - energy(src[1])
            + energy(np.ones(car.n, dtype=int))
        )
        assert abs(connected - row["expectation"]) < 1e-12
        assert (row["operator_norm"] == 0) == (row["distance"] > 1)
    assert len(paths) == 5 and len(edges) == 60


def test_two_physical_clock_types_reject_a_common_lapse():
    data = json.loads((measure.HERE / "clock.json").read_text())["records"]
    for row in data:
        hits = [r["click"] for r in row["counts"]]
        fitted, se = clock.fit(hits)
        assert abs(fitted - row["frequency"]) < 5 * se
        if row["distance"] > 1:
            baseline = next(
                r
                for r in data
                if r["distance"] == row["distance"]
                and r["clock_load"] == row["clock_load"]
                and r["source_gateway_deviation"] == 0
            )
            assert row["frequency"] == baseline["frequency"]
    near = {
        (r["clock_load"], r["source_gateway_deviation"]): r["sampled_frequency"]
        for r in data
        if r["distance"] == 1
    }
    assert abs(near[1, 1] / near[1, 0] - near[-1, 1] / near[-1, 0]) > 0.15


def test_admission_counts_recover_spectral_preparation():
    from . import filter_read

    primary = json.loads((measure.HERE / "results.json").read_text())
    assay = json.loads((measure.HERE / "filter.json").read_text())
    assert filter_read.run(primary) == assay
    for src, read in zip(primary["placements"], assay["placements"], strict=True):
        assert (
            abs(src["filter_gram"][0][1] - read["measured_cross_real"])
            < 5 * read["cross_standard_errors"][0]
        )
        assert abs(read["measured_cross_imag"]) < 5 * read["cross_standard_errors"][1]
