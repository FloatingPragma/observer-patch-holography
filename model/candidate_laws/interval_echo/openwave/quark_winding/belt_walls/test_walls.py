"""Independent sector, receiver and breakup checks for the committed measurements."""

from fractions import Fraction
import json
from math import comb
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.executable_model import carrier
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import (
    belt,
    composite,
    multiplets,
    intact,
    receivers,
    stationary,
)


HERE = Path(__file__).parent


def archive(name):
    return json.loads((HERE / (name + ".json")).read_text())


def compare(a, b):
    if isinstance(b, dict):
        assert a.keys() == b.keys()
        for k in b:
            compare(a[k], b[k])
    elif isinstance(b, list):
        for x, y in zip(a, b, strict=True):
            compare(x, y)
    elif isinstance(b, float):
        assert a == pytest.approx(b, abs=2e-8, rel=0)
    else:
        assert a == b


@pytest.mark.parametrize(
    "name,module",
    [
        ("belt", belt),
        ("receivers", receivers),
        ("composite", composite),
        ("multiplets", multiplets),
        ("intact", intact),
    ],
)
def test_complete_receipt(name, module):
    compare(json.loads(json.dumps(module.run())), archive(name))


def test_charge_winding_and_independent_derivative():
    records = archive("belt")["experiments"]
    for row in records.values():
        length, n = row["length"], row["loads"]
        q = Fraction(3 * n - length, 3)
        sec = BinarySector(belt.ring(length), n)
        basis = [sec.index[tuple(s)] for s in row["wall_addresses"]]
        for i in basis:
            read = belt.static_read(sec.states[i], length)
            assert (
                np.argmax(read["right_clock"]) - np.argmax(read["left_clock"])
            ) % 3 == 1
            assert all(
                Fraction(c["sum_above_bulk"], c["unit_calibration"]) == q
                for c in read["cuts"]
            )
        assert len(row["native_hops"]) == length
        for i, j, dm, dw in row["native_hops"]:
            assert sec.K[j, i] == -0.5 and Fraction(int(dm), int(dw)) == q
        init, other = (sec.index[tuple(s)] for s in row["launch"])
        psi = np.zeros(len(sec.states), complex)
        psi[[init, other]] = [1, 1j]
        psi /= np.linalg.norm(psi)
        derivative = 2 * (psi.conj() * (-1j * sec.hamiltonian(row["kappa"]) @ psi)).real
        assert derivative[other] == pytest.approx(0.5)
        assert derivative[init] == pytest.approx(-0.5)
        assert np.count_nonzero(abs(derivative) > 1e-12) == 2
        p = abs(expm_multiply(-1j * sec.hamiltonian(row["kappa"]), psi)) ** 2
        assert 1 - p[basis].sum() == pytest.approx(
            row["reads"][4]["outside"][0], abs=1e-9
        )
        assert max(p[basis].sum() - 1, 0) < 1e-11


def test_voltage_and_attached_unit():
    j = archive("receivers")
    assert j["unit_attachment"] == {"added": [1], "removed": []}
    for n, expected in [("3", -1), ("4", 2)]:
        row = j["cases"][n]
        for cut in row["contrast_and_voltage"]:
            assert cut["relative_excess"] == pytest.approx(expected / 3)
            assert cut["regional_full_turn_bright_x"] == 1
            assert cut["sharp_fraction_full_turn_bright_x"] == 0.25
        tour = row["written_tour"]
        assert sum(tour["completed_prefix_counts"].values()) == 458
        assert tour["completed_prefix_counts"]["10"] == 0
        assert (
            tour["unit_material_crossings"]
            == tour["independent_head_crossings"]
            == expected
        )
        assert tour["unwritten_control"]["head_crossings"] == 0


def test_three_wall_bounds_and_independent_small_ring_spectrum():
    raw = archive("composite")
    for n, minimum in [("10", 14), ("11", 18)]:
        native = raw["native"]["sectors"][n]
        assert native["placements"] == 3000 and native["min_unique_exits"] == minimum
        for summary in native["orderings"]:
            assert sum(g["count"] for g in summary["placements_by_distance"]) == 1000
    records = archive("stationary")
    for r in records:
        n, l = r["loads"], r["length"]
        dims = r["momentum_dimensions"]
        assert 2 * sum(dims) - dims[0] - (dims[-1] if l % 2 == 0 else 0) == comb(l, n)
        for v in r["maxima"].values():
            assert v["probability"] <= 1 + 1e-10
        if l == 12:
            compare(stationary.scan(l, n - l // 3 + 1), r)
            sec = BinarySector(belt.ring(l), n)
            energies, vectors = np.linalg.eigh(sec.hamiltonian(8))
            ones = n - l // 3 + 1
            arcs = np.asarray(
                [composite.gap_arc(s, l, ones, 3 - ones) for s in sec.states], float
            )
            for key, mask in {
                "minimum_arc": arcs == np.nanmin(arcs),
                "three_walls": np.isfinite(arcs),
            }.items():
                max_weight = max(
                    np.linalg.norm(vectors[mask][:, ix], 2) ** 2
                    for ix in multiplets.groups(energies)
                )
                assert max_weight == pytest.approx(
                    r["maxima"][key]["probability"], abs=2e-9
                )
    # Check the layer error bound with an independent full native cell and projector.
    cell = carrier.cell()
    sec = BinarySector(cell, 3)
    allowed = {composite.mask(f) for f in cell.faces()}
    source = min(allowed)
    flight = composite.closed_ball(source, composite.hop_list(cell.seams), 8, allowed)
    psi = np.array([int(composite.mask(s) == source) for s in sec.states], float)
    outside = [i for i, s in enumerate(sec.states) if composite.mask(s) not in allowed]
    for rd in flight["reads"]:
        prob = (
            abs(expm_multiply(-1j * rd["time"] * sec.hamiltonian(8), psi)[outside]) ** 2
        )
        lo, hi = rd["full_outside_bounds"]
        assert lo - 1e-12 <= sum(prob) <= hi + 1e-12


def test_rotating_wall_records():
    for n, r in archive("intact").items():
        assert r["exhaust_intact_minimum"] and all(
            abs(e["spatial_character_norm"] - 1) < 1e-9 for e in r["levels"]
        )
        certificate = r["exact_colour_certificate"]
        assert (certificate["spatial_irrep_copies"], certificate["krylov_rank"]) == (
            (60, 60) if n == "3" else (132, 131)
        )
        assert max(map(int, r["all_energy_commutant_dimensions"])) < 9
        assert r["reads"][-1]["any_chart"] > 0.97
        assert r["core_pair_multiplicity"] == ({"3": 20} if n == "3" else {"1": 30})
