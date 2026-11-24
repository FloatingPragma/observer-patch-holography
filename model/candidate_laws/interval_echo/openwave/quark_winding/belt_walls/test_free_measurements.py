"""Independent finite checks of free binding, polarization and neutral-pair records."""

from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import (
    bulk_order,
    flavour_charge,
    free_spectrum,
    ground_turn,
    pinned_string,
    routed_order,
)
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls.test_walls import (
    compare,
)
from model.candidate_laws.interval_echo.openwave.quark_winding.written_junction import (
    integer_cut,
)
from model.executable_model import carrier

HERE = Path(__file__).parent


@pytest.mark.parametrize(
    "module",
    [
        bulk_order,
        flavour_charge,
        free_spectrum,
        ground_turn,
        routed_order,
        integer_cut,
        pinned_string,
    ],
)
def test_complete_measurement(module):
    data = json.loads(Path(module.__file__).with_suffix(".json").read_text())
    compare(json.loads(json.dumps(module.run())), data)


def test_trimer_band_against_full_diagonalization():
    car = carrier.cell()
    sector = BinarySector(car, 3)
    w, v = np.linalg.eigh(sector.hamiltonian(8))
    face = sector.agreed == 3
    stats = free_spectrum.face_band(car, 8)
    n = stats["faces"]
    assert n == sum(face) == 20
    assert w[n] - w[n - 1] >= stats["spectral_gap_lower"] - 1e-10
    smallest = np.linalg.svd(v[face, :n], compute_uv=False)[-1] ** 2
    assert smallest >= stats["face_probability_in_low_band_lower"] - 1e-10
    p = np.eye(len(w))[:, np.flatnonzero(face)]
    for t in (0.25, 4, 32):
        matrix = (v @ (np.exp(-1j * t * w)[:, None] * (v.T @ p)))[face]
        assert (
            np.linalg.svd(matrix, compute_uv=False)[-1] ** 2
            >= stats["face_probability_after_arbitrary_wait_lower"]
        )
    assert stats["energy_below_separated_at_least"] > 26
    data = json.loads((HERE / "free_spectrum.json").read_text())
    for r in data["finite_excitation_energy"]:
        assert abs(r["local_slack_minimum"]) < 1e-9
        assert (r["upper_per_material"], r["upper_per_mark"]) == (
            (198, 33) if r["kappa"] == 8 else (6, 1)
        )
    for r in data["boards"]:
        n = sum(r["loads"])
        raw = r["separated_read"]
        assert (
            sum(raw["total_both_odd_other"])
            == raw["shots_per_seam"] * r["cells"] * 63 // 2
        )
        assert abs(raw["energy"] - r["separated_mean"]) < 0.06
        if r["coupling"] == 8:
            assert r["entire_lower_band_wall_upper"] < 0.003
            assert r["compact_single_fragment_line"] >= 3 * n - 1e-9


def test_wall_incidence_bounds_on_actual_cluster():
    board = carrier.federation(0)
    owners = [0] + sorted(
        {
            w // 12
            for edge in board.external
            for v, w in (edge, edge[::-1])
            if v // 12 == 0
        }
    )[:2]

    def words(occupied):
        split = {}
        for j in occupied:
            c, i = divmod(j, 12)
            split[c] = split.get(c, 0) | (1 << i)
        return list(split.values())

    def accepted(occupied):
        groups = words(occupied)
        return len(groups) == 3 and all(w in selected for w in groups)

    def hops(occupied):
        return [
            occupied - {p} | {q}
            for p in occupied
            for q in board.neighbours[p]
            if q not in occupied
        ]

    selected = free_spectrum.local_data(3)[0] | free_spectrum.local_data(4)[0]
    rng = np.random.default_rng(141302)
    for loads, forward in (((3, 3, 4), 47), ((3, 4, 4), 49)):
        for _ in range(10):
            source = set()
            for owner, n in zip(owners, loads, strict=True):
                w = int(rng.choice(sorted(free_spectrum.local_data(n)[0])))
                source.update(owner * 12 + i for i in range(12) if w >> i & 1)
            assert accepted(source)
            out = [w for w in hops(source) if not accepted(w)]
            assert len(out) <= forward
            assert max(sum(accepted(z) for z in hops(w)) for w in out) <= 4


def test_routing_and_twist_control():
    _, ring, meta, graphs = bulk_order.geometry()
    assert len(graphs["completed_square"] - graphs["native_with_chords"]) == 6
    assert all(meta["rotation"][ring[i]] == ring[(i + 4) % 12] for i in range(12))
    for w in combinations(range(12), 4):
        phase = np.exp(2j * np.pi * sum(w) / 12)
        moved = np.exp(2j * np.pi * sum((j + 4) % 12 for j in w) / 12)
        assert abs(moved - phase * np.exp(2j * np.pi / 3)) < 1e-12
    # Verify an unavailable four-load gate directly as a matrix action.
    car = carrier.cell()
    sec = BinarySector(car, 4)
    u, v = min(graphs["completed_square"] - set(car.seams))
    m = min(set(car.neighbours[u]) & set(car.neighbours[v]))
    source = np.random.default_rng(13).normal(size=len(sec.states)).astype(complex)
    source /= np.linalg.norm(source)
    state = source.copy()
    first = routed_order.pairs(sec, u, m)
    routed_order.pulse(state, first, np.pi, 8)
    routed_order.pulse(state, routed_order.pairs(sec, m, v), 0.31, 8)
    routed_order.pulse(state, first, np.pi, 8)
    expected = source.copy()
    gate = expm(-0.31j * np.array([[8.5, -0.5], [-0.5, 8.5]]))
    left, right = routed_order.pairs(sec, u, v)
    expected[left], expected[right] = gate @ np.array([source[left], source[right]])
    np.testing.assert_allclose(state, expected, atol=1e-12)


def test_memory_histories_and_uniform_half_turn():
    data = flavour_charge.run()
    assert max(data["native_charge_commutators"]) == 0
    for panel in data["core_reads"]:
        p = panel["pairs"]
        assert sum(r["probability"] for r in panel["joint"]) == pytest.approx(1)
        assert all(
            r["compensating_marks"] + r["vacancies"] == p for r in panel["joint"]
        )
        assert panel["reverse_after_read"] == pytest.approx(2**-p)
    raw = integer_cut.run()
    assert {x["net"] for x in raw["all_two_end_assignments"]} == {-2, -1, 0, 1, 2}
    assert all(x["imaginary_contrast"] == 0 for x in raw["all_two_end_assignments"])
    for row in raw["panels"]:
        assert sum(row["y_counts"]) == raw["shots_per_panel"]
        true_y = row["y_probabilities"][0] - row["y_probabilities"][1]
        assert abs(true_y - row["measured_y"]) < 0.005
        assert (abs(true_y) > 0.6) == bool(row["width"])


def test_reference_string_against_full_four_load_hamiltonian():
    from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import (
        belt,
    )
    from model.candidate_laws.contrast_memory_spherical_v1 import law

    raw = json.loads((HERE / "pinned_string.json").read_text())
    case = next(
        r
        for r in raw["panels"]
        if (r["length"], r["positive_cores"], r["kappa"], r["rung_duty"])
        == (12, 1, 8, 1)
    )
    c = BinarySector(belt.ring(12), 4)
    h = c.hamiltonian(8)
    flags = np.arange(12) % 3 != 0
    contact = law.seam_term(8)
    potentials = []
    for s in c.states:
        ids = [law.SRC.index((int(i in s), int(r))) for i, r in enumerate(flags)]
        potentials.append(contact.diagonal()[ids].sum() - (8.5) * (4 + sum(flags)))
    h += np.diag(potentials)
    arcs = np.array([pinned_string.core_read(s, 12, 1) for s in c.states])
    low, high = case["arc_range"]
    far = arcs >= (low + high) / 2
    es = np.linalg.eigvalsh(h[np.ix_(far, far)])
    np.testing.assert_allclose(es[[0, -1]], case["stretched_interval"], atol=1e-8)
    e, v = np.linalg.eigh(h)
    tight = arcs == low
    groups = np.split(np.arange(len(e)), np.flatnonzero(np.diff(e) > 1e-8) + 1)
    maxima = max(
        np.linalg.svd(v[tight][:, g], compute_uv=False)[0] ** 2 for g in groups
    )
    assert maxima == pytest.approx(
        case["optimal_stationary_tight"]["all"]["tight_probability"], abs=1e-8
    )
    for entry in raw["panels"]:
        disjoint = entry["disjoint_input_sites"]
        assert disjoint["tight_words"] == 0
        assert disjoint["words"] == (0 if entry["length"] == 12 else 15)
        if entry["length"] == 15:
            assert disjoint["minimum_arc"] == entry["arc_range"][1]
    for panel in raw["domain_energy"]:
        for r in panel["rows"]:
            assert r["energy_change"] == r["predicted"]
            assert abs(r["measured_change"] - r["predicted"]) < 0.04
            assert sum(r["record"]) == panel["shots_per_rung"] * panel["length"]
            if r["shifted_loads"] == panel["length"] // 3:
                assert r["anomalous_gaps"] == 0 and r["energy_change"] < 0
