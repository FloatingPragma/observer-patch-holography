"""Independent checks on the closed source, native propagation and distal instrumentation."""

from itertools import product
import json
import numpy as np
from scipy.linalg import expm
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import federation
from . import experiment as ex, walls, wall_flight as flight


def test_closed_source_disappearance_requires_opposite_adjacent_kinks():
    vacuum = walls.reference(0)
    for state in product((0, 1), repeat=4):
        phase = np.array(state) ^ vacuum
        signs = np.roll(phase, -1) - phase
        assert signs.sum() == 0
        for i in range(4):
            after = phase.copy()
            after[i] ^= 1
            if np.any(after != np.roll(after, -1)) or not np.any(signs):
                continue
            assert sorted(signs.tolist()) == [-1, 0, 0, 1]
            assert signs[i - 1] == -signs[i] != 0


def test_guarded_ring_columns_from_individual_native_seams():
    board = federation(0)
    ring, anchors, wings = flight.extended_geometry()
    sites = ring + wings[0] + wings[1]
    native_edges = list(zip(ring, ring[1:] + ring[:1]))
    for anchor, wing in zip(anchors, wings, strict=True):
        native_edges += [(ring[anchor], wing[0]), *zip(wing, wing[1:])]
    assert len(set(sites)) == len(sites)
    assert all(b in board.neighbours[a] for a, b in native_edges)
    assert all(
        board.degree[a] == board.degree[b] for a, b in zip(sites[::2], sites[1::2])
    )
    states = list(product((0, 1), repeat=len(sites) // 2))
    words = [tuple(walls.loads(row)) for row in states]
    lookup = {w: i for i, w in enumerate(words)}
    native = law.seam_term(8)
    h = walls.operators(8, False, 3)
    for j in (0, 1, 19, 233, 574, 1023):
        w = words[j]
        column = np.zeros(len(words))
        for pa, pb in native_edges:
            a, b = sites.index(pa), sites.index(pb)
            source = law.SRC.index((w[a], w[b]))
            for dest in np.flatnonzero(native[:, source]):
                target = list(w)
                target[a], target[b] = law.SRC[dest]
                allowed = lookup.get(tuple(target))
                if allowed is not None:
                    column[allowed] += native[dest, source]
        np.testing.assert_allclose(column, h[:, j], atol=1e-12)


def test_regional_energy_continuity_and_support():
    h = ex.wing(8, 3, 0)
    rng = np.random.default_rng(19)
    psi = rng.normal(size=8) + 1j * rng.normal(size=8)
    psi /= np.linalg.norm(psi)
    dt = 1e-6
    for start in (1, 2):
        region = flight.regions(8, start)
        current = 1j * (h @ region - region @ h)
        assert np.linalg.norm(current) > 1
        assert np.max(abs(current - current.conj().T)) < 1e-10
        f, b = [expm(s * 1j * dt * h) @ psi for s in (-1, 1)]
        derivative = (np.vdot(f, region @ f) - np.vdot(b, region @ b)) / (2 * dt)
        np.testing.assert_allclose(derivative, np.vdot(psi, current @ psi), atol=1e-7)
    end = flight.regions(8, 2)
    # The farthest region includes half the last bond, and is blind to the first dimer.
    np.testing.assert_allclose(end, np.kron(np.eye(2), end[:4, :4]), atol=1e-12)


def test_exact_current_compression_to_the_wall_free_guide_sector():
    # In a three-dimer wing the two zero-internal-wall records differ by three flips.
    clean = np.array([2, 5])  # binary 010, 101
    for boundary in (0, 2):
        h = ex.wing(8, 3, boundary)
        for start in (1, 2):
            region = flight.regions(8, start)
            current = 1j * (h @ region - region @ h)
            np.testing.assert_allclose(current[np.ix_(clean, clean)], 0, atol=1e-12)


def test_recorded_energy_flow_and_cut_sham():
    short = json.loads((ex.HERE / "walls.json").read_text())
    for run in short["cases"][1]["rows"]:
        assert all(c["wing_energies"]["counts"][3] == 0 for c in run["clocks"])
    d = json.loads((ex.HERE / "wall_flight.json").read_text())
    sham, source = d["cases"]
    assert (
        source["both_wing_excitation"]["exact"][3]
        > 1000 * sham["both_wing_excitation"]["exact"][3]
    )
    for case in (sham, source):
        readings = {}
        for t in case["times"]:
            for i, region in enumerate(t["distal"]):
                for j, spec in enumerate(region["spectra"]):
                    counts = np.asarray(spec["measurement"]["counts"])
                    vals = np.asarray(spec["energies"])
                    assert counts.sum() == d["shots"]
                    readings[t["t"], i, j] = np.average(vals, weights=counts[:-1])
        for j in (0, 1):
            for i, t in enumerate((2, 4)):
                change = readings[t, i, j] - readings[0, i, j]
                assert change > 1.9 if case["wall_pair"] else abs(change) < 0.001


def test_unconditional_bond_budget_from_instrument_records():
    from . import energy_budget

    bank = json.loads((ex.HERE / "energy_budget.json").read_text())
    for case in bank["cases"]:
        for setting, report in zip(case["readings"], case["balance"], strict=True):
            sums = dict.fromkeys(("ring", "left", "right", "attachments"), 0.0)
            for bond in setting["terms"]:
                counts = bond["high_counts"]
                assert all(0 <= c <= bank["shots"] for c in counts)
                sums[bond["region"]] += (
                    (counts[1] - counts[0]) * bond["span"] / bank["shots"]
                )
            total = report["regions"]["total"]
            np.testing.assert_allclose(
                sum(sums.values()), total["measured_surplus"], atol=1e-12
            )
            assert (
                abs(total["measured_surplus"] - (8 * case["contrast"] + 1))
                < total["pooled_99_halfwidth"]
            )
    late = bank["cases"][0]["balance"][-1]["regions"]
    assert (
        1.7
        < 100
        * (late["left"]["measured_surplus"] + late["right"]["measured_surplus"])
        / 65
        < 1.8
    )
    assert energy_budget.SHOTS == bank["shots"]
