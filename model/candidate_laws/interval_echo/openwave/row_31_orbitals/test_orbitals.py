"""Rebuild the small-board spectra and decode the deposited clock and distance counts."""

import json
from itertools import permutations

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from model.experiments.interval_echo_encounter.conjugation import defect_generator

from .measure import HERE, spectral_bounds


def assembled_cell(background, amounts):
    board = cell()
    law = candidate.native_law
    kinetic, cost = law.local_term_parts()
    h = (kinetic + candidate.KAPPA * cost)[: law.N_SRC, : law.N_SRC]
    pairs = list(permutations(range(board.n), 2))
    lookup = {v: i for i, v in enumerate(pairs)}
    out = np.zeros((len(pairs), len(pairs)))
    for index, (x, y) in enumerate(pairs):
        loads = np.full(board.n, background)
        loads[x], loads[y] = amounts
        out[index, index] += (
            candidate.KAPPA * (6 - board.degree) @ (loads**2 - background**2)
        )
        for a, b in board.seams:
            p, q = loads[[a, b]]
            src = law.SRC.index((p, q))
            out[index, index] += h[src, src]
            if p != q:
                swapped = {a: b, b: a}
                dst = lookup[(swapped.get(x, x), swapped.get(y, y))]
                out[dst, index] += h[law.SRC.index((q, p)), src]
    return pairs, out


def test_complete_cell_multiset_spectrum_and_radial_gap():
    board = cell()
    pairs, h = assembled_cell(0, (1, 2))
    e, v = np.linalg.eigh(h - 240 * np.eye(len(h)))
    edges = set(board.seams)
    inside = np.array([tuple(sorted(p)) in edges for p in pairs])
    assert sum(e < -1e-10) == int(inside.sum()) == 60
    cap = spectral_bounds(board)["sectors"][1]["bound_band"][1]
    assert e[59] <= cap + 1e-12 and e[60] >= -1e-12
    record = json.loads((HERE / "results.json").read_text())
    target = record["isolated_cell"]["settings"]["registered"]
    source = inside / np.sqrt(sum(inside))
    weights = (v.T @ source) ** 2
    for i, energy in enumerate(target["energies"]):
        hit = abs(e - energy) < 1e-10
        np.testing.assert_allclose(
            weights[hit].sum(),
            target["radial_probabilities_of_lines"][i][0],
            atol=1e-12,
        )


def test_filled_host_correction_as_complete_matrix():
    board = cell()
    pairs, full = assembled_cell(1, (0, 2))
    idx = [a * board.n + b for a, b in pairs]
    id1 = np.eye(board.n)
    vac = defect_generator(board, -1, candidate.KAPPA, 0).toarray()
    excess = defect_generator(board, 1, candidate.KAPPA, 0).toarray()
    separate = np.kron(vac, id1) + np.kron(id1, excess)
    correction = full - separate[np.ix_(idx, idx)]
    e = np.linalg.eigvalsh(correction)
    values, multiplicities = np.unique(np.round(e, 10), return_counts=True)
    np.testing.assert_array_equal(values, [0, 15, 16])
    np.testing.assert_array_equal(multiplicities, [72, 30, 30])


def test_full_federation_thresholds_and_held_counts():
    j = json.loads((HERE / "results.json").read_text())
    for i, archived in enumerate(j["boards"]):
        assert archived == spectral_bounds(federation(i))
        for sector in archived["sectors"]:
            for low, high in sector["hydrogenic_windows"].values():
                assert sector["bound_band"][1] < low <= high < 0
    energy = np.array(j["relative_bin_energies"])
    held = (-9 <= energy) & (energy <= -2)
    for panel in j["twenty_cells"]["settings"].values():
        counts = np.asarray(panel["counts_by_clock_then_distance"])
        assert counts.sum() == j["shots_per_setting"]
        assert counts[held].sum() == panel["held_window_clicks"]
    native = j["twenty_cells"]["settings"]["registered"]
    assert native["bound_window_clicks"] > 0.98 * j["shots_per_setting"]
    assert native["adjacent_in_bound_window"] / native["bound_window_clicks"] > 0.998
    control = j["isolated_cell"]["settings"]["assigned_detector_control"]
    assert control["held_window_clicks"] / j["shots_per_setting"] > 0.60
    for ref in j["twenty_cells"]["storage_referenced_uniform_single_loads"].values():
        assert sum(ref["clock_bins"]) == j["shots_per_setting"]
        assert 2 * np.argmax(ref["clock_bins"]) == ref["peak_energy"]


def test_induced_energy_on_two_writes_and_independent_receiver():
    from scipy.linalg import eigh

    j = json.loads((HERE / "mediator.json").read_text())
    for level, s in enumerate(j["studies"]):
        board = federation(level)
        first = s["first"]
        row = next(v for v in s["samples"] if v["distance"] == 4)
        applied = [first, row["second"]]
        energy = []
        for writes in (applied[:1], applied[1:], applied):
            # Rebuild the whole one-position graph from its conducting seams.
            mat = np.zeros((board.n, board.n))
            for a, b in board.seams:
                if any({a, b} == set(pair) for pair in writes):
                    mat[a, a] -= candidate.KAPPA
                    mat[b, b] -= candidate.KAPPA
                else:
                    mat[a, a] += 0.5
                    mat[b, b] += 0.5
                    mat[a, b] -= 0.5
                    mat[b, a] -= 0.5
            energy.append(
                float(eigh(mat, subset_by_index=[0, 0], eigvals_only=True)[0])
            )
        np.testing.assert_allclose(min(energy[:2]) - energy[2], row["gain"], atol=1e-12)
        assert (
            s["independent_three_link_calibration"]["coulomb_dark_probability_floor"]
            > 0.999
        )
        assert all(
            v["fresh_meter_calibration"]["dark"] / j["shots"] < 0.004
            for v in s["samples"]
            if v["distance"] >= 4
        )
