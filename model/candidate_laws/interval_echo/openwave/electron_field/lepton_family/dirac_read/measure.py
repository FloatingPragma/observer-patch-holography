"""Compare two incidence Dirac lifts against the actual written-guide end spectrum."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.bands.measure import DEGREE_FIVE_COST, HOP, bloch, verify_route
from model.executable_model.carrier import federation, seam_digest


def register_effect():
    cyclic = np.roll(np.eye(3), 1, axis=0)
    theta = 1/12 + 1/6 - 2/72
    z = np.exp(1j*theta)
    # Compression of the historical single-load register onto its three cyclic addresses.
    amplitude = (32*np.eye(3) + 16*np.sqrt(2)*(z*cyclic + z.conjugate()*cyclic.T))/(96*np.sqrt(3))
    effect = amplitude.conj().T @ amplitude
    vals = np.linalg.eigvalsh(effect)
    assert np.linalg.eigvalsh(amplitude).min() > 0
    np.testing.assert_allclose(vals.sum(), 2/9, atol=1e-13)
    return effect, vals


def incidence(k, cost):
    rows = [[-1, 1, 0], [0, -1, 1], [np.exp(1j*k), 0, -1]]
    bonds = np.sqrt(HOP)*np.array(rows)
    if cost:
        bonds = np.vstack([bonds, [0, np.sqrt(cost), 0]])
    np.testing.assert_allclose(bonds.conj().T @ bonds, bloch(k, cost)[0], atol=1e-13)
    return bonds


def generator(bonds, effects):
    nv, ne, nf = bonds.shape[1], bonds.shape[0], len(effects)
    link = np.kron(bonds, np.eye(nf))
    return np.block([[np.kron(np.eye(nv), effects), link.conj().T],
                     [link, -np.kron(np.eye(ne), effects)]])


def quotient(es):
    return float(sum(es)/sum(np.sqrt(es))**2)


def fixture(n, effects, masses, cost):
    low = bloch(0, cost)
    kinetic_floor = max(0., float(low[1][0]))
    # At zero onsite cost the constant vector is an exact zero of the incidence.
    if cost == 0:
        kinetic_floor = 0.
    curv = float(low[-1][0])
    rest = np.sqrt(kinetic_floor + masses*masses)
    error, dark_weight, dark_count = 0., 0., []
    for i in range(n):
        k = 2*np.pi*i/n
        b = incidence(k, cost)
        d = generator(b, effects)
        levels, vectors = np.linalg.eigh(d)
        kinetic = bloch(k, cost)[1]
        paired = np.sqrt(np.maximum(kinetic[:, None], 0.) + masses**2).ravel()
        extra = -masses if cost else np.array([])
        predicted = np.sort(np.concatenate([paired, -paired, extra]))
        error = max(error, float(np.max(abs(levels-predicted))))
        if cost:
            ids = [np.flatnonzero(abs(levels + m) < 1e-10) for m in masses]
            assert all(len(js) == 1 for js in ids)
            dark_weight = max(dark_weight, float(sum(np.linalg.norm(vectors[:9,js])**2 for js in ids)))
            dark_count.append(sum(len(js) for js in ids))
    assert error < 1e-9 and dark_weight < 1e-20
    _, sample, _, slope, _ = bloch(np.pi/2, cost)
    reverse = bloch(-np.pi/2, cost)
    energy = np.sqrt(sample[0]+masses*masses)
    np.testing.assert_allclose(slope[0], -reverse[3][0], atol=1e-12)
    group_velocity = (slope[0]/(2*energy)).tolist()
    rest_levels, rest_vectors = np.linalg.eigh(generator(incidence(0, cost), effects))
    vertex_weights = []
    for energy, mass in zip(rest, masses):
        indices = np.flatnonzero(abs(rest_levels-energy) < 1e-10)
        assert len(indices) == 1
        weight = float(np.linalg.norm(rest_vectors[:9, indices])**2)
        np.testing.assert_allclose(weight, (1+mass/energy)/2, atol=1e-10)
        vertex_weights.append(weight)
    reverse_velocity = reverse[3][0]/(2*np.sqrt(reverse[1][0]+masses*masses))
    np.testing.assert_allclose(reverse_velocity, -np.array(group_velocity), atol=1e-12)
    gaps = []
    # Two-step finite differences of the squared dispersion avoid catastrophic
    # cancellation for the smallest mass. At rest (E²)''/2 equals E E''.
    for delta in (1e-2, 5e-3):
        plus = np.sqrt(bloch(delta,cost)[1][0]+masses*masses)
        minus = np.sqrt(bloch(-delta,cost)[1][0]+masses*masses)
        gaps.append(((plus**2+minus**2-2*rest**2)/(2*delta**2)).tolist())
    np.testing.assert_allclose(gaps[-1], curv/2, rtol=2e-5)
    return {"degree_five_cost": cost, "physical_momenta": n, "charged_vertex_kinetic_floor": kinetic_floor,
            "lowest_positive_energies_by_character": rest.tolist(), "ratios": (rest/rest[0]).tolist(),
            "Q": quotient(rest), "vertex_weights_at_positive_poles": vertex_weights,
            "common_rest_energy_times_curvature": curv/2, "squared_dispersion_difference_checks": gaps,
            "group_velocities_at_pi_over_two": group_velocity, "group_velocities_at_minus_pi_over_two": reverse_velocity.tolist(),
            "additional_auxiliary_negative_modes_each_momentum": max(dark_count,default=0),
            "largest_summed_vertex_weight_of_extra_modes": dark_weight, "maximum_full_spectrum_error": error}


def longest_unpinned_route():
    candidates = []
    for level in range(3):
        board = federation(level)
        high = [i for i in range(board.n) if board.degree[i] == 6]
        graph = nx.Graph(board.seams).subgraph(high)
        for nodes in nx.connected_components(graph):
            endpoints = sorted(p for p in nodes if graph.degree[p] == 1)
            assert len(endpoints) == 2
            route = tuple(nx.shortest_path(graph, *endpoints))
            candidates.append((-len(route), route, level))
    _, route, level = min(candidates)
    return federation(level), route


def head_paths(route):
    """A set bit points to the lower physical address, as in Patch.charge_moves."""
    edges = [tuple(sorted(e)) for e in zip(route, route[1:])]
    rows = []
    for sign in (-1, 1):
        first = sum(1 << j for j, (a, _) in enumerate(edges) if (route[j] == a) == (sign == 1))
        words = [first ^ ((1 << j)-1) for j in range(len(route))]
        background = None
        enabled = 0
        head_vectors = []

        for pos, word in zip(route, words):
            heads = dict.fromkeys(route, 0)
            for j, (a, b) in enumerate(edges):
                head = a if word & (1 << j) else b
                heads[head] += 1
                if pos in (a, b) and (pos == head) == (sign == 1):
                    other = b if pos == a else a
                    assert word ^ (1 << j) == words[route.index(other)]
                    enabled += 1
            head_vectors.append(list(heads.values()))
            heads[pos] -= sign
            if background is None:
                background = heads
            assert background == heads
        assert enabled == 2*len(edges)
        # The incidence row at edge i is a new grade attached to its left port.
        # Give it exactly that port's charge and native head word. A cross-port
        # incidence entry then flips the native seam; the other entry is onsite.
        locations = np.r_[np.arange(len(route)), np.arange(len(edges))]
        charge = sign*np.eye(len(route), dtype=int)[locations]
        local_flux = np.array(head_vectors, dtype=int)[locations] - np.array(list(background.values()))
        np.testing.assert_array_equal(local_flux, charge)
        for cell in {p//12 for p in route}:
            includes = [p//12 == cell for p in route]
            np.testing.assert_array_equal(local_flux[:, includes].sum(axis=1), charge[:, includes].sum(axis=1))
        rows.append({"sign": sign, "auxiliary_anchor_ports":list(route[:-1]),
                     "max_regional_flux_minus_charge_on_all_grades": 0, "reference_head_counts_minus_charge": list(background.values()),
                     "initial_word": first, "written_words_in_port_order": words,
                     "enabled_directed_hops": enabled})
    return rows


def path_incidence(board, route):
    n = len(route)
    raw = np.zeros((n-1, n))
    for j in range(n-1):
        assert route[j+1] in board.neighbours[route[j]]
        raw[j, j:j+2] = [-np.sqrt(HOP), np.sqrt(HOP)]
    pinned = [i for i, p in enumerate(route) if board.degree[p] == 5]
    return np.vstack([raw, np.sqrt(DEGREE_FIVE_COST)*np.eye(n)[pinned]])


def physical_path(effects, masses):
    board, route = longest_unpinned_route()
    n = len(route)
    incidence_map = path_incidence(board, route)
    assert incidence_map.shape == (n-1, n)
    np.testing.assert_array_equal(incidence_map @ np.ones(n), 0)
    kinetic, spatial = np.linalg.eigh(incidence_map.T @ incidence_map)
    kinetic[0] = 0.
    assert kinetic[1] > 0
    field_compensated = head_paths(route)
    cells = list(dict.fromkeys(p//12 for p in route))
    assert len(cells) > 1
    cut_index = (n-1)//2
    open_h = generator(incidence_map, effects)
    cut_h = generator(np.delete(incidence_map, cut_index, axis=0), effects)
    hvals, hvec = np.linalg.eigh(open_h)
    cvals, cvec = np.linalg.eigh(cut_h)
    _, characters = np.linalg.eigh(effects)
    modes = []
    times = [1, 2, 4, 8, 12]
    for j, mass in enumerate(masses):
        iso = np.kron(np.ones(n)/np.sqrt(n), characters[:, j])
        rest_input = np.zeros(len(hvals), dtype=complex)
        rest_input[:3*n] = iso
        np.testing.assert_allclose(open_h @ rest_input, mass*rest_input, atol=1e-13)
        branch = np.sqrt(mass*mass+kinetic)
        assert max(min(abs(hvals-e)) for e in branch) < 1e-12
        local = np.zeros(len(hvals), dtype=complex)
        local[:3] = characters[:,j]
        coeffs = hvec.conj().T @ local
        cut_coeffs = cvec.conj().T @ local[:len(cvals)]
        readings = []
        for time in times:
            wave = hvec @ (np.exp(-1j*time*hvals)*coeffs)
            cut_wave = cvec @ (np.exp(-1j*time*cvals)*cut_coeffs)
            sites = wave[:3*n].reshape(n,3)
            arrival = sites[-1]
            probability = float(np.vdot(arrival, arrival).real)
            selected = float(abs(np.vdot(characters[:,j], arrival))**2)
            cut_prob = float(np.linalg.norm(cut_wave[3*(n-1):3*n])**2)
            np.testing.assert_allclose(probability, selected, atol=1e-13)
            expected_norm = np.linalg.norm(wave)**2
            assert abs(expected_norm-1) < 1e-12 and cut_prob < 1e-25
            readings.append({"time":time, "source_vertex":float(np.linalg.norm(sites[0])**2),
                             "distal_vertex_same_character":selected,
                             "other_vertices":float(np.linalg.norm(sites[1:-1])**2),
                             "auxiliary":float(np.linalg.norm(wave[3*n:])**2),
                             "cut_distal":cut_prob})
        first_energy = branch[1]
        first_vertex = spatial[:,1]*np.sqrt((1 + mass/first_energy)/2)
        first_aux = incidence_map @ first_vertex/(first_energy+mass)
        first_excited = np.kron(np.r_[first_vertex, first_aux], characters[:,j])
        np.testing.assert_allclose(open_h @ first_excited, first_energy*first_excited, atol=1e-12)
        # Choose the two positive-energy amplitudes to cancel only at the distal
        # port at preparation. This source is distributed on the path.
        weights = np.array([first_vertex[-1], -1/np.sqrt(n)])
        weights /= np.linalg.norm(weights)
        positive_reads = []
        for time in [0, *times]:
            wave = (weights[0]*np.exp(-1j*time*mass)*rest_input +
                    weights[1]*np.exp(-1j*time*first_energy)*first_excited)
            np.testing.assert_allclose(np.linalg.norm(wave), 1., atol=1e-12)
            distal = float(np.linalg.norm(wave[3*(n-1):3*n])**2)
            if time == 0:
                assert distal < 1e-25
            positive_reads.append({"time":time, "distal_charge_same_character":distal})
        modes.append({"rest_charged_vertex_energy":float(mass),
                      "local_source_weight_in_the_rest_mode":1/n,
                      "homologous_positive_spatial_levels":branch.tolist(), "local_source_reads":readings,
                      "two_positive_energy_mode_source_reads":positive_reads})
    extra_port = min(q for q in board.neighbours[route[0]] if board.degree[q] == 5)
    extended = path_incidence(board, (extra_port, *route))
    floor = np.linalg.eigvalsh(extended.T @ extended)[0]
    shifted = np.sqrt(floor + masses*masses)
    return {"board_cells":board.cells, "board_sha256":seam_digest(board), "ordered_ports":list(route),
            "ordered_distinct_cells":cells, "frozen_compensator_outside_active_path":True,
            "native_sign_sectors":field_compensated, "physical_degrees":board.degree[list(route)].tolist(),
            "charged_kinetic_spectrum":kinetic.tolist(), "middle_cut_seam":list(route[cut_index:cut_index+2]),
            "rest_ratios":(masses/masses[0]).tolist(), "rest_Q":quotient(masses), "modes":modes,
            "one_pinned_port_control":{"extra_port":extra_port,"charged_kinetic_floor":float(floor),
                                       "Dirac_threshold_ratios":(shifted/shifted[0]).tolist(),
                                       "Q":quotient(shifted)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    n, route = verify_route()
    effects, masses = register_effect()
    out = {"route": route, "source_phase": "2/9", "effect_eigenvalues_increasing": masses.tolist(),
           "comparison_mass_coupling": 1, "orbital_grading_components": "vertex and incidence row",
           "records": [fixture(n, effects, masses, cost) for cost in (DEGREE_FIVE_COST, 0)],
           "physical_degree_six_path": physical_path(effects, masses),
           "code_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    with args.output.open("x") as f:
        json.dump(out, f, indent=2, allow_nan=False)
        f.write("\n")


if __name__ == "__main__":
    main()
