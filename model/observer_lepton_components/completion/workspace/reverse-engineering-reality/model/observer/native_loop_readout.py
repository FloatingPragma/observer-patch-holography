"""Exact loop test for the coin/current information emitted by native repair.

The fixed observer sees two adaptive seams. An additional, explicitly ideal
protocol measures both coin-conditioned current laws at every labelled seam
on fresh identical total-two preparations. Those laws select integer edge
contrasts. This gate asks whether those selected contrasts contain a
nontrivial closed-loop transport, without supplying a link phase.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from model.common_source import slow_geometry
from model.observer import native
from model.observer.source_discrimination import circuit
from model.observer.source_scalar import compiler


HERE = Path(__file__).resolve().parent
DEFAULT_RECEIPT = HERE / "results" / "native_loop_readout.json"


def carrier() -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int, int], ...]]:
    source = slow_geometry.source()
    seams = tuple(source.seams())
    faces = tuple(source.oriented_faces())
    if (len(seams) != 30 or len(faces) != 20 or len(set(seams)) != 30
            or any(a >= b or not 0 <= a < b < 12 for a, b in seams)
            or any(len(set(face)) != 3 or any(not 0 <= x < 12 for x in face)
                   for face in faces)):
        raise ArithmeticError("native carrier incidence changed")
    edge_set = set(seams)
    if any(tuple(sorted((u, v))) not in edge_set
           for face in faces for u, v in zip(face, face[1:] + face[:1])):
        raise ArithmeticError("native face has an absent seam")
    return seams, faces


def cycle_basis(seams: tuple[tuple[int, int], ...]) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, int], ...]]:
    """A deterministic spanning-tree basis of all independent graph loops."""
    neighbors: dict[int, list[int]] = {i: [] for i in range(12)}
    for a, b in seams:
        neighbors[a].append(b)
        neighbors[b].append(a)
    parent = {0: None}
    queue = [0]
    for a in queue:
        for b in sorted(neighbors[a]):
            if b not in parent:
                parent[b] = a
                queue.append(b)
    if len(parent) != 12:
        raise ArithmeticError("native seam graph disconnected")
    tree = {tuple(sorted((a, b))) for a, b in parent.items() if b is not None}
    chords = tuple(edge for edge in seams if edge not in tree)
    cycles = []
    for a, b in chords:
        ancestors_a = []
        node = a
        while node is not None:
            ancestors_a.append(node)
            node = parent[node]
        ancestors_b = []
        node = b
        while node is not None:
            ancestors_b.append(node)
            node = parent[node]
        meet = next(node for node in ancestors_a if node in ancestors_b)
        path_a = ancestors_a[:ancestors_a.index(meet) + 1]
        path_b = ancestors_b[:ancestors_b.index(meet)]
        cycles.append(tuple(path_a + path_b[::-1] + [a]))
    if len(tree) != 11 or len(chords) != 19 or len(cycles) != 19:
        raise ArithmeticError("native cycle rank changed")
    return tuple(cycles), chords


def read_contrasts(loads: tuple[int, ...], seams: tuple[tuple[int, int], ...]) -> dict[tuple[int, int], int]:
    """Derive the integer edge field from both actual native current records."""
    if loads not in native.states():
        raise ValueError("native total-two load required")
    contrasts = {}
    for seam in seams:
        records = [native.native(loads, seam, coin)[1] for coin in (0, 1)]
        if tuple(record[0] for record in records) != (0, 1):
            raise ArithmeticError("native coin labels changed")
        contrasts[seam] = sum(record[1] for record in records)
    return contrasts


def integer_holonomy(contrasts: dict[tuple[int, int], int], closed: tuple[int, ...]) -> int:
    """Oriented additive loop value; each character maps zero to unit phase."""
    if len(closed) < 4 or closed[0] != closed[-1]:
        raise ValueError("a closed loop with at least three seams is required")
    total = 0
    for u, v in zip(closed, closed[1:]):
        edge = tuple(sorted((u, v)))
        if edge not in contrasts:
            raise ValueError("loop uses a seam without a measured contrast")
        total += contrasts[edge] if u < v else -contrasts[edge]
    return total


def binary_transfer(seams: tuple[tuple[int, int], ...]) -> tuple[list[list[int]], tuple[int, int, int]]:
    """Exact 60-branch transfer counts for the separately uniform-seam protocol."""
    states = native.states(binary=True)
    index = {loads: i for i, loads in enumerate(states)}
    counts = [[0] * len(states) for _ in states]
    for i, loads in enumerate(states):
        for seam in seams:
            for coin in (0, 1):
                counts[i][index[native.native(loads, seam, coin)[0]]] += 1
    if any(sum(row) != 60 for row in counts):
        raise ArithmeticError("native transfer row lost mass")
    if any(counts[i][j] != counts[j][i]
           for i in range(len(states)) for j in range(len(states))):
        raise ArithmeticError("uniform-seam native transfer is not reversible")
    triangle = next(((i, j, k) for i in range(len(states))
                     for j in range(i + 1, len(states)) if counts[i][j]
                     for k in range(j + 1, len(states))
                     if counts[j][k] and counts[k][i]), None)
    if triangle is None:
        raise ArithmeticError("native transfer triangle witness disappeared")
    return counts, triangle


def same_endpoint_phase_readout() -> dict:
    """Source-execute a declared eraser on the fixed observer's two-event tape.

    This is a phase-sensitive *path-record* control, not a carrier face-loop
    holonomy. The first coin is committed in the fixed observer; the eraser
    requires keeping that register coherent until its final source read.
    """
    loads = tuple(int(i in (1, 11)) for i in range(12))
    paths = [(first_weight * second_weight, last)
             for first_weight, first in native.tick(native.State(loads))
             for second_weight, last in native.tick(first)]
    if len(paths) != 4 or any(weight != Q(1, 4) for weight, _ in paths):
        raise ArithmeticError("fixed native path probabilities changed")

    # Little-endian tape: c1, two C current bits, c2, two C current bits.
    # The current-code dictionary is the one used by the source-executed C
    # dilation. The final source state is retained as a separate fibre.
    def code(word: tuple[tuple[int, int], ...]) -> int:
        return sum((coin << (3 * event))
                   + (circuit.CURRENT_CODES[current] << (3 * event + 1))
                   for event, (coin, current) in enumerate(word))

    endpoints: dict[tuple[int, ...], np.ndarray] = {}
    for weight, last in paths:
        assert weight == Q(1, 4)
        vector = endpoints.setdefault(last.loads, np.zeros(64, dtype=complex))
        record = code(last.prefix)
        if vector[record] != 0:
            raise ArithmeticError("native record histories collided before readout")
        vector[record] = .5  # the native isometry's positive square-root weight
    returning = [(code(last.prefix), last.prefix)
                 for _, last in paths if last.loads == loads]
    if returning != [
        (0, ((0, 0), (0, 0))),
        (21, ((1, -1), (0, 1))),
    ]:
        raise ArithmeticError("same-endpoint native record words changed")

    # Recode only the second returning word to |1>, then use the ordinary
    # source-compiled real mixer to read |0> +/- |1>. The scalar first-coin
    # phase is another executable source word. Neither service is selected
    # by the committed native observer.
    phase = [compiler.pulse("Z", math.pi / 2, [0]),
             compiler.pulse("I", -math.pi / 2, [0])]
    eraser = (circuit.basis_swap_word(returning[1][0], 1, 6)
              + [compiler.pulse("Y", math.pi / 4, [0])])
    source_rows = {}
    for label, word in (("base", eraser), ("first_coin_phase", phase + eraser)):
        outputs = {endpoint: compiler.execute(word, 6, vector)
                   for endpoint, vector in endpoints.items()}
        mass = sum(float(np.vdot(v, v).real) for v in outputs.values())
        target = (Q(1, 2), Q(0)) if label == "base" else (Q(0), Q(1, 2))
        returned = tuple(float(abs(outputs[loads][bit]) ** 2) for bit in (0, 1))
        remainder = sum(float(abs(value) ** 2)
                        for endpoint, vector in outputs.items()
                        for bit, value in enumerate(vector)
                        if endpoint != loads or bit not in (0, 1))
        if (abs(mass - 1) > 2e-10 or
                max(abs(returned[bit] - float(target[bit])) for bit in (0, 1)) > 2e-10 or
                abs(remainder - .5) > 2e-10):
            raise ArithmeticError("source pulse word lost a path or the phase contrast")
        source_rows[label] = {
            "exact_returning_plus_minus_probabilities": list(map(str, target)),
            "other_endpoint_probability": "1/2",
            "all_outcome_probability_mass_within_2e_10": abs(mass - 1) <= 2e-10,
            "source_readout_probability_residual_below_2e_10": all(
                abs(returned[bit] - float(target[bit])) <= 2e-10 for bit in (0, 1)),
        }
    # First-coin Z multiplies each complete committed branch by a scalar.
    # This identity holds for arbitrary source/reference input and survives
    # any later adaptive classical control on the first recorded coin.
    committed = {last.prefix: weight for weight, last in paths}
    if len(committed) != 4 or sum(committed.values(), Q(0)) != 1:
        raise ArithmeticError("the complete native current law changed")
    phase_branch_counts = {}
    for law in ("C", "E"):
        signs: dict[tuple[int, ...], set[int]] = {}
        branches = 0
        for initial in native.states():
            for _, first in native.tick(native.State(initial)):
                first_seam = native.E
                first_label = (first.prefix[0][1] if law == "C" else
                               initial[first_seam[0]] - initial[first_seam[1]])
                second_seam = native.selected_seam(first)
                for _, last in native.tick(first):
                    second_label = (last.prefix[1][1] if law == "C" else
                                    first.loads[second_seam[0]] - first.loads[second_seam[1]])
                    key = (first.prefix[0][0], first_label, *second_seam,
                           last.prefix[1][0], second_label)
                    signs.setdefault(key, set()).add(1 - 2 * first.prefix[0][0])
                    branches += 1
        if branches != 312 or any(len(values) != 1 for values in signs.values()):
            raise ArithmeticError("first-coin phase is not scalar on a complete committed branch")
        phase_branch_counts[law] = branches
    return {
        "first_and_second_seams_on_both_returning_paths": [list(native.E), list(native.E)],
        "returning_tape_words": [[list(record) for record in word] for _, word in returning],
        "returning_tape_codes": [value for value, _ in returning],
        "all_four_committed_word_probabilities": {
            str(code(word)): str(weight) for word, weight in sorted(committed.items())},
        "first_coin_scalar_phase_leaves_every_committed_CP_branch_unchanged": True,
        "complete_C_E_adaptive_branches_checked": phase_branch_counts,
        "source_pulse_counts": {"first_coin_phase": len(phase), "record_eraser": len(eraser)},
        "source_executed_eraser": source_rows,
        "read_code_is_a_declared_coherent_tape_analyzer_not_a_native_current": True,
        "phase_selected_by_committed_native_records": False,
        "carrier_face_loop_or_nontrivial_link_holonomy": False,
        "reversed_seam_pair_has_unit_wilson_product_for_any_U1_link_assignment": True,
        "physical_defect_or_mass_identification": False,
    }


FACE_ROUTES = (((0, 1), (1, 5)), ((0, 11), (11, 5)))
FACE_INITIAL = tuple(int(i in (0, 2)) for i in range(12))
FACE_TERMINAL = tuple(int(i in (2, 5)) for i in range(12))


def oriented_link_sum(links: dict[tuple[int, int], int], path: tuple[int, ...]) -> int:
    """Add declared quarter-turn link angles along an oriented carrier path."""
    seams, _ = carrier()
    if len(path) < 2 or any(type(port) is not int or not 0 <= port < 12 for port in path):
        raise ValueError("an oriented native port path is required")
    if any(edge not in seams or type(value) is not int
           for edge, value in links.items()):
        raise ValueError("integer quarter-turn phases on canonical native seams required")
    result = 0
    for a, b in zip(path, path[1:]):
        edge = tuple(sorted((a, b)))
        if edge not in seams:
            raise ValueError("path uses a seam outside the native carrier")
        result += (1 if a < b else -1) * links.get(edge, 0)
    return result


def gauge_shift(links: dict[tuple[int, int], int], potentials: tuple[int, ...]) -> dict[tuple[int, int], int]:
    """An endpoint rephasing changes each declared oriented link by dg."""
    seams, _ = carrier()
    if (len(potentials) != 12 or any(type(value) is not int for value in potentials)):
        raise ValueError("twelve integer port phases required")
    return {(a, b): links.get((a, b), 0) + potentials[b] - potentials[a]
            for a, b in seams}


def registered_face_interferometer() -> dict:
    """Compare two equal-length routes under a declared coherent tape interface.

    Native repair supplies every endpoint and current word. Existing source
    pulses prepare three fair control bits, assign a chosen path phase, erase
    the two successful path words, and read both quadratures. Coherent route
    control, the native-history lift and the link angles remain declared; this
    is an executable interferometer test, not a native holonomy discovery.
    """
    seams, faces = carrier()
    if ((0, 11, 5) not in faces or (0, 5, 1) not in faces
            or any(tuple(sorted(edge)) not in seams for route in FACE_ROUTES for edge in route)):
        raise ArithmeticError("registered two-face route changed")

    # Bits: route, coin_1, current_1[0:2], coin_2, current_2[0:2].
    # The three control bits are prepared with existing proper-code pulses.
    preparation = [compiler.pulse("Y", -math.pi / 4, [bit]) for bit in (0, 1, 4)]
    blank = np.zeros(128)
    blank[0] = 1
    controls = compiler.execute(preparation, 7, blank)
    endpoints: dict[tuple[int, ...], np.ndarray] = {}
    histories = []
    for route, edges in enumerate(FACE_ROUTES):
        for first_coin in (0, 1):
            middle, first = native.native(FACE_INITIAL, tuple(sorted(edges[0])), first_coin)
            for second_coin in (0, 1):
                terminal, second = native.native(middle, tuple(sorted(edges[1])), second_coin)
                input_code = route + 2 * first_coin + 16 * second_coin
                code = (route + 2 * first_coin + 4 * circuit.CURRENT_CODES[first[1]]
                        + 16 * second_coin + 32 * circuit.CURRENT_CODES[second[1]])
                if abs(controls[input_code] - 1 / math.sqrt(8)) > 2e-13:
                    raise ArithmeticError("coherent fair-route preparation changed")
                vector = endpoints.setdefault(terminal, np.zeros(128, complex))
                if abs(vector[code]) > 2e-13:
                    raise ArithmeticError("distinct native histories collided in the tape")
                vector[code] = controls[input_code]
                histories.append({
                    "route": route, "coins": [first_coin, second_coin],
                    "currents": [first[1], second[1]], "record_code": code,
                    "terminal_occupied_ports": [i for i, count in enumerate(terminal) for _ in range(count)],
                    "committed_probability": "1/8",
                })
    if (len(histories) != 8 or len({row["record_code"] for row in histories}) != 8
            or abs(sum(float(np.vdot(v, v).real) for v in endpoints.values()) - 1) > 2e-12):
        raise ArithmeticError("the complete native two-route history law changed")
    successful = [row for row in histories if tuple(row["terminal_occupied_ports"]) == (2, 5)]
    if ([(row["route"], row["coins"], row["currents"], row["record_code"])
         for row in successful] != [(0, [0, 0], [1, 1], 36),
                                    (1, [0, 1], [1, -1], 85)]):
        raise ArithmeticError("the equal-length same-endpoint face paths changed")
    if set(np.flatnonzero(abs(endpoints[FACE_TERMINAL]) > 1e-10)) != {36, 85}:
        raise ArithmeticError("unexpected native path reached the comparison endpoint")

    eraser = (circuit.basis_swap_word(36, 0, 7)
              + circuit.basis_swap_word(85, 1, 7))
    if any(abs(vector[0]) > 1e-12 or abs(vector[1]) > 1e-12
           for vector in endpoints.values()):
        raise ArithmeticError("eraser destination was not blank")
    route_ports = ((0, 1, 5), (0, 11, 5))
    zero: dict[tuple[int, int], int] = {}
    inserted = {(5, 11): 1}  # a declared nonnative quarter-turn control
    shifted = gauge_shift(inserted, tuple(range(12)))
    contrasts = read_contrasts(FACE_INITIAL, seams)
    cases = (("flat", zero), ("native_current_character", contrasts),
             ("inserted_link_phase", inserted), ("port_rephased_inserted_link", shifted))
    rows = {}
    cosine = (1, 0, -1, 0)
    sine = (0, 1, 0, -1)
    for name, links in cases:
        phases = tuple(oriented_link_sum(links, path) for path in route_ports)
        closed = (0, 1, 5, 11, 0)
        wilson_quarter_turns = oriented_link_sum(links, closed)
        if ((phases[0] - phases[1] - wilson_quarter_turns) % 4
                or (oriented_link_sum(links, closed[::-1]) + wilson_quarter_turns) % 4
                or any(oriented_link_sum(links, (a, b, a)) != 0 for a, b in FACE_ROUTES[0])):
            raise ArithmeticError("oriented geometric path composition failed")
        # The diagonal source word assigns exp(i*pi*p_route/2). It collapses
        # declared path angles to one chosen route pulse; it does not derive
        # individual link operators from native repair.
        phase_word = [compiler.pulse("Z", math.pi * (phases[0] - phases[1]) / 4, [0]),
                      compiler.pulse("I", math.pi * (phases[0] + phases[1]) / 4, [0])]
        reads = {}
        for quadrature, read_shift in (("X", 0), ("Y", -math.pi / 2)):
            analyzer = ([compiler.pulse("Z", read_shift / 2, [0]),
                         compiler.pulse("I", -read_shift / 2, [0])]
                        + [compiler.pulse("Y", math.pi / 4, [0])])
            word = phase_word + eraser + analyzer
            output = {terminal: compiler.execute(word, 7, vector)
                      for terminal, vector in endpoints.items()}
            probabilities = {}
            for terminal, vector in output.items():
                occupied = ",".join(str(i) for i, count in enumerate(terminal) for _ in range(count))
                retained = {}
                for code, amplitude in enumerate(vector):
                    value = float(abs(amplitude) ** 2)
                    if value < 1e-11:
                        continue
                    exact = Q(round(16 * value), 16)
                    if abs(value - float(exact)) > 2e-10:
                        raise ArithmeticError("source readout disagrees with exact quarter-turn law")
                    retained[str(code)] = str(exact)
                probabilities[occupied] = retained
            selected = output[FACE_TERMINAL]
            observed = tuple(float(abs(selected[bit]) ** 2) for bit in (0, 1))
            expected = ((Q(1 + (cosine if quadrature == "X" else sine)[wilson_quarter_turns % 4], 8),
                         Q(1 - (cosine if quadrature == "X" else sine)[wilson_quarter_turns % 4], 8)))
            if (any(abs(observed[bit] - float(expected[bit])) > 2e-10 for bit in (0, 1))
                    or any(abs(selected[bit]) > 2e-10 for bit in range(2, 128))
                    or sum(Q(value) for bins in probabilities.values() for value in bins.values()) != 1
                    or sum(Q(value) for endpoint, bins in probabilities.items()
                           if endpoint != "2,5" for value in bins.values()) != Q(3, 4)):
                raise ArithmeticError("source face-loop quadrature or complete outcome mass changed")
            reads[quadrature] = {
                "exact_success_endpoint_plus_minus": [str(x) for x in expected],
                "all_endpoint_and_read_code_probabilities": probabilities,
                "source_pulses": len(word),
            }
        rows[name] = {"route_quarter_turn_sums": list(phases),
                      "closed_loop_quarter_turns_mod_4": wilson_quarter_turns % 4,
                      "quadratures": reads}
    if (rows["flat"]["closed_loop_quarter_turns_mod_4"] != 0
            or rows["native_current_character"]["closed_loop_quarter_turns_mod_4"] != 0
            or rows["inserted_link_phase"]["closed_loop_quarter_turns_mod_4"] != 1
            or rows["port_rephased_inserted_link"]["closed_loop_quarter_turns_mod_4"] != 1
            or rows["inserted_link_phase"]["quadratures"]["Y"]["exact_success_endpoint_plus_minus"]
            != ["1/4", "0"]):
        raise ArithmeticError("flat or nonnative link control did not separate")
    return {
        "closed_registered_route": [0, 1, 5, 11, 0],
        "adjacent_registered_faces": [[0, 5, 1], [0, 11, 5]],
        "two_equal_length_routes": [list(path) for path in route_ports],
        "initial_occupied_ports": [0, 2], "success_terminal_occupied_ports": [2, 5],
        "coherent_route_and_two_coin_preparation_pulses": len(preparation),
        "all_eight_native_histories": histories,
        "successful_history_codes": [36, 85],
        "success_probability_before_eraser": "1/4",
        "other_terminal_probability": "3/4",
        "source_eraser_pulses": len(eraser),
        "readout_cases": rows,
        "same_committed_native_history_probabilities_under_every_phase_case": True,
        "link_angles_and_coherent_route_interface_declared": True,
        "route_phase_pulse_is_not_an_individual_edge_transport": True,
        "native_repair_source_executes_individual_link_unitaries": False,
        "native_current_character_has_nontrivial_holonomy": False,
        "source_selected_defect_holonomy": False,
    }


def registered_triangle_control() -> dict:
    """A one-face phase discriminator with an explicit empty-seam idle arm.

    The moving arms have lengths two and one. Both scheduled arms execute two
    native attempts by adding an empty (2,3) seam attempt after 0->5. Only
    one of its two coin-labelled idle histories is erased into the chosen
    comparison bin; the other remains a retained same-endpoint outcome.
    The arms have unequal counts of moving events despite equal attempt
    counts. Nothing in the committed native records fixes the idle or
    move-count dynamical phases, so the uniform minus sign is a supplied
    coherent control rather than a calibrated geometric holonomy.
    """
    seams, faces = carrier()
    if (0, 5, 1) not in faces or (2, 3) not in seams:
        raise ArithmeticError("registered triangle or empty idle seam changed")
    initial = tuple(int(i in (0, 7)) for i in range(12))
    terminal = tuple(int(i in (5, 7)) for i in range(12))
    route_edges = (((0, 1), (1, 5)), ((0, 5), (2, 3)))
    route_ports = ((0, 1, 5), (0, 5))
    blank = np.zeros(128)
    blank[0] = 1
    prepare = [compiler.pulse("Y", -math.pi / 4, [bit]) for bit in (0, 1, 4)]
    controls = compiler.execute(prepare, 7, blank)
    endpoints: dict[tuple[int, ...], np.ndarray] = {}
    histories = []
    for route, edges in enumerate(route_edges):
        for first_coin in (0, 1):
            middle, first = native.native(initial, edges[0], first_coin)
            for second_coin in (0, 1):
                final, second = native.native(middle, edges[1], second_coin)
                code = (route + 2 * first_coin + 4 * circuit.CURRENT_CODES[first[1]]
                        + 16 * second_coin + 32 * circuit.CURRENT_CODES[second[1]])
                control_code = route + 2 * first_coin + 16 * second_coin
                if abs(controls[control_code] - 1 / math.sqrt(8)) > 2e-13:
                    raise ArithmeticError("triangle route/coin preparation changed")
                vector = endpoints.setdefault(final, np.zeros(128, complex))
                if abs(vector[code]) > 2e-13:
                    raise ArithmeticError("triangle native histories collided in the tape")
                vector[code] = controls[control_code]
                histories.append({
                    "route": route, "coins": [first_coin, second_coin],
                    "currents": [first[1], second[1]], "record_code": code,
                    "terminal_occupied_ports": [i for i, count in enumerate(final) for _ in range(count)],
                    "committed_probability": "1/8",
                })
    if (len(histories) != 8 or len({row["record_code"] for row in histories}) != 8
            or abs(sum(float(np.vdot(v, v).real) for v in endpoints.values()) - 1) > 2e-12
            or set(np.flatnonzero(abs(endpoints[terminal]) > 1e-10)) != {5, 21, 36}):
        raise ArithmeticError("one-face native branch or idle law changed")
    eraser = circuit.basis_swap_word(36, 0, 7) + circuit.basis_swap_word(5, 1, 7)
    if any(abs(vector[0]) > 1e-12 or abs(vector[1]) > 1e-12 for vector in endpoints.values()):
        raise ArithmeticError("triangle eraser destination was not blank")
    uniform_minus = {edge: 2 for edge in seams}  # exp(i*pi) on every moving edge
    cases = (("flat", {}), ("uniform_minus_moving_link_control", uniform_minus))
    readouts = {}
    for name, links in cases:
        phases = tuple(oriented_link_sum(links, path) for path in route_ports)
        closed = oriented_link_sum(links, (0, 1, 5, 0)) % 4
        if ((phases[0] - phases[1] - closed) % 4
                or oriented_link_sum(links, (0, 5, 1, 0)) % 4 != (-closed) % 4
                or oriented_link_sum(links, (2, 3, 2)) != 0):
            raise ArithmeticError("triangle loop, reversal or idle composition changed")
        word_phase = [compiler.pulse("Z", math.pi * (phases[0] - phases[1]) / 4, [0]),
                      compiler.pulse("I", math.pi * (phases[0] + phases[1]) / 4, [0])]
        reads = {}
        for quadrature, shift in (("X", 0), ("Y", -math.pi / 2)):
            analyzer = [compiler.pulse("Z", shift / 2, [0]),
                        compiler.pulse("I", -shift / 2, [0]),
                        compiler.pulse("Y", math.pi / 4, [0])]
            output = {final: compiler.execute(word_phase + eraser + analyzer, 7, vector)
                      for final, vector in endpoints.items()}
            probabilities = {}
            for final, vector in output.items():
                occupied = ",".join(str(i) for i, count in enumerate(final) for _ in range(count))
                bins = {}
                for code, amplitude in enumerate(vector):
                    value = float(abs(amplitude) ** 2)
                    if value < 1e-11:
                        continue
                    exact = Q(round(16 * value), 16)
                    if abs(value - float(exact)) > 2e-10:
                        raise ArithmeticError("triangle source readout differs from exact branch law")
                    bins[str(code)] = str(exact)
                probabilities[occupied] = bins
            observed = tuple(float(abs(output[terminal][bit]) ** 2) for bit in (0, 1))
            expected = ((Q(1, 4), Q(0)) if quadrature == "X" and closed == 0 else
                        (Q(0), Q(1, 4)) if quadrature == "X" and closed == 2 else
                        (Q(1, 8), Q(1, 8)))
            if (any(abs(observed[bit] - float(expected[bit])) > 2e-10 for bit in (0, 1))
                    or sum(Q(value) for bins in probabilities.values() for value in bins.values()) != 1
                    or sum(Q(value) for code, value in probabilities["5,7"].items()
                           if code not in ("0", "1")) != Q(1, 8)):
                raise ArithmeticError("triangle source quadrature or retained mass changed")
            reads[quadrature] = {
                "exact_selected_plus_minus": [str(value) for value in expected],
                "all_endpoint_and_read_code_probabilities": probabilities,
            }
        readouts[name] = {"route_quarter_turn_sums": list(phases),
                          "one_face_wilson_quarter_turns_mod_4": closed,
                          "quadratures": reads}
    if (readouts["flat"]["one_face_wilson_quarter_turns_mod_4"] != 0
            or readouts["uniform_minus_moving_link_control"]
            ["one_face_wilson_quarter_turns_mod_4"] != 2):
        raise ArithmeticError("odd-face uniform-sign control did not separate")
    return {
        "registered_face": [0, 5, 1],
        "comparison_loop": [0, 1, 5, 0],
        "moving_routes": [list(path) for path in route_ports],
        "native_attempt_counts": [2, 2],
        "moving_event_counts_on_selected_arms": [2, 1],
        "second_attempt_on_short_route": [2, 3],
        "empty_idle_current_on_both_coins": 0,
        "initial_occupied_ports": [0, 7],
        "comparison_terminal_occupied_ports": [5, 7],
        "all_eight_native_histories": histories,
        "selected_path_codes": [36, 5],
        "unselected_idle_coin_code_at_same_terminal": 21,
        "selected_interfering_mass": "1/4",
        "unselected_same_terminal_mass": "1/8",
        "all_other_terminal_mass": "5/8",
        "readout_cases": readouts,
        "uniform_sign_gives_two_face_loop_unit_but_one_face_minus_unit": True,
        "idle_phase_assumed_zero_in_declared_coherent_lift": True,
        "idle_and_move_count_dynamical_phases_independently_calibrated": False,
        "phase_difference_identified_as_geometric_holonomy": False,
        "source_selected_one_face_holonomy": False,
    }


def compact_interferometer_receipt(row: dict) -> dict:
    """Bind every computed outcome without repeating complete tables in JSON."""
    histories = row.pop("all_eight_native_histories")
    if (len(histories) != 8 or len({history["record_code"] for history in histories}) != 8
            or sum((Q(history["committed_probability"]) for history in histories), Q(0)) != 1):
        raise ArithmeticError("interferometer receipt lost a native history")
    row["complete_native_history_distribution_sha256"] = hashlib.sha256(
        canonical(histories).encode("utf-8")).hexdigest()
    row["native_history_count"] = 8
    row["complete_committed_probability_mass"] = "1"
    for case in row["readout_cases"].values():
        for read in case["quadratures"].values():
            probabilities = read.pop("all_endpoint_and_read_code_probabilities")
            total = sum((Q(value) for bins in probabilities.values()
                         for value in bins.values()), Q(0))
            if total != 1:
                raise ArithmeticError("interferometer receipt discarded an outcome")
            read["complete_outcome_distribution_sha256"] = hashlib.sha256(
                canonical(probabilities).encode("utf-8")).hexdigest()
            read["nonzero_endpoint_and_read_code_outcomes"] = sum(
                len(bins) for bins in probabilities.values())
            read["complete_probability_mass"] = "1"
    return row


def run() -> dict:
    seams, faces = carrier()
    cycles, chords = cycle_basis(seams)
    binding = slow_geometry.registered_binding()
    if not binding["all_pass"] or binding["registered_gram_entries_checked_exactly"] != 144:
        raise ArithmeticError("native port chart lost its registered carrier binding")
    all_states = native.states()
    for loads in all_states:
        contrasts = read_contrasts(loads, seams)
        if any(contrasts[a, b] != loads[a] - loads[b] for a, b in seams):
            raise ArithmeticError("native branch records are not a port-potential gradient")
        if any(integer_holonomy(contrasts, loop) != 0 for loop in cycles):
            raise ArithmeticError("native edge contrast has a nonzero fundamental loop")
        if any(integer_holonomy(contrasts, face + (face[0],)) != 0 for face in faces):
            raise ArithmeticError("native edge contrast has a nonzero face loop")

    witness_loads = tuple(int(i in (0, 1)) for i in range(12))
    witness_contrasts = read_contrasts(witness_loads, seams)
    corrupted = dict(witness_contrasts)
    corrupted[chords[0]] += 1
    changed_cycles = [i for i, loop in enumerate(cycles)
                      if integer_holonomy(corrupted, loop) != 0]
    if len(changed_cycles) != 1 or changed_cycles[0] != 0:
        raise ArithmeticError("single-chord flux control was not detected")

    counts, triangle = binary_transfer(seams)
    i, j, k = triangle
    forward = counts[i][j] * counts[j][k] * counts[k][i]
    reverse = counts[i][k] * counts[k][j] * counts[j][i]
    if not forward or forward != reverse:
        raise ArithmeticError("native transition-loop balance changed")
    # This stronger binary witness returns to its own source configuration
    # along distinct native tape paths, one with two nonzero currents.
    collision_loads = tuple(int(i in (1, 11)) for i in range(12))
    collision_paths = {
        last.prefix: first_weight * second_weight
        for first_weight, first in native.tick(native.State(collision_loads))
        for second_weight, last in native.tick(first)
        if last.loads == collision_loads
    }
    expected_paths = (((0, 0), (0, 0)), ((1, -1), (0, 1)))
    if tuple(sorted(collision_paths)) != expected_paths or any(
        weight != Q(1, 4) for weight in collision_paths.values()
    ):
        raise ArithmeticError("same-endpoint record-labelled path witness changed")
    return {
        "schema": "oph.native_loop_readout.v2",
        "status": "NATIVE_CURRENT_CONTRAST_FLAT__COHERENT_HOLONOMY_OPEN",
        "source": {
            "native_rule": "model/observer/native.py: native",
            "carrier_incidence": "model/common_source/slow_geometry.py: source().seams(), oriented_faces()",
            "registered_carrier_binding": binding,
            "mass_or_particle_target_data_consumed": False,
        },
        "ideal_all_seam_readout": {
            "protocol": "both exact coin-conditioned current laws on each labelled seam, with fresh identical preparations",
            "native_states_exhausted": len(all_states),
            "branch_records_checked": len(all_states) * len(seams) * 2,
            "edge_contrast_identity": "j_0(a,b)+j_1(a,b)=n_a-n_b",
            "independent_fundamental_loops": len(cycles),
            "oriented_face_loops": len(faces),
            "all_integer_loop_sums": "0",
            "every_character_holonomy": "1 for every U(1) character of the measured additive integer contrast",
            "formal_connection": (
                "The measured gradient is a pure gauge under every character; "
                "Lean/Screen/SeamU1HolonomyClassification.lean proves the abstract "
                "coboundary_range_eq_chordHolonomyHom_ker criterion on the "
                "registered seam carrier. The Python record-to-gradient check is "
                "not itself a Lean theorem."
            ),
            "single_chord_plus_one_control": {
                "chord": list(chords[0]),
                "nonzero_fundamental_loop_indices": changed_cycles,
                "native_branch_packet": False,
            },
        },
        "uniform_seam_binary_transfer": {
            "protocol": "independent uniform choice of 30 seams and fair native coin",
            "state_count": len(counts),
            "nonzero_unordered_transition_pairs": sum(bool(counts[a][b]) for a in range(len(counts))
                                                      for b in range(a + 1, len(counts))),
            "all_transition_counts_symmetric": True,
            "example_closed_transition_triangle": list(triangle),
            "forward_reverse_count_products": [forward, reverse],
            "all_positive_transition_cycle_affinities": "0 because P(x,y)=P(y,x)",
        },
        "fixed_observer": {
            "first_seam": list(native.E),
            "possible_second_seams": [list(native.E), list(native.F)],
            "distinct_seam_count": 2,
            "enclosed_carrier_face_count": 0,
            "same_endpoint_path_control": {
                "initial_and_final_occupied_ports": [1, 11],
                "separate_retained_record_words": [[list(record) for record in path]
                                                   for path in expected_paths],
                "probability_each": "1/4",
                "relative_phase_readout_in_committed_tape": False,
            },
            "coherent_path_record_control": same_endpoint_phase_readout(),
        },
        "registered_face_interferometer": compact_interferometer_receipt(
            registered_face_interferometer()),
        "registered_triangle_control": compact_interferometer_receipt(
            registered_triangle_control()),
        "scope": (
            "Exact native current records select a flat additive port contrast, not a "
            "U(1) connection with nontrivial flux. The two-event observer does not "
            "traverse a carrier face. The ideal all-seam protocol and uniform-seam "
            "transfer are separate declared protocols; their flat/balanced outputs "
            "do not rule out a different coherent source-derived holonomy. A "
            "record-accessible relative phase in the fixed committed interface, "
            "source-selected link transport, protected flux sector and continuum "
            "spectral pole remain unconstructed. A declared coherent tape "
            "eraser distinguishes scalar branch phases on the two returning "
            "E/E paths; this is not a carrier-loop holonomy. On a separate "
            "registered two-face loop, exact native branch endpoints and a "
            "source-pulse eraser resolve both quadratures of a supplied "
            "quarter-turn route phase. The native current character is flat; "
            "the nonzero control phase is supplied and not a discovered defect. "
            "The one-face idle-arm control has two versus one moving events; "
            "its dynamical phase remains uncalibrated."
        ),
        "all_pass": True,
    }


def canonical(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", type=Path, default=DEFAULT_RECEIPT)
    group.add_argument("--write", type=Path)
    args = parser.parse_args()
    encoded = canonical(run())
    if args.write is not None:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(encoded)
    elif args.check.read_text(encoding="utf-8") != encoded:
        raise SystemExit("Native loop-readout receipt mismatch")
    print(json.dumps({"all_pass": True, "sha256": hashlib.sha256(encoded.encode()).hexdigest()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
