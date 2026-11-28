"""Compare two written charge paths with distinct fan readers and an intactness guard."""

import argparse
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations, pairwise, product
import json
from pathlib import Path
from runpy import run_path
import sys

import networkx as nx
import numpy as np

HERE = Path(__file__).resolve().parent
FIELD = HERE.parents[1]


def outer_sites(board, cell):
    return [p for p in range(12 * cell, 12 * (cell + 1)) if board.degree[p] == 6]


def geometry(board):
    pairs = {}
    triple_centers = []
    for c in range(board.cells):
        ports = outer_sites(board, c)
        assert len(ports) == 3
        adjacent = [
            pair
            for pair in combinations(ports, 2)
            if pair[1] in board.neighbours[pair[0]]
        ]
        assert len(adjacent) <= 1
        if adjacent:
            pairs[c] = adjacent[0]
        for center in range(12 * c, 12 * (c + 1)):
            if set(ports) <= {center, *board.neighbours[center]}:
                assert center not in ports
                triple_centers.append([center, *ports])

    external_pairs = nx.Graph()
    external_pairs.add_nodes_from(pairs)
    for a, b in board.external:
        if a in pairs.get(a // 12, ()) and b in pairs.get(b // 12, ()):
            external_pairs.add_edge(a // 12, b // 12)
    assert nx.is_forest(external_pairs)
    components = Counter(map(len, nx.connected_components(external_pairs)))
    return {
        "board_cells": board.cells,
        "cells_with_adjacent_external_pair": len(pairs),
        "links_for_external_hub_transit_charts": len(external_pairs.edges),
        "external_hub_transit_component_sizes": dict(components),
        "external_hub_transit_cycles": 0,
        "any_hub_four_face_support_can_enroll_three_external_ports": triple_centers,
    }


def ordering_gauge(path):
    """Prove the sign comparator is a diagonal gauge on every guarded site move."""
    cells = dict((port, c) for c, fan in enumerate(path.cells) for port in fan)
    edges = {tuple(sorted(e)) for e in path.bridges}
    for fan in path.cells:
        edges.update(tuple(sorted((fan[a], fan[b]))) for a, b in path.field.edges)
    adjacency = nx.Graph(edges)
    pairs = [
        pair
        for pair in combinations(sorted(cells), 2)
        if cells[pair[0]] != cells[pair[1]]
    ]
    tested = 0
    for pair in pairs:
        initial = 1 if cells[pair[0]] < cells[pair[1]] else -1
        for source, spectator in (pair, pair[::-1]):
            for target in adjacency[source]:
                if cells[target] == cells[spectator]:
                    continue
                assert (cells[source] < cells[spectator]) == (
                    cells[target] < cells[spectator]
                )
                final = (
                    1
                    if cells[min(target, spectator)] < cells[max(target, spectator)]
                    else -1
                )
                sign = (
                    -1 if min(source, target) < spectator < max(source, target) else 1
                )
                assert sign == initial * final
                tested += 1
    return {
        "fan_hubs": [fan[0] for fan in path.cells],
        "distinct_cell_position_pairs": len(pairs),
        "tested_guarded_directed_material_hops": tested,
    }


def cycle_word_control(board, terminal_fan):
    """Use a physical cycle that leaves the external-hub chart family."""
    owners = [132, 52, 40, 144, 296]
    owner_pairs = list(zip(owners, owners[1:] + owners[:1], strict=True))
    external = {
        (a // 12, b // 12): (a, b)
        for seam in board.external
        for a, b in (seam, seam[::-1])
    }
    bridges = [external[pair] for pair in owner_pairs]
    fans = []
    for i, c in enumerate(owners):
        needed = {bridges[i - 1][1], bridges[i][0]}
        choices = [
            terminal_fan(board, h)
            for h in range(12 * c, 12 * (c + 1))
            if needed <= {h, *board.neighbours[h]}
        ]
        chart = min(choices, key=lambda fan: (fan[0] not in needed, fan))
        assert needed <= set(chart)
        fans.append(list(chart))

    # The topology control follows a common circulation direction and checks
    # the bridge word. No internal field amplitude is assigned to these moves.
    start = (0, 1)
    positions = list(start)
    word = 0
    directions = [(1, 1, 2), (1, 2, 3), (1, 3, 4), (0, 0, 1), (1, 4, 0)]
    for passenger, source, destination in directions:
        assert positions[passenger] == source and destination not in positions
        assert destination == (source + 1) % len(owners)
        positions[passenger] = destination
        word ^= 1 << source
    assert tuple(positions) == start[::-1] and word == (1 << len(owners)) - 1
    # Explore a superset of native enabled moves. The actual head can only
    # restrict these paths. Ring moves do not touch the bridge word.
    first = (0, 1, 0)
    queue = deque([first])
    reached = {first}
    while queue:
        a, b, bits = queue.popleft()
        for i, (p, spectator) in enumerate(((a, b), (b, a))):
            for sense in (1, -1):
                q = (p + sense) % len(owners)
                if q == spectator:
                    continue
                crossed = p if sense == 1 else q
                target = (
                    (q, b, bits ^ (1 << crossed))
                    if i == 0
                    else (a, q, bits ^ (1 << crossed))
                )
                if target not in reached:
                    reached.add(target)
                    queue.append(target)
    recorded = {(min(a, b), max(a, b), word) for a, b, word in reached}
    assert len(recorded) == len(reached)
    return {
        "physical_fans": fans,
        "physical_bridge_order": [list(e) for e in bridges],
        "geometric_guarded_crossings": [
            [owners[src], owners[dst]] for _, src, dst in directions
        ],
        "exchange_and_park_final_bridge_word_xor": word,
        "complete_position_and_bridge_outputs_in_exhaustive_geometric_lift": len(
            recorded
        ),
        "outputs_with_both_assignments_of_material_origins": 0,
        "topology_and_head_tag_control_only": True,
    }


def junction_control(board, terminal_fan, archive):
    """Show what a third incident bridge changes at a real 20-cell junction."""
    hub = 18
    centers = (14, 20, 22)
    central = terminal_fan(board, hub)
    assert set(centers) <= set(central)
    exterior = {a: b for e in board.external for a, b in (e, e[::-1])}
    leaves = [exterior[p] for p in centers]
    a, b, empty = range(3)
    junction = 3
    swap = [(a, junction, empty), (b, junction, a), (empty, junction, b)]
    park = [(a, junction, empty, junction, a), (b, junction, b)]
    counts = []
    for paths in (swap, park):
        pos = [a, b]
        uses = Counter()
        for path in paths:
            passenger = pos.index(path[0])
            for source, dest in pairwise(path):
                assert pos[passenger] == source and dest not in pos
                assert junction in (source, dest)
                uses[min(source, dest)] += 1
                pos[passenger] = dest
        assert uses == {a: 2, b: 2, empty: 2}
        counts.append(pos)
    assert counts == [[b, a], [a, b]]
    asymmetric = [
        p
        for p in leaves
        if set(board.neighbours[p]).intersection(outer_sites(board, p // 12))
    ]
    symmetric = [terminal_fan(board, p) for p in leaves if p not in asymmetric]
    for fan in symmetric:
        assert tuple(board.degree[list(fan)]) == (6, 5, 5, 5, 5, 5)
    prove_leaf_reflection(archive)
    return {
        "central_ordinary_hub_fan": list(central),
        "external_contacts": [list(pair) for pair in zip(centers, leaves, strict=True)],
        "geometric_exchange_and_reference_return_every_bridge_head": True,
        "terminal_gateways_with_adjacent_second_gateway": asymmetric,
        "symmetric_terminal_fans": [list(fan) for fan in symmetric],
        "symmetric_leaf_vacuum_reflection_and_doublet_action_both_charge_signs": [
            1,
            [[0, -1], [-1, 0]],
        ],
        "topology_only_field_code_transport_untested": True,
    }


def prove_leaf_reflection(archive):
    search, witness = archive
    field = search.fan(4)
    heads = tuple(witness["release"]["head_sectors"][0])
    reflected = (0, 5, 4, 3, 2, 1)
    degrees = [sum(i in pair for pair in field.edges) for i in range(6)]
    assert {tuple(sorted((reflected[i], reflected[j]))) for i, j in field.edges} == set(
        field.edges
    )
    assert {frozenset(reflected[i] for i in f) for f in field.rings} == {
        frozenset(f) for f in field.rings
    }
    for background in (
        heads,
        tuple(deg - h for deg, h in zip(degrees, heads, strict=True)),
    ):
        assert tuple(background[i] for i in reflected) == background
        words = field.words[background]
        order = []
        for word in words:
            transformed = 0
            for (i, j), bit in field.bits.items():
                edge = tuple(sorted((reflected[i], reflected[j])))
                old_head = i if word & bit else j
                if reflected[old_head] == edge[0]:
                    transformed |= field.bits[edge]
            order.append(words.index(transformed))
        vac, code = search.vacuum_code(field, background)
        np.testing.assert_array_equal(vac[order], vac)
        np.testing.assert_array_equal(code[order], -code[:, [1, 0]])


def pulse(state, contact, fixture):
    """Use the native isolated contact at exposure pi/8."""
    first, second, word = state
    a, b = contact
    if (a in (first, second)) == (b in (first, second)):
        return state
    source, dest = (a, b) if a in (first, second) else (b, a)
    bit = fixture.bits[contact]
    if (fixture.sign == 1) != ((source == a) == bool(bit & word)):
        return state
    return (
        dest if first == source else first,
        dest if second == source else second,
        word ^ bit,
    )


def collision_control(fixture, left_head):
    if left_head:
        contacts = [(24, 35), (24, 34), (14, 24), (24, 29)]
    else:
        contacts = [tuple(sorted(e)) for e in ((14, 23), (14, 22), (14, 24), (14, 16))]
    exchange = contacts * 2
    return_without_swap = [edge for edge in contacts for _ in range(2)]
    assert Counter(exchange) == Counter(return_without_swap)
    a, b = fixture.bridge
    retained = np.array([[(-1) ** u, (-1) ** v] for u, v in product((0, 1), repeat=2)])
    retained = np.kron(retained, retained)  # source applies an exact factor 1/4
    swap = np.eye(4, dtype=int)[[0, 2, 1, 3]]
    records = []
    for should_swap, steps in ((False, return_without_swap), (True, exchange)):
        collisions = []
        words = []
        odd = []
        current = [(a, b, word) for _, word in fixture.seeds]
        for e in steps:
            current = [pulse(state, e, fixture) for state in current]
            collisions.append(sum(i // 12 == j // 12 for i, j, _ in current))
        for (i, j, w), (_, initial_w) in zip(current, fixture.seeds, strict=True):
            assert {i, j} == {a, b} and w == initial_w
            words.append(1)
            odd.append(1 if i < j else -1)
        assert odd == [1 - 2 * should_swap] * 16
        native = retained.T @ np.diag(words) @ retained // 16
        comparator = retained.T @ np.diag(odd) @ retained // 16
        np.testing.assert_array_equal(native, np.eye(4, dtype=int))
        if should_swap:
            assert max(collisions) == 16 and min(collisions) == 0
        else:
            assert max(collisions) == 0
        records.append(
            {
                "material_origins_exchanged": should_swap,
                "physical_pulses": [list(e) for e in steps],
                "collocated_charge_word_count_out_of_16_by_prefix": collisions,
                "both_gateways_and_complete_word_return_probability": "1",
                "native_code_amplitude": native.tolist(),
                "odd_worldline_comparator_code_amplitude": comparator.tolist(),
                "two_original_reference_bell_pass_probability": str(
                    Fraction(int(np.trace(native)) ** 2, 16)
                ),
                "two_crossed_reference_bell_pass_probability": str(
                    Fraction(int(np.trace(swap @ native)) ** 2, 16)
                ),
            }
        )

    identity, move = records
    fringe = {}
    for tag in ("native_code_amplitude", "odd_worldline_comparator_code_amplitude"):
        ref = np.array(identity[tag])
        arm = np.array(move[tag])
        bright, dark = (arm + ref) // 2, (arm - ref) // 2
        np.testing.assert_array_equal(bright.T @ bright + dark.T @ dark, np.eye(4))
        fringe[tag] = {
            "bright_probability_all_inputs": int((bright.T @ bright)[0, 0]),
            "dark_probability_all_inputs": int((dark.T @ dark)[0, 0]),
            "both_doublets_refused_probability": 0,
        }
    return {
        "sign": fixture.sign,
        "bridge_initial_head_at_left_before_conjugation": left_head,
        "source_refusal_in_ideal_prepared_ensemble": 0,
        "records": records,
        "matched_fringe": fringe,
    }


def run():
    sys.path.insert(0, str(HERE.parents[6]))
    from model.executable_model import federation
    from model.candidate_laws.interval_echo.openwave.electron_field.dressed_exchange.search import (
        Fixture,
    )
    from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.fixture import (
        FanPath,
        archived_fixture,
        terminal_fan,
    )

    boards = [federation(i) for i in range(3)]
    topology = [geometry(board) for board in boards]
    assert [
        len(row["any_hub_four_face_support_can_enroll_three_external_ports"])
        for row in topology
    ] == [8, 0, 0]
    with archived_fixture() as archive:
        gauges = [
            ordering_gauge(FanPath(archive, n, 1, asymmetric=True)) for n in (3, 4)
        ]
        junction = junction_control(boards[0], terminal_fan, archive)
    return {
        "geometry": topology,
        "ordering_gauge_on_actual_320_routes": gauges,
        "physical_cycle_bridge_which_way_control": cycle_word_control(
            boards[-1], terminal_fan
        ),
        "physical_three_arm_junction_geometry_control": junction,
        "matched_two_fan_collision_controls": [
            collision_control(Fixture(sign, head), head)
            for sign, head in product((1, -1), (True, False))
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = json.loads(json.dumps(run()))
    if args.output:
        with args.output.open("x") as f:
            f.write(json.dumps(record, indent=2, allow_nan=False) + "\n")
    else:
        compare = run_path(str(FIELD / "replay.py"))["compare"]
        compare(
            record,
            json.loads((HERE / "results.json").read_text()),
            "intact_two_core_path",
        )
        print(
            "Verified intact path restrictions and the matched written two-fan control."
        )


if __name__ == "__main__":
    main()
