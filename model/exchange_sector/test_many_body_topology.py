"""Independent integral, path-lift and hostile-replay controls for the catalogue."""
from collections import defaultdict
from copy import deepcopy
from itertools import combinations
from math import comb, gcd

import pytest
import sympy as sp

from . import many_body_topology as topology
from ..clock_defect_energy import defect_history
from ..observer import native_loop_readout


@pytest.fixture(scope="module")
def certificates():
    return {n: topology.reduction_certificate(n) for n in range(13)}


@pytest.fixture(scope="module")
def cocycles():
    return {n: topology.edge_cocycle(n) for n in range(2, 11)}


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def independent_chain(g, polygon, by_subset=None):
    """Build signed edges from occupied subsets rather than producer addresses."""
    if by_subset is None:
        by_subset = {(g.states[a], g.states[b]): i for i, (a, b) in enumerate(g.edges)}
    result = defaultdict(int)
    for a, b in zip(polygon, polygon[1:]+polygon[:1]):
        left, right = g.states[a], g.states[b]
        if left == right:
            continue
        direction = 1 if left < right else -1
        result[by_subset[tuple(sorted((left, right)))]] += direction
    return {k: v for k, v in result.items() if v}


@pytest.mark.parametrize("n", range(13))
def test_complete_cell_census_and_integer_boundaries(n):
    g = topology.geometry(n)
    assert len(g.states) == comb(12, n)
    assert len(g.edges) == 30*choose(10, n-1)
    assert sum(k == "particle_triangle" for k, _ in g.cells) == 20*choose(9, n-1)
    assert sum(k == "disjoint_square" for k, _ in g.cells) == 315*choose(8, n-2)
    assert len(g.tree) == len(g.states)-1
    assert set(g.tree).isdisjoint(g.chords)
    assert set(g.tree)|set(g.chords) == set(range(len(g.edges)))
    for a, b in g.edges:
        assert tuple(sorted(set(g.states[a]) ^ set(g.states[b]))) in g.seams
    positions = {(g.states[a], g.states[b]): i for i, (a, b) in enumerate(g.edges)}
    for kind, polygon in g.cells:
        supports = [set(g.states[x]) for x in polygon]
        fixed = set.intersection(*supports)
        varying = set.union(*supports)-fixed
        assert (len(polygon), len(fixed), len(varying)) == (
            (3, n-1, 3) if kind == "particle_triangle" else (4, n-2, 4))
        incidence = defaultdict(int)
        for e, value in independent_chain(g, polygon, positions).items():
            a, b = g.edges[e]
            incidence[a] -= value
            incidence[b] += value
        assert not any(incidence.values())


@pytest.mark.parametrize("n", range(13))
def test_integer_replay_exhausts_all_occupation_sectors(n, certificates):
    result = topology.verify_reduction(n, certificates[n])
    assert result["free_rank"] == (19 if n == 11 else 0)
    assert result["torsion"] == ([2*(n-1)] if 2 <= n <= 10 else [])
    assert result["unit_pivots"] == len(topology.geometry(n).chords)-(1 if 2 <= n <= 10 else 19 if n == 11 else 0)


def test_two_particle_domain_matches_frozen_predecessor():
    g = topology.geometry(2)
    states, edges, cells, _, _, _ = defect_history.geometry()
    old_edges = {frozenset((states[a], states[b])) for a, b in edges}
    new_edges = {frozenset((g.states[a], g.states[b])) for a, b in g.edges}
    assert old_edges == new_edges
    assert {frozenset(states[i] for i in p) for _, p in cells} == {
        frozenset(g.states[i] for i in p) for _, p in g.cells}


def test_complete_native_endpoint_graph_matches_all_thirteen_complexes():
    from . import many_body
    movements = {n: set() for n in range(13)}
    for mask in range(4096):
        for seam in many_body.SEAMS:
            for coin in (0, 1):
                _, output, _, _ = many_body.column(mask, "C", seam, coin)
                if output != mask:
                    movements[mask.bit_count()].add(tuple(sorted((mask, output))))
    for n, actual in movements.items():
        g = topology.geometry(n)
        addresses = tuple(sum(1 << p for p in state) for state in g.states)
        assert actual == {tuple(sorted((addresses[a], addresses[b]))) for a, b in g.edges}


def test_independent_sparse_integer_substitution_replays_two_particle_trace(certificates):
    """No producer/replay matrix helper is used in this alternate oracle."""
    g = topology.geometry(2)
    variable = {edge: k for k, edge in enumerate(g.chords)}
    relations = [{variable[e]: x for e, x in independent_chain(g, polygon).items()
                  if e in variable} for _, polygon in g.cells]
    active = set(range(len(g.chords)))
    for row, col in certificates[2]["eliminations"]:
        equation = relations[row]
        assert equation[col] in (-1, 1) and col in active
        for j, relation in enumerate(relations):
            if j == row or col not in relation:
                continue
            factor = sp.Rational(relation[col], equation[col])
            updated = {c: relation.get(c, 0)-factor*equation.get(c, 0)
                       for c in set(relation)|set(equation)}
            assert all(x.q == 1 for x in updated.values())
            relations[j] = {c: int(x) for c, x in updated.items() if x}
        relations[row] = {}
        active.remove(col)
    assert len(active) == 1
    divisor = 0
    for relation in relations:
        for value in relation.values():
            divisor = gcd(divisor, value)
    assert divisor == 2


def test_verifier_accepts_an_independently_chosen_valid_elimination_order(certificates):
    g = topology.geometry(2)
    variables = {e: i for i, e in enumerate(g.chords)}
    rows = [{variables[e]: x for e, x in independent_chain(g, p).items() if e in variables}
            for _, p in g.cells]
    remaining, steps = set(range(len(g.chords))), []
    while True:
        choices = [(len(r), i, -c) for i, r in enumerate(rows)
                   for c, x in r.items() if abs(x) == 1]
        if not choices:
            break
        _, i, minus_c = min(choices)
        c, equation = -minus_c, rows[i]
        steps.append([i, c])
        for j, row in enumerate(rows):
            if j == i or c not in row:
                continue
            factor = row[c]*equation[c]
            changed = {d: row.get(d, 0)-factor*equation.get(d, 0)
                       for d in set(row)|set(equation)}
            rows[j] = {d: x for d, x in changed.items() if x}
        rows[i] = {}
        remaining.remove(c)
    divisor = 0
    for row in rows:
        for x in row.values(): divisor = gcd(divisor, x)
    packet = deepcopy(certificates[2])
    packet["eliminations"] = steps
    packet["result"] = {"unit_pivots": len(steps), "remaining_generators": sorted(remaining),
        "nonzero_residual_rows": sum(bool(row) for row in rows),
        "residual_gcd": divisor, "free_rank": 0, "torsion": [divisor]}
    assert steps != certificates[2]["eliminations"]
    assert topology.verify_reduction(2, packet)["torsion"] == [2]


@pytest.mark.parametrize("n", range(2, 11))
def test_every_exchange_is_primitive_and_every_cell_is_flat(n, cocycles):
    g, values = topology.geometry(n), cocycles[n]
    modulus = 2*(n-1)
    positions = {(g.states[a], g.states[b]): i for i, (a, b) in enumerate(g.edges)}
    for _, polygon in g.cells:
        assert sum(values[e]*x for e, x in independent_chain(g, polygon, positions).items()) % modulus == 0
    loops = topology.exchange_loops(n)
    assert len(loops) == 20*comb(9, n-2)
    for face, spectators, path in loops:
        assert set(face).isdisjoint(spectators)
        assert path[0] == path[-1]
        chain = independent_chain(g, path[:-1], positions)
        assert sum(values[e]*x for e, x in chain.items()) % modulus == 1
        backwards = independent_chain(g, tuple(reversed(path[:-1])), positions)
        assert sum(values[e]*x for e, x in backwards.items()) % modulus == modulus-1


@pytest.mark.parametrize("n", range(2, 11))
def test_cyclic_generator_mod_two_equals_actual_ordered_token_transport(n, cocycles):
    """Compare native permutation parity up to its necessary vertex chart gauge."""
    g, values = topology.geometry(n), cocycles[n]
    parity = []
    for a, b in g.edges:
        left, right = g.states[a], g.states[b]
        origin = next(iter(set(left)-set(right)))
        target = next(iter(set(right)-set(left)))
        transported = tuple(target if p == origin else p for p in left)
        parity.append(sum(transported[i] > transported[j]
                          for i in range(n) for j in range(i+1, n)) % 2)
    adjacency = defaultdict(list)
    for e in g.tree:
        a, b = g.edges[e]
        adjacency[a].append((b, e)); adjacency[b].append((a, e))
    gauge, queue = {0: 0}, [0]
    for a in queue:
        for b, e in adjacency[a]:
            if b not in gauge:
                gauge[b] = gauge[a] ^ parity[e] ^ (values[e] % 2)
                queue.append(b)
    assert len(gauge) == len(g.states)
    assert all((gauge[a] ^ gauge[b] ^ parity[e]) == values[e] % 2
               for e, (a, b) in enumerate(g.edges))


def test_one_hole_twenty_face_loops_have_exactly_one_relation():
    g = topology.geometry(11)
    assert not g.cells
    rows = [tuple(independent_chain(g, p[:-1]).get(e, 0) for e in g.chords)
            for _, _, p in topology.exchange_loops(11)]
    matrix = sp.Matrix(rows)
    assert list(sp.ones(1, 20)*matrix) == [0]*19
    assert abs(matrix[:-1, :].det(method="domain-ge")) == 1
    assert topology.hole_face_certificate()["each_nineteen_face_basis_determinant_abs"] == 1
    assert topology.loop_class(11, topology.exchange_loops(11)[0][2]) == rows[0]


@pytest.mark.parametrize("n", [1, 2, 5, 10])
def test_particle_hole_complement_does_not_preserve_filled_triangles(n, cocycles):
    old, new = topology.geometry(n), topology.geometry(12-n)
    index = {state: i for i, state in enumerate(new.states)}
    triangle = next(p for k, p in old.cells if k == "particle_triangle")
    path = tuple(index[tuple(sorted(set(range(12))-set(old.states[i])))] for i in triangle)
    chain = independent_chain(new, path)
    assert chain
    assert not any(frozenset(path) == frozenset(poly) for _, poly in new.cells)
    if new.number == 11:
        assert any(chain.get(e, 0) for e in new.chords)
    else:
        assert sum(cocycles[new.number][e]*x for e, x in chain.items()) % (2*new.number-2) in (1, 2*new.number-3)


@pytest.mark.parametrize("corruption", ["number_bool", "number_float", "domain", "extra", "step_bool",
    "step_float", "step_tuple", "bad_row", "bad_col", "repeat", "missing_step", "order", "free_bool", "torsion_float", "nan"])
def test_replay_rejects_forged_or_malformed_certificates(corruption, certificates):
    packet = deepcopy(certificates[2])
    if corruption == "number_bool": packet["number"] = True
    elif corruption == "number_float": packet["number"] = 2.0
    elif corruption == "domain": packet["domain_sha256"] = "0"*64
    elif corruption == "extra": packet["discarded_cell"] = True
    elif corruption == "step_bool": packet["eliminations"][0][0] = True
    elif corruption == "step_float": packet["eliminations"][0][0] = float(packet["eliminations"][0][0])
    elif corruption == "step_tuple": packet["eliminations"][0] = tuple(packet["eliminations"][0])
    elif corruption == "bad_row": packet["eliminations"][0][0] = 10**6
    elif corruption == "bad_col": packet["eliminations"][0][1] = -1
    elif corruption == "repeat": packet["eliminations"][1] = packet["eliminations"][0][:]
    elif corruption == "missing_step": packet["eliminations"].pop()
    elif corruption == "order": packet["result"]["torsion"] = [4]
    elif corruption == "free_bool": packet["result"]["free_rank"] = False
    elif corruption == "torsion_float": packet["result"]["torsion"] = [2.0]
    elif corruption == "nan": packet["result"]["residual_gcd"] = float("nan")
    with pytest.raises((ValueError, ArithmeticError)):
        topology.verify_reduction(2, packet)


@pytest.mark.parametrize("number", [True, False, -1, 13, 2.0, "2", None])
def test_strict_particle_number(number):
    with pytest.raises(ValueError): topology.geometry(number)
    with pytest.raises(ValueError): topology.reduction_certificate(number)


def test_source_geometry_mutation_invalidates_warm_certificate(monkeypatch, certificates):
    seams, faces = native_loop_readout.carrier()
    monkeypatch.setattr(native_loop_readout, "carrier", lambda: (seams, faces[:-1]))
    with pytest.raises(ArithmeticError): topology.verify_reduction(2, certificates[2])


def test_source_orientation_mutation_is_not_silently_accepted(monkeypatch):
    seams, faces = native_loop_readout.carrier()
    changed = (tuple(reversed(faces[0])), *faces[1:])
    monkeypatch.setattr(native_loop_readout, "carrier", lambda: (seams, changed))
    with pytest.raises(ArithmeticError): topology.geometry(2)


def test_closed_path_validation():
    for bad in ((), (0, 1), (True, True), (0, 10000, 0), [0.0, 0.0]):
        with pytest.raises(ValueError): topology.loop_class(2, bad)
    g = topology.geometry(2)
    a, b = next((a, b) for a, b in combinations(range(len(g.states)), 2) if (a, b) not in g.edges)
    with pytest.raises(ArithmeticError): topology.loop_class(2, (a, b, a))
    assert topology.loop_class(0, (0, 0)) == 0
    with pytest.raises(ValueError): topology.edge_cocycle(11)


def test_compact_report_retains_homology_not_fundamental_group_scope():
    report = topology.certificate()
    assert len(report["sectors"]) == 13
    assert report["fundamental_group_computed_for_N_at_least3"] is False
    assert report["higher_rank_unitary_representations_classified"] is False
    assert report["physical_statistics_or_particle_identity_selected"] is False
