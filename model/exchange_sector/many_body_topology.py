"""Integral first homology of every finite native binary occupation sector.

The declared cells are one particle crossing a triangular face with all other
particles outside it, and two moves on vertex-disjoint seams.  A spanning tree
identifies the integral cycle lattice with its non-tree edge coordinates.
Eliminating a relation with coefficient +1 or -1 is an invertible integer
presentation operation.  The replay below checks each such operation.  At
N=2,...,10 the remaining presentation is one generator modulo the gcd of its
relations.  Every local oriented triangular exchange is that primitive
generator.  Thus this is an exhaustive U(1) flat-character catalogue, not a
fundamental-group computation or a classification of higher-rank transports.

At N=11 there is just one vacancy: no declared two-cell is available and the
movement graph has free first homology of rank19.  Its twenty oriented hole
face loops generate with their single sum-zero relation.  Complementing all
occupations is a graph map but takes filled particle triangles to unfilled
exchange loops, so it is not a map of the declared two-complexes.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
from heapq import heappop, heappush
from itertools import combinations
from math import comb, gcd
import json

from ..observer import native_loop_readout


SCHEMA = "oph.native_many_body_homology.v1"


def _number(value):
    if type(value) is not int or not 0 <= value <= 12:
        raise ValueError("an exact particle number from0 through12 is required")
    return value


def _hash(value):
    return sha256(json.dumps(value, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


def _source():
    seams, faces = native_loop_readout.carrier()
    seams, faces = tuple(map(tuple, seams)), tuple(map(tuple, faces))
    if (len(seams) != 30 or len(set(seams)) != 30
            or any(len(e) != 2 or any(type(v) is not int or not 0 <= v < 12 for v in e)
                   or e[0] >= e[1] for e in seams)
            or len(faces) != 20 or len({tuple(sorted(f)) for f in faces}) != 20
            or any(len(f) != 3 or len(set(f)) != 3
                   or any(type(v) is not int or not 0 <= v < 12 for v in f) for f in faces)):
        raise ArithmeticError("the actual native seam/face inventory is malformed")
    incidence = defaultdict(list)
    for face in faces:
        for a, b in zip(face, face[1:]+face[:1]):
            edge = tuple(sorted((a, b)))
            if edge not in seams:
                raise ArithmeticError("a face has a nonnative edge")
            incidence[edge].append(1 if a < b else -1)
    if set(incidence) != set(seams) or any(sorted(v) != [-1, 1] for v in incidence.values()):
        raise ArithmeticError("the native faces must coherently cover every seam twice")
    return seams, faces


@dataclass(frozen=True)
class Geometry:
    number: int
    states: tuple
    edges: tuple
    cells: tuple
    tree: tuple
    chords: tuple
    seams: tuple
    faces: tuple


def geometry(number):
    """The actual finite complex; state addresses are lexicographic subsets."""
    n = _number(number)
    seams, faces = _source()
    states = tuple(combinations(range(12), n))
    index = {frozenset(s): i for i, s in enumerate(states)}
    edges = set()
    for i, state in enumerate(states):
        occupied = set(state)
        for a, b in seams:
            if (a in occupied) != (b in occupied):
                edges.add(tuple(sorted((i, index[frozenset(occupied ^ {a, b})]))))
    edges = tuple(sorted(edges))
    cells = []

    def at(s):
        return index[frozenset(s)]

    if 1 <= n <= 10:
        for face in faces:
            for spectators in combinations(sorted(set(range(12))-set(face)), n-1):
                cells.append(("particle_triangle", tuple(at((*spectators, p)) for p in face)))
    if 2 <= n <= 10:
        for (a, b), (c, d) in combinations(seams, 2):
            if len({a, b, c, d}) != 4:
                continue
            for s in combinations(sorted(set(range(12))-{a, b, c, d}), n-2):
                cells.append(("disjoint_square", tuple(at(p) for p in
                    ((*s, a, c), (*s, b, c), (*s, b, d), (*s, a, d)))))
    parents = list(range(len(states)))

    def find(a):
        while parents[a] != a:
            parents[a] = parents[parents[a]]
            a = parents[a]
        return a

    tree, chords = [], []
    for i, (a, b) in enumerate(edges):
        aa, bb = find(a), find(b)
        if aa == bb:
            chords.append(i)
        else:
            parents[aa] = bb
            tree.append(i)
    if len(tree) != len(states)-1:
        raise ArithmeticError("the binary movement graph is disconnected")
    return Geometry(n, states, edges, tuple(cells), tuple(tree), tuple(chords), seams, faces)


def _boundary(path, positions):
    result = {}
    for a, b in zip(path, path[1:]+path[:1]):
        if a == b:
            continue
        edge = positions.get(tuple(sorted((a, b))))
        if edge is None:
            raise ArithmeticError("a cell/path uses a nonnative movement")
        result[edge] = result.get(edge, 0)+(1 if a < b else -1)
    return {i: x for i, x in result.items() if x}


def _cycle_rows(g):
    positions = {edge: i for i, edge in enumerate(g.edges)}
    variables = {edge: c for c, edge in enumerate(g.chords)}
    rows = []
    for _, polygon in g.cells:
        boundary = _boundary(polygon, positions)
        # Check d1*d2=0 on the full edge chain, before dropping tree edges.
        vertices = defaultdict(int)
        for e, coefficient in boundary.items():
            a, b = g.edges[e]
            vertices[a] -= coefficient
            vertices[b] += coefficient
        if any(vertices.values()):
            raise ArithmeticError("an integral cell boundary is not a cycle")
        rows.append({variables[e]: x for e, x in boundary.items() if e in variables})
    return rows


def _domain_hash(g):
    return _hash([g.number, g.states, g.edges, g.cells, g.tree, g.chords])


def _produce_steps(original):
    """Greedy sparse unit pivots; this search is not the verifier."""
    rows = {i: row.copy() for i, row in enumerate(original) if row}
    columns, version, queue = defaultdict(set), defaultdict(int), []
    for i, row in rows.items():
        for c in row:
            columns[c].add(i)
        heappush(queue, (len(row), i, 0))
    steps = []
    while queue:
        _, i, v = heappop(queue)
        if i not in rows or version[i] != v:
            continue
        row = rows[i]
        units = [c for c, x in row.items() if abs(x) == 1]
        if not units:
            continue
        pivot = min(units, key=lambda c: (len(columns[c]), c))
        affected = tuple(sorted(columns[pivot]-{i}))
        steps.append([i, pivot])
        del rows[i]
        for c in row:
            columns[c].remove(i)
        for j in affected:
            other = rows[j]
            multiplier = other[pivot]*row[pivot]
            for c, x in row.items():
                old = other.get(c, 0)
                new = old-multiplier*x
                if new:
                    other[c] = new
                    if not old:
                        columns[c].add(j)
                elif old:
                    del other[c]
                    columns[c].remove(j)
            version[j] += 1
            if other:
                heappush(queue, (len(other), j, version[j]))
            else:
                del rows[j]
    return steps


def _replay(original, column_count, steps):
    """Check integer presentation operations independently of pivot search."""
    if type(steps) is not list or len(steps) > column_count:
        raise ValueError("a bounded list of unit-pivot operations is required")
    rows = [dict(row) for row in original]
    occurrences = defaultdict(set)
    for i, row in enumerate(rows):
        for c in row:
            occurrences[c].add(i)
    remaining = set(range(column_count))
    replacements = []
    for step in steps:
        if (type(step) is not list or len(step) != 2
                or any(type(x) is not int for x in step)):
            raise ValueError("each pivot is a pair of exact integer row/column addresses")
        i, c = step
        if not 0 <= i < len(rows) or c not in remaining or abs(rows[i].get(c, 0)) != 1:
            raise ValueError("pivot must be an available relation with coefficient plus/minus one")
        relation = rows[i]
        sign = relation[c]
        replacement = {d: -x*sign for d, x in relation.items() if d != c}
        replacements.append((c, replacement))
        for j in tuple(occurrences[c]-{i}):
            old = rows[j]
            coefficient = old[c]
            new = dict(old)
            del new[c]
            for d, x in replacement.items():
                new[d] = new.get(d, 0)+coefficient*x
                if not new[d]:
                    del new[d]
            for d in set(old)-set(new):
                occurrences[d].remove(j)
            for d in set(new)-set(old):
                occurrences[d].add(j)
            rows[j] = new
        for d in relation:
            occurrences[d].remove(i)
        rows[i] = {}
        remaining.remove(c)
    active = set().union(*(set(row) for row in rows)) if rows else set()
    if any(c not in remaining for c in active):
        raise ArithmeticError("a removed generator survived replay")
    if len(active) > 1:
        raise ValueError("the supplied trace did not reach the declared one-column Smith form")
    divisor = 0
    for row in rows:
        for x in row.values():
            divisor = gcd(divisor, x)
    free_rank = len(remaining)-len(active)
    torsion = [divisor] if divisor > 1 else []
    summary = {"unit_pivots": len(steps), "remaining_generators": sorted(remaining),
               "nonzero_residual_rows": sum(bool(row) for row in rows),
               "residual_gcd": divisor, "free_rank": free_rank, "torsion": torsion}
    return summary, replacements, active


def reduction_certificate(number):
    """Return a replayable sparse integer trace; the retained report uses its hash."""
    g = geometry(number)
    original = _cycle_rows(g)
    steps = _produce_steps(original)
    summary, _, _ = _replay(original, len(g.chords), steps)
    return {"schema": SCHEMA, "number": number, "domain_sha256": _domain_hash(g),
            "eliminations": steps, "result": summary}


def verify_reduction(number, packet):
    """Verify every operation against fresh actual cells, with strict JSON types."""
    g = geometry(number)
    if (type(packet) is not dict
            or set(packet) != {"schema", "number", "domain_sha256", "eliminations", "result"}
            or packet["schema"] != SCHEMA or type(packet["number"]) is not int
            or packet["number"] != number or packet["domain_sha256"] != _domain_hash(g)):
        raise ValueError("certificate does not describe the complete declared native complex")
    result, _, _ = _replay(_cycle_rows(g), len(g.chords), packet["eliminations"])
    supplied = packet["result"]
    if (type(supplied) is not dict or set(supplied) != set(result)
            or any(type(supplied[k]) is not int for k in
                   ("unit_pivots", "nonzero_residual_rows", "residual_gcd", "free_rank"))
            or any(type(supplied[k]) is not list
                   or any(type(x) is not int for x in supplied[k])
                   for k in ("remaining_generators", "torsion"))):
        raise ValueError("the Smith summary requires exact integer scalars and lists")
    # Type-sensitive JSON comparison rejects bool/int and int/float aliases.
    if json.dumps(supplied, sort_keys=True, allow_nan=False) != json.dumps(result, sort_keys=True):
        raise ValueError("retained Smith summary differs from exact integer replay")
    expected = ([2*(number-1)] if 2 <= number <= 10 else [])
    if result["torsion"] != expected or result["free_rank"] != (19 if number == 11 else 0):
        raise ArithmeticError("the actual native integral homology changed")
    return result


def exchange_loops(number):
    """All oriented triangular exchanges, with fixed spectators outside the face."""
    g = geometry(number)
    if not 2 <= number <= 11:
        return ()
    index = {frozenset(s): i for i, s in enumerate(g.states)}
    result = []
    for a, b, c in g.faces:
        for spectators in combinations(sorted(set(range(12))-{a, b, c}), number-2):
            path = tuple(index[frozenset((*spectators, *pair))]
                         for pair in ((a, b), (b, c), (a, c), (a, b)))
            result.append(((a, b, c), spectators, path))
    return tuple(result)


def _checked_path(g, path, positions=None):
    if (type(path) not in (list, tuple) or not path or path[0] != path[-1]
            or any(type(x) is not int or not 0 <= x < len(g.states) for x in path)):
        raise ValueError("a closed path of exact native configuration addresses is required")
    if positions is None:
        positions = {edge: i for i, edge in enumerate(g.edges)}
    return _boundary(tuple(path[:-1]), positions)


def edge_cocycle(number):
    """Primitive cyclic H1 coordinate, normalized to +1 on every face exchange."""
    _number(number)
    if not 2 <= number <= 10:
        raise ValueError("a cyclic occupation sector from2 through10 is required")
    g = geometry(number)
    original = _cycle_rows(g)
    steps = _produce_steps(original)
    result, replacements, active = _replay(original, len(g.chords), steps)
    modulus = 2*(number-1)
    if result["torsion"] != [modulus] or result["free_rank"] or len(active) != 1:
        raise ArithmeticError("the cyclic integer presentation changed")
    values = {next(iter(active)): 1}
    for c, replacement in reversed(replacements):
        values[c] = sum(x*values[d] for d, x in replacement.items()) % modulus
    edges = [0]*len(g.edges)
    for c, e in enumerate(g.chords):
        edges[e] = values[c]
    loops = exchange_loops(number)
    positions = {edge: i for i, edge in enumerate(g.edges)}
    first = sum(edges[e]*x for e, x in _checked_path(g, loops[0][2], positions).items()) % modulus
    if gcd(first, modulus) != 1:
        raise ArithmeticError("the local exchange is not a primitive homology generator")
    scale = pow(first, -1, modulus)
    edges = tuple((x*scale) % modulus for x in edges)
    if any(sum(edges[e]*x for e, x in _boundary(poly, positions).items()) % modulus
           for _, poly in g.cells):
        raise ArithmeticError("the integral character coordinate is not flat on an admitted cell")
    if any(sum(edges[e]*x for e, x in _checked_path(g, path, positions).items()) % modulus != 1
           for _, _, path in loops):
        raise ArithmeticError("the oriented local exchanges do not share one primitive class")
    return edges


def loop_class(number, path):
    """Cyclic residue, rank19 integer coordinates, or the zero trivial class."""
    g = geometry(number)
    boundary = _checked_path(g, path)
    if number == 11:
        return tuple(boundary.get(e, 0) for e in g.chords)
    if 2 <= number <= 10:
        values = edge_cocycle(number)
        return sum(values[e]*x for e, x in boundary.items()) % (2*(number-1))
    return 0


def _determinant(matrix):
    """Fraction-free determinant, used only on the nineteen hole-face cycles."""
    a = [list(row) for row in matrix]
    n, sign, previous = len(a), 1, 1
    for k in range(n-1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        denominator = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = a[i][j]*denominator-a[i][k]*a[k][j]
                quotient, residue = divmod(numerator, previous)
                if residue:
                    raise ArithmeticError("Bareiss elimination lost exact divisibility")
                a[i][j] = quotient
            a[i][k] = 0
        previous = denominator
    return sign*a[-1][-1] if n else 1


def hole_face_certificate():
    """Twenty face loops freely generate H1 subject to their one sum relation."""
    g = geometry(11)
    vectors = tuple(tuple(_checked_path(g, path).get(e, 0) for e in g.chords)
                    for _, _, path in exchange_loops(11))
    if len(vectors) != 20 or any(sum(v[j] for v in vectors) for j in range(19)):
        raise ArithmeticError("the oriented twenty-face relation changed")
    determinants = tuple(_determinant(vectors[:i]+vectors[i+1:]) for i in range(20))
    if any(abs(d) != 1 for d in determinants):
        raise ArithmeticError("nineteen hole-face cycles fail to be an integral basis")
    return {"face_cycles": 20, "free_rank": 19,
            "single_relation": "sum of the twenty coherently oriented hole-face loops is zero in H1",
            "each_nineteen_face_basis_determinant_abs": 1,
            "cycle_coordinates_sha256": _hash(vectors)}


def certificate():
    """Compact deterministic report; full pivot traces are generated on demand."""
    sectors = []
    for number in range(13):
        g = geometry(number)
        packet = reduction_certificate(number)
        result = verify_reduction(number, packet)
        triangles = sum(kind == "particle_triangle" for kind, _ in g.cells)
        row = {"number": number, "vertices": len(g.states), "movement_edges": len(g.edges),
               "particle_triangles": triangles, "disjoint_squares": len(g.cells)-triangles,
               "tree_edges": len(g.tree), "cycle_rank": len(g.chords),
               "free_rank": result["free_rank"], "torsion": result["torsion"],
               "unit_pivots": result["unit_pivots"],
               "residual_gcd": result["residual_gcd"],
               "domain_sha256": packet["domain_sha256"],
               "integer_replay_sha256": _hash(packet)}
        if 2 <= number <= 10:
            values = edge_cocycle(number)
            row.update({"local_exchange_loops": len(exchange_loops(number)),
                        "every_oriented_exchange_class": 1,
                        "exchange_generator_order": 2*(number-1),
                        "edge_cocycle_sha256": _hash(values)})
        sectors.append(row)
    return {"schema": SCHEMA, "sectors": sectors, "one_hole_faces": hole_face_certificate(),
            "coefficient_domain": "integers; exact unit pivots and final one-column gcd",
            "flat_scalar_characters": "N2..10: exp(2*pi*i*k/(2*N-2)), k=0..2*N-3; N11: U(1)^19; N0,1,12: trivial",
            "deformation_relations": "one-particle triangles and disjoint-edge squares only",
            "particle_hole_complement": "movement-graph isomorphism, not a map of the declared two-cells",
            "fundamental_group_computed_for_N_at_least3": False,
            "higher_rank_unitary_representations_classified": False,
            "physical_statistics_or_particle_identity_selected": False}
