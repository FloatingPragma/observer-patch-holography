"""Two spectator-stable scalar flat families on the declared native complexes.

This adds an explicit compatibility requirement to the unrestricted H1
catalogue.  It does not derive that requirement from native repair, choose a
sign, classify nonabelian transport, or identify physical particles.

For every coherently oriented native triangular face, exchanging two occupied
vertices around the empty third vertex must have the same U(1) holonomy,
independent of fixed spectators outside that face.  The exact H1 calculation
makes these exchanges generators.  At N=2 their order is two, hence their
common phase is +1 or -1.  At N=11 the twenty face cycles have their sole
sum-zero relation, so the common phase also obeys z**20=1.  Both signs satisfy
every relation.  Connected graphs and tree gauge identify scalar flat
connections modulo vertex gauges with characters of integral H1; therefore
exactly two compatible families remain.

Existence is checked from actual ordered paths.  No creation algebra is an
input: ordering the occupied ports gives movement parity and, on the
subcomplex avoiding p, the insertion gauge (-1)**#{s<p}.  Adding a particle
does not define a map of the whole complex.  In particular a filling which
visits p is outside that subcomplex even if its boundary avoids p.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from math import comb, gcd

from . import many_body_topology as topology


SCHEMA = "oph.native_spectator_compatible_characters.v1"


def _character(value):
    if type(value) is not int or value not in (-1, 1):
        raise ValueError("an exact integer character +1 or -1 is required")
    return value


def _state(value):
    if (type(value) is not tuple or len(value) > 12
            or any(type(p) is not int or not 0 <= p < 12 for p in value)
            or tuple(sorted(set(value))) != value):
        raise ValueError("a sorted tuple of distinct exact native port integers is required")
    return value


def _port(value):
    if type(value) is not int or not 0 <= value < 12:
        raise ValueError("an exact native port integer is required")
    return value


def _turn(value):
    # Fraction's private slots can be changed with object.__setattr__.  Check
    # the exact canonical representation before invoking rational arithmetic.
    if type(value) is not F:
        raise ValueError("an exact canonical Fraction turn is required")
    numerator, denominator = value.numerator, value.denominator
    if (type(numerator) is not int or type(denominator) is not int
            or denominator <= 0 or gcd(numerator, denominator) != 1
            or not 0 <= numerator < denominator):
        raise ValueError("a normalized canonical Fraction turn is required")
    return value


def _ordered_parity(source, target):
    """Parity of the actual transported labels, assuming one legal movement."""
    if source == target:
        return 0
    a = next(iter(set(source)-set(target)))
    b = next(iter(set(target)-set(source)))
    lifted = tuple(b if p == a else p for p in source)
    return sum(x > y for i, x in enumerate(lifted) for y in lifted[i+1:]) % 2


def movement_sign(source, target, character=-1):
    """Sign from the ordered lift, with identity holds and no CAR assumption."""
    source, target = _state(source), _state(target)
    character = _character(character)
    if len(source) != len(target):
        raise ValueError("a number-preserving movement is required")
    if source != target:
        removed, added = set(source)-set(target), set(target)-set(source)
        if len(removed) != 1 or len(added) != 1:
            raise ValueError("one collision-free movement is required")
        if tuple(sorted((*removed, *added))) not in topology._source()[0]:
            raise ValueError("an actual native seam is required")
    return character**_ordered_parity(source, target)


def insertion_phase(state, port, character=-1):
    """Endpoint rechart for insertion at an empty port; not a unitary creator."""
    state, port, character = _state(state), _port(port), _character(character)
    if port in state:
        raise ValueError("the inserted spectator port must be empty")
    return character**sum(p < port for p in state)


def insert_path(number, path, port):
    """Embed a closed path in the p-empty subcomplex, with the spectator fixed."""
    number, port = topology._number(number), _port(port)
    if number == 12:
        raise ValueError("a full configuration has no insertion domain")
    source, target = topology.geometry(number), topology.geometry(number+1)
    topology._checked_path(source, path)
    if any(port in source.states[i] for i in path):
        raise ValueError("the entire path must avoid the inserted spectator port")
    index = {tuple(s): i for i, s in enumerate(target.states)}
    return tuple(index[tuple(sorted((*source.states[i], port)))] for i in path)


def _tree_gauge(g, parity):
    neighbors = [[] for _ in g.states]
    for e in g.tree:
        a, b = g.edges[e]
        neighbors[a].append((b, parity[e]))
        neighbors[b].append((a, parity[e]))
    values = {0: 0}
    queue = [0]
    for a in queue:
        for b, change in neighbors[a]:
            if b not in values:
                values[b] = values[a] ^ change
                queue.append(b)
    if len(values) != len(g.states):
        raise ArithmeticError("the declared tree does not reach every configuration")
    return values


def character_family(character):
    """Exact turn coordinates in each sector's checked H1 presentation.

    Trivial sectors have (), cyclic sectors one rational turn, and the
    one-hole sector nineteen turns in the native spanning-tree chord basis.
    """
    character = _character(character)
    g = topology.geometry(11)
    parity = tuple(_ordered_parity(g.states[a], g.states[b]) for a, b in g.edges)
    gauge = _tree_gauge(g, parity)
    hole = tuple(F((parity[e] ^ gauge[g.edges[e][0]] ^ gauge[g.edges[e][1]])
                   if character == -1 else 0, 2) for e in g.chords)
    return tuple(() if n in (0, 1, 12) else hole if n == 11
                 else (F(character == -1, 2),) for n in range(13))


def uniform_hole_character(turn):
    """The unique N11 character with one common oriented face phase.

    There are twenty such characters within N11 alone.  Their reduction to
    two is a cross-sector condition, not a deletion of the other eighteen.
    """
    turn = _turn(turn)
    if (20*turn).denominator != 1:
        raise ValueError("a normalized exact twentieth-root turn is required")
    g = topology.geometry(11)
    rows = tuple(tuple(topology._checked_path(g, path).get(e, 0) for e in g.chords)
                 for _, _, path in topology.exchange_loops(11))
    if len(rows) != 20 or any(sum(row[j] for row in rows) for j in range(19)):
        raise ArithmeticError("the twenty face cycles lost their sum-zero relation")
    if abs(topology._determinant(rows[:-1])) != 1:
        raise ArithmeticError("nineteen face cycles must be an integral basis")
    matrix = [[F(x) for x in row]+[turn] for row in rows[:-1]]
    for j in range(19):
        pivot = next(i for i in range(j, 19) if matrix[i][j])
        matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
        scale = matrix[j][j]
        matrix[j] = [x/scale for x in matrix[j]]
        for i in range(19):
            if i != j:
                scale = matrix[i][j]
                matrix[i] = [x-scale*y for x, y in zip(matrix[i], matrix[j])]
    coordinates = tuple(row[-1] % 1 for row in matrix)
    if any(sum((x*y for x, y in zip(row, coordinates)), F(0)) % 1 != turn for row in rows):
        raise ArithmeticError("the solved face phases violate the final closure relation")
    return coordinates


def validate_family(family):
    """Reject incomplete, nonflat, or spectator-dependent scalar characters.

    Fractions are normalized turns in [0,1); floats and integer/bool aliases
    are rejected.  This is an exact admission gate, not a phase tolerance.
    Completeness of these coordinates is certified by topology's H1 replay.
    """
    if type(family) is not tuple or len(family) != 13:
        raise ValueError("all thirteen occupation-sector characters are required")
    for n, coordinates in enumerate(family):
        count = 0 if n in (0, 1, 12) else 19 if n == 11 else 1
        if type(coordinates) is not tuple or len(coordinates) != count:
            raise ValueError("complete normalized exact rational turn coordinates are required")
        for turn in coordinates:
            _turn(turn)
        if 2 <= n <= 10 and (2*(n-1)*coordinates[0]).denominator != 1:
            raise ValueError("the phase violates its sector's integral torsion relation")
    common = family[2][0]
    if common not in (F(0), F(1, 2)):
        raise ArithmeticError("the two-particle order-two relation was lost")
    if any(family[n][0] != common for n in range(3, 11)):
        raise ValueError("local exchange holonomy depends on the added spectators")
    g = topology.geometry(11)
    for _, _, path in topology.exchange_loops(11):
        boundary = topology._checked_path(g, path)
        phase = sum((family[11][j]*boundary.get(e, 0)
                     for j, e in enumerate(g.chords)), F(0)) % 1
        if phase != common:
            raise ValueError("a one-hole face exchange differs from the common phase")
    if (20*common).denominator != 1:
        raise ArithmeticError("the twenty oriented hole-face relation was lost")
    return 1 if common == 0 else -1


def _cell_key(kind, polygon):
    # A geometric filled triangle/square is independent of start and orientation.
    rotations = tuple(polygon[i:]+polygon[:i] for i in range(len(polygon)))
    reverse = polygon[::-1]
    rotations += tuple(reverse[i:]+reverse[:i] for i in range(len(reverse)))
    return kind, min(rotations)


def spectator_certificate():
    """Every permitted vertex, edge and two-cell insertion, with exact parity."""
    geometries = tuple(topology.geometry(n) for n in range(13))
    counts = Counter()
    for n, source in enumerate(geometries[:-1]):
        target = geometries[n+1]
        target_index = {s: i for i, s in enumerate(target.states)}
        target_edges = set(target.edges)
        target_cells = {_cell_key(kind, polygon) for kind, polygon in target.cells}
        for p in range(12):
            inclusion = {i: target_index[tuple(sorted((*s, p)))]
                         for i, s in enumerate(source.states) if p not in s}
            counts["vertex_inclusions"] += len(inclusion)
            for a, b in source.edges:
                if a not in inclusion or b not in inclusion:
                    continue
                aa, bb = inclusion[a], inclusion[b]
                if tuple(sorted((aa, bb))) not in target_edges:
                    raise ArithmeticError("spectator insertion lost a native movement")
                before = _ordered_parity(source.states[a], source.states[b])
                after = _ordered_parity(target.states[aa], target.states[bb])
                change = (sum(s < p for s in source.states[a])
                          + sum(s < p for s in source.states[b])) % 2
                if before ^ change != after:
                    raise ArithmeticError("the ordered-lift spectator rechart failed")
                counts["edge_inclusions"] += 1
            for kind, polygon in source.cells:
                if all(i in inclusion for i in polygon):
                    image = tuple(inclusion[i] for i in polygon)
                    if _cell_key(kind, image) not in target_cells:
                        raise ArithmeticError("an avoiding-subcomplex cell failed to embed")
                    counts[kind+"_inclusions"] += 1
    expected = {"vertex_inclusions": 12*2**11,
                "edge_inclusions": 30*10*2**9,
                "particle_triangle_inclusions": 20*9*2**8,
                "disjoint_square_inclusions": 315*8*2**7}
    if dict(counts) != expected:
        raise ArithmeticError("the complete spectator-inclusion inventory changed")
    return {**dict(counts), "both_characters_intertwine_every_allowed_edge": True,
            "domain": "vertices, edges and fillings all avoid the inserted port",
            "whole_complex_creation_map_claimed": False}


def particle_hole_obstruction():
    """A filled N1 face becomes a nonzero unfilled one-hole exchange loop."""
    one, hole = topology.geometry(1), topology.geometry(11)
    one_index = {s: i for i, s in enumerate(one.states)}
    hole_index = {s: i for i, s in enumerate(hole.states)}
    face = one.faces[0]
    path = tuple(one_index[(p,)] for p in (*face, face[0]))
    image = tuple(hole_index[tuple(q for q in range(12) if q != p)]
                  for p in (*face, face[0]))
    if _cell_key("particle_triangle", path[:-1]) not in {
            _cell_key(kind, polygon) for kind, polygon in one.cells}:
        raise ArithmeticError("the original one-particle face is not filled")
    coordinates = topology.loop_class(11, image)
    if hole.cells or not any(coordinates):
        raise ArithmeticError("the one-hole noncell obstruction disappeared")
    return {"face": face, "one_particle_boundary": path,
            "complemented_one_hole_boundary": image,
            "one_hole_homology_coordinates": coordinates,
            "one_particle_face_is_filled": True,
            "one_hole_image_is_filled": False,
            "complement_is_a_two_complex_map": False}


def certificate():
    """Replay the integral classification and check both compatible realizations."""
    sectors = []
    total_cells, total_loops = 0, 0
    for n in range(13):
        g = topology.geometry(n)
        packet = topology.reduction_certificate(n)
        result = topology.verify_reduction(n, packet)
        parity = tuple(_ordered_parity(g.states[a], g.states[b]) for a, b in g.edges)
        positions = {edge: i for i, edge in enumerate(g.edges)}
        for _, polygon in g.cells:
            if sum(parity[e]*v for e, v in topology._boundary(polygon, positions).items()) % 2:
                raise ArithmeticError("the actual ordered lift is not flat on a declared cell")
        loops = topology.exchange_loops(n)
        cyclic = topology.edge_cocycle(n) if 2 <= n <= 10 else None
        for _, spectators, path in loops:
            boundary = topology._checked_path(g, path)
            if sum(parity[e]*v for e, v in boundary.items()) % 2 != 1:
                raise ArithmeticError("an actual local exchange does not swap ordered labels")
            if cyclic is not None and sum(cyclic[e]*v for e, v in boundary.items()) % (2*n-2) != 1:
                raise ArithmeticError("an exchange fails to generate its cyclic homology")
            if len(spectators) != n-2:
                raise ArithmeticError("a spectator was omitted from the exchange inventory")
        if len(loops) != (20*comb(9, n-2) if 2 <= n <= 11 else 0):
            raise ArithmeticError("an oriented exchange loop is missing")
        total_cells += len(g.cells)
        total_loops += len(loops)
        sectors.append({"number": n, "free_rank": result["free_rank"],
                        "torsion": result["torsion"], "local_exchange_loops": len(loops),
                        "negative_lift_flat_on_every_declared_cell": True})
    for character in (1, -1):
        family = character_family(character)
        if validate_family(family) != character:
            raise ArithmeticError("a compatible character family failed admission")
        if family[11] != uniform_hole_character(F(character == -1, 2)):
            raise ArithmeticError("the ordered lift differs from the unique common-phase hole character")
    return {"schema": SCHEMA, "sectors": sectors,
            "ordered_lift_flat_cells": total_cells,
            "ordered_lift_odd_exchange_loops": total_loops,
            "one_hole_face_relation": topology.hole_face_certificate(),
            "one_hole_common_phase_constraint": "z^20=1",
            "uniform_face_phase_characters_within_N11_alone": 20,
            "two_particle_common_phase_constraint": "z^2=1",
            "compatible_scalar_families_up_to_vertex_gauge": 2,
            "common_exchange_characters": [1, -1],
            "spectator_inclusions": spectator_certificate(),
            "particle_hole_obstruction": particle_hole_obstruction(),
            "compatibility_is_an_explicit_additional_requirement": True,
            "nonabelian_transport_classified": False,
            "creation_algebra_assumed_to_select_a_sign": False,
            "physical_statistics_or_particle_identity_selected": False}
