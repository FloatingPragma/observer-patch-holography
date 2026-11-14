"""Local face attachment for the actual cubic-routed carrier graph.

Faces are declared geometric data. The carrier graph and its transport law
are unchanged. No magnetic dynamics or continuum field is inferred here.
"""
from collections import deque
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sys

import sympy as s
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'reverse-engineering-reality'))
from model.native_spatial_refinement import source as geometry


def build(size=2):
    g = geometry.graph(size)
    edges = tuple(g.seams)
    edge_id = {e: i for i, e in enumerate(edges)}
    cells = {c: i for i, c in enumerate(g.cell_coordinates)}
    ico = tuple(geometry.slow_geometry.source().seams())
    adjacency = [[] for _ in range(12)]
    for a, b in ico:
        adjacency[a].append(b)
        adjacency[b].append(a)
    previous = {0: None}
    todo = deque([0])
    while todo:
        x = todo.popleft()
        for y in sorted(adjacency[x]):
            if y not in previous:
                previous[y] = x
                todo.append(y)
    paths = {}
    for v in range(12):
        p = [v]
        while p[-1]:
            p.append(previous[p[-1]])
        paths[v] = tuple(reversed(p))

    def chain(vertices):
        out = s.zeros(len(edges), 1)
        for a, b in zip(vertices, vertices[1:]):
            out[edge_id[tuple(sorted((a, b)))]] += 1 if a < b else -1
        return out

    def shift(c, axis):
        p = list(c)
        p[axis] = (p[axis] + 1) % size
        return tuple(p)

    B = s.zeros(g.ports, len(edges))
    for e, (a, b) in enumerate(edges):
        B[a, e], B[b, e] = 1, -1
    face_columns, labels = [], []
    ico_set = set(ico)
    triangles = [t for t in combinations(range(12), 3)
                 if all(e in ico_set for e in combinations(t, 2))]
    assert len(triangles) == 20
    for cell in range(g.cells):
        for a, b, c in triangles:
            face_columns.append(chain([12 * cell + x for x in (a, b, c, a)]))
            labels.append(('internal_triangle', cell, a, b, c))
    macro = {}
    periods = s.zeros(3, len(edges))
    for c, cell in cells.items():
        for axis, row in enumerate(geometry.links()):
            nxt = shift(c, axis)
            other = cells[nxt]
            for channel, (a, b) in enumerate(row):
                vertices = ([12 * cell + x for x in paths[a]]
                            + [12 * other + x for x in reversed(paths[b])])
                macro[c, axis, channel] = chain(vertices)
                left, right = 12 * cell + a, 12 * other + b
                e = edge_id[tuple(sorted((left, right)))]
                periods[axis, e] = s.Rational(1 if left < right else -1, size)
            face_columns.append(macro[c, axis, 1] - macro[c, axis, 0])
            labels.append(('intercell_ribbon', cell, axis))
    for c, cell in cells.items():
        for a, b in combinations(range(3), 2):
            face_columns.append(macro[c, a, 0] + macro[shift(c, a), b, 0]
                                - macro[shift(c, b), a, 0] - macro[c, b, 0])
            labels.append(('macro_plaquette', cell, a, b))
    A = s.Matrix.hstack(*face_columns)
    assert B * A == s.zeros(g.ports, len(face_columns))
    assert periods * A == s.zeros(3, len(face_columns))
    windings = []
    for axis in range(3):
        c = (0, 0, 0)
        winding = s.zeros(len(edges), 1)
        for _ in range(size):
            winding += macro[c, axis, 0]
            c = shift(c, axis)
        windings.append(winding)
    Z = s.Matrix.hstack(*windings)
    assert B * Z == s.zeros(g.ports, 3) and periods * Z == s.eye(3)
    return g, B, A, periods, Z, labels


def run():
    source_paths = [Path(__file__), Path(geometry.__file__)]
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in source_paths}
    g, B, A, periods, winding, labels = build(2)
    # Cotree coordinates identify the integral graph cycle lattice with Z^b1.
    parent = list(range(g.ports))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    tree, cotree = [], []
    for e, (a, b) in enumerate(g.seams):
        a, b = find(a), find(b)
        if a != b:
            parent[a] = b
            tree.append(e)
        else:
            cotree.append(e)
    assert len(tree) == g.ports - 1 and len(cotree) == 193
    cycle_faces = A[cotree, :]
    normal = smith_normal_form(cycle_faces, domain=s.ZZ)
    factors = [abs(int(normal[i, i])) for i in range(min(normal.shape)) if normal[i, i]]
    assert factors == [1] * 190
    # Three explicit unit periods split the free quotient; there is no torsion.
    assert winding[cotree, :].row_join(cycle_faces).rank() == len(cotree)
    max_cells = 0
    max_edges = 0
    for col in range(A.cols):
        support = [e for e in range(A.rows) if A[e, col]]
        cells = {x // 12 for e in support for x in g.seams[e]}
        max_cells = max(max_cells, len(cells))
        max_edges = max(max_edges, len(support))
        assert len(cells) <= 4
    # Internal triangles alone omit all41 intercell cycle classes at size2.
    internal = A[:, [i for i, label in enumerate(labels) if label[0] == 'internal_triangle']]
    assert internal.rank() == 19 * g.cells
    result = {'schema': 'oph.local12.local_cell_complex.v1',
              'candidate_id': 'native-local12-current-v0-proposed',
              'checks': {'actual_graph_unchanged': True, 'all_faces_have_zero_divergence': True,
                         'all_faces_have_zero_period': True, 'three_unit_period_cycles': True,
                         'integral_cycle_quotient_is_free_rank_three': True,
                         'bounded_local_face_support': True, 'omitted_glue_faces_control': True},
              'vertices': g.ports, 'edges': len(g.seams), 'faces': A.cols,
              'face_counts': {kind: sum(label[0] == kind for label in labels)
                             for kind in ('internal_triangle', 'intercell_ribbon', 'macro_plaquette')},
              'graph_cycle_rank': len(cotree), 'face_boundary_rank': len(factors),
              'cycle_quotient_Smith_nonzero_factors': factors,
              'harmonic_cycle_rank': 3, 'internal_only_remaining_cycle_rank': len(cotree) - internal.rank(),
              'maximum_cells_per_face': max_cells, 'maximum_edges_per_face': max_edges,
              'orientation': 'B has +1 at lower-index source and -1 at target; n+B F is conserved by a positive source-to-target move.',
              'scope': 'Declared local two-cells on the existing carrier graph. The complete integral homology check is size2 development. No field Hamiltonian, magnetic update, radiation or physical continuum is claimed.',
              'source_sha256': pins, 'official_rows_changed': []}
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in source_paths}
    (HERE / 'CELL_COMPLEX.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'vertices', 'edges', 'faces',
                      'face_counts', 'graph_cycle_rank', 'face_boundary_rank',
                      'harmonic_cycle_rank', 'internal_only_remaining_cycle_rank',
                      'maximum_cells_per_face', 'maximum_edges_per_face')}))


if __name__ == '__main__':
    run()
