"""Literal finite source expressions with an explicit routing certificate."""
from fractions import Fraction as F
from functools import lru_cache
import importlib
from itertools import permutations
import json
import struct

import numpy as np

from model.federated_observer import source
from . import refinement as r

# Every source mesh edge key low*vertex_count+high fits signed int64 here.
# This is a local expression domain, not a claim of full graph allocation.
MAX_LOCAL_LEVEL = 14
SCORE_ERROR = F(1, 2**40)


def contraction_error_bound():
    """Two length-three double contractions, two legal accumulation orders."""
    unit = F(1, 2**53)
    gamma = 5*unit/(1-5*unit)
    result = 288*gamma+144*gamma**2+F(1, 2**900)
    if result >= SCORE_ERROR:
        raise ArithmeticError('the exact double-contraction allowance is insufficient')
    return result


def modules():
    source.verify_source_pins()
    source._source()
    return (importlib.import_module('oph_fpe.core.icosahedral'),
            importlib.import_module('oph_fpe.core.screen_ports'))


def _float32(bits):
    return struct.unpack('<f', struct.pack('<I', bits))[0]


def round32(value):
    """Exact nearest-even binary32 rounding of a bounded rational."""
    if type(value) is not F or abs(value) > 2:
        raise ValueError('an exact rational score in [-2,2] required')
    if value < 0:
        return -round32(-value)
    if value == 0:
        return F()
    bits = struct.unpack('<I', struct.pack('<f', float(value)))[0]
    for _ in range(3):
        point = F(_float32(bits))
        low = (F(_float32(bits-1))+point)/2 if bits else F()
        high = (point+F(_float32(bits+1)))/2
        if value < low or (value == low and bits % 2):
            bits -= 1
        elif value > high or (value == high and bits % 2):
            bits += 1
        else:
            return point
    raise ArithmeticError('binary32 rounding did not reach its exact nearest-even cell')


def certified_assignment(values, error=SCORE_ERROR):
    """Certify the literal first-argmax/collision-repair decision.

    The bounded double contraction variation is enclosed before binary32
    rounding. Collision repair requires a unique maximum-weight assignment,
    so its result cannot depend on the global row order or a solver tie.
    """
    if type(error) is not F or error < 0 or error > F(1, 1024):
        raise ValueError('a small exact nonnegative score allowance required')
    values = np.asarray(values)
    if values.shape != (3, 12) or values.dtype != np.float64 or not np.isfinite(values).all():
        raise ValueError('three complete finite binary64 score rows required')
    intervals = tuple(tuple((round32(F(float(x))-error),
                            round32(F(float(x))+error)) for x in row) for row in values)
    selected = tuple(int(i) for i in np.argmax(values.astype(np.float32), axis=1))
    for row, chosen in zip(intervals, selected, strict=True):
        for other, (low, high) in enumerate(row):
            if other == chosen:
                continue
            if (other < chosen and row[chosen][0] <= high) or (
                    other > chosen and row[chosen][0] < high):
                raise ValueError('routing score interval does not fix the source first argmax')
    collided = len(set(selected)) != 3
    if collided:
        candidates = tuple(permutations(range(12), 3))
        midpoint_scores = tuple(sum(F(float(np.float32(values[i, p])))
                                     for i, p in enumerate(choice)) for choice in candidates)
        selected = candidates[max(range(len(candidates)), key=midpoint_scores.__getitem__)]
        lower = sum(intervals[i][p][0] for i, p in enumerate(selected))
        alternatives = [sum(intervals[i][p][1] for i, p in enumerate(choice))
                        for choice in candidates if choice != selected]
        gap = lower-max(alternatives)
        if gap <= 0:
            raise ValueError('collision assignment lacks a unique certified source optimum')
    else:
        gap = min(intervals[i][p][0]-intervals[i][j][1]
                  for i, p in enumerate(selected) for j in range(12) if j != p)
    return {'ports': selected, 'collision_repaired': collided,
            'certified_score_gap': gap, 'score_variation_allowance': error}


class LocalSource(r.Topology):
    """Compute the source's local midpoint and routing expressions.

    Cell addresses and every normalization/coordinate-cleaning step are
    inherited. No source graph, clock or endpoint law is replaced.
    """
    def __init__(self):
        super().__init__()
        self.ico, self.ports_module = modules()
        base, self.base_faces = self.ico._base_icosahedron()
        self.base = self.ico._clean_coordinates(base)

    @lru_cache(None)
    def point(self, level, vertex):
        birth, parents = self.vertices[vertex]
        if level == 0:
            return self.base[vertex].copy()
        if birth < level:
            point = self.point(level-1, vertex).copy()
        else:
            a, b = parents
            point = self.point(level-1, a)+self.point(level-1, b)
            point /= np.linalg.norm(point)
        # _refine_level and _make_level each perform this cleaning.
        return self.ico._clean_coordinates(
            self.ico._clean_coordinates(np.asarray([point])))[0]

    @lru_cache(None)
    def cell(self, level, cell):
        r.integer(level, upper=MAX_LOCAL_LEVEL)
        r.integer(cell, upper=20*4**level-1)
        if level == 0:
            return tuple(map(int, self.base_faces[cell]))
        raw = self.raw_child(level, cell)
        points = np.asarray([self.point(level, v) for v in raw])
        face = self.ico._orient_faces_outward(points, np.asarray([[0, 1, 2]]))[0]
        return tuple(raw[int(i)] for i in face)

    @lru_cache(None)
    def center(self, level, cell):
        points = np.asarray([[self.point(level, v) for v in self.cell(level, cell)]])
        center = np.sum(points, axis=1)
        center /= np.linalg.norm(center, axis=1, keepdims=True)
        return center[0]

    def score_rows(self, level, cell):
        neighbours = tuple(self.neighbour(level, cell, side)[0] for side in range(3))
        points = np.asarray([self.center(level, cell),
                             *(self.center(level, n) for n in neighbours)])
        frames = self.ports_module._local_surface_frames(points)
        vectors = points[1:]-points[0]
        lengths = np.linalg.norm(vectors, axis=1, keepdims=True)
        if np.any(lengths <= 1e-15):
            raise ValueError('routing reaches the source normalization floor')
        vectors /= np.maximum(lengths, 1e-15)
        local = np.einsum('ni,nij->nj', vectors, frames[[0]*3], optimize=True)
        directions = self.ico.build_geodesic_icosahedral_tower(0).levels[0].vertices
        if any(not np.isfinite(array).all() or np.max(np.abs(array)) > 2
               for array in (vectors, frames, directions)):
            raise ValueError('source direction/frame component bound failed')
        contraction_error_bound()
        return neighbours, np.asarray(local@directions.T, dtype=np.float64)

    @lru_cache(None)
    def route(self, level, cell):
        neighbours, values = self.score_rows(level, cell)
        certificate = certified_assignment(values)
        return tuple(zip(neighbours, certificate['ports'], strict=True)), certificate


def routed_window(level, radius):
    r.integer(level, upper=MAX_LOCAL_LEVEL); r.integer(radius)
    if level < r.minimum_level(radius):
        raise ValueError('a proved corner collar is required')
    local = LocalSource()
    cells = local.ball(level, radius)
    expected = r.signature(level, radius)
    actual = tuple((r.address(level, c), tuple(r.address(level,
        local.neighbour(level, c, k)[0]) for k in range(3))) for c in cells)
    if actual != expected:
        raise ValueError('floating orientation differs from the exact corner topology')
    routing = {c: dict(local.route(level, c)[0]) for c in cells}
    names = tuple((r.address(level, c), p) for c in cells for p in range(12))
    name_to_index = {name: i for i, name in enumerate(names)}
    edges = []
    local_seams = source.transfer.native_counts().seams
    for c in cells:
        for a, b in local_seams:
            edges.append((name_to_index[(r.address(level, c), a)],
                          name_to_index[(r.address(level, c), b)]))
        for other, port in routing[c].items():
            if other in routing and c < other:
                a = name_to_index[(r.address(level, c), port)]
                b = name_to_index[(r.address(level, other), routing[other][c])]
                edges.append((a, b))
    source_ports = tuple(12*(face*4**level+offset)+port
                        for (face, offset), port in names)
    if tuple(sorted(source_ports)) != source_ports:
        raise ArithmeticError('the source address transport reversed endpoint orientation')
    return {'level': level, 'cell_radius': radius,
            'canonical_ports': names, 'source_ports': source_ports,
            'internal_seams': tuple(sorted(edges)),
            'root_ports': tuple(name_to_index[((0, 0), p)] for p in range(12)),
            'routing': tuple((r.address(level, c), tuple(sorted(
                (r.address(level, n), p) for n, p in routing[c].items()))) for c in cells),
            'routing_certificates': tuple(local.route(level, c)[1] for c in cells),
            'orientation_preserved': True,
            'global_constructor_or_native_history_executed': False}


def same_window(first, second):
    authenticate_window(first); authenticate_window(second)
    keys = ('cell_radius', 'canonical_ports', 'internal_seams', 'root_ports', 'routing')
    return all(first[k] == second[k] for k in keys)


def typed(value):
    if type(value) is F:
        return {'rational': str(value)}
    if type(value) is tuple:
        return {'tuple': [typed(v) for v in value]}
    if type(value) is dict and all(type(k) is str for k in value):
        return {k: typed(v) for k, v in value.items()}
    if value is None or type(value) in (str, int, bool):
        return value
    raise ValueError('an exact typed window required')


def authenticate_window(value):
    if type(value) is not dict or 'level' not in value or 'cell_radius' not in value:
        raise ValueError('a complete rooted window required')
    expected = routed_window(value['level'], value['cell_radius'])
    if json.dumps(typed(value), sort_keys=True) != json.dumps(typed(expected), sort_keys=True):
        raise ValueError('rooted source window differs from authenticated replay')
    return True
