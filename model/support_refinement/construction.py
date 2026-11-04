"""Rational subdivision retracts of the actual finite spherical support."""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
from types import MappingProxyType

ROOT = Path(__file__).resolve().parents[2]
SOURCE_SHA256 = MappingProxyType({
    'code/source_selection_model/geometry.py': 'de47ae19dd4ae1677d0cefb0163513f2f1f6d24dfb4dd0264506d0b17fceb843',
    'Lean/ObserverPatchHolography/CoreAxioms.lean': 'addcea08fe55c0ccc3edbeb90ebd9df4b640c96bd20e9755417bcdb84f093fe4',
    'code/a5_closure/verify_port_gram_support_attachment_independent.py': 'c2432831d20029e8eff4a2b015733075a751e263fd9fbad19f2b287f04aa73ac',
    'code/m1_source_realization/receipt.json': 'cef14e5c14bdf02e14f21d6f470bcfc6c6f396dd301e8067ca6454194f25ab49',
})
GEOMETRY, CORE, SEARCH, RECEIPT = SOURCE_SHA256


def need(condition, message):
    if not condition:
        raise ValueError(message)


def basis(vertices):
    """An oriented simplex; repeated vertices represent zero in an exterior power."""
    need(type(vertices) is tuple and 1 <= len(vertices) <= 3 and
         all(type(v) is int and v >= 0 for v in vertices), 'exact simplex vertices')
    if len(set(vertices)) != len(vertices):
        return {}
    sign = (-1) ** sum(a > b for i, a in enumerate(vertices) for b in vertices[i+1:])
    return {tuple(sorted(vertices)): Q(sign)}


def combine(*terms):
    """Add rational sparse chains, each supplied as (coefficient, chain)."""
    result = defaultdict(Q)
    for scale, chain in terms:
        need(type(scale) in (int, Q), 'exact rational coefficient')
        for simplex, value in chain.items():
            need(type(value) in (int, Q), 'exact rational coefficient')
            result[simplex] += scale * value
    return {simplex: value for simplex, value in result.items() if value}


def boundary(chain):
    return combine(*((value * (-1)**i, basis(simplex[:i]+simplex[i+1:]))
                     for simplex, value in chain.items() if len(simplex) > 1
                     for i in range(len(simplex))))


def push(chain, permutation):
    return combine(*((value, basis(tuple(permutation[v] for v in simplex)))
                     for simplex, value in chain.items()))


class Level:
    """A checked oriented level of the four-stage support, with canonical bases."""
    def __init__(self, row):
        need(type(row) is dict and set(row) == {'level', 'vertices', 'faces', 'coarsen'},
             'support row schema')
        n, level, faces = row['vertices'], row['level'], row['faces']
        need(type(level) is int and 0 <= level <= 3 and type(n) is int and
             n == 10 * 4**level + 2, 'retained support size')
        need(type(faces) is list and len(faces) == 20 * 4**level, 'retained face census')
        self.index, self.vertices = level, n
        self.faces = tuple(tuple(f) for f in faces)
        face_chain, edges = {}, Counter()
        for face in self.faces:
            need(len(face) == 3 and len(set(face)) == 3 and
                 all(type(v) is int and 0 <= v < n for v in face), 'support face')
            term = basis(face)
            key = next(iter(term))
            need(key not in face_chain, 'duplicate support face')
            face_chain.update(term)
            a, b, c = face
            edges.update(((a, b), (b, c), (c, a)))
        need(all(count == 1 and edges[b, a] == 1 for (a, b), count in edges.items()),
             'oriented closed surface')
        self.simplices = tuple(sorted({s for face in self.faces for k in (1, 2, 3)
                                      for s in combinations(sorted(face), k)}))
        self.members = frozenset(self.simplices)
        need({v for s in self.simplices for v in s} == set(range(n)) and
             len(edges) == 60 * 4**level, 'complete support incidence')
        self.fundamental = MappingProxyType(face_chain)
        parent = row['coarsen']
        if level == 0:
            need(parent is None, 'seed has no coarsening')
        else:
            need(type(parent) is list and len(parent) == n and
                 all(type(v) is int and 0 <= v < 10*4**(level-1)+2 for v in parent),
                 'endpoint coarsening type')
        self.endpoint = None if parent is None else tuple(parent)

    def check(self, chain):
        need(isinstance(chain, dict) or isinstance(chain, MappingProxyType), 'sparse chain mapping')
        for simplex, coefficient in chain.items():
            need(type(simplex) is tuple and simplex in self.members and
                 all(type(v) is int for v in simplex), 'chain outside the declared complex')
            need(type(coefficient) in (int, Q), 'exact rational coefficient')


class Refinement:
    """Local maps C, S, H with dH+Hd=id-SC and CS=id."""
    def __init__(self, coarse, fine):
        need(isinstance(coarse, Level) and isinstance(fine, Level) and
             fine.index == coarse.index+1, 'consecutive support levels')
        self.coarse, self.fine = coarse, fine
        n = coarse.vertices
        adjacency = defaultdict(set)
        for simplex in fine.simplices:
            if len(simplex) == 2:
                a, b = simplex
                adjacency[a].add(b)
                adjacency[b].add(a)
        endpoints = {v: tuple(sorted(x for x in adjacency[v] if x < n))
                     for v in range(n, fine.vertices)}
        old_edges = {s for s in coarse.simplices if len(s) == 2}
        need(all(len(edge) == 2 and edge in old_edges for edge in endpoints.values()) and
             len(set(endpoints.values())) == len(endpoints) and
             set(endpoints.values()) == old_edges, 'unique complete midpoint incidence')
        self.midpoint = {edge: v for v, edge in endpoints.items()}
        self.endpoints = endpoints
        need(fine.endpoint[:n] == tuple(range(n)) and
             all(fine.endpoint[v] in edge for v, edge in endpoints.items()),
             'endpoint coarsening changes ancestry')
        expected = combine(*((1, self._subdivide(face)) for face in coarse.faces))
        need(expected == fine.fundamental, 'oriented midpoint refinement')
        self._C = {s: self._coarsen(s) for s in fine.simplices}
        self._S = {s: self._subdivide(s) for s in coarse.simplices}
        self._H = {s: self._homotopy(s) for s in fine.simplices}
        for mapping, target in ((self._C, coarse), (self._S, fine), (self._H, fine)):
            for image in mapping.values():
                target.check(image)

    def _subdivide(self, simplex):
        if len(simplex) == 1:
            return basis(simplex)
        a, b = simplex[:2]
        ab = self.midpoint[tuple(sorted((a, b)))]
        if len(simplex) == 2:
            return combine((1, basis((a, ab))), (1, basis((ab, b))))
        c = simplex[2]
        bc, ca = (self.midpoint[tuple(sorted(edge))] for edge in ((b, c), (c, a)))
        return combine(*((1, basis(face)) for face in
                         ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))))

    def _coarsen(self, simplex):
        choices = [(v,) if v < self.coarse.vertices else self.endpoints[v] for v in simplex]
        coefficient = Q(1, 2**sum(len(xs) == 2 for xs in choices))
        return combine(*((coefficient, basis(term)) for term in product(*choices)))

    def _homotopy(self, simplex):
        n = self.coarse.vertices
        if len(simplex) == 1:
            v = simplex[0]
            return {} if v < n else combine(*((Q(1, 2), basis((a, v)))
                                               for a in self.endpoints[v]))
        if len(simplex) == 3 or any(v < n for v in simplex):
            return {}
        u, v = simplex
        shared = set(self.endpoints[u]) & set(self.endpoints[v])
        need(len(shared) == 1, 'interior fine edge has one corner')
        a = shared.pop()
        b = next(x for x in self.endpoints[u] if x != a)
        c = next(x for x in self.endpoints[v] if x != a)
        return combine((1, basis((a, u, v))), (Q(-1, 4), self._subdivide((a, b, c))))

    @staticmethod
    def _apply(chain, mapping):
        return combine(*((value, mapping[s]) for s, value in chain.items()))

    def coarsen(self, chain):
        self.fine.check(chain)
        return self._apply(chain, self._C)

    def subdivide(self, chain):
        self.coarse.check(chain)
        return self._apply(chain, self._S)

    def homotopy(self, chain):
        self.fine.check(chain)
        return self._apply(chain, self._H)

    def lift(self, permutation):
        need(type(permutation) is tuple and all(type(v) is int for v in permutation) and
             sorted(permutation) == list(range(self.coarse.vertices)) and
             push(self.coarse.fundamental, permutation) == self.coarse.fundamental,
             'proper coarse support action')
        result = permutation + tuple(self.midpoint[tuple(sorted(permutation[a] for a in self.endpoints[v]))]
                                     for v in range(self.coarse.vertices, self.fine.vertices))
        need(push(self.fine.fundamental, result) == self.fine.fundamental,
             'proper lifted support action')
        return result


class Tower:
    def __init__(self, rows, seed_actions):
        need(type(rows) is list and len(rows) == 4, 'complete retained tower')
        self.levels = tuple(Level(row) for row in rows)
        self.steps = tuple(Refinement(a, b) for a, b in zip(self.levels, self.levels[1:]))
        actions = tuple(tuple(p) for p in seed_actions)
        need(len(actions) == 60 and len(set(actions)) == 60 and
             all(sorted(p) == list(range(12)) and all(type(v) is int for v in p) and
                 push(self.levels[0].fundamental, p) == self.levels[0].fundamental for p in actions),
             'complete proper seed action')
        self.actions = (actions,)
        for step in self.steps:
            self.actions += (tuple(step.lift(p) for p in self.actions[-1]),)

    def _range(self, lower, upper):
        need(type(lower) is int and type(upper) is int and 0 <= lower <= upper <= 3,
             'declared tower interval')

    def coarsen(self, chain, lower=0, upper=3):
        self._range(lower, upper)
        self.levels[upper].check(chain)
        for i in range(upper-1, lower-1, -1):
            chain = self.steps[i].coarsen(chain)
        return dict(chain)

    def subdivide(self, chain, lower=0, upper=3):
        self._range(lower, upper)
        self.levels[lower].check(chain)
        for i in range(lower, upper):
            chain = self.steps[i].subdivide(chain)
        return dict(chain)

    def homotopy(self, chain, lower=0, upper=3):
        self._range(lower, upper)
        self.levels[upper].check(chain)
        if lower == upper:
            return {}
        step = self.steps[upper-1]
        inner = self.homotopy(step.coarsen(chain), lower, upper-1)
        return combine((1, step.homotopy(chain)), (1, step.subdivide(inner)))


def load_support(root=None):
    """Re-read and validate source bytes on every call; execute no cached module."""
    root = ROOT if root is None else Path(root)
    source = {}
    for path, expected in SOURCE_SHA256.items():
        location = (root/path).resolve()
        need(location.is_relative_to(root.resolve()), 'source path outside repository')
        data = location.read_bytes()
        need(hashlib.sha256(data).hexdigest() == expected, 'support source changed: '+path)
        source[path] = data
    namespaces = []
    for path in (GEOMETRY, SEARCH):
        namespace = {'__name__': '_support_refinement_source', '__file__': str(root/path)}
        exec(compile(source[path], str(root/path), 'exec'), namespace)
        namespaces.append(namespace)
    geometry, search = namespaces
    rows = geometry['tower']()
    geometry['check_tower'](rows)
    supports = [row['support'] for row in json.loads(source[RECEIPT])['evidence']['topology']]
    need(json.loads(json.dumps(rows)) == supports, 'source-bound support receipt differs')
    local = search['read_local']()
    need(local == geometry['seed'](), 'local ordered incidence differs')
    return Tower(rows, search['automorphisms'](local))


def _check_retract(C, S, H, coarse, fine):
    for simplex in coarse.simplices:
        x = basis(simplex)
        need(C(S(x)) == x and boundary(S(x)) == S(boundary(x)) and H(S(x)) == {},
             'subdivision chain/retraction identities')
    for simplex in fine.simplices:
        x = basis(simplex)
        need(boundary(C(x)) == C(boundary(x)), 'coarsening chain identity')
        need(combine((1, boundary(H(x))), (1, H(boundary(x)))) == combine((1, x), (-1, S(C(x)))),
             'chain homotopy identity')
        need(C(H(x)) == {} and H(H(x)) == {}, 'strong contraction side conditions')


def certificate():
    tower = load_support()
    transitions = []
    for index, step in enumerate(tower.steps):
        _check_retract(step.coarsen, step.subdivide, step.homotopy, step.coarse, step.fine)
        endpoint_symmetries = 0
        for before, after in zip(tower.actions[index], tower.actions[index+1]):
            for simplex in step.fine.simplices:
                x = basis(simplex)
                need(step.coarsen(push(x, after)) == push(step.coarsen(x), before), 'C equivariance')
                need(step.homotopy(push(x, after)) == push(step.homotopy(x), after), 'H equivariance')
            for simplex in step.coarse.simplices:
                x = basis(simplex)
                need(step.subdivide(push(x, before)) == push(step.subdivide(x), after), 'S equivariance')
            c = step.fine.endpoint
            endpoint_symmetries += all(c[after[v]] == before[c[v]] for v in range(step.fine.vertices))
        need(endpoint_symmetries == 1, 'retained endpoint negative control')
        transitions.append(dict(lower=index, upper=index+1,
            coarse_basis_count=len(step.coarse.simplices), fine_basis_count=len(step.fine.simplices),
            chain_retraction=True, equivariant_homotopy=True, strong_side_conditions=True,
            proper_actions_checked=60, endpoint_coarsening_symmetries=endpoint_symmetries))
    compositions = []
    for lower in range(3):
        for upper in range(lower+1, 4):
            C = lambda x: tower.coarsen(x, lower, upper)
            S = lambda x: tower.subdivide(x, lower, upper)
            H = lambda x: tower.homotopy(x, lower, upper)
            _check_retract(C, S, H, tower.levels[lower], tower.levels[upper])
            compositions.append(dict(lower=lower, upper=upper, strong_chain_retract=True))
    swaps = []
    for edge in (s for s in tower.levels[0].simplices if len(s) == 2):
        a, b = edge
        p = next(p for p in tower.actions[0] if p[a] == b and p[b] == a and
                 all(p[p[v]] == v for v in range(12)))
        swaps.append(dict(edge=list(edge), permutation=list(p)))
    files = ('__init__.py', '__main__.py', 'construction.py', 'README.md', 'test_refinement.py')
    return dict(schema='oph.support_refinement.v1', coefficient_field='rational',
        source_sha256=dict(SOURCE_SHA256),
        implementation_sha256={name: hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                               for name in files},
        support_counts=[dict(vertices=level.vertices, edges=30*4**level.index,
                             faces=len(level.faces)) for level in tower.levels],
        transitions=transitions, tower_compositions=compositions,
        endpoint_obstruction=dict(edge_involutions_checked=len(swaps), example=swaps[0]),
        physical_attachment_selected=False, responsive_geometry_constructed=False,
        changes_source_endpoint_process=False)
