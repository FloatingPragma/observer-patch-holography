"""Independent minors, local fillings, tower controls and damaged-input checks."""
from fractions import Fraction as Q
from itertools import combinations, permutations
import json
from pathlib import Path
import sys
import types

import pytest

from . import construction as c


@pytest.fixture(scope='module')
def tower():
    return c.load_support()


@pytest.fixture(scope='module')
def checked():
    return c.certificate()


def determinant(matrix):
    n = len(matrix)
    result = Q(0)
    for p in permutations(range(n)):
        term = Q((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n)))
        for i, j in enumerate(p):
            term *= matrix[i][j]
        result += term
    return result


def minor_oracle(simplex, vertex_images):
    vertices = sorted({v for source in simplex for (v,) in vertex_images[source]})
    result = {}
    for target in combinations(vertices, len(simplex)):
        value = determinant([[vertex_images[source].get((v,), Q(0)) for source in simplex]
                             for v in target])
        if value:
            result[target] = value
    return result


def test_complete_certificate_matches_retained_result(checked):
    path = Path(__file__).with_name('results.json')
    assert checked == json.loads(path.read_text())
    assert [x['fine_basis_count'] for x in checked['transitions']] == [242, 962, 3842]
    assert len(checked['tower_compositions']) == 6
    assert not checked['physical_attachment_selected']
    assert not checked['responsive_geometry_constructed']
    assert not checked['changes_source_endpoint_process']


def test_every_extensor_column_matches_independent_determinants(tower):
    for step in tower.steps:
        vertex_images = {v: ({(v,): Q(1)} if v < step.coarse.vertices else
                            {(a,): Q(1, 2) for a in step.endpoints[v]})
                         for v in range(step.fine.vertices)}
        for simplex in step.fine.simplices:
            assert step.coarsen(c.basis(simplex)) == minor_oracle(simplex, vertex_images)


def test_three_level_composition_matches_direct_vertex_minors(tower):
    images = {v: tower.coarsen(c.basis((v,))) for v in range(tower.levels[-1].vertices)}
    for simplex in tower.levels[-1].simplices:
        assert tower.coarsen(c.basis(simplex)) == minor_oracle(simplex, images)
    for level in range(1, 4):
        assert tower.coarsen(tower.levels[level].fundamental, 0, level) == tower.levels[0].fundamental
        assert tower.subdivide(tower.levels[0].fundamental, 0, level) == tower.levels[level].fundamental


def solve_full_column_rank(matrix, rhs):
    """Independent rational elimination, including all overdetermined rows."""
    rows = [list(map(Q, row))+[Q(value)] for row, value in zip(matrix, rhs)]
    columns = len(matrix[0])
    for j in range(columns):
        pivot = next(i for i in range(j, len(rows)) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        divisor = rows[j][j]
        rows[j] = [x/divisor for x in rows[j]]
        for i in range(len(rows)):
            if i != j:
                scale = rows[i][j]
                rows[i] = [a-scale*b for a, b in zip(rows[i], rows[j])]
    assert all(not row[-1] for row in rows[columns:])
    return [rows[j][-1] for j in range(columns)]


def test_local_homotopy_is_unique_rational_filling(tower):
    step = tower.steps[0]
    parent = step.coarse.faces[0]
    children = step.subdivide(c.basis(parent))
    fine_edges = sorted({edge for child in children for edge in combinations(child, 2)})
    boundaries = [c.boundary(c.basis(child)) for child in children]
    matrix = [[column.get(edge, 0) for column in boundaries] for edge in fine_edges]
    for edge in fine_edges:
        x = c.basis(edge)
        cycle = c.combine((1, x), (-1, step.subdivide(step.coarsen(x))),
                          (-1, step.homotopy(c.boundary(x))))
        coefficients = solve_full_column_rank(matrix, [cycle.get(e, 0) for e in fine_edges])
        independent = {face: value for face, value in zip(children, coefficients) if value}
        assert step.homotopy(x) == independent


def test_actual_endpoint_rule_and_all_endpoint_choices_have_obstruction(tower):
    step = tower.steps[0]
    a, b = next(iter(step.midpoint))
    middle = step.midpoint[a, b]
    p = next(p for p in tower.actions[0] if p[a] == b and p[b] == a)
    lifted = step.lift(p)
    assert lifted[middle] == middle
    assert all(p[v] != v for v in (a, b))
    assert step.fine.endpoint[lifted[middle]] != p[step.fine.endpoint[middle]]
    assert step.coarsen(c.basis((middle,))) == {(a,): Q(1, 2), (b,): Q(1, 2)}
    assert c.push(step.coarsen(c.basis((middle,))), p) == step.coarsen(c.basis((middle,)))


def test_nonzero_homotopy_cannot_be_replaced_by_identity_or_zero(tower):
    step = tower.steps[0]
    middle = next(iter(step.endpoints))
    x = c.basis((middle,))
    assert step.subdivide(step.coarsen(x)) != x
    assert step.homotopy(x)
    assert c.boundary(step.homotopy(x)) == c.combine((1, x), (-1, step.subdivide(step.coarsen(x))))


@pytest.mark.parametrize('damage', ['face_orientation', 'duplicate_face', 'wrong_ancestry', 'missing_level', 'wrong_count'])
def test_changed_incidence_and_refinement_fail_closed(tower, damage):
    rows = [dict(level=x.index, vertices=x.vertices, faces=[list(f) for f in x.faces],
                 coarsen=None if x.endpoint is None else list(x.endpoint)) for x in tower.levels]
    if damage == 'face_orientation':
        rows[1]['faces'][0][1:] = reversed(rows[1]['faces'][0][1:])
    elif damage == 'duplicate_face':
        rows[1]['faces'][0] = rows[1]['faces'][1]
    elif damage == 'wrong_ancestry':
        middle, endpoints = next(iter(tower.steps[0].endpoints.items()))
        rows[1]['coarsen'][middle] = next(v for v in range(12) if v not in endpoints)
    elif damage == 'missing_level':
        rows.pop(1)
    else:
        rows[1]['vertices'] = True
    with pytest.raises(ValueError):
        c.Tower(rows, tower.actions[0])


@pytest.mark.parametrize('chain', [{(False,): Q(1)}, {(0,): True}, {(0,): 0.5},
                                   {(1, 0): Q(1)}, {(999,): Q(1)}, {(0, 0): Q(1)}])
def test_undeclared_chain_inputs_are_rejected(tower, chain):
    with pytest.raises(ValueError):
        tower.coarsen(chain)


@pytest.mark.parametrize('interval', [(False, 1), (0, True), (-1, 3), (2, 1), (0, 4)])
def test_tower_bounds_are_exact(tower, interval):
    with pytest.raises(ValueError):
        tower.homotopy({}, *interval)


@pytest.mark.parametrize('path', tuple(c.SOURCE_SHA256))
def test_source_changes_rejected_even_after_prior_load(tmp_path, path):
    for source in c.SOURCE_SHA256:
        destination = tmp_path/source
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes((c.ROOT/source).read_bytes())
    c.load_support(tmp_path)
    changed = tmp_path/path
    changed.write_bytes(changed.read_bytes()+b'\n')
    with pytest.raises(ValueError, match='support source changed'):
        c.load_support(tmp_path)


def test_cached_geometry_module_cannot_replace_checked_source(monkeypatch):
    module = types.ModuleType('source_selection_model.geometry')
    module.tower = lambda: []
    monkeypatch.setitem(sys.modules, 'source_selection_model.geometry', module)
    assert c.load_support().levels[-1].vertices == 642


def test_missing_proper_action_rejected(tower):
    rows = [dict(level=x.index, vertices=x.vertices, faces=[list(f) for f in x.faces],
                 coarsen=None if x.endpoint is None else list(x.endpoint)) for x in tower.levels]
    with pytest.raises(ValueError, match='complete proper seed action'):
        c.Tower(rows, tower.actions[0][:-1])


def test_certificate_catches_broken_homotopy(monkeypatch):
    monkeypatch.setattr(c.Refinement, '_homotopy', lambda self, simplex: {})
    with pytest.raises(ValueError, match='chain homotopy identity'):
        c.certificate()
