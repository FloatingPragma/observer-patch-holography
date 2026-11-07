"""Independent endpoint, spectral and source-binding controls."""
from fractions import Fraction as F
from itertools import combinations, product
from math import exp, log

import numpy as np
import pytest
import sympy as sp

from . import exchange as e, packet as p


def complete_family(ports):
    edges = tuple(combinations(range(ports), 2))
    return e.sector.source.Family((1,)*ports, edges, (1,)*len(edges))


def independent_native_matrix(family, occupied):
    ports = len(family.weights)
    states = tuple(combinations(range(ports), occupied))
    index = {state: i for i, state in enumerate(states)}
    result = sp.zeros(len(states))
    for col, state in enumerate(states):
        mask = sum(1 << i for i in state)
        for edge, coin in product(family.seams, (0, 1)):
            target, _, _ = e.sector.source.image(mask, edge, coin, ports)
            out = tuple(i for i in range(ports) if target >> i & 1)
            result[index[out], col] += F(1, 2*len(family.seams))
    return states, result


@pytest.mark.parametrize('ports', (4, 5, 6))
def test_every_native_column_matches_local_generator_on_every_number_sector(ports):
    family = complete_family(ports)
    for q in range(ports + 1):
        states, native = independent_native_matrix(family, q)
        g = 2*len(family.seams)*(sp.eye(len(states))-native)
        for col, state in enumerate(states):
            actual = e.generator(family, {state: 1}, q)
            assert [actual.get(y, 0) for y in states] == list(g[:, col])


@pytest.mark.parametrize('ports', (4, 6))
def test_full_sector_endpoint_bound_and_generator_counterbound(ports):
    family = complete_family(ports)
    for q in range(ports + 1):
        _, native = independent_native_matrix(family, q)
        vals = np.linalg.eigvalsh(np.array(native, dtype=float))
        assert vals.min() >= 10/21 - 1e-13
        delta = 1-vals
        for n in (1, 2, 10, 100):
            actual = np.max(np.abs(vals**n-np.exp(-n*delta)))
            bound = e.response_budget(len(family.seams), n)['best_response_error_upper']
            assert actual <= float(bound)+1e-13
        if q not in (0, ports):
            actual = np.max(-2*len(family.seams)*np.log(vals)-2*len(family.seams)*delta)
            bound = e.full_generator_obstruction(ports, len(family.seams), q)
            assert actual >= float(bound['rescaled_log_minus_generator_norm_lower'])-1e-12


def test_scalar_full_sector_and_floor_time_inequalities_on_independent_grid():
    for n in (1, 2, 3, 10, 100, 1000):
        for probability in np.linspace(10/21, 1, 301):
            exact = exp(-n*(1-probability))-probability**n
            assert -1e-14 <= exact <= 3/(5*n)+1e-14
            for kappa in (F(0), F(1, 7), F(99, 100)):
                target = exp(-(n+float(kappa))*(1-probability))
                assert abs(target-probability**n) <= float((F(3, 5)+kappa/2)/n)+1e-14


def test_bounded_primitive_degree_generator_bound_on_larger_occupation():
    family = complete_family(6)
    states2, native2 = independent_native_matrix(family, 2)
    states3, native3 = independent_native_matrix(family, 3)
    lowering = sp.Matrix([[int(i in state) for state in states2] for i in range(6)])
    primitive = sp.Matrix.hstack(*lowering.nullspace())
    inclusion = sp.Matrix([[int(set(a) <= set(b)) for a in states2] for b in states3])
    image = inclusion*primitive
    assert native3*image == inclusion*native2*primitive
    gram = np.array(image.T*image, dtype=float)
    norm_vals, norm_vecs = np.linalg.eigh(gram)
    ortho = np.array(image, dtype=float) @ (norm_vecs/np.sqrt(norm_vals))
    restricted = ortho.T @ np.array(native3, dtype=float) @ ortho
    vals = np.linalg.eigvalsh(restricted)
    g = 30*(1-vals)
    error = np.max(-30*np.log(vals)-g)
    bound = e.degree_budget(15, 5, 2)
    assert g.max() <= bound['spectral_ceiling']+1e-12
    assert error <= float(bound['generator_error_upper'])+1e-12


def test_exact_moments_and_intervals_against_independent_dense_transfer():
    family = complete_family(6)
    h = e.sector.contrast((0, 1), (2, 3), 6)
    states, native = independent_native_matrix(family, 2)
    vector = sp.Matrix([h.get(x, 0) for x in states])
    generator = 30*(sp.eye(len(states))-native)
    expected = tuple(F((vector.T*(generator**k)*vector)[0], 4) for k in range(9))
    assert e.moments(family, h, 2, 8) == expected
    vals, vecs = np.linalg.eigh(np.array(generator, dtype=float))
    v = np.array(vector, dtype=float).reshape(-1)
    weights = (vecs.T @ v)**2/4
    for n in (0, 1, 2, 7, 15):
        row = e.response_intervals(family, h, 2, n, 8)
        actual = F((vector.T*(native**n)*vector)[0], 4)
        assert row['native_response_interval'][0] <= actual <= row['native_response_interval'][1]
        exact_exp = float(weights @ np.exp(-float(row['native_time'])*vals))
        lo, hi = row['local_exchange_response_interval']
        assert float(lo)-1e-13 <= exact_exp <= float(hi)+1e-13
        assert hi-lo <= row['exponential_interval_width_upper']


def test_small_complete_backacting_experiment_matches_response_not_naive_score():
    family = complete_family(6)
    h = e.sector.contrast((0, 1), (2, 3), 6)
    for q in (2, 3, 4):
        variance = 4*e.sector.norm_factor(6, 2, q)
        for n in (0, 1, 2):
            observed = p.records.direct_probe(family, q, n, (0, 1), (2, 3))/variance
            lo, hi = e.response_intervals(family, h, 2, n)['native_response_interval']
            assert lo <= observed <= hi


def test_time_zero_floor_and_rational_budget_composition():
    assert e.response_budget(630, 0, 6, 2)['best_response_error_upper'] == 0
    assert e.rounded_time_budget(630, F(0))['best_response_error_upper'] == 0
    assert e.rounded_time_budget(630, F(1, 10000))['attempts'] == 0
    row = e.rounded_time_budget(630, F(1, 32), 6, 2)
    assert row['attempts'] == 39 and row['fractional_attempt'] == F(3, 8)
    assert e.preparation_composition(F(1, 20), F(1, 1000), F(1, 500)) == F(11, 500)
    assert e.full_generator_obstruction(240, 630, 120)['rescaled_log_minus_generator_norm_lower'] >= F(630, 16)


@pytest.mark.parametrize('bad', (True, 1.0, -1, '1'))
def test_exact_domains_reject_aliased_or_invalid_numbers(bad):
    with pytest.raises(ValueError): e.native_time(630, bad)
    with pytest.raises(ValueError): e.degree_budget(630, 6, bad)
    with pytest.raises(ValueError): e.rounded_time_budget(630, bad)
    with pytest.raises(ValueError): e.full_generator_obstruction(240, 630, bad)


def test_source_and_degree_guards_are_not_silent_extensions():
    with pytest.raises(ValueError): e.generator(e.sector.source.local_family(), {(0, 1): 1}, 2)
    with pytest.raises(ValueError): e.degree_budget(12, 6, 2)
    with pytest.raises(ValueError): e.degree_budget(630, True, 2)
    with pytest.raises(ValueError): e.preparation_composition(F(0), F(0), F(0))
    family = complete_family(6)
    with pytest.raises(ValueError): e.response_intervals(family, {(0, 1): 1}, 2, 1)
    with pytest.raises(ValueError): e.moments(family, {}, 2, 2)
    with pytest.raises(ValueError): e.moments(family, {(0, 1): 1}, 2, 13)
    with pytest.raises(ValueError): e.response_intervals(family, {(0, 1): 1}, 2, 1, 3)


def test_complete_record_length_control_is_not_endpoint_equivalence():
    # The count of an unread Poisson representation is itself observable if
    # exported. It cannot be identified with a deterministic native count.
    assert 1-exp(-1) > 1/2
    assert 1-exp(-1) > float(e.response_budget(630, 100)['best_response_error_upper'])


def test_native_count_mutation_is_rejected(monkeypatch):
    old = e.sector.counts_step
    def changed(*args, **kwargs):
        result = old(*args, **kwargs)
        key = next(iter(result))
        result[key] += 1
        return result
    monkeypatch.setattr(e.sector, 'counts_step', changed)
    with pytest.raises(ArithmeticError, match='owning transfer counts'):
        p.check_source_generator()


@pytest.fixture(scope='module')
def produced():
    return p.produce()


def test_retained_receipt_has_typed_exact_replay_and_source_closure(produced):
    assert p.canonical(p.strict_json.load(p.OUTPUT)) == p.canonical(produced)
    assert produced['schema'] == 'oph.native_exchange_clock.v1'
    assert len(produced['evaluated_federations']) == 3
    assert all(name in produced['source_bindings'] for name in p.SOURCE_FILES)
    assert not any(name.startswith('docs/') for name in produced['source_bindings'])


@pytest.mark.parametrize('field,value', (('schema', True),
    ('native_attachment', {'fair_H66_columns': True}), ('source_bindings', {})))
def test_receipt_type_or_evidence_mutation_fails_closed(produced, field, value, tmp_path, monkeypatch):
    altered = dict(produced); altered[field] = value
    path = tmp_path/'changed.json'; path.write_text(p.canonical(altered))
    monkeypatch.setattr(p, 'produce', lambda: produced)
    with pytest.raises(ValueError, match='exact replay'):
        p.verify(path)


def test_bool_integer_alias_is_rejected_by_canonical_receipt(produced, tmp_path, monkeypatch):
    altered = p.strict_json.loads(p.canonical(produced))
    altered['evaluated_federations'][0]['level'] = False
    path = tmp_path/'changed.json'; path.write_text(p.canonical(altered))
    monkeypatch.setattr(p, 'produce', lambda: produced)
    with pytest.raises(ValueError, match='exact replay'):
        p.verify(path)
