"""Independent complete Kraus effects, slice kernels and preparation bounds."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb

import pytest
import sympy as sp

from . import preparation as p


def graph(ports=6):
    seams = tuple(sorted({(i, i+1) for i in range(ports-1)} | {(0, ports-1)}))
    return p.s.source.Family((1,)*ports, seams, (1,)*len(seams))


def direct_effect(family, q, warmup, lag, past, future, law):
    """Actual source label/current maps on all matrix units; no dual formula."""
    ports = len(family.weights)
    states = tuple(combinations(range(ports), q)); index = {x: i for i, x in enumerate(states)}
    commands = [('uniform', None)]*warmup+[('out', e) for e in past]
    commands += [('uniform', None)]*lag+[('in', e) for e in future]
    effect = sp.eye(len(states))
    for kind, fixed in reversed(commands):
        result = sp.zeros(len(states)); edges = family.seams if fixed is None else (fixed,)
        laws = ('C', 'E') if fixed is None and law == 'mixture' else (law if fixed is None else 'E',)
        for edge, tag, coin in product(edges, laws, (0, 1)):
            decoding = p.s.source.interpreter.circuit.record_values(tag, coin)
            for raw in range(1 << p.s.source.interpreter.circuit.RECORD_BITS[tag]):
                label = decoding.get(raw)
                current = (label if tag == 'C' else
                           (label-(2*coin-1)*(abs(label) % 2))//2) if label is not None else 0
                score = (label-2*current if kind == 'out' else label if kind == 'in' else 1) if label is not None else 0
                matrix = sp.zeros(len(states))
                for column, state in enumerate(states):
                    target, j, d = p.r.native_step(state, edge, coin)
                    if (j if tag == 'C' else d) == label:
                        matrix[index[target], column] = 1
                result += sp.Rational(score, 2*len(edges)*len(laws))*matrix.T*effect*matrix
        effect = result
    return states, effect


@pytest.mark.parametrize('law', ('C', 'E', 'mixture'))
@pytest.mark.parametrize('q', (1, 2, 3, 4, 5))
def test_complete_quantum_effect_equals_degree_four_dual(law, q):
    family = graph(); past = ((0, 1), (2, 3)); future = ((0, 5), (3, 4))
    states, actual = direct_effect(family, q, 1, 2, past, future, law)
    expected = p.evolve(family, p.cross_effect(family, past, future, 2), 1)
    diagonal = [p.evaluate(expected, x, 6) for x in states]
    assert actual == sp.diag(*diagonal)
    assert all(-1 <= value <= 1 for value in diagonal)
    if q in (1, 5):
        assert not any(diagonal)


@pytest.mark.parametrize('q', (2, 3, 4))
def test_reproducing_kernel_has_constant_diagonal_without_graph_symmetry(q):
    ports = 6; degree = min(4, q, ports-q)
    states = tuple(combinations(range(ports), q)); size = len(states)
    labels = tuple(a for d in range(degree+1) for a in combinations(range(ports), d))
    evaluation = sp.Matrix([[int(set(a) <= set(x)) for a in labels] for x in states])
    basis = sp.Matrix.hstack(*evaluation.columnspace())
    projection = basis*(basis.T*basis).inv()*basis.T-sp.ones(size)/size
    rank = comb(ports, degree)-1
    assert projection.trace() == rank
    # L2(tau) reproducing kernel is size times the Euclidean projector.
    assert all(size*projection[i, i] == rank for i in range(size))
    # A path graph has no transitive vertex symmetry. The kernel result
    # concerns its full polynomial space, not automorphisms of that graph.
    assert len({1 if i in (0, ports-1) else 2 for i in range(ports)}) == 2


def test_degree_four_kernel_is_proper_subspace_on_large_slice():
    ports = 10; q = 5; degree = 4
    states = tuple(combinations(range(ports), q)); labels = tuple(combinations(range(ports), degree))
    # Exact inclusion rank is tested separately by the owning decomposition;
    # this independent numerical projection checks the proper non-full case.
    import numpy as np
    matrix = np.array([[set(a) <= set(x) for a in labels] for x in states], dtype=float)
    left, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    assert np.count_nonzero(singular > 1e-10) == comb(ports, degree)
    projection = left@left.T-np.ones((len(states), len(states)))/len(states)
    expected = F(comb(ports, degree)-1, len(states))
    assert np.max(np.abs(np.diag(projection)-float(expected))) < 2e-14


def test_preparation_is_needed_and_warmup_reduces_exact_effect_error():
    family = p.s.source.Family((1,)*4, ((0, 1), (1, 2), (2, 3)), (1,)*3)
    value = p.score_effect(family, (0, 1), (2, 3), 0)
    states = tuple(combinations(range(4), 2)); center = F(2, 3)
    errors = []
    for n in range(7):
        advanced = p.evolve(family, value, n)
        errors.append(max(abs(p.evaluate(advanced, x, 4)-center) for x in states))
        assert p.uniform_mean(advanced, 4, 2) == center
    assert errors[0] == F(2, 3) and errors[-1] < errors[0]
    assert all(a >= b for a, b in zip(errors, errors[1:]))


def test_reference_valued_score_uses_only_initial_diagonal_blocks():
    family = graph(4); past = ((0, 1), (2, 3))
    states, effect = direct_effect(family, 2, 2, 1, past, past, 'mixture')
    # An exact entangled vector with nonzero off-diagonal source blocks.
    amplitudes = {(0, 0): sp.Rational(1, 2), (1, 1): sp.Rational(1, 2),
                  (4, 0): sp.Rational(1, 2), (5, 1): sp.Rational(1, 2)}
    joint = sp.zeros(2*len(states))
    for (i, a), left in amplitudes.items():
        for (j, b), right in amplitudes.items():
            joint[2*i+a, 2*j+b] = left*right
    weighted = sp.kronecker_product(effect, sp.eye(2))*joint
    reduced = sp.Matrix(2, 2, lambda a, b: sum(weighted[2*i+a, 2*i+b] for i in range(len(states))))
    direct = sp.Matrix(2, 2, lambda a, b: sum(effect[i, i]*joint[2*i+a, 2*i+b]
                                             for i in range(len(states))))
    assert reduced == direct
    assert any(joint[i, j] for i in range(2*len(states)) for j in range(2*len(states)) if i//2 != j//2)


@pytest.mark.parametrize('ports,q', ((4, 2), (12, 2), (240, 80), (960, 480), (10**122, 2), (10**122, 10**121)))
def test_integer_budget_certifies_normalized_error_without_sector_enumeration(ports, q):
    result = p.preparation_budget(ports, q, F(1, 1000), 10)
    v = F(result['stationary_variance']); c = result['normalized_prefactor_integer_upper']
    assert (c-1)**2*v < result['mean_zero_polynomial_rank'] <= c*c*v
    assert F(result['normalized_score_error_upper']) <= F(1, 1024)
    assert result['warmup_attempts'] == 1000*result['halving_blocks']
    assert F(result['unnormalized_reference_score_trace_norm_upper']) == v*F(result['normalized_score_error_upper'])


@pytest.mark.parametrize('law', ('C', 'E', 'mixture'))
def test_actual_h66_interpreter_records_match_sparse_effect(law):
    result = p.source_effect_control(1, 1, law)
    assert result['source_columns'] == 66 and result['checked_raw_slots'] > 100
    assert result['zero_operator_slots'] > 0 and result['degree_upper'] <= 4


def test_all_input_h66_warmup_is_propagated_with_exact_source_counts():
    result = p.prepared_h66_control()
    assert result['budget']['warmup_attempts'] == 352
    assert result['complete_attempts'] == 358
    assert F(result['all_input_normalized_error_upper_80bit']) <= F(result['budget']['normalized_score_error_upper'])
    assert result['expanded_history_executed'] is False


def test_generic_polynomial_moments_and_snapshot_bound_are_exact():
    family = graph(); q = 3
    left = {(0, 1): F(3, 2), (2,): F(-1), (): F(1, 3)}
    right = {(1, 4): F(2), (3, 5): F(-2), (4,): F(1, 5)}
    effect = p.snapshot_effect(family, left, right, 2)
    values = [p.evaluate(effect, x, 6) for x in combinations(range(6), q)]
    moments = p.polynomial_moments(effect, 6, q)
    mean = sum(values, F())/len(values)
    assert moments['mean'] == mean
    assert moments['second_moment'] == sum((x*x for x in values), F())/len(values)
    assert moments['centered_l2_squared'] == sum(((x-mean)**2 for x in values), F())/len(values)
    budget = p.polynomial_budget(6, q, F(1, 18), moments['centered_l2_squared'], F(2), 8)
    assert F(budget['normalized_score_error_upper']) <= F(1, 256)
    zero = p.polynomial_budget(6, q, F(1, 18), F(0), F(1), 8)
    assert zero['warmup_attempts'] == 0 and zero['normalized_score_error_upper'] == '0'
    for norm, scale in ((0.0, F(1)), (F(-1), F(1)), (F(1), 0.5), (F(1), F(0))):
        with pytest.raises(ValueError): p.polynomial_budget(6, q, F(1, 18), norm, scale)


@pytest.mark.parametrize('bad', (True, 1.0, -1))
def test_budget_and_effect_reject_ambiguous_inputs(bad):
    with pytest.raises(ValueError): p.preparation_budget(bad, 2, F(1, 10))
    with pytest.raises(ValueError): p.preparation_budget(12, bad, F(1, 10))
    with pytest.raises(ValueError): p.evolve(graph(), {}, bad)


def test_boundary_domains_and_exact_source_limits_are_preserved():
    for q in (0, 1, 11, 12):
        with pytest.raises(ValueError): p.preparation_budget(12, q, F(1, 10))
    for gap in (0, 0.1, F(0), F(2)):
        with pytest.raises(ValueError): p.preparation_budget(12, 2, gap)
    with pytest.raises(ValueError): p.preparation_budget(12, 2, F(1, 10), True)
    with pytest.raises(ValueError): p.multiply({(0, 1, 2): 1}, {(3, 4, 5): 1}, 6)
    with pytest.raises(ValueError): p.polynomial({(0,): True}, 6)
    with pytest.raises(ValueError): p.cross_effect(graph(), ((0, 1), (1, 2)), ((0, 1), (2, 3)), 0)
    with pytest.raises(ValueError): p.prepared_program(0, 2, 2**32-3, 0)
    with pytest.raises(ValueError): p.prepared_program(True, 2, 0, 0)
    with pytest.raises(ValueError): p.source_budgets(3)
    with pytest.raises(ValueError): p.source_budgets(0, lag=2**40)
    with pytest.raises(ValueError): p.source_budgets(0, lag=1.5)
    program = p.prepared_program(0, 2, 5, 2)
    assert program['declared_resources']['attempts_upper'] == 11
    assert program['memory']['overflow'] == 'export_oldest'
    initial = {(3, 0, 3, 0): F(1)}
    stopped = p.s.source.interpreter.execute(program, initial, max_branches=1)
    assert stopped.status == 'RESOURCE_STOP' and stopped.attempted_steps == 0
    assert stopped.branches[0].density == initial
