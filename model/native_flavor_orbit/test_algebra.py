"""Independent full-current, complete-record and complex-basis controls."""
import numpy as np
import pytest
import sympy as s

from ..common_source import native_repair_mass as native
from ..common_source import native_family_record_memory as source_memory
from ..common_source import native_matter_attachment as matter
from . import algebra as a


def test_complete_current_commutant_includes_all_five_charge_blocks():
    report = a.commutant_certificate()
    assert [(r['q6Y'], r['dimension'], r['exact_rank']) for r in report['charge_blocks']] == [
        (-4, 3, 8), (-3, 2, 3), (1, 6, 35), (2, 3, 8), (6, 1, 0)]
    assert report['family_commutant_complex_dimension'] == 45


def exact_unitaries():
    rotation = s.Matrix([[s.Rational(3, 5), s.Rational(4, 5), 0],
                         [-s.Rational(4, 5), s.Rational(3, 5), 0], [0, 0, 1]])
    cycle = s.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    return {charge: rotation**index*s.diag(1, s.I**index, (-1)**index)*cycle**index
            for index, charge in enumerate(a.current_blocks())}


def test_all_45_rows_and_both_Higgs_components_transform_by_the_same_current_frame():
    units = exact_unitaries()
    frame = a.embedding(units)
    assert a.zero(frame.H*frame-s.eye(45))
    for current in matter.representation()['matter_currents']:
        current45 = s.kronecker_product(current, s.eye(3))
        assert a.zero(current45*frame-frame*current45)
    coefficients = {
        'up': s.Matrix([[1, s.I, 0], [0, 2, 1], [1-s.I, 0, 3]]),
        'down': s.Matrix([[2, 1, 0], [s.I, 3, 1], [0, 0, 4]]),
        'charged_lepton': s.Matrix([[1, 0, 1], [0, 3, 0], [s.I, 0, 5]])}
    source = a.source_gates()
    for component in range(2):
        before, after = s.zeros(45), s.zeros(45)
        for name, row in source.items():
            left, right = row['left_projector'], row['right_projector']
            directed = left*row['components'][component]*right
            lc = 1 if name in ('up', 'down') else -3
            rc = {'up': -4, 'down': 2, 'charged_lepton': 6}[name]
            Y = coefficients[name]
            transformed = units[lc].H*Y*units[rc].conjugate()
            term = s.kronecker_product(directed, Y.conjugate())
            term_after = s.kronecker_product(directed, transformed.conjugate())
            before += term+term.T
            after += term_after+term_after.T
        assert a.zero(frame.T*before*frame-after)


def test_rotating_only_one_quark_weak_component_changes_the_actual_current():
    data = matter.representation()
    selected = list(data['selected'])
    t3 = data['operators']['T3'].extract(selected, selected)
    bad = s.eye(45)
    for i in a.current_blocks()[1]:
        if t3[i, i] > 0:
            bad[3*i:3*i+3, 3*i:3*i+3] = s.diag(1, -1, 1)
    assert any(not a.zero(s.kronecker_product(g, s.eye(3))*bad
                          -bad*s.kronecker_product(g, s.eye(3)))
               for g in data['matter_currents'])


def test_independent_6144_record_sum_reconstructs_every_duty_coefficient():
    left, right = s.Rational(1, 3), s.Rational(2, 5)
    p, q = left/12, right/6
    source = source_memory.architecture_from_factors(p, q)
    weights = {policy: s.zeros(12, 1) for policy in a.POLICIES}
    total = s.Integer(0)
    count = 0
    for cc, cr, ec, er, advance, native_weight in native.native_records():
        for cl, el in ((0, 0), (0, 1), (1, 0), (1, 1)):
            weight = native_weight*(left if cl else 1-left)*(right if el else 1-right)
            c, e = advance and bool(cl), ec == 0 and er in (4, 5) and bool(el)
            for old in range(12):
                count += 1
                total += source['stationary'][old]*weight
                destination = source['gC'][old] if c else 2
                if c and e:
                    destination = source['gE'][destination]
                if destination in (0, 1):
                    destination = 2
                for policy in a.POLICIES:
                    if policy == 'always' or policy == 'C' and c or policy == 'E' and e or policy == 'CE' and c and e:
                        weights[policy][destination] += source['stationary'][old]*weight
    assert count == 6144 and total == 1
    for policy, values in weights.items():
        expected = sum((weight*P for weight, P in zip(values, source['projectors'], strict=True)), s.zeros(3))
        assert a.zero(expected-a.coefficient(p, q, policy))
    assert any(w == 0 for *_, w in native.native_records())


@pytest.mark.parametrize('p,q', [(s.Rational(1, 12), s.Rational(1, 6)),
                               (s.Rational(1, 8), s.Rational(3, 4)),
                               (0, 0), (s.Rational(1, 4), 0), (s.Rational(1, 4), 1)])
def test_stationary_and_duty_formulas_cover_boundaries_without_conditioning(p, q):
    source = source_memory.architecture_from_factors(p, q)
    assert a.zero(a.coefficient(p, q)-source['Y'])
    for policy, trace in zip(a.POLICIES, (1, p, q, p*q), strict=True):
        result = a.coefficient(p, q, policy)
        assert s.simplify(s.trace(result)-trace) == 0
        assert np.linalg.eigvalsh(np.asarray(result, float)).min() > -1e-12


def test_real_source_mixing_and_CP_are_distinct_basis_invariant_statements():
    up = a.coefficient(s.Rational(1, 12), s.Rational(1, 7))
    down = a.coefficient(s.Rational(1, 12), s.Rational(1, 8))
    assert not a.zero(up*down-down*up)
    assert not a.zero(up**2*down**2-down**2*up**2)
    assert a.cp_invariant(up, down) == 0
    units = exact_unitaries()
    units[1] = units[1]*s.diag(1, s.I, 1)
    transformed_up = units[1].H*up*units[-4].conjugate()
    transformed_down = units[1].H*down*units[2].conjugate()
    assert any(s.im(x) != 0 for x in transformed_up)
    assert a.cp_invariant(transformed_up, transformed_down) == 0
    assert a.zero(transformed_up*transformed_up.H-units[1].H*up**2*units[1])
    common_down = units[1].H*up*units[2].conjugate()
    assert a.zero(transformed_up*transformed_up.H-common_down*common_down.H)


def test_complex_flavor_control_has_nonzero_CP_invariant():
    up = s.diag(1, 2, 3)
    down = s.Matrix([[1, 1, 0], [0, 2, s.I], [1, 0, 3]])
    assert a.cp_invariant(up, down) != 0


@pytest.mark.parametrize('p,q', [(-1, 0), (1, 0), (0, -1), (0, 2), (s.I, 0), (0.1, 0)])
def test_invalid_memory_probabilities_are_rejected(p, q):
    with pytest.raises(ValueError):
        a.coefficient(p, q)


def test_incomplete_or_nonunitary_family_embeddings_are_rejected():
    values = exact_unitaries()
    values.pop(6)
    with pytest.raises(ValueError):
        a.embedding(values)
    values = exact_unitaries()
    values[6] = 2*s.eye(3)
    with pytest.raises(ValueError):
        a.embedding(values)
