"""Exact native record reversal, reflection and transfer-Hilbert controls."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import product

import pytest
import sympy as s

from . import response as r, transfer as t
from ..observer import native


@pytest.fixture(scope='module')
def data():
    value = r.data();r.validate(value)
    return value


@pytest.fixture(scope='module')
def record_words():
    a = ('repair','E',(0,5),0,5)
    b = ('repair','C',(1,7),0,1)
    c = ('repair','C',(0,1),1,0)
    return ((), (a,), (t.reverse_label(a),), (c,), (a,b),
            (('hold',None,None,None,None),),
            (('repair','E',(0,5),0,3),))


def test_all_721_actual_raw_labels_have_weighted_adjoint(data):
    value = t.branch_reversal_certificate(data)
    assert value['complete_raw_slots'] == 721
    assert value['weighted_source_columns_checked'] == 7920
    assert value['zero_slots_including_hold'] == 421
    assert value['reserved_abort_slots'] == 240
    assert value['moving_slots'] == 120


def test_all_empty_native_slots_have_literal_zero_score(data):
    empty = [row for row in data['branches'] if not row.pairs]
    assert len(empty) == 421
    assert all(row.feature(edge) == 0 for row in empty for edge in r.local.SEAMS)
    assert all(abs(row.feature(edge)) <= 60 for row in data['branches'] for edge in r.local.SEAMS)


def test_adjoint_rejects_changed_weight_and_missing_zero_label(data):
    changed = dict(data);rows = list(data['branches'])
    at = next(i for i,row in enumerate(rows) if any(a!=b for a,b in row.pairs))
    rows[at] = replace(rows[at],weight=rows[at].weight*2)
    changed['branches'] = tuple(rows)
    with pytest.raises(ArithmeticError,match='weighted'):
        t.branch_reversal_certificate(changed)
    changed['branches'] = data['branches'][1:]
    with pytest.raises(ArithmeticError,match='721'):
        t.branch_reversal_certificate(changed)


def actual_word_probability(data, word):
    """Unnormalized indicator probability using each complete source kernel."""
    value = r.source();density = {(i,0,i,0): F(p) for i,p in enumerate(data['pi'])}
    memory = archive = ()
    for label in word:
        all_outputs = r.interpreter._outcomes(value,value['commands'][0],density,memory,archive)
        assert len(all_outputs) == 721
        matched = [out for out in all_outputs if
                   (out.record.kind,out.record.law,out.record.seam,out.record.coin,out.record.raw)==label]
        assert len(matched) == 1
        outcome, = matched
        density,memory,archive = outcome.density,outcome.memory,outcome.archive
    return r.interpreter.trace(density)


def test_complete_word_probabilities_and_true_reversal_match_native_kernel(data,record_words):
    for word in record_words:
        actual = actual_word_probability(data,word)
        assert actual == sum(t.forward_word(data,word))
        assert actual == (data['pi'].T*t.conditional_word(data,word))[0]
        assert actual == actual_word_probability(data,t.reverse_word(word))
        assert t.reverse_word(t.reverse_word(word)) == word


@pytest.mark.parametrize('gap',(0,1,2))
def test_full_reflected_word_gram_matches_chronological_process(data,record_words,gap):
    words = record_words if gap < 2 else record_words[:3]
    gram = t.reflected_gram(data,words,gap)
    assert gram == gram.T
    for i,j in product(range(len(words)),repeat=2):
        assert gram[i,j] == t.reflected_joint(data,words[i],words[j],gap)
    if gap == 0:
        for i,j in ((1,2),(3,4),(4,4),(6,2)):
            assert gram[i,j] == actual_word_probability(data,t.reverse_word(words[i])+words[j])
    coefficients = s.Matrix([s.Rational(i+1,7)+s.I*s.Rational(2-i,5) for i in range(len(words))])
    norm = s.expand((coefficients.H*gram*coefficients)[0])
    assert s.im(norm) == 0 and s.re(norm) >= 0


def independent_two_record_form(data,reflected):
    """Use only endpoint arithmetic for the two selected signed contrasts."""
    edges = ((0,1),(0,5));signs = (-1,1)
    theta = (25,25)+(8,)*10
    result = F()
    for state,pi in zip(data['states'],data['pi'],strict=True):
        for edge,sign in zip(edges,signs,strict=True):
            a,b = edge
            for coin in (0,1):
                after,(_,current) = native.native(state,edge,coin)
                first = state[a]-state[b]-(2*current if reflected else 0)
                future = sum(sgn*(after[u]-after[v]) for (u,v),sgn in zip(edges,signs,strict=True))
                # The first record score's 60 cancels its tagged-seam probability.
                result += F(pi)*F(theta[a if coin else b],theta[a]+theta[b])*sign*first*future
    mean = data['mean'][1]-data['mean'][0]
    return result-mean*mean


def test_time_oriented_contrast_has_negative_naive_and_positive_true_form(data):
    vector = s.Matrix([-1,1]+[0]*28)
    ordinary = (vector.T*r.covariance(data,1)*vector)[0]
    reflected = (vector.T*t.reflected_response(data)*vector)[0]
    assert ordinary == independent_two_record_form(data,False) == -F(83481968,1858725825)
    assert reflected == independent_two_record_form(data,True) == F(21304304,56325025)
    for gap in (0,1,2):
        gram = t.reflected_response(data,gap)
        assert gram == gram.T
        assert (vector.T*gram*vector)[0] > 0


def test_positive_transfer_bound_uses_weighted_hilbert_norm(data):
    pi,Q = s.diag(*data['pi']),data['P'].T
    assert Q.T*pi == pi*Q
    assert Q*s.ones(66,1) == s.ones(66,1)
    minimum = min(Q[i,i] for i in range(66))
    assert minimum == F(791,990)
    assert 2*minimum-1 == t.LOWER
    f = data['D'][0,:].T+s.I*data['D'][1,:].T
    energy = s.expand((f.H*pi*(s.eye(66)-Q)*f)[0])
    norm = s.expand((f.H*pi*f)[0])
    assert 0 <= energy <= (1-t.LOWER)*norm
    assert t.logarithm_action(data,s.ones(66,1),8) == s.zeros(66,1)
    action = t.logarithm_action(data,f,4)
    assert s.expand((f.H*pi*action)[0]) > 0
    assert 0 < t.logarithm_tail(16) < F(1,10**7)


def test_actual_record_cylinders_generate_full_reflection_quotient(data):
    ranks = r.response_ranks(data)
    assert ranks['observable_rank_by_depth'][-1] == 66
    # Literal one-step score functions, no initial hidden-state measurement.
    conditional = [s.zeros(66,1) for _ in range(30)]
    for row in data['branches']:
        if row.law == 'E':
            e = r.local.SEAMS.index(row.seam)
            for i,_ in row.pairs:
                conditional[e][i] += row.weight*row.feature(row.seam)
    assert all(conditional[e] == data['D'][e,:].T for e in range(30))


def test_native_working_modular_clock_and_transfer_clock_are_inequivalent(data):
    value = t.clock_distinction(data)
    x,y = value['source_indices'];pi = s.diag(*data['pi'])
    assert (data['P']*pi-pi*data['P'])[y,x] == F(68,148599)
    assert value['transfer_H_commutes_with_native_modular_K'] is False
    assert value['transfer_H_is_scalar_rescaling_and_shift_of_native_modular_K'] is False


@pytest.mark.parametrize('value', (None,(),('repair','E',(0,5),2,0),('RESOURCE_STOP',None,None,None,None)))
def test_reject_outside_stationary_raw_label_domain(value):
    with pytest.raises(ValueError):
        t.reverse_label(value)


def test_transfer_certificate_does_not_promote_scalar_reconstruction(data):
    value = t.certificate(data)
    assert value['dimension'] == 66
    scope = value['scope']
    assert scope['same_stationary_scalar_complete_record_law']
    assert scope['raw_attempts_operational_only_at_integer_spacings']
    assert not scope['full_quantum_operator_history_positivity']
    assert not scope['unitary_analytic_continuation_is_an_executed_native_service']
    assert not scope['physical_clock_spatial_field_or_particle_masses_identified']
