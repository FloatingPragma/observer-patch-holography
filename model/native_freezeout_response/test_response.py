from fractions import Fraction as F

import pytest

from .response import (hitting_data,terminal_law,limiting_state,statistics,reset,
                       boundary_payoff,universal_count_bound,tv,trajectory)
from model.native_abundance_epochs import model as native


def test_all_hitting_equations_and_exact_independent_native_aggregate():
    data=hitting_data();old=native.persistent_probability()
    for s in native.STATES:
        occupied=sum(v for part,v in zip(data['classes'],data['absorption'][s]) if part[0][3])
        assert occupied==old[s]
    p=native.initial_distribution()
    assert sum(p[s]*data['mean'][s] for s in p)==F(5557709,1125000)
    assert max(data['mean'].values())==15
    expected=[F(12817,13824),F(499,27648),F(25,1728),F(1,128),F(25,1728),F(499,27648)]
    assert statistics(limiting_state(p))['oriented_face_class_shares']==expected


def test_finite_class_distance_is_bounded_without_postselection():
    p=native.initial_distribution();limit=statistics(limiting_state(p))['oriented_face_class_shares']
    for n in range(17):
        obs=statistics(p)
        distance=sum(abs(a-b) for a,b in zip(obs['oriented_face_class_shares'],limit))/2
        assert distance<=obs['unresolved_probability']
        assert sum(p.values())==1
        p=native.evolve(p)


def test_asymmetric_input_response_and_boundary_payoff_are_exact():
    p=native.initial_distribution();q={(0,1,5,1):F(1)};eps=F(1,10)
    mixed=reset(p,next(iter(q)),eps)
    assert terminal_law(mixed)==[(1-eps)*a+eps*b for a,b in zip(terminal_law(p),terminal_law(q))]
    assert statistics(limiting_state(q))['oriented_1_minus_5']==-F(91,384)
    boundary=[F(i%3,2) for i in range(46)];h=boundary_payoff(boundary)
    assert all(0<=v<=1 for v in h.values())
    for s,row in native.kernel().items():assert h[s]==sum(v*h[t] for t,v in row.items())
    # Distinct frozen classes attain the universal contraction coefficient one.
    a,b={(2,2,2,1):F(1)},{(2,2,2,5):F(1)}
    assert tv(a,b)==1 and sum(abs(x-y) for x,y in zip(terminal_law(a),terminal_law(b)))/2==1


def test_count_budget_and_boundary_domain_controls():
    assert universal_count_bound(F(1,1000))['native_attempts']==96
    assert universal_count_bound(F(1,1000))['unresolved_bound']==F(1,15625)
    with pytest.raises(ValueError):universal_count_bound(0)
    with pytest.raises(ValueError):universal_count_bound(0.001)
    with pytest.raises(ValueError):terminal_law({(0,0,0,0):F(1,2)})
    with pytest.raises(ValueError):reset(native.initial_distribution(),(0,1,5,1),F(11,10))


def test_original_reset_bit_joint_law_is_retained_at_every_post_reset_cut():
    case=trajectory({(0,1,5,0):F(1)},boundary_reset=(8,(0,1,5,1),F(1,10)))
    for cut in case['cuts']:
        branches=cut['complete_boundary_record_state_laws']
        if cut['attempts']<8:
            assert branches==[]
            continue
        assert [b['record'] for b in branches]==['no_reset','reset']
        assert [b['original_probability'] for b in branches]==[F(9,10),F(1,10)]
        assert all(sum(v)==p for v,p in zip(zip(*(b['unnormalized_state_law'] for b in branches)),cut['complete_state_law']))
