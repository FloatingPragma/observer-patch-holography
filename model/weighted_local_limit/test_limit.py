"""Independent finite enumeration, nondivisor windows and large-capacity bounds."""
from collections import defaultdict
from fractions import Fraction as F
from math import comb, log

import pytest
import sympy as s

from . import limit as l


def test_fugacity_is_the_unique_positive_mean_root():
    z = (1+s.sqrt(6401))/1600
    assert s.simplify(800*z*z-z-2) == 0
    assert s.simplify(2*25*z/(1+25*z)+10*8*z/(1+8*z)-4) == 0
    assert (1-s.sqrt(6401))/1600 < 0
    bounds = l.parameters()
    for key, exact in [('fugacity', z), ('high',25*z/(1+25*z)), ('low',8*z/(1+8*z))]:
        lo, hi = bounds[key]
        assert s.Rational(lo.numerator,lo.denominator) < exact < s.Rational(hi.numerator,hi.denominator)


def test_native_weight_word_is_inherited_and_drift_rejected(monkeypatch):
    assert l.native_word() == (25,25)+(8,)*10
    class Changed:
        weights = (1,)*12
    monkeypatch.setattr(l.owning, 'local_family', lambda: Changed())
    with pytest.raises(ArithmeticError, match='weight word'):
        l.native_word()


@pytest.mark.parametrize('mean',[2,4])
@pytest.mark.parametrize('cells', [1, 2, 3, 5, 12, 20])
def test_conditioning_bound_has_exact_finite_binomial_witness(cells,mean):
    row = l.adjacent_count_control(cells,mean=mean)
    assert row['high_count']+row['low_count'] == mean*cells
    assert row['event_probability_lower'] > F(1,240*cells)
    condition = l.conditioning_interval(cells,mean=mean)
    assert condition[0] > row['event_probability_lower']


def test_wrong_total_does_not_inherit_the_mean_conditioning_bound():
    # The empty sector is legitimate, but exponentially unlikely in this product.
    assert l.conditioning_interval(20, occupation=0)[1] < F(1,240*20)


@pytest.mark.parametrize('high,low', [(0,0),(2,3),(2,10),(4,20)])
def test_type_coefficient_matches_existing_weighted_polynomial(high,low):
    weights = (25,)*high+(8,)*low
    for degree in range(-1,high+low+2):
        assert l.coefficient(high,low,degree) == l.owning.elementary(weights,degree)


def test_complete_one_cell_counts_match_independent_configuration_enumeration():
    canonical = l.owning.canonical(l.WORD,4)
    counts = defaultdict(F)
    for mask,probability in canonical.items():
        a=(mask&3).bit_count(); b=(mask>>2).bit_count()
        counts[a,b] += probability
    reported = l.canonical_counts(1,1)
    assert len(reported)==33
    assert all(reported[key]==counts[key] for key in reported)
    assert sum(reported.values())==1
    assert all(not value for (a,b),value in reported.items() if a+b!=4)


@pytest.mark.parametrize('cells,retained', [(1,0),(1,1),(2,1),(3,2),(5,2),(20,1)])
@pytest.mark.parametrize('mean',[2,4])
def test_full_state_distance_is_inside_the_all_domain_bound(cells,retained,mean):
    lower,upper=l.distance_interval(cells,retained,mean=mean)
    assert 0<=lower<=upper<=l.local_bound(cells,retained)['half_trace_distance_upper']
    assert upper-lower<F(1,10**40)
    assert sum(l.canonical_counts(cells,retained,mean=mean).values())==1


def test_full_window_distance_is_the_complement_conditioning_probability():
    low,high=l.distance_interval(2,2)
    qlo,qhi=l.conditioning_interval(2)
    assert low<=1-qlo and 1-qhi<=high
    assert max(high-(1-qhi),(1-qlo)-low)<F(1,10**40)


def test_cell_restriction_retains_every_number_sector_at_positive_density():
    values=l.canonical_counts(20,1)
    sectors=[sum((p for (a,b),p in values.items() if a+b==n),F()) for n in range(13)]
    assert all(p>0 for p in sectors)
    assert sum(sectors)==1
    # Replacing the full restriction by its two-load conditional state loses mass.
    assert 1-sectors[2]>F(1,2)


def _kl(cells,retained):
    canonical=l.canonical_counts(cells,retained)
    intervals=l.product_counts(retained)
    return sum(float(p)*log(float(p/((intervals[key][0]+intervals[key][1])/2)))
               for key,p in canonical.items() if p)


def test_entropy_bound_for_nondivisor_window_and_changed_marginal_mean():
    # This diagnostic checks the entropy argument at k=2, B=5, without a
    # disjoint-block divisibility assumption. The all-domain proof is in PROOF.md.
    local,global_=_kl(5,2),_kl(5,5)
    assert 0<local<(F(2,5)*global_)-0.01
    counts=l.canonical_counts(1,1)
    canonical_high=sum(F(a)*p for (a,_),p in counts.items())/2
    lo,hi=l.parameters()['high']
    assert canonical_high<lo or canonical_high>hi


def test_any_subset_of_ports_contracts_the_complete_count_distance():
    # One high-weight port is a stochastic projection of the full cell state.
    counts=l.canonical_counts(5,1)
    probability=sum(F(a,2)*p for (a,_),p in counts.items())
    lo,hi=l.parameters()['high']
    projected=max(abs(probability-lo),abs(probability-hi))
    assert projected<=l.distance_interval(5,1)[0]


def test_postselection_cannot_be_silently_added_to_the_complete_bound():
    epsilon=F(1,1000)
    p=(epsilon,F(),1-epsilon); q=(F(),epsilon,1-epsilon)
    assert sum(abs(a-b) for a,b in zip(p,q,strict=True))/2==epsilon
    # The successful event {0,1}, after renormalization, gives orthogonal laws.
    assert abs(p[0]/epsilon-q[0]/epsilon)==1


def test_bounded_gain_map_does_not_inherit_trace_preserving_contractivity():
    # Phi=G*identity is CP. Its two orthogonal-state outputs have half norm G.
    gain=F(257,256)
    assert l.instrument_error(F(1),gain=gain)==gain>1
    epsilon,eta=F(1,10**6),F(1,10**7)
    assert l.instrument_error(epsilon,gain=gain,implementation_error=eta)==gain*epsilon+eta
    assert l.instrument_error(epsilon)==epsilon
    with pytest.raises(ValueError):l.instrument_error(F(2),gain=gain)


def test_huge_levels_do_not_call_any_coefficient_or_graph_enumerator(monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError('global enumeration was attempted')
    monkeypatch.setattr(l,'coefficient',forbidden)
    monkeypatch.setattr(l.owning.geometry,'federation',forbidden)
    finite=l.level_bound(200)
    enormous=l.level_bound(1000)
    assert finite['ports']>10**122
    assert finite['half_trace_distance_upper']<F(2,10**60)
    assert enormous['half_trace_distance_upper']<F(1,10**298)
    assert enormous['half_trace_distance_upper']<finite['half_trace_distance_upper']
    assert l.two_size_bound(20*4**200,20*4**201)<F(3,10**60)


@pytest.mark.parametrize('value', [F(),F(1,4),F(2),F(1,10**1000)])
def test_square_root_certificate_is_outward_at_arbitrary_scale(value):
    upper=l.sqrt_upper(value)
    assert upper*upper>=value
    if value:
        assert upper*upper/value<F(1000001,1000000)


@pytest.mark.parametrize('cells,retained', [(True,1),(0,0),(-1,0),(2,True),(2,3),(2,-1),(F(2),1)])
def test_limit_domain_rejects_implicit_coercion_or_invalid_windows(cells,retained):
    with pytest.raises(ValueError):
        l.local_bound(cells,retained)


def test_empty_window_has_exactly_zero_error():
    assert l.local_bound(20,0)['half_trace_distance_upper']==0
    assert l.distance_interval(3,0)==(0,0)


def test_small_binomial_modes_and_neighbor_factors_independently():
    # Rational points exercise all integer and neighboring-mode edge cases.
    for n in range(1,12):
        for p in (F(1,4),F(1,2),F(3,5),F(3,4)):
            values=[F(comb(n,j))*p**j*(1-p)**(n-j) for j in range(n+1)]
            mode=max(values); mean=n*p
            below=mean.numerator//mean.denominator
            above=-(-mean.numerator//mean.denominator)
            assert values[below]>=n*(1-p)/(n*(1-p)+1)*mode
            assert values[above]>=n*p/(n*p+1)*mode
            assert mode*mode>=F(1,16*n)


def test_mean_two_root_and_exact_single_cell_H66_attachment():
    z=(s.sqrt(1256)-16)/1000
    assert s.simplify(1000*z*z+32*z-1)==0
    assert s.simplify(2*25*z/(1+25*z)+10*8*z/(1+8*z)-2)==0
    native=l.owning.local.reference(l.owning.local.Parameters(port_weights=l.WORD))
    states=l.owning.local.basis()
    expected=defaultdict(F)
    for state,probability in zip(states,native,strict=True):
        expected[sum(state[:2]),sum(state[2:])]+=probability
    reported=l.canonical_counts(1,1,mean=2)
    assert all(value==expected[key] for key,value in reported.items())
    # Its growing-family restriction has additional sectors despite mean two.
    grown=l.canonical_counts(20,1,mean=2)
    assert sum(p for (a,b),p in grown.items() if a+b!=2)>F(1,2)


def test_declared_means_are_distinct_and_neither_is_selected():
    assert l.parameters(mean=2)['high'][1]<l.parameters(mean=4)['high'][0]
    assert l.local_bound(20,mean=2)['total_occupation']==40
    assert l.local_bound(20,mean=4)['total_occupation']==80
    assert l.local_bound(20,mean=2)['half_trace_distance_upper']==l.local_bound(20,mean=4)['half_trace_distance_upper']
    for value in (True,2.0,3,0,F(2)):
        with pytest.raises(ValueError):l.parameters(mean=value)
        with pytest.raises(ValueError):l.local_bound(20,mean=value)
