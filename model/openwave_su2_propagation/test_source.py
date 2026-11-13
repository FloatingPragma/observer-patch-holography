from fractions import Fraction as F
import pytest
from .source import plaquettes,graph_counts,exponential_tail_upper,commutator_bound,boundary_bound,run
from .source import compact_group_gap_map,su_fundamental_gap

def test_compact_group_stability_normalization_and_SU_minimum():
    from model.openwave_su2_volume_gap.source import premise_map
    actual=compact_group_gap_map(F(3,4))
    parent=premise_map(F(32))
    assert actual['local_perturbation_norm_upper_exact']==parent['local_perturbation_operator_norm_exact']
    assert actual['theorem_to_source_energy_scale_exact']==parent['theorem_to_source_energy_scale_exact']
    for n in range(2,33):
        assert su_fundamental_gap(n)==F(n*n-1,2*n)
        assert all(F(k*(n-k)*(n+1),2*n)>=su_fundamental_gap(n) for k in range(1,n))


def test_incidence_counts_and_actual_chain_enumeration():
    faces=plaquettes(3)
    for side in (3,4,5):
        r=graph_counts(side)
        assert r['all_links']==r['all_plaquettes']==3*side**3
        assert r['max_faces_per_link']==4
    x=next(iter(faces[0]));chains=[p for p in faces if x in p]
    assert len(chains)==4
    count=sum(len([q for q in faces if p&q]) for p in chains)
    assert count<=4*16
    with pytest.raises(ValueError):plaquettes(2)


def test_tail_rational_bound_against_high_precision():
    import mpmath as mp
    with mp.workdps(100):
        for x in (F(1,2),F(2),F(8)):
            xx=mp.mpf(x.numerator)/x.denominator
            for d in (0,8,16,32):
                bound=exponential_tail_upper(x,d)
                actual=mp.fsum(xx**n/mp.factorial(n) for n in range(d,400))
                assert actual<=mp.mpf(bound.numerator)/bound.denominator
                assert abs(actual-mp.mpf(bound.numerator)/bound.denominator)<mp.mpf('1e-60')


def test_volume_error_decreases_without_changing_flux_domain():
    bounds=[boundary_bound(32,1,1,r,4) for r in (8,16,24,32)]
    assert all(a>b for a,b in zip(bounds,bounds[1:]))
    assert bounds[-1]<F(1,10**18)
    assert commutator_bound(32,1,0,8,4)==0
    r=run();assert r['proved']['no_Peter_Weyl_cutoff']
    assert not any(r['scope'].values())
    assert all(not x['numeric_coupling_gap_admitted'] for x in r['full_flux_propagation_bounds'])
