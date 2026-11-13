from fractions import Fraction as F
import pytest
from .source import solve,pi_interval,run


def test_length_elimination_for_independent_exact_cases():
    for g in (F(1,7),F(37)):
        for a in (F(1,13),F(81)):
            r=solve(g,F(3),a,F(11),F(1,19))
            assert F(r['strength'])/F(r['speed'])==g/9
            assert not r['length_alone_can_match_both']
    assert solve(F(9,19),3,4,11,F(1,19))['length_alone_can_match_both']
    with pytest.raises(ValueError):solve(1,0,1,1,1)


def test_machin_enclosure_independent_high_precision():
    import mpmath as mp
    with mp.workdps(120):
        lo,hi=pi_interval()
        m=lambda x:mp.mpf(x.numerator)/x.denominator
        assert m(lo)<mp.pi<m(hi)
        assert hi-lo<F(1,10**80)


def test_both_unit_equations_after_explicit_solution():
    r=run();x=r['fixed_source'];u=F(r['external_reference_inputs']['hbar_c_GeV_fm'])/F(x['light_matched_length_fm'])
    g=F(x['invariant'])*u*u
    hc=F(r['external_reference_inputs']['hbar_c_GeV_fm']);target=F(r['external_reference_inputs']['Newton_G_natural_GeV_minus2'])
    assert u*F(x['light_matched_length_fm'])/hc==1
    assert g*F(x['gravity_matched_length_fm'])/(u*hc)==target
    assert F(x['required_stiffness_ratio'])>10**42
    assert not r['source_stiffness_calibrated_to_G']
