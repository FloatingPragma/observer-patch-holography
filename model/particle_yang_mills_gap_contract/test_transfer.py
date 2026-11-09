"""Exact transfer matrices distinguish spectral decay, clocks and probe coverage."""
from fractions import Fraction as F
import json
import pytest
import sympy as s
from . import transfer as g


def two_state(ratio):
    omega=s.Matrix([3,4])/5; p=omega*omega.T
    return p+ratio*(s.eye(2)-p),omega,s.diag(0,1)


def test_nondiagonal_transfer_keeps_actual_variance_at_every_time():
    t,omega,a=two_state(s.Rational(1,2))
    receipt=g.finite_transfer(t,omega,[a])
    assert receipt['full_excitation_span']
    assert receipt['probes'][0]['variance']=='144/625'
    for k in range(9):
        assert g.correlation(t,omega,a,k)==s.Rational(144,625)*s.Rational(1,2)**k


def test_incomplete_fast_probe_does_not_certify_global_gap_of_coupled_source():
    lap=s.Matrix([[1,-1,0],[-1,2,-1],[0,-1,1]])
    t=s.eye(3)-lap/4; omega=s.ones(3,1)/s.sqrt(3)
    fast=s.diag(1,-2,1);slow=s.diag(-1,0,1)
    result=g.finite_transfer(t,omega,[fast],ratio='1/4')
    assert not result['full_excitation_span']
    assert g.correlation(t,omega,fast,5)==2*s.Rational(1,4)**5
    assert g.correlation(t,omega,slow,5)==s.Rational(2,3)*s.Rational(3,4)**5
    with pytest.raises(ValueError,match='retained spectral weight'):
        g.finite_transfer(t,omega,[fast,slow],ratio='1/4')
    assert g.finite_transfer(t,omega,[fast,slow],ratio='3/4')['full_excitation_span']


def test_positive_finite_gaps_need_not_be_uniform():
    for n in (1,2,8,32):
        ratio=1-s.Rational(1,n+2); t,omega,a=two_state(ratio)
        assert g.correlation(t,omega,a,0)==s.Rational(144,625)
        with pytest.raises(ValueError,match='retained spectral weight'):
            g.finite_transfer(t,omega,[a],ratio='1/2')
    assert s.limit(1-1/(s.Symbol('n',positive=True)+2),s.Symbol('n',positive=True),s.oo)==1


def test_relabeling_source_ticks_changes_the_physical_decay():
    # Same source eigenvalue 1/2. One physical unit equals half a source unit.
    t,omega,a=two_state(1/s.sqrt(2))
    assert g.correlation(t,omega,a,2)==s.Rational(72,625)
    with pytest.raises(ValueError,match='retained spectral weight'):
        g.finite_transfer(t,omega,[a],ratio='1/2',time_step='1')


def test_divergent_gap_has_no_strongly_continuous_nontrivial_limit():
    variance=s.Rational(144,625)
    for n in (2,4,8):
        t,omega,a=two_state(s.Rational(1,2)**n)
        assert g.correlation(t,omega,a,1)==variance/s.Integer(2)**n
        assert g.correlation(t,omega,a,0)==variance
    # The limiting correlator is zero at every positive integer step, while
    # the continuous-time form exp(-n log(2)t) has a jump at zero. Its
    # original vector energy is n log(2) Var(A), so no finite uniform M exists.
    n=s.Symbol('n',positive=True);time=s.Symbol('time',positive=True)
    assert s.limit(s.exp(-n*s.log(2)*time),n,s.oo)==0
    assert s.limit(n*s.log(2)*variance,n,s.oo)==s.oo


def test_vanishing_original_residue_can_hide_a_closing_prelimit_gap():
    # Every prelimit has a slow mode; its original residue tends to zero.
    # Fixed-time errors must vanish before the long-time limit is taken.
    for n in (2,4,8):
        epsilon=s.Rational(1,n);qslow=1-s.Rational(1,n+1)
        t=s.diag(1,qslow,s.Rational(1,2));omega=s.Matrix([1,0,0])
        v=s.Matrix([0,s.sqrt(epsilon),s.sqrt(1-epsilon)])
        a=v*omega.T+omega*v.T
        for k in range(8):
            actual=g.correlation(t,omega,a,k)
            assert actual==epsilon*qslow**k+(1-epsilon)*s.Rational(1,2)**k
            assert actual<=s.Rational(1,2)**k+epsilon
        with pytest.raises(ValueError,match='retained spectral weight'):
            g.finite_transfer(t,omega,[a],ratio='1/2')
        bound=g.low_energy_weight_bound(1,'1/2',str(qslow),5,str(epsilon))
        assert F(epsilon)<=bound


def test_noncollapse_rejects_incompatible_physical_gap_and_energy():
    with pytest.raises(ValueError,match='contradicts'):
        g.lattice_gap_window('1/8','2','1','1')


def test_nonpositive_transfer_cannot_define_the_declared_finite_hamiltonian():
    t,omega,a=two_state(0)
    with pytest.raises(ValueError,match='spectrum'):
        g.finite_transfer(t,omega,[a])


def test_stationary_low_energy_density_is_not_a_ground_correlator():
    h=s.diag(0,1);a=s.Matrix([[0,1],[0,0]])
    for n in (2,8,32):
        epsilon=s.Rational(1,n);rho=s.diag(1-epsilon,epsilon)
        assert rho*h==h*rho
        assert s.trace(rho*h)==epsilon
        for k in range(5):
            # Exact imaginary times t=k log(2), with every thermal weight.
            decay=s.diag(1,s.Rational(1,2)**k)
            growth=s.diag(1,s.Integer(2)**k)
            value=s.trace(rho*a.H*decay*a*growth)
            assert value==epsilon*s.Integer(2)**k


def test_retained_native_diagnostics_use_the_declared_clock_not_a_mass_fit():
    record=g.native_clock_diagnostics();rows=record['rows']
    assert [row['ports'] for row in rows]==[240,960,3840]
    assert all(row['exact_normalization_identity'] for row in rows)
    for before,after in zip(rows,rows[1:]):
        ratio=F(before['per_attempt_gap_diagnostic'])/F(after['per_attempt_gap_diagnostic'])
        assert 15<ratio<18
        assert after['attempts_per_comparison_time']==16*before['attempts_per_comparison_time']
    assert not record['certified_eigenvalue_enclosures']
    assert not record['same_source_as_selected_Wilson_Gauss_program']


def test_actual_primary_action_rejects_a_term_without_the_common_clock_scale():
    retained=json.loads((g.ROOT/'model/particle_common_contract/contract.json').read_text())
    source=retained['fixed_primary_member']['source']
    g.verify_primary_scaling_source(source)
    changed={**source,'onsite_one_body':'3 beta_D'}
    with pytest.raises(ValueError,match='scaling theorem premises'):
        g.verify_primary_scaling_source(changed)
