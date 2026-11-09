from fractions import Fraction as F
from pathlib import Path
from itertools import permutations

import pytest
import sympy as s

from . import readout as m


def creator(n, j):
    out=s.MutableSparseMatrix.zeros(1<<n)
    for mask in range(1<<n):
        if not (mask>>j)&1:
            out[mask|(1<<j),mask]=(-1)**((mask&((1<<j)-1)).bit_count())
    return s.ImmutableSparseMatrix(out)


def test_full_exterior_translation_retains_CAR_signs_and_group_law():
    for p in permutations(range(3)):
        t=m.signed_fock_permutation(p)
        assert t.H*t==s.eye(8)
        for j in range(3):assert t*creator(3,j)*t.H==creator(3,p[j])
    t=m.signed_fock_permutation((1,2,0))
    assert t**3==s.eye(8)
    unsigned=abs(m.signed_fock_permutation((1,0,2)))
    assert unsigned*creator(3,0)*unsigned.H != creator(3,1)


def test_position_effects_are_complete_on_vacuum_and_every_number_sector():
    e=m.position_effects(3,2)
    assert sum(e,s.zeros(64))==s.eye(64)
    assert all(x>=0 for effect in e for x in effect.diagonal())
    assert sum(e[:-1],s.zeros(64))!=s.eye(64)
    assert e[-1].trace()==1
    # It is an unsharp random occupied-mode detector, not a site PVM.
    assert e[0]**2 != e[0]


def test_cyclic_boundary_wrap_and_number_dependent_Weyl_relation():
    n=4;t=m.signed_fock_permutation((1,2,3,0))
    e=m.position_effects(n,1); q=s.eye(16)-e[-1]
    x=sum((j*e[j] for j in range(n)),s.zeros(16))
    assert t.H*x*t==x+q-n*e[-2]
    number=s.diag(*[mask.bit_count() for mask in range(16)])
    dipoles=[sum(j for j in range(n) if mask>>j&1) for mask in range(16)]
    z=s.diag(*[s.I**j for j in dipoles])
    twist=s.diag(*[s.I**number[j,j] for j in range(16)])
    assert t.H*z*t==twist*z
    assert t.H*z*t!=s.I*z  # N=1 phase cannot be used on the whole Fock space.


def test_full_Fourier_projectors_zero_labels_and_actual_log_branch():
    t=m.signed_fock_permutation((1,2,3,0)); ps=m.fourier_projectors(t,4)
    assert sum(ps,s.zeros(16))==s.eye(16)
    for k,p in enumerate(ps):
        assert p.H==p and p*p==p
        assert t*p==s.I**k*p
        for j,q in enumerate(ps):
            if j!=k:assert p*q==s.zeros(16)
    # A shorter actual period keeps the unused Fourier labels at zero.
    short=m.fourier_projectors(s.eye(2),4)
    assert short==(s.eye(2),s.zeros(2),s.zeros(2),s.zeros(2))


def test_reference_dressing_preserves_joint_Gauss_domain_not_ambient_commutation():
    c=m.small_reference_control();t=c['joint_translation'];gs=c['joint_Gauss_generators']
    assert t*gs[0]*t.H==gs[1] and t*gs[1]*t.H==gs[0]
    physical=(s.eye(64)+gs[0])*(s.eye(64)+gs[1])/4
    assert t*physical==physical*t
    assert t*gs[0]!=gs[0]*t
    assert c['source_only_translation_Gauss_residual_squared']==2
    assert t*c['canonical_state']==c['canonical_state']


def test_complete_position_then_momentum_matches_the_binomial_formula():
    c=m.small_reference_control(); C=m.initial_overlap(2,1);A=s.Rational(3,8)
    assert C==s.Rational(1,8)
    rows={(x['position'],x['momentum_index']):x['original_probability']
          for x in c['joint_position_momentum_records']}
    assert rows['site0',0]==rows['site1',0]==(A+C)/2
    assert rows['site0',1]==rows['site1',1]==(A-C)/2
    assert rows['vacuum',0]==s.Rational(1,4) and rows['vacuum',1]==0
    assert sum(rows.values())==1


@pytest.mark.parametrize('v,d',[(2,2),(3,1),(4,2)])
def test_binomial_overlap_and_strict_momentum_weight_from_all_Fock_inputs(v,d):
    e=m.position_effects(v,d)
    direct=sum(s.sqrt(e[0][i,i]*e[1][i,i]) for i in range(1<<(v*d)))/2**(v*d)
    C=m.initial_overlap(v,d)
    A=(1-s.Rational(1,2**(v*d)))/v
    assert s.simplify(C-direct)==0
    assert bool(s.simplify(A-C-s.Rational(d,2**(v*d)))>=0)
    assert s.simplify(v*((A+(v-1)*C)/v+(v-1)*(A-C)/v)+s.Rational(1,2**(v*d)))==1


def test_actual_Robertson_input_and_unsharp_noise_without_finite_CCR():
    c=m.small_reference_control()
    assert c['X_variance']==s.Rational(11,64)
    assert c['record_coordinate_variance']==s.Rational(15,64)
    assert c['unsharp_coordinate_noise']==s.Rational(1,16)
    assert c['P_variance']==3*s.pi**2/16
    assert c['commutator_expectation']==s.I*s.pi/4
    assert s.simplify(c['Robertson_lhs']-c['Robertson_rhs'])>0
    t=c['joint_translation']; p=s.pi*(s.eye(64)-t)/2
    x=s.diag(*range(64))
    assert s.trace(x*p-p*x)==0
    assert s.trace(s.I*s.eye(64))!=0


def test_native_link_and_matter_translation_on_every_one_body_column():
    c=m.native_covariance_control()
    assert c['one_body_dimension']==480 and c['native_internal_rows']==15
    assert len(c['all_axes'])==3
    assert all(row['complete_columns']==480 and row['joint_covariance_residual_zero']
               and row['matter_only_translation_with_fixed_links_fails'] for row in c['all_axes'])


def test_native_law_and_source_domains_keep_all_records_and_references():
    c=m.certificate(); law=c['native_complete_law']
    assert law['cells']==8 and law['modes_per_cell']==60
    assert law['all_joint_labels']==72
    assert law['vacuum_probability']['inverse_exponent']==480
    assert law['A_minus_C_positive_lower']['numerator']==60
    assert not law['full_ambient_gauge_commutation_claimed']
    assert not c['joint_crystal_momentum_is_physical_particle_momentum']
    root=Path(__file__).resolve().parents[2]
    assert all((root/p).is_file() for p in c['source_inputs'])


def test_invalid_or_intractable_explicit_inputs_are_not_silently_truncated():
    with pytest.raises(ValueError):m.position_effects(8,60)
    with pytest.raises(ValueError):m.signed_fock_permutation((0,0))
    with pytest.raises(ValueError):m.native_law(0)
    with pytest.raises(ValueError):m.fourier_projectors(s.Matrix([[0,1],[1,0]]),3)
    with pytest.raises(ValueError):m.initial_overlap(8,60)
