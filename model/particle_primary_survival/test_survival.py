"""Primary native contractions and independent inference countercontrols."""
from copy import deepcopy
from fractions import Fraction as F
from math import sin

import numpy as np
import pytest
import sympy as s

from . import survival as v


@pytest.fixture(scope='module')
def member():
    return v.contract.member(1)


@pytest.fixture(scope='module')
def witness(member):
    return v.selected_witness(member)


def test_actual_native_metric_and_all_dirac_modes_fix_the_loop_constant():
    value=v.native_witness_control()
    assert value['native_generator_count']==12
    assert value['Dirac_CAR_modes_per_cell']==60
    assert {r['integer_Q6_charge']:F(r['trace_metric_Casimir']) for r in value['all_native_blocks']} == {
        -4:F(4,5),-3:F(9,20),1:F(21,20),2:F(7,10),6:F(3,10)}
    assert F(value['sum_squared_block_Casimirs'])==F(101,40)
    assert value['trial_occupied_modes_per_cell']==30
    assert F(value['trial_onsite_energy_times_spacing'])==-90
    assert not value['loop_character_span_claimed_invariant_under_H']


def test_compressed_multiplications_commute_on_retained_blank_not_generally():
    # Independent finite two-coordinate Fourier multiplication control.
    shift=s.zeros(5)
    for j in range(4):shift[j+1,j]=1
    cosine=(shift+shift.T)/2
    a=s.kronecker_product(cosine,s.eye(5));b=s.kronecker_product(s.eye(5),cosine)
    assert a*b==b*a
    indices=[5*(x+2)+(y+2) for x in range(-2,3) for y in range(-2,3) if abs(x)+abs(y)<=1]
    p=s.diag(*[int(i in indices) for i in range(25)])
    blank=s.eye(25)[:,12]
    assert p*a*blank==a*blank and p*b*blank==b*blank
    compressed=p*a*p*b*p-p*b*p*a*p
    assert compressed*blank==s.zeros(25,1)
    assert compressed!=s.zeros(25)


def test_full_fock_extremes_and_actual_negative_beta_slater_trial():
    # All sixty modes per cell are present in this two-site one-body identity.
    beta=v.density.gamma_data()[3];gamma=v.density.gamma_data()[0]
    onsite=s.kronecker_product(3*beta,s.eye(15))
    hop=s.kronecker_product((-beta+beta*gamma)/2,s.eye(15))
    h=onsite.row_join(hop).col_join(hop.H.row_join(onsite))
    assert h==h.H and s.trace(h)==0
    occupied=s.kronecker_product((s.eye(4)-beta)/2,s.eye(15))
    slater=s.diag(occupied,occupied)
    assert slater*slater==slater and s.trace(slater)==60
    assert s.trace(h*slater)==-180
    # Exterior CAR action of each off-diagonal source monomial on N=0,120.
    full=(1<<120)-1
    def apply_number_monomial(mask,i,j):
        if not mask&(1<<j):return None
        reduced=mask^(1<<j)
        if reduced&(1<<i):return None
        return reduced|(1<<i)
    for i,j in h.todok():
        assert apply_number_monomial(0,i,j) is None
        assert apply_number_monomial(full,i,j)==(full if i==j else None)
    # A nontrivial actual native group element has determinant one.
    link=v.density.native_weak_link(1,s.Rational(3,5),s.Rational(4,5))
    assert link.det()==1


def test_actual_polynomial_form_gives_positive_negative_pair_dominance():
    energies=[F(-2),F(0),F(3)]
    p=v.passive_polynomial_probabilities(energies)
    assert sum(p)==1 and p[0]>p[1]>p[2]>0
    a=s.Matrix([[0,s.Rational(1,2),s.Rational(1,3)],
                [s.Rational(1,2),0,s.Rational(1,4)],
                [s.Rational(1,3),s.Rational(1,4),0]])
    h=s.diag(*map(s.Rational,energies));rho=s.diag(*map(s.Rational,p))
    comm=h*a-a*h
    moment=sum(s.Rational(p[i])*abs(a[j,i])**2*(energies[j]-energies[i])**2
               for i in range(3) for j in range(3))
    assert s.trace(rho*comm.H*comm)==moment
    for i in range(3):
        for j in range(i+1,3):
            assert p[i]*abs(a[j,i])**2 >= p[j]*abs(a[i,j])**2
    inverted=tuple(reversed(p))
    assert inverted[0]*abs(a[1,0])**2 < inverted[1]*abs(a[0,1])**2
    # Independent spectral integration checks the passivity band inequality.
    atoms=[]
    for i in range(3):
        for j in range(3):
            energy=float(energies[j]-energies[i])
            multiplier=1 if energy==0 else (sin(energy/10)/(energy/10))**4
            atoms.append((energy,float(p[i]*abs(a[j,i])**2)*multiplier))
    second=sum(e*e*w for e,w in atoms)
    lower=max(0,(second-.5**2)/(2*5**2))
    assert 0<lower<=sum(w for e,w in atoms if .5<=e<=5)


def test_primary_contrast_has_its_actual_level_dependent_domain(member,witness):
    assert witness['probe']['active'] and witness['probe']['kind']=='plaquette'
    assert witness['clock_id']==member['ids']['clock_id']
    assert witness['reference_id']==member['ids']['reference_id']
    assert F(witness['blank_commutator_norm_squared'])==F(4949,20)/F(member['source']['spacing'])**2
    assert F(witness['filter_half_width'])*F(witness['source_spectral_width_upper'])==1
    assert F(witness['continuous_original_contrast_lower']['coefficient'])>0
    assert F(witness['original_positive_band_weight_lower']['coefficient'])>0
    actual=F(witness['finite_apparatus_original_positive_band_weight_lower']['coefficient'])
    ideal=F(witness['original_positive_band_weight_lower']['coefficient'])
    contrast=F(witness['continuous_original_contrast_lower']['coefficient'])
    assert actual==ideal-contrast/4>0
    assert witness['actual_primary_contrast_lower_bound_derived']
    assert not witness['fixed_positive_energy_band_survival_proved']
    assert not witness['postselection_normalization_applied']
    assert witness['physical_mass_or_GeV_assignment'] is None
    damaged=deepcopy(member);damaged['clock']['GeV_or_proper_time_conversion']='fitted clock'
    with pytest.raises(ValueError):v.selected_witness(damaged)
    with pytest.raises(ValueError):v.selected_witness(v.contract.member(0))


def test_extreme_sector_bound_uses_complete_weight_and_actual_plaquette_count(member):
    result=v.extreme_sector_bound(member)
    volume=member['source']['side']**3
    a=F(member['source']['spacing']);eps=F(member['state']['target_ground_excess'])
    assert F(result['global_ground_energy_upper'])==-90*volume/a
    assert F(result['combined_original_sector_probability_upper'])*F(result['extreme_sector_ground_relative_energy_lower'])==eps
    assert F(result['combined_original_sector_probability_upper'])<F(7,627*2**15)
    assert not result['conditional_sector_spectrum_transferred_to_primary']


def test_factored_counts_close_inherited_complete_response_errors_without_allocation():
    for coefficient in (F(1,7),F(101,13)):
        for exponent in (4,25):
            plan=v.factored_apparatus(coefficient,exponent)
            cells=2**plan['quadrature_cells']['exponent']
            steps=2**plan['Cayley_subdivisions']['exponent']
            exact=v.apparatus.error_budget(-2,3,F(1,5),cells,F(1,5),steps,'0')
            assert exact['original_contrast_error_to_continuous_triangle'] <= coefficient/F(2**(exponent+1))
    huge=v.factored_apparatus(F(1,17),10**12)
    assert huge['quadrature_cells']['exponent']>10**12
    assert not huge['counts_or_spatial_matrices_materialized']
    assert v.dyadic_leq(v.dyadic(1,10**12),v.dyadic(1,10**12-1))
    assert not v.dyadic_leq(v.dyadic(3,10**12),v.dyadic(1,10**12))


def test_transfer_band_stability_is_not_a_final_energy_projector():
    energies=np.array([0.,1.,3.]);probabilities=np.array([.6,.3,.1])
    b=np.array([[.25,1/3,0],[1/3,0,.25],[0,.25,-.25]])
    changed=b+np.array([[0,.01,0],[.01,0,0],[0,0,0]])
    assert max(np.linalg.norm(b,2),np.linalg.norm(changed,2))<1
    transfer=energies[:,None]-energies[None,:]
    mask=(transfer>=1)&(transfer<=2)
    def weight(operator):return ((operator*np.sqrt(probabilities))**2)[mask].sum()
    original=weight(b);perturbed=weight(changed)
    assert 0<abs(original-perturbed)<=2*np.linalg.norm(changed-b,2)
    final_energy=((b*np.sqrt(probabilities))**2)[1,:].sum()
    assert abs(final_energy-original)>.01


def test_cooling_faithfulness_and_shrinking_filter_do_not_prevent_uv_escape():
    # Logical counterfamily only. It is not the primary spatial Hamiltonian.
    for n in (0,2,5,10):
        energy=F(16*2**n);excited=F(1,2**(4*n+12));eps=F(1,2**(n+4))
        assert 0<excited<1-excited and energy*excited<=eps
        weights=[(-energy,excited),(energy,1-excited)]
        assert sum(weight*frequency**2 for frequency,weight in weights)==energy**2
        mass=sum(float(weight)*(sin(float(frequency)/2)/(float(frequency)/2))**4
                 for frequency,weight in weights)
        assert mass<=float(16/energy**4)
        shrinking_mass=sum(float(weight)*(sin(.5)/.5)**4 for _,weight in weights)
        assert shrinking_mass>.84
        assert sum(weight for frequency,weight in weights if 1<=frequency<=8)==0
        conditioned=[float(weight)*(sin(float(frequency)/2)/(float(frequency)/2))**4/mass
                     for frequency,weight in weights]
        assert mass>0 and abs(sum(conditioned)-1)<1e-14
        assert sum(w for (frequency,_),w in zip(weights,conditioned) if 1<=frequency<=8)==0


def test_declared_sources_exist_and_dyadic_types_are_exact():
    assert all((v.contract.ROOT/path).is_file() for path in v.SOURCE_INPUTS)
    for coefficient,exponent in [(True,1),(0.1,1),(1,True),(1,-1)]:
        with pytest.raises(ValueError):v.dyadic(coefficient,exponent)
