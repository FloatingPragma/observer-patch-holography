"""Exact polynomial ratios, complete label refinement and physical multiplicity."""
from copy import deepcopy
from fractions import Fraction as F
from math import exp,factorial
import json

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_dyadic_suppression_uses_exact_outward_arithmetic_without_underflow():
    row=m.dyadic_exponential_upper(F(22271,2))
    assert row['certified_negative_power_of_two']==11135
    assert row['arithmetic_upper']==F(1,65536)
    assert not row['floating_exponential_evaluated']
    assert m.dyadic_exponential_upper(F(7,3),16)['arithmetic_upper']==F(1,4)
    assert m.dyadic_exponential_upper(0)['arithmetic_upper']==1
    with pytest.raises(ValueError):m.dyadic_exponential_upper(-1)
    assert m.paired_window_lower(F(3,5),F(1,4))==F(12,25)


def test_actual_polynomial_bracket_is_not_an_assumed_exact_Gibbs_identity():
    beta=F(2);lower=F(-1);upper=F(1);order=2;scale=4
    bound=m.modular.polynomial_slope(beta,lower,upper,order,scale)
    energies=(-1,0,1)
    def weight(E):
        y=beta*(upper-E)/(2*scale)
        return sum((y**k/F(factorial(k)) for k in range(order+1)),F())**(2*scale)
    weights=tuple(weight(E) for E in energies)
    operator=((0,F(1,3),0),(F(1,2),0,F(2,5)),(F(1,7),0,0))
    seen=0
    for i,Ei in enumerate(energies):
        for j,Ej in enumerate(energies):
            if Ej<=Ei or not operator[j][i]:continue
            plus=weights[i]*operator[j][i]**2
            reverse=weights[j]*operator[j][i]**2
            ratio=reverse/plus;omega=Ej-Ei;seen+=1
            assert exp(-float(bound['modular_slope_upper'])*omega)<=float(ratio)
            assert float(ratio)<=exp(-float(bound['modular_slope_lower'])*omega)
            assert abs(float(ratio)-exp(-float(beta)*omega))>1e-5
    assert seen==2 and bound['relative_modular_slope_defect_upper']>0


def test_actual_inversion_is_distinct_from_the_central_two_pi_spin_lift():
    control=m.inversion_control()
    assert control['proper_binary_cube_elements']==48
    assert control['full_rotation_inversion_group_order']==96
    assert control['beta_reverses_actual_Wilson_hop']
    assert not control['central_two_pi_CAR_parity_is_spatial_inversion']


def test_literal_fusion_tables_are_complete_and_flip_only_the_spin_center():
    c=m.fusion_control();assert len(c['complete_total_cube_irrep_labels'])==16
    for label,total in [('G1',28),('H',44)]:
        row=c['probe_tables'][label];assert row['allowed_initial_final_pairs']==total
        assert row['allowed_transition_multiplicities']['A1g']=={label+'g':1}
        assert row['allowed_transition_multiplicities']['A1u']=={label+'u':1}
        assert all(a[-1]==b[-1] for a,outputs in row['allowed_transition_multiplicities'].items() for b in outputs)
    assert c['probe_tables']['G1']['allowed_transition_multiplicities']['Hg']=={'Eg':1,'T1g':1,'T2g':1}
    assert c['probe_tables']['H']['allowed_transition_multiplicities']['Hg']['T1g']==2


def test_complete_initial_and_final_labels_preserve_energy_law_but_not_reference_state():
    # A stationary, label-invariant source with coherent purification.
    rho=np.diag([.4,.3,.2,.1]);H=np.diag([0,1,2,3]);projections=[np.diag([1,0,1,0]),np.diag([0,1,0,1])]
    A=np.array([[0,.2,0,0],[.3,0,0,0],[0,0,0,.4],[0,0,.5,0]],complex)
    failure=np.diag(np.sqrt(1-np.diag(A.conj().T@A).real));K=(A,failure)
    energy=[np.diag([int(i==j) for i in range(4)]) for j in range(4)]
    for initial in energy:
        for final in energy:
            for operator in K:
                plain=np.trace(final@operator@initial@rho@initial@operator.conj().T)
                labelled=sum(np.trace(beta@final@operator@alpha@initial@rho@initial@alpha@operator.conj().T@final@beta)
                             for alpha in projections for beta in projections)
                assert abs(plain-labelled)<1e-15
    psi=np.diag(np.sqrt(np.diag(rho))).reshape(-1)
    pure=np.outer(psi,psi.conj())
    dephased=sum(np.kron(P,np.eye(4))@pure@np.kron(P,np.eye(4)) for P in projections)
    assert np.linalg.norm(pure-dephased)>.1
    assert np.linalg.norm(np.einsum('irjr->ij',dephased.reshape(4,4,4,4))-rho)<1e-15
    # Source invariance is needed for an unchanged marginal.
    initial=np.ones((2,2))/2;probe=np.array([[1,1],[0,0]])/np.sqrt(2)
    pinched=np.diag(np.diag(initial))
    assert abs(np.trace(probe@initial@probe.T)-np.trace(probe@pinched@probe.T))>.4


def test_actual_480_mode_witnesses_have_same_labels_but_distinct_energy_means():
    control=m.physical_multiplicity_witness()
    assert control['actual_orthonormal_source_vectors']==12
    assert control['uniform_three_defect_sites']==8 and control['uniform_seven_defect_site_pairs']==56
    assert control['source_anchor_native_occupations']==(48,32,96,48,16)
    assert control['neutral_native_increment']==(0,0,0,0,0)
    assert control['neutral_operator_cube_irrep']=='A1g'
    assert control['seven_defect_source_energy_mean']-control['three_defect_source_energy_mean']==12
    assert not control['distinct_means_are_claimed_eigenvalues']
    assert not control['isolated_lowest_branch_is_ruled_out']
    for label,component,three,seven in m.source_multiplicity_vectors():
        assert label in ('G1','H') and not set(three)&set(seven)
        assert all(isinstance(mask,int) and mask.bit_length()<=480 for mask in (*three,*seven))
        assert s.simplify(sum(s.conjugate(v)*v for v in three.values()))==1
        assert s.simplify(sum(s.conjugate(v)*v for v in seven.values()))==1


def test_actual_first_residual_action_proves_original_measure_is_not_a_single_line():
    c=m.original_transfer_spread_control()
    assert c['canonical_six_component_residual_Gram_diagonal']==3
    assert c['surviving_positive_beta_hopping_image_Gram']==F(3,4)
    assert c['actual_each_component_residual_anchor_norm_squared']==F(9,4)
    assert c['actual_oriented_cube_residual_profiles']==24
    assert c['all_first_action_link_words_admitted_by_degree_one_cutoff']
    assert c['Gauss_projection_preserves_the_first_residual_action']
    assert c['every_component_original_transfer_measure_is_non_single_line']
    assert not c['nonzero_variance_excludes_an_isolated_pole']
    assert not c['actual_original_mean_transfer_computed']
    assert not c['current_clock_error_resolves_this_variance_floor']


@pytest.fixture(scope='module')
def actual():return m.bind_member(m.contract.member(0))


def test_same_original_source_has_stronger_complete_forward_windows(actual):
    r=actual['result'];assert F(r['reported_original_total_selected_phase_weight_lower'])>F(1,3)
    assert F(r['reported_original_total_complete_Walsh_row_weight_lower'])>F(1,24)
    assert r['paired_adjoint_ratio_upper']['certified_negative_power_of_two']==11135
    for channel,simple in zip(r['spin_resolved_records'],(F(5,48),F(11,48))):
        assert F(channel['reported_original_selected_phase_weight_lower'])>simple
        assert F(channel['reported_original_complete_Walsh_row_weight_lower'])>simple/8
        assert not channel['particular_robustly_populated_initial_final_pair_identified']
    assert not r['actual_polynomial_state_replaced_by_Gibbs']
    assert r['negative_probe_is_adjoint_not_source_C_conjugated_creator']
    assert r['physical_mass_table_GeV']==[]
    spread=r['original_transfer_spread']
    assert spread['physical_density_floor_factored']['negative_exponent']==582407133
    assert F(spread['normalized_original_transfer_variance_lower_factored']['multiplier_of_physical_density_floor'])==F(9,4)


def test_explicit_original_total_cube_branches_retain_tiny_weight_and_clock_scope(actual):
    r=actual['result']
    for branch,coefficient in zip(r['explicit_original_final_energy_branches'],(F(85,4096),F(85,2048))):
        factored=branch['complete_uniform_Walsh_original_weight_lower_factored']
        assert factored['physical_density_floor']['negative_exponent']==582407133
        assert F(factored['multiplier'])>coefficient
        assert branch['initial_total_cube_irrep']=='A1g'
        assert branch['final_native_occupations']==[48,34,96,48,17]
        assert F(branch['inferred_excitation_interval'][0])>0
        assert not branch['initial_energy_record_applied'] and not branch['energy_transfer_law_claimed']
        assert not branch['source_input_replaced_by_witness']
    assert r['original_energy_law_is_marginal_of_refined_complete_records']
    assert not r['initial_label_measurement_preserves_source_reference_joint_state']


def test_bound_rejects_source_changes_or_pole_promotion(actual):
    m.verify_bound(m.contract.member(0),actual)
    wrong=deepcopy(actual);wrong['result']['isolated_original_particle_branch_proved']=True
    with pytest.raises(ValueError):m.verify_bound(m.contract.member(0),wrong)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))


def test_certificate_is_canonical_and_pins_the_modular_and_original_source_proofs():
    c=m.certificate();assert json.loads(json.dumps(c,sort_keys=True))==c
    assert 'model/particle_primary_gravity_premises/PROOF.md' in c['source_sha256']
    assert 'model/particle_primary_coherent_odd_response/KERNEL_PROOF.md' in c['source_sha256']
    assert len(c['members'])==1
