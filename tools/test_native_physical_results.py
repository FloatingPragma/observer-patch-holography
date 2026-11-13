from fractions import Fraction as F
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import json
import pytest
from tools import native_physical_results as summary

ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def authenticated_result():
    # Authenticate the immutable receipt DAG once, then isolate every projection
    # mutation with a deep copy. Separate custody mutations still call collect.
    return summary.collect(ROOT)


def test_projection_preserves_calibrations_failed_comparisons_and_domains(authenticated_result):
    result=deepcopy(authenticated_result)
    assert len(result['receipts_sha256'])==len(summary.PACKAGES)
    assert not result['all_sources_identified_as_one_physical_theory']
    assert {r['name'] for r in result['strange_masses'] if r['comparison_role']=='calibration'}=={'Lambda','K+'}
    assert not result['nuclear_masses']['helion']['matches']
    assert not result['nuclear_masses']['alpha']['matches']
    assert not result['independent_unit_attachment']['physical_G_strength_match']
    assert all(r['rows']==15 and r['fitted_chronometer_parameters']==0
               for r in result['conditional_expansion_comparisons'])
    text='\n'.join(summary.render(result))
    assert 'Both unit comparisons fail' in text and 'virtual response' in text
    assert 'not independent likelihoods' in text
    # Measured input ancestry is projected from its source, not a second literal.
    result['conditional_expansion_comparisons'][0]['anchors']['H0_km_s_Mpc']='68.125'
    result['conditional_expansion_comparisons'][0]['anchors']['Omega_m']='0.3125'
    changed='\n'.join(summary.render(result))
    assert 'H0=68.125 and Omega_m=0.3125' in changed
    assert 'H0=67.4' not in changed.split('### Resolved charge profiles')[0]
    for name in ('neutron','proton'):
        line=next(x for x in text.splitlines() if x.startswith(f'| {name} band |'))
        display=list(map(F,line.split('[',1)[1].split(']',1)[0].split(', ')))
        source=list(map(F,result['bulk_baryons'][name+'_GeV_interval']))
        assert display[0]<=source[0]<=source[1]<=display[1]
    gravity=result['transverse_photon_and_clocks']
    for label,centre_key,error_key in (
        ('mean arrival parameter','conditional_arrival_delay','arrival_delay_error_upper'),
        ('transverse position / cell','transverse_position_difference','position_difference_error_upper'),
        ('Fourier direction / rad','separate_Fourier_direction_deflection','direction_error_upper'),
    ):
        line=next(x for x in text.splitlines() if x.startswith(f'| {label} |'))
        centre,error=map(F,line.split('|')[2:4])
        source_centre,source_error=F(gravity[centre_key]),F(gravity[error_key])
        assert centre-error<=source_centre-source_error
        assert source_centre+source_error<=centre+error


def test_projected_errors_round_outwards():
    for x in (0.0005076474164068462,F(1,3),F(-1,3)):
        assert F(summary.decimal(x))<=F(x)<=F(summary.decimal(x,upper=True))


@pytest.mark.parametrize('value,error',[
    (F(1,3),F(1,10**14)),(-F(1,3),F(1,10**14)),
    (F(1,10**200),F(1,10**250)),(-F(1,10**200),F(1,10**250)),
    (F(1,10**400),F(0)),(-F(1,10**400),F(0)),(F(0),F(0)),
])
def test_printed_centre_error_contains_exact_signed_and_tiny_intervals(value,error):
    centre,radius=map(F,summary.difference_with_error(value,error))
    assert centre-radius<=value-error<=value+error<=centre+radius
    with pytest.raises(ValueError,match='nonnegative'):
        summary.difference_with_error(value,-F(1,10**200))


def test_incomplete_custody_inventory_fails_before_projection(tmp_path,monkeypatch):
    monkeypatch.setattr(summary,'PACKAGES',('probe',))
    monkeypatch.setattr(summary.importlib,'import_module',lambda _:SimpleNamespace(pins=lambda:{'a':'one','b':'two'}))
    p=tmp_path/'model/openwave_probe/results.json';p.parent.mkdir(parents=True)
    p.write_text(json.dumps({'source_sha256':{'a':'one'},'verdict':'PASS'}))
    with pytest.raises(ValueError,match='inventory changed'):
        summary.collect(tmp_path)


def test_scale_projection_keeps_member_choices_failures_and_exact_intervals(authenticated_result):
    value=deepcopy(authenticated_result)
    scale=value['native_scale_and_structure']
    assert scale['core']['selection']['selection_uses_resolution_outcome']
    assert not scale['core']['selection']['selection_uses_empirical_radius_or_running_data']
    assert not scale['core']['scope']['reference_bound_energies_are_observed_particle_poles']
    assert not scale['unit_constraint']['scope']['physical_G_and_c_match']
    assert not scale['capacity_histories']['same_data_alternatives_are_independent']
    assert len(scale['capacity_histories']['comparisons'])==8
    assert all(x['rows']==15 and x['fitted_chronometer_parameters']==0
               for x in scale['capacity_histories']['comparisons'])
    assert not scale['full_Haar_dynamics']['scope']['Clay_problem_solved']
    text='\n'.join(summary.render_scale_and_structure(scale))
    for row in scale['core']['members']:
        line=next(x for x in text.splitlines() if x.startswith('| '+row['family']+' |'))
        for display, key in zip(line.split('|')[2:5],
                ('reference_bound_energy_GeV_interval','radius_cells_interval','maximum_suppression_interval')):
            lo,hi=map(F,display.strip().strip('[]').split(', '))
            actual_lo,actual_hi=map(F,row[key])
            assert lo<=actual_lo<=actual_hi<=hi
    scale['capacity_histories']['anchors']['H0_km_s_Mpc']='68.125'
    assert 'H0=68.125' in '\n'.join(summary.render_scale_and_structure(scale))


def test_propagation_projection_distinguishes_shell_from_decay_and_nuclear_members(authenticated_result):
    value=deepcopy(authenticated_result)
    assert len(value['receipts_sha256'])==56
    group=value['propagation_and_nuclei']
    b=group['baryon']['positive_correlated_propagation']
    assert b['k0_complete_N3_operator_change_exact']=='0'
    assert not b['other_full_Fock_spectra_transferred']
    assert not b['exact_Lorentz_dispersion_or_scattering_pole']
    assert group['strange']['source']['complete_two_cell_operator_preserved']
    assert not group['strange']['native_weak_scope']['new_on_shell_weak_matrix_element_certified']
    rows=group['nuclei']['complete_finite_spectra']
    assert [r['target_in_interval'] for r in rows]==[True,False,False]
    assert not group['nuclei']['declared_electrostatic_member']['seven_parent_mass_calibrations_transferred']
    text='\n'.join(summary.render_propagation_and_nuclei(group))
    assert 'strange decay rate is not supplied' in text
    for label,pair in [('Lambda',group['strange']['Lambda_rest_GeV_interval']),
        ('charged pion',group['strange']['pion_rest']['full_ground_GeV_interval']),
        ('charged kaon',group['strange']['kaon_rest']['full_ground_GeV_interval'])]:
        line=next(x for x in text.splitlines() if x.startswith('| '+label+' |'))
        lo,hi=map(F,line.split('|')[2].strip().strip('[]').split(', '))
        a,b=map(F,pair)
        assert lo<=a<=b<=hi


def new_source_data():
    names=('baryon_poles','pole_weak_current','strange_shell_current','coupled_haar_source',
           'coupled_haar_spectrum','euclidean_refinement','source_ward','einstein_error','radiation_capacity')
    return {name:json.loads((ROOT/f'model/openwave_{name}/results.json').read_text())for name in names}


def test_pole_gauge_and_radiation_projection_keeps_exact_domains_and_no_score(authenticated_result):
    result=deepcopy(authenticated_result);group=result['pole_currents_and_full_field_records']
    assert group['physical_comparison_count_added']==0
    assert len(result['receipts_sha256'])==56
    assert group['baryon_poles']['momentum_fibre_count']==12
    assert group['strange_current']['complete_component_probe']['probabilities_are_fixed_total_momentum_fibre_quantities']
    assert not group['strange_current']['complete_component_probe']['normalizable_wave_packet_transfer_executed']
    assert len(group['coupled_Haar_source']['graph']['edges'])==7
    assert len(group['coupled_Haar_spectrum']['records'])==6
    assert not group['Euclidean_response']['continuum_acceptance']['fixed_graph_fixed_coupling_trajectory_admitted']
    assert not group['radiation_capacity']['scope']['independent_observational_validation']
    text='\n'.join(summary.render_pole_currents_and_full_fields(group))
    for boundary in ('not percentage agreements with Nature','total-momentum-zero relative-fibre',
                     'normalizable center-of-mass wave packet is not executed','no Lambda decay lifetime',
                     'without a physical GeV glueball','zero Nature-validation rows'):
        assert boundary in text
    for row in group['coupled_Haar_spectrum']['records']:
        if not row['shared_edge']:continue
        line=next(x for x in text.splitlines()if x.startswith(f"| {row['g_squared_exact']} | {row['spacing_exact']} | {row['diagnostic_max_twice_loop_spin']} |"))
        for shown,key in zip(line.split('|')[4:6],('full_Haar_gap_interval_exact','connected_Wilson_covariance_interval_exact')):
            lo,hi=map(F,shown.strip().strip('[]').split(', '));a,b=map(F,row[key]);assert lo<=a<=b<=hi


@pytest.mark.parametrize('package,path',[
    ('baryon_poles',('physical_3plus1_Lorentz_particle_pole_identification',)),
    ('pole_weak_current',('autonomous_weak_decay_or_neutron_lifetime_proved',)),
    ('strange_shell_current',('physical_on_shell_decay_amplitude_proved',)),
    ('strange_shell_current',('complete_component_probe','normalizable_wave_packet_transfer_executed')),
    ('coupled_haar_source',('scope','physical_GeV_mass')),
    ('coupled_haar_spectrum',('scope','Clay_problem_solved')),
    ('euclidean_refinement',('scope','four_dimensional_interacting_continuum_constructed')),
    ('source_ward',('scope','physical_Newton_strength_matched')),
    ('einstein_error',('scope','cosmological_constant_fixed_by_null_records')),
    ('radiation_capacity',('scope','independent_observational_validation')),
])
def test_unproved_physical_promotion_rejected(package,path):
    data=new_source_data();row=data[package]
    for key in path[:-1]:row=row[key]
    row[path[-1]]=True
    with pytest.raises(ValueError,match='unproved physical validation'):
        summary.pole_currents_and_full_fields(data)


def test_missing_fixed_fibre_caveat_rejected():
    data=new_source_data();data['strange_shell_current']['complete_component_probe']['probabilities_are_fixed_total_momentum_fibre_quantities']=False
    with pytest.raises(ValueError,match='fixed-fibre domain omitted'):
        summary.pole_currents_and_full_fields(data)



def latest_source_data():
    names=('uniform_pole_bands','normalizable_particle_records','scattering_current_support',
        'haar_coupling_flow','cubic_haar_source','gauge_refinement_clock','staggered_ward',
        'replicated_geometry','capacity_perturbations')
    return {name:json.loads((ROOT/f'model/openwave_{name}/results.json').read_text())for name in names}


def test_normalizable_spatial_projection_keeps_calibration_and_full_domains(authenticated_result):
    result=deepcopy(authenticated_result);group=result['normalizable_particles_and_spatial_fields']
    assert len(result['receipts_sha256'])==56 and group['physical_comparison_count_added']==0
    assert group['normalizable_packets']['normalizable_packet']['carrier_is_infinite_chain']
    assert not group['normalizable_packets']['complete_protocol']['finite_local_compilation_of_exact_band_filter_proved']
    assert group['strange_composite_current']['channel']['independent_composite_rank']==23
    assert not group['strange_composite_current']['complete_Lambda_pole_transfer']['normalizable_Lambda_packet_executed']
    assert not group['strange_composite_current']['actual_kinetic_substep']['complete_autonomous_strong_Hamiltonian_used']
    assert group['spatial_Haar_blocking']['geometry']['fine_links']==54
    assert not group['replicated_geometry']['scope']['physical_Newton_strength_matched']
    text='\n'.join(summary.render_normalizable_particles_and_spatial_fields(group))
    for boundary in ('not new rest masses','declared nonlocal spectral service','normalizable Lambda packet',
                     'not the complete autonomous strong evolution','magnetic-factor blocking',
                     'failed G/c strength is unchanged','established holographic-quintessence','zero Nature-validation rows'):
        assert boundary in text
    for row in group['cubic_Haar_source']['records']:
        line=next(x for x in text.splitlines()if x.startswith('| '+row['g_squared_exact']+' |'))
        lo,hi=map(F,line.split('|')[2].strip().strip('[]').split(', '))
        a,b=map(F,row['full_source_gap_interval_exact']);assert lo<=a<=b<=hi
    for row in group['normalizable_packets']['members']:
        line=next(x for x in text.splitlines()if x.startswith('| '+row['assignment']+' |'))
        lo,hi=map(F,line.split('|')[4].strip().strip('[]').split(', '))
        a,b=map(F,row['complete_packet_trial_mean_energy_GeV_interval']);assert lo<=a<=b<=hi


@pytest.mark.parametrize('package,path',[
    ('uniform_pole_bands',('physical_four_dimensional_particle_dictionary_proved',)),
    ('normalizable_particle_records',('complete_protocol','finite_local_compilation_of_exact_band_filter_proved')),
    ('scattering_current_support',('complete_Lambda_pole_transfer','normalizable_Lambda_packet_executed')),
    ('scattering_current_support',('actual_kinetic_substep','complete_autonomous_strong_Hamiltonian_used')),
    ('scattering_current_support',('interacting_spectral_test','Moller_wave_operator_constructed')),
    ('haar_coupling_flow',('Yang_Mills_continuum_constructed',)),
    ('cubic_haar_source',('physical_gap_GeV_or_continuum_identification',)),
    ('gauge_refinement_clock',('scope','autonomous_full_Hamiltonian_RG_derived')),
    ('staggered_ward',('scope','physical_metric_or_stress_derived')),
    ('replicated_geometry',('scope','physical_Newton_strength_matched')),
    ('capacity_perturbations',('canonical_scalar_completion','source_scalar_transfer_proved')),
])
def test_latest_unproved_promotion_rejected(package,path):
    data=latest_source_data();row=data[package]
    for key in path[:-1]:row=row[key]
    row[path[-1]]=True
    with pytest.raises(ValueError,match='unproved spatial-field physical promotion'):
        summary.normalizable_particles_and_spatial_fields(data)


def test_latest_carrier_and_clock_boundary_rejected():
    data=latest_source_data();data['gauge_refinement_clock']['scope']['higher_coarse_operator_remainder_retained']=False
    with pytest.raises(ValueError,match='complete carrier or retained-clock domain omitted'):
        summary.normalizable_particles_and_spatial_fields(data)


def complete_source_data():
    names=('full_strong_current','lambda_wavepacket','cluster_weyl_sequence',
        'refined_haar_source','refined_haar_spectrum','haar_memory_closure',
        'tensor_record_transport','collective_gravity_strength','capacity_transfer')
    return {name:json.loads((ROOT/f'model/openwave_{name}/results.json').read_text())for name in names}


def test_complete_source_projection_keeps_rest_energy_and_calibration_separate(authenticated_result):
    group=deepcopy(authenticated_result['complete_dynamics_and_constraints'])
    assert group['physical_comparison_count_added']==0
    text='\n'.join(summary.render_complete_dynamics_and_constraints(group))
    energy=group['source_packages']['lambda_wavepacket']['calibrated_energy']
    for label,key in [('Lambda rest pole','source_rest_pole_GeV_interval'),
                      ('Lambda packet mean, including kinetic energy','trial_packet_mean_GeV_interval')]:
        line=next(x for x in text.splitlines()if x.startswith('| '+label+' |'))
        lo,hi=map(F,line.split('|')[2].strip().strip('[]').split(', '))
        a,b=map(F,energy[key]);assert lo<=a<=b<=hi
    for boundary in ('not independent mass predictions','on-shell weak spectral weight',
                     'specified Haar-vacuum defect bound','explicit calibration',
                     'complete readout work','not observed matches','zero Nature-validation rows'):
        assert boundary in text


@pytest.mark.parametrize('package,path',[
    ('full_strong_current',('physical_3plus1_decay_or_lifetime_proved',)),
    ('lambda_wavepacket',('calibrated_energy','mean_energy_is_rest_mass')),
    ('cluster_weyl_sequence',('nonzero_native_weak_on_shell_spectral_density_proved',)),
    ('refined_haar_source',('uniform_repeated_refinement_or_4D_limit_proved',)),
    ('refined_haar_spectrum',('derived_running_coupling_or_continuum_limit',)),
    ('haar_memory_closure',('retained_201_modes_exactly_autonomous',)),
    ('tensor_record_transport',('scope','quantum_curvature_acquisition_from_native_source_proved')),
    ('capacity_transfer',('scope','CMB_large_scale_structure_likelihood')),
])
def test_complete_source_projection_rejects_unproved_physical_promotions(package,path):
    data=complete_source_data();row=data[package]
    for key in path[:-1]:row=row[key]
    row[path[-1]]=True
    with pytest.raises(ValueError,match='unproved complete-dynamics physical promotion'):
        summary.complete_dynamics_and_constraints(data)


def test_observed_G_selection_and_uncompiled_spectral_preparation_cannot_be_hidden():
    data=complete_source_data()
    data['collective_gravity_strength']['calibrated_member']['calibration_input']='no observed inputs'
    with pytest.raises(ValueError,match='calibration or carrier boundary omitted'):
        summary.complete_dynamics_and_constraints(data)
    data=complete_source_data()
    data['cluster_weyl_sequence']['constructive_Weyl_family']['finite_local_observer_preparation_of_exact_spectral_vectors_compiled']=True
    with pytest.raises(ValueError,match='source domain or preparation boundary omitted'):
        summary.complete_dynamics_and_constraints(data)
