"""Independent replay and false-promotion controls for native repair constraints."""
from __future__ import annotations

import builtins
import copy
import importlib.util
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import verify_native_repair_flavor_constraints as verifier


@pytest.fixture(scope='module')
def receipt():
    return verifier.strict_load(verifier.OUTPUT)


def test_exact_replay_with_native_and_producer_imports_blocked(receipt,monkeypatch):
    original = builtins.__import__
    forbidden = {'native_repair_flavor_constraints','record_counting_mechanism_certificate',
                 'source_repair_generator_certificate','paper_math',
                 'conditional_quark_mass_replay','quark_rscc_completion_candidate'}
    def guarded(name,*args,**kwargs):
        if name.split('.')[-1] in forbidden:
            raise AssertionError('forbidden producer or target import: '+name)
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    verifier.independent_expectation.cache_clear()
    result = verifier.verify_payload(receipt)
    assert result['independent_states'] == 73789
    assert result['independent_settled_states'] == 303
    assert result['physical_mass_identification'] is False
    global_menu = receipt['global_stationary_menu']
    assert global_menu['all_uniform_shift_classes'] == 6077
    assert global_menu['individual_quadrupole_spectral_types'] == 31
    assert len(global_menu['chi_spectral_type_counts']) == 11
    assert global_menu['conditional_positive_log_gap_ratios'] == 19
    assert global_menu['nondegenerate_spectral_types'] == 22
    assert sum(global_menu['spectral_types_by_total_residue']['0'].values()) == 303
    assert len(global_menu['spectral_types_by_total_residue']['0']) == 4


def test_deterministic_exact_producer_rebuild(receipt):
    path = HERE/'native_repair_flavor_constraints.py'
    spec = importlib.util.spec_from_file_location('native_repair_producer_test',path)
    producer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(producer)
    # Integer/rational/algebraic outputs carry no platform-dependent ODE
    # floats, so byte-equivalent canonical JSON is the appropriate contract.
    assert verifier.serialized(producer.build_payload()) == verifier.serialized(receipt)


@pytest.mark.parametrize('path,value',[
    (('scope','source_selected_vacuum'),True),
    (('scope','physical_quark_operator_derived'),True),
    (('scope','quark_masses_emitted'),True),
    (('scope','quark_targets_loaded'),True),
    (('scope','quark_fit_performed'),True),
    (('scope','stationary_law_unique'),True),
    (('scope','accepted_covariance_equals_proposal_covariance'),True),
    (('scope','physical_clock_derived'),True),
    (('scope','continuum_hessian_claimed'),True),
    (('scope','positive_physical_evidence_claimed'),True),
    (('scope','exact_finite_count_constraints'),1),
    (('scope','individual_spectrum_menu_applies_to_arbitrary_mixtures'),True),
    (('specification','one_record_sector','alphabet'),[-1,0,1,2]),
    (('specification','one_record_sector','P_loaded'),True),
    (('specification','one_record_sector','meaning'),'the unique physical vacuum'),
    (('specification','finite_fiber','clock'),'source-derived Poisson physical clock'),
    (('specification','finite_fiber','proposal_hypothesis'),'no additional hypothesis'),
    (('specification','finite_fiber','history_hypothesis'),'marginal seam symmetry proves IID'),
    (('specification','finite_fiber','conditional_linear_gap_mass_bound'),'all quark mass ratios must be at most30'),
    (('specification','conditional_affine_log_readout','physical_law_derived'),True),
    (('specification','conditional_affine_log_readout','empirical_comparison'),True),
    (('specification','global_stationary_classification','coordinate_bound'),'N_i in {0,1,2} at every total'),
    (('specification','global_stationary_classification','integer_mean_corollary'),'Every signed state is ternary before repair'),
    (('specification','global_stationary_classification','termination'),'Every scheduler terminates, including one that never proposes a legal seam'),
    (('specification','global_stationary_classification','individual_state_scope'),'All averaged mixture responses have one of31spectra'),
    (('specification','global_stationary_classification','physical_boundary'),'The total residue is P-selected'),
    (('specification','global_stationary_classification','sharp_quadrupole_bound','domain'),'All real native settled counts'),
    (('specification','conditional_affine_log_readout','all_totals_extension'),'All physical quark masses obey19ratios unconditionally'),
    (('specification','conditional_affine_log_readout','gap_invariant_equation'),'chi(g)=(g-1)^2*(g+2)^2*(2*g+1)^2/(24*(g^2+g+1)^3)+1/10^100'),
    (('enumeration','states'),73788),
    (('enumeration','settled_states'),302),
    (('enumeration','settled_A5_orbits'),1),
    (('enumeration','A5_invariant_stationary_simplex_dimension'),0),
    (('enumeration','complete_transition_sha256'),'0'*64),
    (('spectrum_K_multiplicities','29/30'),1321),
    (('matching','maximal_matchings_checked'),13132),
    (('matching','complete_matching_sha256'),'0'*64),
    (('matching','largest_minimum_cost'),6),
    (('matching','largest_maximum_steps'),4),
    (('minimal133','settled_states'),1),
    (('minimal133','spectrum_K_multiplicities','1'),1),
    (('nonunique_absorption_witness','terminal_orbits',0,'probability'),'1'),
    (('settled_orbits',1,'trace_Q2'),'16/5'),
    (('settled_orbits',1,'trace_Q3'),'1/1000000000000000000000000000000000'),
    (('settled_orbits',3,'chi'),'1/6'),
    (('settled_orbits',1,'port_covariance_P3'),'1/6'),
    (('quadrupole_type_multiplicities','8/5;0'),179),
    (('accepted_response','full_pair_increment_Gram_coefficient'),'1/5'),
    (('accepted_response','quadrupole_map_squared_norm'),'8/5'),
    (('accepted_response','quadrupole_covariance_upper_coefficient'),'16/25'),
    (('accepted_response','stationary_accepted_moments'),'same as proposal moments'),
    (('conditional_log_gap_ratios',1),'(3/sqrt(5)-1)/2+1/10^100'),
    (('global_stationary_menu','sum_distances'),17),
    (('global_stationary_menu','all_uniform_shift_classes'),6076),
    (('global_stationary_menu','range_three_classes'),0),
    (('global_stationary_menu','individual_quadrupole_spectral_types'),30),
    (('global_stationary_menu','complete_height_sha256'),'0'*64),
    (('global_stationary_menu','nondegenerate_spectral_types'),23),
    (('global_stationary_menu','conditional_positive_log_gap_ratios'),18),
    (('global_stationary_menu','sharp_quadrupole_bound','pair_matchings_checked'),14),
    (('global_stationary_menu','sharp_quadrupole_bound','trace_Q2_upper_bound'),'24/5+1/10^100'),
    (('global_stationary_menu','signed_outside_sector_witness','initial',0),0),
    (('global_stationary_menu','signed_outside_sector_witness','endpoint',0),3),
])
def test_scientific_and_numerical_corruption_rejected(receipt,path,value):
    mutant = copy.deepcopy(receipt)
    cursor = mutant
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    with pytest.raises(ValueError):
        verifier.verify_payload(mutant)


def test_extra_physical_promotion_cannot_be_appended(receipt):
    mutant = copy.deepcopy(receipt)
    mutant['quark_mass_agreement'] = {'derived':True}
    with pytest.raises(ValueError):
        verifier.verify_payload(mutant)


def test_source_pins_cannot_be_rewritten(receipt):
    mutant = copy.deepcopy(receipt)
    mutant['source_pins'][verifier.PINS[0]] = '0'*64
    with pytest.raises(ValueError):
        verifier.verify_payload(mutant)


def test_cached_verification_still_checks_live_source_bytes(receipt,monkeypatch):
    original = Path.read_bytes
    source = verifier.ROOT/verifier.PINS[0]
    def changed(path):
        content = original(path)
        return content+b' ' if path == source else content
    monkeypatch.setattr(Path,'read_bytes',changed)
    with pytest.raises(ValueError,match='source pin'):
        verifier.verify_payload(receipt)


def test_orbit_omission_and_false_normalization_rejected(receipt):
    mutant = copy.deepcopy(receipt)
    mutant['settled_orbits'].pop()
    mutant['enumeration']['settled_A5_orbits'] = 8
    mutant['enumeration']['A5_invariant_stationary_simplex_dimension'] = 7
    with pytest.raises(ValueError):
        verifier.verify_payload(mutant)


@pytest.mark.parametrize('text',['{"x":0,"x":1}','{"x":NaN}','{"x":Infinity}'])
def test_malformed_json_rejected(tmp_path,text):
    path = tmp_path/'receipt.json'
    path.write_text(text)
    with pytest.raises(ValueError):
        verifier.strict_load(path)
