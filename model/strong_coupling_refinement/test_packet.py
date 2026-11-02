from copy import deepcopy
import pytest
from . import packet


@pytest.fixture(scope='module')
def receipt():
    return packet._produce()


def test_current_source_cold_replay_is_deterministic(receipt):
    assert packet.report() == receipt


@pytest.mark.parametrize('key', [
    'interacting_relativistic_continuum_constructed','all_growing_coupling_trajectories_excluded',
    'physical_particle_species_identified','Einstein_dynamics_derived',
    'old_c1_over8_U1_cube_in_new_global_expansion_domain','global_bound_failure_is_a_model_no_go'])
def test_finite_interaction_cannot_be_promoted_to_physics_or_global_no_go(receipt,key):
    row=deepcopy(receipt);row['theorem'][key]=True
    with pytest.raises(ValueError):packet.validate(row)


def test_unbounded_U_compactness_and_both_charge_branches_required(receipt):
    for key in ('new_unbounded_U_compactness_argument','both_charge_branches_retained','outside_cyclicity_promotes_to_full_hierarchy'):
        row=deepcopy(receipt);row['theorem'][key]=False
        with pytest.raises(ValueError):packet.validate(row)


def test_inherited_custody_status_cannot_be_hidden(receipt):
    row=deepcopy(receipt)
    prior=row['prior_custody'][packet.PRIORS[0]]
    prior['current_source_mismatches']=['deliberate/changed/source.py']
    prior['current_source_verified']=True
    with pytest.raises(ValueError):packet.validate(row)


def test_complete_intermediate_car_carrier_is_retained(receipt):
    row=deepcopy(receipt);row['controls']['basis_truncation_before_hopping']=True
    with pytest.raises(ValueError):packet.validate(row)


@pytest.mark.parametrize('section,key,replacement', [
    ('theorem','interacting_relativistic_continuum_constructed',0),
    ('theorem','both_charge_branches_retained',1),
    ('theorem','free_limit_valleys',8.0),
    ('bounds','finite_spin_limit_is_interacting',1),
    ('bounds','global_norm_failure_is_a_model_no_go',0),
])
def test_scope_and_bound_flags_reject_numeric_aliases(receipt,section,key,replacement):
    row=deepcopy(receipt);row[section][key]=replacement
    with pytest.raises(ValueError):packet.validate(row)


def test_exact_control_payload_rejects_boolean_and_integer_aliases(receipt):
    row=deepcopy(receipt)
    row['controls']['interaction']['entanglement_witness_survives_full_source_error']=1
    with pytest.raises(ValueError):packet.validate(row)
    row=deepcopy(receipt)
    row['controls']['exact_graph_audits'][0]['half_filled_dimension']=6.0
    with pytest.raises(ValueError):packet.validate(row)


def test_historical_custody_rejects_numeric_flags_and_fabricated_freshness(receipt):
    row=deepcopy(receipt)
    row['prior_custody'][packet.PRIORS[0]]['historical_all_pass']=1
    with pytest.raises(ValueError):packet.validate(row)
    row=deepcopy(receipt)
    for prior in row['prior_custody'].values():
        prior['current_source_mismatches']=[]
        prior['current_source_verified']=True
    row['dependency_custody_pass']=row['all_pass']=True
    # Deliberately invalidate one claimed historical state even if a future
    # checkout happens to have refreshed every inherited receipt.
    row['prior_custody'][packet.PRIORS[0]]['historical_all_pass']=False
    with pytest.raises(ValueError):packet.validate(row)
