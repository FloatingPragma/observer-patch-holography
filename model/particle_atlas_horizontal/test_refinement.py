"""The illustrated thermal/strain join rejects changed source domains."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path

import pytest

from .refinement_readouts import validate_inputs


@pytest.fixture(scope='module')
def report():
    path=Path(__file__).resolve().parents[1]/'particle_observer_outputs/results.json'
    data=json.loads(path.read_text())
    return {key:data[key] for key in ('refinement_uniform_causality',
            'interacting_field_thermal_records','declared_spatial_geometry_response')}


def test_frozen_shared_thermal_and_strain_contract(report):
    validate_inputs(report)


@pytest.mark.parametrize('path,value',[
    (('declared_spatial_geometry_response','shared_quantum_rotor_work','magnetic_coefficient'),'1/7'),
    (('declared_spatial_geometry_response','shared_quantum_rotor_work','source'),'another Hamiltonian'),
    (('interacting_field_thermal_records','finite_state_preparation','dimension'),207),
    (('interacting_field_thermal_records','complete_flux_records','source_dimension'),207),
    (('declared_spatial_geometry_response','shared_quantum_rotor_work','inverse_temperature_target'),'1/32'),
    (('interacting_field_thermal_records','finite_state_preparation','positive_polynomial','beta'),'1/32'),
    (('interacting_field_thermal_records','finite_state_preparation','positive_polynomial','order'),13),
    (('interacting_field_thermal_records','complete_flux_records','measurement_times'),['0','1/32']),
    (('interacting_field_thermal_records','finite_state_preparation','joint_purification_stationary'),True),
    (('interacting_field_thermal_records','source_thermal_identities','polynomial_preparation_exactly_Gibbs'),True),
    (('interacting_field_thermal_records','source_thermal_identities','arbitrary_reference_operations_passive'),True),
    (('interacting_field_thermal_records','physical_vacuum_or_continuum_derived'),True),
    (('refinement_uniform_causality','Lorentz_covariance_continuum_vacuum_or_chirality_proved'),True),
    (('declared_spatial_geometry_response','metric_identification_boundary','proper_geometric_stress_or_Einstein_limit_claimed'),True),
    (('declared_spatial_geometry_response','compact_coframe_causality','compact_family','proper_geometry_Lorentz_covariance_or_limiting_state_proved'),True),
    (('declared_spatial_geometry_response','compact_coframe_causality','shared_rotor_thermal_state_transferred_to_spatial_refinement'),True),
    (('declared_spatial_geometry_response','shared_quantum_rotor_work','rotor_is_a_normalizable_reducing_subspace_of_full3D_quantum_field'),True),
    (('declared_spatial_geometry_response','shared_quantum_rotor_work','pressure_is_a_physical_cosmological_equation_of_state'),True),
    (('declared_spatial_geometry_response','shared_quantum_rotor_work','global_Gibbs_entropy_slope_is_regional_geometric_entropy'),True),
])
def test_changed_source_state_clock_or_interpretation_is_rejected(report,path,value):
    changed=deepcopy(report); target=changed
    for key in path[:-1]:target=target[key]
    target[path[-1]]=value
    with pytest.raises(ValueError):validate_inputs(changed)


def test_preparation_error_cannot_be_transferred_from_a_different_state(report):
    changed=deepcopy(report)
    work=changed['declared_spatial_geometry_response']['shared_quantum_rotor_work']
    work['target_state_trace_error_upper']=str(2*Fraction(work['target_state_trace_error_upper']))
    with pytest.raises(ValueError,match='preparation errors disagree'):validate_inputs(changed)
