"""Reject domain transfer and hidden premises in the unified requirement map."""
import copy

import pytest

from model.particle_joint_reconstruction import run
from tools import strict_json


@pytest.fixture
def graph():
    return strict_json.load(run.HERE / 'dependencies.json')


def test_all_component_requirements_and_conditional_assumptions_have_destinations(graph):
    run.validate_dependencies(graph)
    omitted = copy.deepcopy(graph)
    for row in omitted['nodes']:
        row['conditional_assumptions'] = [item for item in row['conditional_assumptions']
                                          if item != 'chiral_measure_and_mirror_control']
    with pytest.raises(ValueError, match='omitted'):
        run.validate_dependencies(omitted)


@pytest.mark.parametrize('node', [3, 4, 5, 6, 7, 8, 9])
def test_declared_status_cannot_prove_physics(graph, node):
    graph['nodes'][node]['status'] = 'finite_constructed'
    with pytest.raises(ValueError, match='promoted'):
        run.validate_dependencies(graph)


def test_circular_physical_justification_rejected(graph):
    graph['nodes'][0]['depends_on'] = ['particle_observables']
    with pytest.raises(ValueError, match='cyclic'):
        run.validate_dependencies(graph)


def test_missing_evidence_and_unknown_original_requirement_rejected(graph):
    damaged = copy.deepcopy(graph)
    damaged['nodes'][2]['evidence'][0]['path'] = '../unproved.md'
    with pytest.raises(ValueError, match='existing repository'):
        run.validate_dependencies(damaged)
    graph['requirement_maps']['yang_mills_Y1_Y6']['Y6'] = ['conjecture_as_theorem']
    with pytest.raises(ValueError, match='destination'):
        run.validate_dependencies(graph)


@pytest.mark.parametrize('field', ['discharges_primary_physical_dependencies',
                                  'calibrated_tree_coefficients_are_independent_predictions',
                                  'full_theory_positive_global_gap_claimed'])
def test_calibration_or_assumed_qft_cannot_discharge_source_output(graph, field):
    graph['conditional_candidate'][field] = True
    with pytest.raises(ValueError, match='conditional extension'):
        run.validate_dependencies(graph)


def test_pure_gauge_must_supply_its_own_shared_nodes(graph):
    graph['nodes'][8]['required_shared_nodes'] = ['relativistic_field']
    with pytest.raises(ValueError, match='own source'):
        run.validate_dependencies(graph)


def test_finite_current_cannot_be_promoted_to_einstein_dynamics(graph):
    graph['nodes'][9]['proper_geometry_or_Einstein_limit_proved'] = True
    with pytest.raises(ValueError, match='Einstein'):
        run.validate_dependencies(graph)


def test_missing_component_lists_and_erased_dependencies_rejected(graph):
    damaged = copy.deepcopy(graph)
    damaged['requirement_maps'] = {}
    with pytest.raises(ValueError, match='coverage'):
        run.validate_dependencies(damaged)
    damaged = copy.deepcopy(graph)
    del damaged['requirement_maps']['native_full_spectrum_ten']['6']
    with pytest.raises(ValueError, match='coverage'):
        run.validate_dependencies(damaged)
    for row in graph['nodes']:
        row['depends_on'] = []
    with pytest.raises(ValueError, match='dependency was removed'):
        run.validate_dependencies(graph)


@pytest.mark.parametrize('node', [0, 1, 2])
def test_finite_domain_cannot_be_relabelled_as_physical_continuum(graph, node):
    graph['nodes'][node]['domain'] = 'All physical particles of an interacting continuum'
    with pytest.raises(ValueError, match='actual experiment'):
        run.validate_dependencies(graph)
