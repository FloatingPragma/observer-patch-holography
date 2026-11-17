"""Canonical dependency coverage and finite/physical claim boundaries."""
from copy import deepcopy
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from . import drawing, source_connections


@pytest.fixture(scope='module')
def graph():
    root = Path(__file__).resolve().parents[2]
    return json.loads((root/'model/particle_joint_reconstruction/dependencies.json').read_text())


def test_dependency_panel_retains_graph_edges_domains_and_height(graph):
    canvas = drawing.Canvas()
    end = source_connections.render(canvas, {'joint_dependencies': graph}, 100)
    root = ET.fromstring('<svg xmlns="http://www.w3.org/2000/svg">'+''.join(canvas.items)+'</svg>')
    text = ' '.join(' '.join(root.itertext()).split())
    assert end == 4530
    assert root.find('{http://www.w3.org/2000/svg}g').get('id') == 'joint-field-particle-dependencies'
    for node in graph['nodes']:
        assert node['title'] in text
    assert '1 → 2 → 3 → 4 → 5 → 6 → 7 → 8' in text
    assert 'Deferred pure gauge: 5 → 9' in text
    assert 'Deferred gravity: 5 → 10' in text
    assert text.count('FINITE CONSTRUCTED') == 3
    assert text.count('PHYSICAL PROOF REQUIRED') == 6
    assert text.count('CONDITIONAL GAP CRITERION') == 1
    assert 'massless photon excludes a positive global vacuum gap' in text
    assert 'not physical pole masses' in text
    hrefs = [node.get('href') for node in root.iter('{http://www.w3.org/2000/svg}a')]
    assert '../../model/ccf_lt2/PROOF.md' in hrefs
    assert '../../model/ccf_lt2/JOINT_GEOMETRY_PROOF.md' in hrefs
    assert '../../model/ccf_lt2/EINSTEIN_BRIDGE.md' in hrefs
    assert '../../model/particle_joint_reconstruction/CANDIDATE.md' in hrefs
    official, marker, comparison = text.partition('INDEPENDENT WILSON COMPARISON: FOUR FINITE TESTS')
    assert marker and 'CCF-LT-2' in official
    assert '96 native fermion modes per vertex' in official
    assert '4567 and 4564' in official and '−9/32, +9/32 and 0' in official
    assert 'uniform physical cone are separate conditions' in official
    assert 'One open Wilson program fixes' not in official
    assert 'ΔE ∈ [63−59/1351350, 63+1/45675]' in comparison
    assert 'original weight >c_g/13' in comparison
    assert 'Each of four raw spin rows has original pair-energy weight >1/16' in comparison
    assert 'M₂ ≤ 4,000,000/a' in comparison
    assert 'electric energy ≤91V/a' in comparison
    assert 'A=h′₀ • B=i[H,A] • applied work retained' in comparison
    assert 'Controlled energy conservation includes the work supplied by the pulse' in comparison
    for owner in ('particle_primary_reduced_sector','particle_primary_regional_fermion',
                  'particle_primary_weighted_spectral_tail','particle_primary_energy_response'):
        assert '../../model/'+owner+'/PROOF.md' in hrefs


@pytest.mark.parametrize('mutation', [
    'missing_node', 'physical_promotion', 'finite_promotion', 'branch_edge',
    'global_gap', 'candidate_discharge', 'calibration_prediction', 'gravity_promotion',
    'target_change','deferred_branch_erased','wrong_program','wrong_capacity_law',
])
def test_graph_mutations_cannot_silently_promote_the_display(graph, mutation):
    value = deepcopy(graph)
    if mutation == 'missing_node': value['nodes'].pop()
    elif mutation == 'physical_promotion': value['nodes'][3]['status'] = 'finite_constructed'
    elif mutation == 'finite_promotion': value['nodes'][0]['status'] = 'physical_theory_proved'
    elif mutation == 'branch_edge': value['nodes'][8]['depends_on'] = ['particle_observables']
    elif mutation == 'global_gap': value['nodes'][8]['implies_full_theory_positive_global_gap'] = True
    elif mutation == 'candidate_discharge': value['conditional_candidate']['discharges_primary_physical_dependencies'] = True
    elif mutation == 'gravity_promotion': value['nodes'][9]['proper_geometry_or_Einstein_limit_proved'] = True
    elif mutation == 'target_change': value['active_target']='Einstein dynamics required before mass comparison'
    elif mutation == 'deferred_branch_erased': value['deferred_branches']=[]
    elif mutation == 'wrong_program': value['primary_program']='independent-Wilson'
    elif mutation == 'wrong_capacity_law': value['candidate']='shared-capacity'
    else: value['conditional_candidate']['calibrated_tree_coefficients_are_independent_predictions'] = True
    with pytest.raises(ValueError): source_connections.validate(value)
