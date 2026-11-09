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
    text = ' '.join(''.join(root.itertext()).split())
    assert end == 3730
    assert root.find('{http://www.w3.org/2000/svg}g').get('id') == 'joint-field-particle-dependencies'
    for node in graph['nodes']:
        assert node['title'] in text
    assert '1 → 2 → 3 → 4 → 5 → 6 → 7 → 8' in text
    assert 'Separate pure-gauge branch: 5 → 9' in text
    assert text.count('FINITE CONSTRUCTED') == 3
    assert text.count('PHYSICAL PROOF REQUIRED') == 5
    assert text.count('CONDITIONAL GAP CRITERION') == 1
    assert 'massless photon excludes a positive global vacuum gap' in text
    assert 'not physical pole masses' in text
    hrefs = [node.get('href') for node in root.iter('{http://www.w3.org/2000/svg}a')]
    assert '../../model/particle_primary_locality/PROOF.md' in hrefs
    assert '../../model/particle_joint_reconstruction/CANDIDATE.md' in hrefs
    assert '../../model/particle_joint_reconstruction/POSITIVE_ENERGY.md' in hrefs
    assert 'continuous subsequential responses, possibly zero' in text


@pytest.mark.parametrize('mutation', [
    'missing_node', 'physical_promotion', 'finite_promotion', 'branch_edge',
    'global_gap', 'candidate_discharge', 'calibration_prediction',
])
def test_graph_mutations_cannot_silently_promote_the_display(graph, mutation):
    value = deepcopy(graph)
    if mutation == 'missing_node': value['nodes'].pop()
    elif mutation == 'physical_promotion': value['nodes'][3]['status'] = 'finite_constructed'
    elif mutation == 'finite_promotion': value['nodes'][0]['status'] = 'physical_theory_proved'
    elif mutation == 'branch_edge': value['nodes'][8]['depends_on'] = ['particle_observables']
    elif mutation == 'global_gap': value['nodes'][8]['implies_full_theory_positive_global_gap'] = True
    elif mutation == 'candidate_discharge': value['conditional_candidate']['discharges_primary_physical_dependencies'] = True
    else: value['conditional_candidate']['calibrated_tree_coefficients_are_independent_predictions'] = True
    with pytest.raises(ValueError): source_connections.validate(value)
