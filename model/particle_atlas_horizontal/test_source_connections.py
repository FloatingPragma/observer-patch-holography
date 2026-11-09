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
    assert end == 3730
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
    assert '../../model/particle_primary_locality/PROOF.md' in hrefs
    assert '../../model/particle_joint_reconstruction/CANDIDATE.md' in hrefs
    assert '../../model/particle_primary_spatial_algebra/PROOF.md' in hrefs
    assert '../../model/particle_primary_response_moments/PROOF.md' in hrefs
    assert '../../model/particle_primary_mass_readout/PROOF.md' in hrefs
    assert '../../model/particle_primary_capacity/PROOF.md' in hrefs
    assert '../../model/particle_primary_field_algebra/PROOF.md' in hrefs
    assert '../../model/particle_primary_stress/PROOF.md' in hrefs
    assert '../../model/particle_positive_window/PROOF.md' in hrefs
    for owner in ('particle_primary_momentum', 'particle_primary_pair_channels',
                  'particle_primary_refinement_response', 'particle_primary_gravity_premises',
                  'particle_primary_residual_refinement', 'particle_primary_occupation_bound',
                  'particle_primary_variational_readout', 'particle_primary_joint_spectroscopy',
                  'particle_primary_branch_resolution', 'particle_primary_energy_tightness',
                  'particle_primary_origin_momentum', 'particle_primary_open_coframe',
                  'particle_primary_charged_string', 'particle_primary_odd_residual',
                  'particle_primary_momentum_transfer', 'particle_primary_local_fluctuation',
                  'particle_primary_spectral_resolution', 'particle_primary_charged_readout',
                  'particle_primary_staggered_response', 'particle_primary_current_tail'):
        assert '../../model/'+owner+'/PROOF.md' in hrefs
    assert 'common real-time histories with positive GNS energy' in text
    assert 'coin and reference remain in the complete output' in text
    assert 'primary physical-mass table in GeV is empty' in text
    assert '[1/(100β),46/a]' in text
    assert 'Unfiltered spatial momentum and fixed-origin invariance are unproved' in text
    assert 'relative-observation sector has commuting H≥0 and P' in text
    assert 'HΩ=PΩ=0' in text
    assert 'joint sinc² energy/momentum kernels; their weight may vanish' in text
    assert 'limit may be purely elastic' in text
    assert 'CAR-even records' in text
    assert 'finite Fourier labels, not exact momentum projections' in text
    assert 'ground energy lies in' in text
    assert 'No level is isolated' in text
    assert 'Source-clock energy transfer on the original selected state: local pairs ≥3/20 in [4,8]' in text
    assert 'opposite Q6=±6 transfers' in text
    assert 'charged endpoint pairs >3/10 in [9.75,33.75]' in text
    assert 'excited energy in [1.06662597,14.69223697]' in text
    assert 'residual ≤4.947618257' in text
    assert 'original prepared weight >1/4' in text
    assert 'Idle memory retains the displaced source state and reference' in text
    assert 'coherent A2g pairs >17/1280 in [3,9]' in text
    assert 'radix-five coherent meters' in text
    assert 'original joint success weight >1/4' in text
    assert 'without reading initial absolute charges' in text
    assert 'Absolute terminal charge is a separate disturbing readout' in text
    assert 'identify no particle species or rest mass' in text
    assert '104.806' not in text and '5.913111' not in text
    assert 'Fourier-to-momentum map require proof' not in text
    assert 'H(I)=H_primary exactly' in text
    assert 'Proper geometry, Lorentz stress' in text
    assert 'source-marginal relative entropy D/β' in text
    assert 'entropy comparison traces out the quantum reference' in text
    assert 'Finite mass comparison needs a resolved branch' in text


@pytest.mark.parametrize('mutation', [
    'missing_node', 'physical_promotion', 'finite_promotion', 'branch_edge',
    'global_gap', 'candidate_discharge', 'calibration_prediction', 'gravity_promotion',
    'target_change','deferred_branch_erased',
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
    else: value['conditional_candidate']['calibrated_tree_coefficients_are_independent_predictions'] = True
    with pytest.raises(ValueError): source_connections.validate(value)
