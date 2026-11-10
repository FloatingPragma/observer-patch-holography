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
                  'particle_primary_staggered_response', 'particle_primary_current_tail',
                  'particle_primary_sector_threshold', 'particle_primary_populated_momentum',
                  'particle_primary_weighted_survival', 'particle_primary_charge_conjugation',
                  'particle_primary_composite_spin', 'particle_primary_composite_transport',
                  'particle_primary_selected_odd_response', 'particle_primary_covariant_stress',
                  'particle_primary_selected_spin_response', 'particle_primary_odd_resolvent',
                  'particle_primary_odd_low_energy', 'particle_primary_geometric_cut',
                  'particle_primary_coherent_odd_response',
                  'particle_primary_original_shell_moments', 'particle_primary_original_branch', 'particle_primary_critical_refinement', 'particle_primary_boundary_information',
                  'particle_primary_odd_sum_rules', 'particle_primary_branch_certificate', 'particle_primary_work_spectrum', 'particle_primary_boundary_correlation',
                  'particle_primary_reduced_sector', 'particle_primary_regional_fermion', 'particle_primary_weighted_spectral_tail', 'particle_primary_energy_response'):
        assert '../../model/'+owner+'/PROOF.md' in hrefs
    assert '../../model/particle_primary_coherent_odd_response/KERNEL_PROOF.md' in hrefs
    assert 'Filtered and unfiltered probes share positive time dynamics' in text
    assert 'primary physical-mass table in GeV is empty' in text
    assert '[1/(100β),46/a]' in text
    assert 'fixed-band survival and spatial momentum require proof' in text
    assert 'relative-observation sector has commuting H≥0 and P' in text
    assert 'HΩ=PΩ=0' in text
    assert 'sinc² energy/momentum weight may vanish' in text
    assert 'weight escapes every fixed upper energy band' in text
    assert 'G1g ⊕ Hg' in text and 'G1u ⊕ Hu' in text
    assert 'disjoint-site fermionic exchange is exact' in text
    assert 'Finite Fourier labels are not exact momentum projections' in text
    assert 'Radix-five charge and radix-three spin counters' in text
    assert 'source spin increments ±1/2' in text
    assert 'each of four complete spin rows has original energy-band weight >1/16' in text
    assert 'Terminal spin/charge readout is a separate disturbing operation' in text
    assert '104.806' not in text and '5.913111' not in text
    assert 'Fourier-to-momentum map require proof' not in text
    assert 'Proper geometry, Lorentz stress' in text
    assert 'Finite mass comparison needs a resolved branch' in text
    assert 'Five coherent qutrit meters' in text
    assert 'Each complete Walsh spatial row has original positive-energy' in text
    assert 'weight >1/24 in [0.25,17.75]' in text
    assert 'G1/H parts exceed 5/384 and 11/384' in text
    assert 'All 96 successes and 8 failures remain' in text
    assert 'separately chosen normalized phase has weight >1/3' in text
    assert 'not total-spin or conserved-momentum eigenvalues' in text
    assert 'Four uncomputed spectral measures determine each summed multiplet' in text
    assert 'Complete selector sums preserve the original energy law' in text
    assert 'blank centered at 9/a retains weight >1/80 at energy scale 1/a' in text
    assert 'finite universal area coefficient' in text
    assert 'same gauge representation and identical reflected records' in text
    assert 'Native Γ5 is not established Lorentz chirality' in text
    assert 'equivariant field/current dictionary' in text
    assert 'one-handed sector' in text
    assert 'true finite gauge transition is isolated near 63 source units' in text
    assert 'ΔE ∈ [63−59/1351350, 63+1/45675]' in text
    assert 'original weight >c_g/13' in text
    assert 'sector weight vanishes under the selected refinement' in text
    assert 'neither a surviving particle nor a GeV mass' in text
    assert 'The global ground is below −723' in text
    assert 'Each of four raw spin rows has original pair-energy weight >1/16' in text
    assert 'A separate A1/T1 coefficient readout has weights >5/64 and >15/64' in text
    assert 'Pair transfer energy cannot be divided by two' in text
    assert 'no electron or Lorentz spin is identified' in text
    assert 'M₂ ≤ 4,000,000/a' in text
    assert 'electric energy ≤91V/a' in text
    assert 'floor of order a and an upper endpoint of order 1/a' in text
    assert 'does not give a nonzero fixed-energy limit' in text
    assert 'Controlled energy conservation includes the work supplied by the pulse' in text
    assert 'A=h′₀ • B=i[H,A] • applied work retained' in text
    assert 'J+C_center≤2β‖V_cut‖+r' in text
    assert 'All 96 successes and 8 failures remain' in text


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
