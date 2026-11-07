from itertools import product
import json

import pytest
import sympy as s

from . import certificate
from . import strings as g


def test_native_representation_identity_uses_actual_currents_not_dimension_counts():
    result = g.representation_certificate()
    assert result['actual_native_currents'] == 12
    assert result['link_dimension'] == 60
    assert [row['dimension'] for row in result['blocks']] == [1, 3, 6, 2, 3]
    assert [row['commutator_constraint_rank'] for row in result['blocks']] == [0, 8, 35, 3, 8]
    assert result['Hom_dimensions'][1][4] == 0  # Equal dimensions, inequivalent representations.
    assert all(row['nonvacuum_single_endpoint_invariant_dimension'] == 0 for row in result['blocks'])


def test_full_one_edge_action_matches_all_owning_matter_link_columns():
    result = g.one_edge_identity()
    assert result['all_matter_and_link_columns_compared'] == 1800
    assert result['electric_projector_and_literal_raising_amplitudes_identical']


def test_every_matter_and_link_block_assignment_is_included_in_the_census():
    result = g.central_sector_census()
    assert result['complete_block_assignments_checked'] == 21600
    assert len(result['admitted_assignments']) == 20
    assert result['no_extra_off_path_flux_or_other_matter_blocks']
    for row in result['admitted_assignments']:
        if row['matter_position'] == 0:
            assert row['link_blocks'] == [None, None, None]


def test_full_tree_operator_and_all_nonabelian_gauss_equations():
    result = g.exact_tree_certificate()
    assert result['exact_Gauss_equations_checked'] == 960
    assert [row['complete_Gauss_kernel_dimension'] for row in result['blocks']] == [4]*5
    assert result['blocks'][2]['ambient_one_fermion_reference_link_dimension'] == 77760000
    assert result['blocks'][2]['string_support_sizes'] == [6, 36, 36, 216]
    assert result['blocks'][2]['relative_Gram_diagonal'] == [1, 6, 6, 36]


def test_independent_full_kernel_rank_inside_every_surviving_weak_block_assignment():
    # Exhaust the raw tensor basis admitted by the central equations for
    # the actual weak block on the branched tree, without using strings.
    tree, charge = g.Tree(), '3/2'
    indices = dict(g.native.native_blocks())[charge]
    _, address = g.link_dictionary()
    basis = []
    for site in range(tree.size):
        path_edges = tree.path(site)[1:]
        for matter, reference in product(indices, repeat=2):
            for endpoint_labels in product(indices, repeat=2*len(path_edges)):
                links = [0]*len(tree.edges)
                for k, child in enumerate(path_edges):
                    links[child-1] = address[charge, endpoint_labels[2*k], endpoint_labels[2*k+1]]
                basis.append((site, matter, reference, tuple(links)))
    assert len(basis) == 100
    output_rows, entries = {}, {}
    for column, key in enumerate(basis):
        for current, generator in enumerate(g.native.native_currents()):
            for vertex in range(tree.size):
                for output, value in g.gauss_action(tree, {key: s.Integer(1)}, generator, vertex).items():
                    row_key = current, vertex, output
                    row = output_rows.setdefault(row_key, len(output_rows))
                    entries[row, column] = value
    constraints = s.SparseMatrix(len(output_rows), len(basis), entries)
    assert constraints.rank() == 96
    assert len(basis)-constraints.rank() == tree.size


def test_root_singleton_and_deeper_chain_are_included():
    for parents in [(-1,), (-1, 0, 1, 2)]:
        tree = g.Tree(parents)
        result = g.exact_tree_certificate(tree)
        assert all(row['complete_Gauss_kernel_dimension'] == len(parents) for row in result['blocks'])


def test_reference_and_oriented_link_are_necessary_on_the_actual_carrier():
    tree = g.Tree()
    center = g.native.native_currents()[-1]
    root = g.string_state(tree, '-1/2', 0)
    away = g.string_state(tree, '-1/2', 3)
    assert g.gauss_action(tree, root, center, 0, include_reference=False)
    assert g.gauss_action(tree, away, center, 3, link_sign=-1)
    assert not g.gauss_action(tree, root, center, 0)
    assert not g.gauss_action(tree, away, center, 3)


def test_nonorthonormal_string_coordinates_require_the_gram_metric():
    h, gram, normalized = g.reduced_model(g.Tree(), '-1/2')
    assert h != h.H
    assert gram*h == h.H*gram
    root = s.diag(*(s.sqrt(x) for x in gram.diagonal()))
    assert root*h*root.inv() == normalized
    assert normalized[1, 0] == s.sqrt(6)
    assert tuple(normalized.diagonal()) == (0, 1, 1, 2)


def test_complete_positions_holds_zero_outcomes_and_detector_work():
    result = g.complete_history(g.Tree(), '-1/2')
    assert len(result['outcomes']) == 64
    assert sum(s.Rational(row['probability']) for row in result['outcomes']) == 1
    assert any(row['word'] == [0, 0, 0] and s.Rational(row['probability']) > 0 for row in result['outcomes'])
    assert sum(s.Rational(row['detector_work']) for row in result['energy_ledger']) == s.Rational(result['total_detector_work'])
    assert s.Rational(result['total_detector_work']) > 0
    stopped = g.complete_history(g.Tree(), '-1/2', hopping=0)
    assert len(stopped['outcomes']) == 64
    assert sum(s.Rational(row['probability']) == 0 for row in stopped['outcomes']) == 63
    assert stopped['total_detector_work'] == '0'


def test_exact_finite_volume_tail_and_stopped_hopping():
    value = g.tail_bound(12, 0, 3, 6, s.Rational(1, 8))
    assert value['first_omitted_order'] == 13
    assert value['z_upper'] == '9/8'
    assert 0 < s.Rational(value['complete_record_total_variation_upper']) < s.Rational(1, 10**8)
    assert g.tail_bound(12, 0, 3, 6, s.Rational(1, 8), hopping=0)['complete_record_total_variation_upper'] == '0'
    assert g.tail_bound(0, 0, 3, 6, s.Integer(1))['complete_record_total_variation_upper'] == '1'


@pytest.mark.parametrize('parents', [(), (0,), (-1, 1), (-1, 0, 2), (-1, True), [-1, 0]])
def test_tree_validation_rejects_cycles_and_ambiguous_vertex_labels(parents):
    with pytest.raises(ValueError):
        g.Tree(parents)


@pytest.mark.parametrize('kwargs', [{'duration': .1}, {'link_energy': -1}, {'rounds': True}, {'rounds': 0}])
def test_exact_control_domain(kwargs):
    with pytest.raises(ValueError):
        g.complete_history(g.Tree(), '-1/2', **kwargs)


def test_receipt_retains_scope_and_rejects_mutated_mass_claim(tmp_path):
    value = certificate.produce()
    assert value['scope']['all_one_fermion_Gauss_states_on_the_declared_tree_classified']
    assert not value['scope']['physical_or_asymptotic_particle_mass_identified']
    assert not value['scope']['cycles_or_additional_matter_included_in_kernel_classification']
    value['scope']['physical_or_asymptotic_particle_mass_identified'] = True
    output = tmp_path/'mutated.json'
    output.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='differs'):
        certificate.verify(output)
