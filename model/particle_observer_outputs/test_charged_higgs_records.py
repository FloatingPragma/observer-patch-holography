"""Literal full-CAR charged sectors and complete gauge-compensated records."""
import pytest
import sympy as s

from model.particle_observer_outputs import charged_higgs_records as c
from model.particle_observer_outputs import higgs_records as n


@pytest.fixture(scope='module')
def data():
    return c.charged_data(s.diag(1, 2, 3), 1, 0)


def test_actual_sector_dimensions_and_every_vertex_input(data):
    assert len(data['plus_basis']) == 1650
    assert len(data['minus_basis']) == 36
    assert data['W_one'].shape == (1632, 18)
    assert c.verify_full_vertex(data) == {'complete_addition_input_columns': 1650,
                                          'complete_removal_input_columns': 36}
    assert data['plus_H'] == data['plus_H'].H
    assert data['minus_energy'] == data['plus_dark_energy'] == 1


def test_literal_gram_and_all_dark_dimensions(data):
    expected = []
    for sigma2 in (1, 4, 9):
        expected.extend([56-sigma2]*4+[56-2*sigma2]*2)
    assert data['gram'] == s.diag(*expected)
    assert min(expected) >= 28
    assert s.trace(data['pair_density']) == 112
    assert data['plus_dark_rank'] == 1632-18 == 1614
    assert data['plus_bright_rank']+data['plus_dark_rank'] == 1650


def test_complex_right_family_transpose_is_necessary():
    y = s.Matrix([[1, s.I, 2], [3, 2, 0], [1-s.I, 1, 4]])
    d = c.charged_data(y)
    e_indices = [4, 10, 16]
    correct = n._matrix(2*(y.H*y).T)
    wrong = n._matrix(2*y.H*y)
    assert d['pair_density'].extract(e_indices, e_indices) == correct
    assert not n._zero(correct-wrong)
    assert n._zero(d['gram']+d['pair_density']-d['neutral']['g_squared']*s.eye(18))


def test_complete_probe_effects_on_all307_inputs(data):
    additions, removals = c.creation_removal_maps(data)
    assert len(additions) == len(removals) == 18
    add = sum((a.H*a for a in additions), s.zeros(307))/18
    remove = sum((a.H*a for a in removals), s.zeros(307))/18
    assert n._zero(add+remove-s.eye(307))
    assert remove == data['neutral']['occupation']/9
    assert add == s.eye(307)-remove


def test_actual_ground_addition_and_removal_vectors(data):
    vectors = c.ground_probe_vectors(data)
    alpha, beta = s.Rational(8, 15), s.Rational(7, 15)
    add, remove = vectors['addition'], vectors['removal']
    assert n._zero(add[:18, :]-s.sqrt(alpha)*s.eye(18))
    assert n._zero(add[18:, :]+s.sqrt(beta/56)*data['W_one'])
    assert n._zero(add.H*add+(remove.H*remove).T-s.eye(18))
    assert s.simplify(s.trace(remove.H*remove)/18-s.Rational(7, 135)) == 0
    assert s.simplify(s.trace(add.H*add)/18-s.Rational(128, 135)) == 0


def test_schur_resolvent_solves_the_actual1650_operator(data):
    vectors = c.ground_probe_vectors(data)
    full = vectors['addition']
    z = s.I
    x = z-data['fermion_energy']
    denom = x*(x-data['pair_cost'])*s.eye(18)-data['gram']
    top = denom.inv()*((x-data['pair_cost'])*full[:18, :]
                      +data['W_one'].H*full[18:, :])
    bottom = (full[18:, :]+data['W_one']*top)/(x-data['pair_cost'])
    solution = s.ImmutableSparseMatrix.vstack(top, bottom)
    assert n._zero((z*s.eye(1650)-data['plus_H'])*solution-full)
    assert n._zero(full.H*solution-c.addition_resolvent(data, z))


def test_complete_lines_have_original_weights_and_common_number_offset(data):
    rows = c.response_lines(data)
    assert sum(row['multiplicity'] for row in rows) == 18
    assert len(rows) == 6
    assert all(row['removal_ground_referenced_energy'] == 8 for row in rows)
    total = 0
    for row in rows:
        per_mode = row['removal_weight_per_normalized_mode']
        for line in row['addition_lines']:
            per_mode += line['weight_per_normalized_mode']
            assert line['weight_per_normalized_mode'].is_nonnegative
            total += line['unconditional_probe_probability']
        assert s.simplify(per_mode-1) == 0
        total += row['removal_unconditional_probe_probability']
    assert s.simplify(total-1) == 0
    shifted = c.charged_data(data['Y'], 1, s.Rational(1, 3))
    assert shifted['neutral']['H'] == data['neutral']['H']
    shifted_rows = c.response_lines(shifted)
    for base, changed in zip(rows, shifted_rows, strict=True):
        assert changed['removal_ground_referenced_energy'] == base['removal_ground_referenced_energy']-s.Rational(1, 3)
        for old, new in zip(base['addition_lines'], changed['addition_lines'], strict=True):
            assert s.simplify(new['ground_referenced_energy']-old['ground_referenced_energy']) == s.Rational(1, 3)
            assert new['weight_per_normalized_mode'] == old['weight_per_normalized_mode']


def test_dark_states_are_unseen_in_ground_and_seen_by_other_admitted_input(data):
    ground = c.ground_probe_vectors(data)['ground']
    assert c.addition_dark_probability(data, ground) == 0
    dark = n.gauge_invariant_dark_probe(data['neutral'])['vector']
    weight = c.addition_dark_probability(data, dark)
    assert weight.is_positive
    assert weight < s.Rational(8, 9)
    # Deleting the1614-dimensional outcome changes the complete law.
    assert weight != 0


def test_zero_singular_family_is_retained_with_its_zero_addition_line():
    data = c.charged_data(s.diag(1, 0, 0), 1, 0)
    zero_rows = [row for row in c.response_lines(data) if row['sigma_squared'] == 0]
    assert len(zero_rows) == 2
    for row in zero_rows:
        assert row['addition_lines'][0]['ground_referenced_energy'] == 0
        assert row['addition_lines'][0]['weight_per_normalized_mode'] == 1
        assert row['addition_lines'][1]['weight_per_normalized_mode'] == 0
        assert row['removal_weight_per_normalized_mode'] == 0


def test_source_coherent_control_keeps_its_own_coefficients():
    data = c.coherent_charged_data()
    rows = c.response_lines(data)
    assert sum(row['multiplicity'] for row in rows) == 18
    assert len(rows) == 4
    assert s.simplify(data['neutral']['g_squared']-(425+100*s.sqrt(2))/1944) == 0
    assert all(c._rational_simplify(row['removal_ground_referenced_energy']
                                  -(1+data['neutral']['gap'])/2) == 0 for row in rows)
