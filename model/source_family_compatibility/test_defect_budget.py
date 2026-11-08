from fractions import Fraction as F

import pytest
import sympy as s

from . import defect_budget as d


@pytest.mark.parametrize('theta', [F(0), F(1, 3), F(1)])
def test_complete_Gram_on_every_Fock_matrix_unit(theta):
    gram, q = d.seam_gram(2, (0, 1), theta)
    assert gram[0, 0] == 0  # Empty endpoints also hold without a distinguishing record.
    assert gram[1, 1] == 1+theta
    assert gram[3, 3] == 2  # Collision remains in the complete archive.
    assert gram[5, 5] == 0
    for a in range(16):
        for b in range(16):
            assert gram[a, b] == (gram[a, a] if a == b else 0)
    assert max(gram[j, j]-2*q[j] for j in range(16)) == 0


def test_coherent_reference_cannot_be_replaced_by_defect_postselection():
    gram, _ = d.seam_gram(2, (0, 1), F(1, 2))
    # Exact source/reference vector with coherent empty, single and double occupations.
    psi = s.zeros(32, 1)
    psi[0], psi[2], psi[7], psi[11] = s.Rational(1, 2), s.I/2, -s.Rational(1, 2), s.I/2
    value = (psi.H*s.kronecker_product(gram, s.eye(2))*psi)[0]
    counts = d.defect_counts(2)
    budget = (psi.H*s.kronecker_product(s.diag(*counts), s.eye(2))*psi)[0]
    assert value == s.Rational(7, 8)
    assert 0 < value <= 2*budget


def test_exact_source_action_and_interacting_leakage():
    control = d.finite_controls()
    assert [x['whole_Fock_dimension'] for x in control['whole_source_controls']] == [16, 64]
    row = control['actual_interacting_time_control']
    assert row['defect_expectation'] == '9/25'
    assert row['doublon_hole_probability'] == '9/50'
    assert row['complete_dilation_difference_squared'] == ['9/25']*3
    assert not row['single_occupation_invariance_assumed']


def test_local_hopping_constant_is_actual_c_not_half_c():
    _, hopping, _ = d.exact_source(2, ((0, 1, -1),), F(3, 2), F(4))
    assert {F(x) for x in hopping.eigenvals()} == {F(-3, 2), F(-3, 4), 0, F(3, 4), F(3, 2)}


def test_degree_bound_does_not_assume_global_strong_projection_regime():
    row = d.energy_budget(100, 100, 2, 1, 16, 0)
    assert 16 < 4*100
    assert row['local_form_bound'] == F(200, 9)
    assert row['defect_expectation_upper'] == F(25, 2)
    assert not row['ground_sector_or_projection_assumed']


def test_sharp_degree_budget_improves_volume_norm_and_controls_full_attempts():
    row = d.uniform_history(8, 12, 3, 1, 1024, 0, 4)
    assert row['local_form_bound'] < row['global_hopping_bound']
    assert row['complete_trace_bound_squared'] < F(1, 50)
    assert row['bound_controls_entire_source_reference_archive']
    assert row['ideal_energy_preservation_or_work_budget_required']


@pytest.mark.parametrize('args', [(8,12,3,1,0,0), (8,12,3,1,16,-1),
                                 (8,12,8,1,16,0), (8,13,3,1,16,0)])
def test_invalid_budget_domains_rejected(args):
    with pytest.raises(ValueError): d.energy_budget(*args)


def test_increased_work_budget_cannot_improve_defect_bound():
    assert d.energy_budget(8,12,3,1,256,0)['defect_expectation_upper'] < d.energy_budget(8,12,3,1,256,1)['defect_expectation_upper']


def test_certificate_keeps_physical_selection_and_limit_unproved():
    result = d.certificate()
    assert not result['capacity_action_or_physical_clock_selected']
    assert not result['interacting_relativistic_particle_limit_proved']
