from fractions import Fraction as F
import pytest
from .source import periodic_graph, premise_map, run


def test_periodic_reassembly_and_Gauss_are_complete():
    for side in (3, 4, 5):
        row = periodic_graph(side)
        assert row['full_Haar_links'] == row['closed_plaquettes'] == 3*side**3
        assert row['all_Gauss_and_star_support_columns_checked'] == 3*side**3
        assert row['every_electric_term_occurs_exactly'] == 4


def test_local_norm_has_no_volume_or_spacing_penalty():
    assert premise_map(F(8))['local_perturbation_operator_norm_exact'] == '1'
    assert premise_map(F(16))['local_perturbation_operator_norm_exact'] == '1/4'
    assert premise_map(F(32))['local_perturbation_operator_norm_exact'] == '1/16'
    assert premise_map(F(32), F(1, 32))['local_perturbation_operator_norm_exact'] == '1/16'
    with pytest.raises(ValueError):
        premise_map(F(0))
    with pytest.raises(ValueError):
        periodic_graph(2)


def test_existential_threshold_does_not_admit_numeric_or_Clay_claims():
    result = run()
    assert result['proved_conclusions']['above_threshold_finite_volume_gap_uniform_in_side']
    assert not any(result['scope_limits'].values())
    assert not any(row['numerical_coupling_admitted_by_uniform_gap_theorem']
                   for row in result['coupling_norms'])
    assert not result['separate_open_cube_control'][-1]['strict_positive_gap_certificate']
