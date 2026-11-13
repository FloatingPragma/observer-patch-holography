from fractions import Fraction as F
import pytest
from .response import run, reflection_form, clock_scaling, dyadic_outer, exponential_decay_interval
from model.openwave_coupled_haar_spectrum.certificate import calculate


def test_full_source_response_and_clock_obstruction():
    result = run()
    assert result['proved']['nonzero_positive_time_Wilson_response']
    for record in result['representation_refinement_records']:
        assert not record['source_Haar_representations_truncated']
        assert len(record['all_finite_Galerkin_spectral_weights_exact']) in (55, 140, 285)
        assert F(record['finite_positive_reflection_form_interval_exact'][0]) >= 0
    rows = result['fixed_coupling_spacing_trajectory']
    gaps = [tuple(map(F, row['physical_clock_gap_interval_exact'])) for row in rows]
    assert all(a[0] <= b[0]/2 <= a[1] for a,b in zip(gaps,gaps[1:]))
    assert all(F(row['zero_time_variance_lower_exact']) > 0 for row in rows)
    assert F(rows[-1]['connected_correlation_upper_exact']) < F(1,10**10)
    assert not result['continuum_acceptance']['fixed_graph_fixed_coupling_trajectory_admitted']
    assert not result['scope']['Clay_problem_solved']


def test_reflection_contract_rejects_missing_hypotheses():
    for weights, energies, times, coefficients in (([-1],[1],[0],[1]),
            ([1],[-1],[0],[1]), ([1],[1],[-1],[1]), ([1],[1],[0],[]),
            ([1,1],[1],[0],[1])):
        with pytest.raises(ValueError):
            reflection_form(weights,energies,times,coefficients)
    lo, hi = reflection_form([F(1,2),F(1,2)],[1,2],[0,F(1,4)],[1,-1])
    assert 0 < lo <= hi
    with pytest.raises(ValueError):
        clock_scaling(calculate('shared_g2_2_cap4')[0],[0])


def test_exact_outer_rounding_and_exponential_tail():
    value = F(-123456789,987654321)
    lo, hi = dyadic_outer((value,value))
    assert lo <= value <= hi and hi-lo <= F(1,2**80)
    assert exponential_decay_interval(F(0)) == (F(1), F(1))
    lo, hi = exponential_decay_interval(F(65))
    assert lo == 0 and 0 < hi < F(1,10**27)
