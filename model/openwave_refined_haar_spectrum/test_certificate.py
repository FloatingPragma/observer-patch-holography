from fractions import Fraction as F
import pytest
from .certificate import calculate,data,run
from model.openwave_colour_string.fixed_witness import spectral_enclosure


def test_full_refinement_gaps_and_all_local_correlations_are_positive():
    result=run()
    for record in result['records']:
        assert F(record['full_gap_interval_exact'][0])>0
        assert record['all_omitted_action_states']>0
        for field in record['observables']:
            assert F(field['variance_interval_exact'][0])>0
            for c in field['correlations']:
                lo,hi=map(F,c['centered_Euclidean_correlation_interval_exact']);assert 0<lo<=hi
    assert result['fine_interior_modes_are_retained']and not result['derived_running_coupling_or_continuum_limit']


def test_exact_matrix_mutation_cannot_reuse_a_complete_frame_certificate():
    record,_=calculate('coarse');H=data('coarse')['H'];mutated=H.copy();mutated[0,0]+=1
    old=record['complete_spectral_certificate'];new=spectral_enclosure(mutated,old)
    assert F(new['all_eigenvalue_error_upper'])>F(old['all_eigenvalue_error_upper'])+F(1,2)


def test_embedded_field_is_kept_outside_the_low_energy_frame():
    record,_=calculate('refined');field=record['observables'][1]
    assert field['loop_length']==8 and field['full_Wilson_trial_output_states']>201
    fine=record['observables'][0]
    assert F(field['normalized_first_spectral_moment_interval_exact'][0])>F(fine['normalized_first_spectral_moment_interval_exact'][1])
    with pytest.raises(ValueError):data('undeclared')
