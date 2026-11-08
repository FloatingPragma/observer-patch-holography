from fractions import Fraction as F
from itertools import islice

import numpy as np
import pytest

from model.native_ordered_ground import source
from . import field_gap


def test_actual_full_source_Hadamard_maps_and_normalization():
    value = field_gap.source_certificate()
    assert [row['all_matrix_rows'] for row in value['finite_source_matrix_controls']] == [2,4,8,16,256]
    assert field_gap.gap(F(1,4)) == F(1,32)
    assert value['source_couplings_after_global_Hadamard']['h'] == 'B/2'
    assert value['gap_theorem_is_an_explicit_external_mathematical_input']
    assert not value['zero_field_gap_or_ground_marginal_modulus_claimed']


def test_true_small_source_gaps_dominate_the_uniform_bound_along_path():
    # Independent numerical controls of the source; not a proof of the gap.
    for sites,edges in ((1,()),(2,((0,1),)),(3,((0,1),(1,2)))):
        for b in (F(1,100),F(1,4),F(3)):
            for j in (F(1,100),F(1,2),F(2)):
                h,scale = source.matrix(sites,edges,j,b)
                spectrum = np.linalg.eigvalsh(h.astype(float)/scale)
                assert spectrum[1]-spectrum[0] >= float(field_gap.gap(b))-1e-12


def test_adiabatic_source_budget_uses_no_ground_state_or_energy_oracle():
    value = field_gap.resource_plan(2,((0,1),))
    assert value['adiabatic_coefficient'] == 231424
    assert value['source_clock_duration'] == 69427200
    assert value['total_prepared_output_half_trace_error_upper'] <= F(1,100)
    assert value['all_first_failure_records_retained']
    assert not value['ground_postselection_used']
    assert not value['infinite_volume_or_zero_field_accuracy_included']
    assert not value['large_ramp_executed']
    tighter = field_gap.resource_plan(2,((0,1),),error=F(1,1000))
    assert tighter['source_clock_duration'] == 10*value['source_clock_duration']
    assert tighter['native_pulses_upper'] > value['native_pulses_upper']


def test_no_bond_product_start_and_invalid_field_are_handled():
    value = field_gap.resource_plan(1,())
    assert value['source_clock_duration'] == value['piecewise_constant_slices'] == 0
    assert value['native_pulses_upper'] == 16
    assert len(list(field_gap.preparation_word(1,(),F(1),F(1,4),F(),0))) == 16
    with pytest.raises(ValueError): field_gap.gap(0)
    with pytest.raises(ValueError): field_gap.resource_plan(3,((0,1),(1,2),(0,2)))
    with pytest.raises(ValueError): field_gap.resource_plan(2,((0,1),),error=F(1))
    with pytest.raises(ValueError): list(field_gap.preparation_word(2,((0,1),),F(1),F(1,4),F(1),0))


def test_literal_ramp_gate_count_and_refined_source_commands():
    word = list(field_gap.preparation_word(2,((0,1),),F(1),F(1,4),F(1,16),2))
    assert len(word) == 16*2+2+62+2*11*2
    commands = list(islice(field_gap.dyadic_preparation_word(
        2,((0,1),),F(1),F(1,4),F(1,16),2,24),2))
    assert len(commands) == 2
    for row in commands:
        command = row['port_command']
        assert len(command['integrated_port_coefficients']) == 12
        assert F(command['generator_and_complete_output_error']) <= 29*F(1,1<<24)


def test_short_native_ramp_keeps_full_rows_and_all_failures():
    value = field_gap.small_native_execution()
    assert value['all_work_rows'] == 16
    assert value['pulses'] == 140
    assert value['all_failure_records_retained']
    assert not value['short_run_is_a_certified_ground_preparation']
    assert not value['long_certified_adiabatic_budget_executed']
