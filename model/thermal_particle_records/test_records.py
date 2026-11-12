from fractions import Fraction as F
import math

from .records import (inputs,boson_rotation,gauge_block_control,color_reference_control,
                      binomial_records,count_law)


def test_exact_partition_and_complete_actual_fock_rotation():
    groups=inputs()[3]
    f=[g for g in groups if g['source_sector']=='fermion'];b=[g for g in groups if g['source_sector']=='boson']
    assert len(f)==84 and len(b)==28
    assert math.prod(g['projector_rank']+1 for g in f)==2**132
    assert math.prod(sum(math.comb(g['projector_rank'],n) for n in range(g['projector_rank']+1)) for g in f)==2**180
    rotation=boson_rotation()
    assert rotation['matrix'].shape==(1653,1653) and rotation['unitary_error']<2e-12
    assert len(gauge_block_control())==52


def test_binomial_original_records_include_numerical_underflow():
    rows=binomial_records(3,170.,.1)
    assert len(rows)==4 and sum(r['projector_rank'] for r in rows)==8
    assert rows[1]['positive_probability_underflowed']
    assert all(math.isfinite(r['log_probability']) for r in rows)
    assert [r['count'] for r in rows]==[0,1,2,3]


def test_count_measurement_preserves_entangled_color_coherence():
    for rank,expected in [(3,F(4,3)),(8,F(7,4))]:
        result=color_reference_control(rank)
        assert F(result['fine_color_dephasing_trace_distance_exact'])==expected
        assert result['group_count_one_Lueders_output_unchanged']


def test_full_finite_thermal_counts_and_initial_detector_complement():
    cold=count_law(.1);warm=count_law(10.)
    for row in (cold,warm):
        records=row['boson_joint_group_count_records']
        assert len(records)==435 and sum(r['projector_rank'] for r in records)==1653
        assert sorted(i for r in records for i in r['original_occupation_indices'])==list(range(1653))
        assert abs(row['boson_probability_sum']-1)<2e-14
        assert abs(sum(x['fraction_of_expected_counts'] for x in row['species_expected_counts'])-1)<2e-12
        assert row['pre_probe_one_excitation_detector']['evaluated_before_any_probe']
        assert not row['pre_probe_one_excitation_detector']['claimed_to_be_the570_post_probe_record_law']
    assert warm['pre_probe_one_excitation_detector']['many_excitation_complement_probability']>cold['pre_probe_one_excitation_detector']['many_excitation_complement_probability']
