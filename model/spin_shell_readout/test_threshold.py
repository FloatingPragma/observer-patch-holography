from copy import deepcopy
from fractions import Fraction as F

import pytest

from . import threshold as t


def test_nonatomic_edge_converges_with_vanishing_original_weight():
    widths, weights = [], []
    for k in range(1,10):
        row=t.non_atomic_control(k); b=row['bracket']
        assert b['threshold_lower_strict'] < row['threshold'] <= b['threshold_upper']
        assert row['atom_weight'] == 0
        assert row['original_success_weight']+row['complete_failure_weight'] == 1
        widths.append(b['threshold_upper']-b['threshold_lower_strict'])
        weights.append(row['original_weight_below_upper_probe'])
    assert all(a>b for a,b in zip(widths,widths[1:]))
    assert all(a>b for a,b in zip(weights,weights[1:]))


def test_actual_source_threshold_retains_every_failure_record():
    row=t.native_control(); source=row['actual_source_control']
    assert row['source_prepared_response_edge'] == F(9,8)
    assert sum(map(F,source['complete_record_probabilities'])) == 1
    assert source['physical_particle_mass'] is None
    assert row['complete_original_probability_bracket']['threshold_lower_strict'] == F(1123,1000)


def test_complete_scan_counts_and_error_are_never_conditioned_on_success():
    row=t.sampled_cdf(500000, 1000, [0,10], F(1,100), F(1,100), F(1,1000))
    assert row['required_trials'] == 500000
    assert row['preparation_failure_count'] == 499000
    assert row['success_above_last_threshold_count'] == 990
    assert row['original_CDF_intervals'][1][1] == F(10,500000)+F(11,1000)
    with pytest.raises(ValueError,match='too few'):
        t.sampled_cdf(499999,1000,[0,10],F(1,100),F(1,100))


def test_strict_tail_evidence_is_needed_even_when_lower_bin_is_empty():
    with pytest.raises(ValueError,match='insufficient'):
        t.bracket(0,(0,F(1,10)),1,(F(1,10),F(1,5)),1,F(1,10),1)
    with pytest.raises(ValueError,match='insufficient'):
        t.bracket(0,(0,0),1,(0,F(1,5)),1,F(1,10),1)


@pytest.mark.parametrize('power',[1,2,3,5])
def test_power_cdf_certifies_true_edge_without_atom(power):
    theta, weight, eps=F(-1,3),F(1,17),F(1,32)
    error=weight*eps**power/3
    row=t.bracket(theta-eps,(0,error),theta+eps,
        (weight*eps**power-error,weight*eps**power+error),eps,weight,power)
    assert row['threshold_lower_strict'] < theta <= row['threshold_upper']


def test_gapless_source_event_obstructs_only_complete_gapped_charge_map():
    for gap,conversion in [(F(1),F(1)),(F(1,10**30),F(7)),(F(10**30),F(1,3))]:
        row=t.gapless_charge_obstruction(gap,conversion)
        assert 0 < row['target_energy_upper_of_source_event'] < gap
        assert 0 < row['complete_dictionary_half_trace_error_lower'] < 1
        assert row['electromagnetic_charge_identification_proved'] is False
        assert row['massive_subbranches_or_different_sources_excluded'] is False


@pytest.mark.parametrize('call',[
    lambda:t.bracket(0,(0,0),1,(F(1,2),1),1,2,1),
    lambda:t.bracket(0,(0,0),1,(F(1,2),1),1,1,True),
    lambda:t.bracket(0,(0,0),1,(F(1,2),1),.1,1,1),
    lambda:t.sampled_cdf(100,2,[2,1],F(1,2),F(1,2)),
    lambda:t.sampled_cdf(100,2,[3],F(1,2),F(1,2)),
    lambda:t.sampled_cdf(True,1,[1],F(1,2),F(1,2)),
    lambda:t.gapless_charge_obstruction(0),
    lambda:t.gapless_charge_obstruction(1,0),
])
def test_invalid_and_scientifically_insufficient_inputs_fail(call):
    with pytest.raises(ValueError):call()


def test_threshold_receipt_promotion_is_rejected():
    from . import packet
    expected=packet.produce(); changed=deepcopy(expected)
    changed['infrared_threshold_readout']['continuum_infraparticle_or_measured_mass_identified']=True
    with pytest.raises(ArithmeticError):packet.validate(changed,expected)
