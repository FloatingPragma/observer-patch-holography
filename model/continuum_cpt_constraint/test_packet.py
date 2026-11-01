from copy import deepcopy

import pytest

from . import packet


def test_cold_replay_equals_direct_proof():
    assert packet.report()==packet._produce()


def test_CPT_dictionary_exclusion_does_not_promote_global_no_go():
    value=packet._produce();packet.validate(value)
    for key in packet.BOUNDARIES:
        wrong=deepcopy(value);wrong['constraints'][key]=True
        with pytest.raises(ValueError):packet.validate(wrong)


def test_zero_mean_does_not_promote_ground_positivity():
    wrong=packet._produce()
    wrong['controls']['centering_failure']['same_reference_is_centered_ground']=True
    with pytest.raises(ValueError):packet.validate(wrong)
