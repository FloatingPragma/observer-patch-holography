from copy import deepcopy

import pytest

from . import packet


def test_complete_field_net_preserves_its_boundaries():
    value=packet._produce();packet.validate(value)
    for key in ('fixed_positive_U_interaction_survival_proved','Lorentz_covariance_proved',
                'physical_species_or_masses_identified','Einstein_G1_to_G6_completed'):
        changed=deepcopy(value);changed['theorem'][key]=True
        with pytest.raises(ValueError):packet.validate(changed)
    changed=deepcopy(value);changed['valleys']['valley_count']=1
    with pytest.raises(ValueError):packet.validate(changed)
    changed=deepcopy(value);changed['bounds']['inverse_weak_limit_interchanged']=True
    with pytest.raises(ValueError):packet.validate(changed)


def test_fresh_worker_and_in_process_payload_agree():
    assert packet.report()==packet._produce()
