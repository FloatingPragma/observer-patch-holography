from copy import deepcopy
import pytest
from . import packet


def test_theorem_preserves_its_scope_boundaries():
    value=packet._produce();packet.validate(value)
    for key in ('positive_energy_time_group_on_entire_enlarged_GNS_proved',
                'composite_limiting_time_continuity_proved','bounded_frequency_memory_survival_proved',
                'Hartree_free_closure_proved','statistical_independence_or_Gaussianity_proved',
                'physical_species_or_masses_identified','Einstein_G1_to_G6_completed'):
        changed=deepcopy(value);changed['theorem'][key]=True
        with pytest.raises(ValueError):packet.validate(changed)
    changed=deepcopy(value);changed['theorem']['valleys_in_each_family']=1
    with pytest.raises(ValueError):packet.validate(changed)


def test_cold_worker_reproduces_exact_payload():
    assert packet.report()==packet._produce()
