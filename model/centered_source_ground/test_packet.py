from copy import deepcopy
import pytest
from . import packet


def test_cold_replay_matches_direct_source():
    assert packet.report()==packet._produce()


def test_centered_reference_cannot_be_reported_as_the_old_one():
    value=packet._produce();packet.validate(value)
    for key in ('old_ground_reused_without_check','interacting_Lorentz_continuum_proved','physical_mass_or_species_identified'):
        wrong=deepcopy(value);wrong['bounds'][key]=True
        with pytest.raises(ValueError):packet.validate(wrong)
    wrong=deepcopy(value);wrong['cube']['bounds']['quantified_neutral_gap_proved']=True
    with pytest.raises(ValueError):packet.validate(wrong)
    for key,changed in (('ground_multiplicity',5),('bipartition_sizes',[3,5]),
                        ('proof_imported_not_reproved',False)):
        wrong=deepcopy(value);wrong['cube']['imported_Lieb_theorem'][key]=changed
        with pytest.raises(ValueError):packet.validate(wrong)
    for key in ('infinite_ground_local_faithfulness_inferred','geometric_modular_charge_identified'):
        wrong=deepcopy(value);wrong['cube_local_state'][key]=True
        with pytest.raises(ValueError):packet.validate(wrong)
