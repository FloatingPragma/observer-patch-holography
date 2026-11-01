from copy import deepcopy

import pytest

from . import packet


def test_cold_replay_equals_in_process_proof():
    assert packet.report()==packet._produce()


def test_microscopic_ground_limit_does_not_promote_continuum_or_masses():
    value=packet._produce();packet.validate(value)
    for key in ('second_moment_equality_passed_to_limit','macroscopic_smeared_field_Gaussianity_inferred',
                'chemical_potential_shift_added','physical_mass_or_species_identified'):
        wrong=deepcopy(value);wrong['constraints'][key]=True
        with pytest.raises(ValueError):packet.validate(wrong)
