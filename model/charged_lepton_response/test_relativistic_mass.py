from fractions import Fraction as F
from functools import lru_cache
import json

import numpy as np
import pytest

from . import relativistic_mass as rm


def test_exact_all_left_clifford_pairing():
    value = rm.clifford_certificate()
    assert value['spatial_Clifford_checks'] == 9
    assert value['spin_pairing_cross_term_checks'] == 3
    assert value['massive_Dirac_spacetime_Clifford_checks'] == 16


@pytest.mark.parametrize('args', [(0, 1, 1), (1, -1, 1), (1, 1, 0), (True, 1, 1)])
def test_invalid_continuum_domains(args):
    with pytest.raises(ValueError): rm.continuum_bound(*args)


def test_compact_momentum_budget_controls_actual_dispersion():
    a, p, c = F(1, 16), F(2), F(3, 2)
    bound = F(rm.continuum_bound(a, p, c, 2, F(1, 100))['symbol_norm_error_upper'])
    for vec in [(2, -2, 2), (.1, .3, -1.7), (0, 0, 0)]:
        actual = float(c)*np.linalg.norm(np.sin(float(a)*np.array(vec))/float(a)-vec)
        assert actual <= float(bound)
    fine = rm.continuum_bound(a/2, p, c)
    assert F(fine['symbol_norm_error_upper']) == bound/4
    with pytest.raises(ValueError): rm.continuum_bound(a, p, c, mass_floor=0)


def test_positive_gap_example_is_scoped_to_its_actual_phase():
    value = rm.phase_zero_gap()
    assert F(value['smallest_charged_action_mass_interval'][0]) > F(1, 1000)
    assert value['program_phase'] == '0'
    assert value['applies_to_other_phases_without_a_separate_check'] is False


@lru_cache(maxsize=1)
def certified():
    return rm.certificate()


def test_same_native_field_spectrum_and_no_mirror_removal():
    value = certified()
    for row in value['finite_controls']:
        assert row['full_native_symbol_control']['complete_Nambu_dimension'] == 180
        assert row['full_native_symbol_control']['zero_action_components'] == 39
        assert row['all_eight_valley_energy_controls'] == 8
    assert value['scope']['controlled_free_relativistic_valley_limit'] is True
    assert value['scope']['interacting_gauge_field_continuum_proved'] is False
    assert value['scope']['physical_measured_mass_ratios'] == []
    json.dumps(value, allow_nan=False)
