"""Independent conversion identities, domain rejection and scale controls."""
from dataclasses import replace
from fractions import Fraction as F

import pytest

from . import units


def clock(gap=(F(1, 3), F(1, 2)), generator='fixture-H'):
    return units.CsClockInput(
        generator_id=generator, angular_gap_interval=gap,
        source_gap_provenance=('conditional exact-gap test premise',),
        physical_process_provenance=('conditional physical-process test premise',),
        same_time_provenance=('test premise: both gaps use the same H and tau',),
        physical_process=units.CESIUM_PROCESS)


def fractions(row, key):
    return tuple(map(F, row[key]))


def test_si_constants_and_exact_cesium_energy_are_independent_of_pi():
    expected = F(662607015, 10**42)*9192631770/F(1602176634, 10**19)
    assert units.CESIUM_TRANSITION_ENERGY_GEV == expected
    result = units.energy_conversion((1, 1), generator_id='fixture-H', clock=clock((1, 1)))
    assert fractions(result, 'energy_GeV_interval') == (expected, expected)
    assert result['absolute_scale_derived_by_this_conversion'] is False


def test_pi_enclosure_against_different_exact_machin_identity():
    # pi=8 atan(1/3)+4 atan(1/7), independent of the production identity.
    def atan(x, n=100):
        s = sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(n)), F())
        adjacent = s+(-1)**n*x**(2*n+1)/(2*n+1)
        return min(s, adjacent), max(s, adjacent)
    a, b = atan(F(1, 3)), atan(F(1, 7))
    reference = 8*a[0]+4*b[0], 8*a[1]+4*b[1]
    lo, hi = units.pi_interval()
    assert lo < reference[0] <= reference[1] < hi
    assert hi-lo < F(1, 10**34)
    hbar = units.hbar_gev_seconds_interval()
    assert F('6.582119569509065e-25') < hbar[0] <= hbar[1] < F('6.582119569509067e-25')


def test_interval_division_attains_both_extreme_energy_endpoints():
    result = units.energy_conversion((F(2, 5), F(7, 9)),
                                     generator_id='fixture-H', clock=clock())
    scale = units.CESIUM_TRANSITION_ENERGY_GEV
    assert fractions(result, 'energy_GeV_interval') == (scale*F(4, 5), scale*F(7, 3))
    for gap in (F(2, 5), F(7, 9)):
        for epsilon in (F(1, 3), F(1, 2)):
            assert F(result['energy_GeV_interval'][0]) <= scale*gap/epsilon <= F(result['energy_GeV_interval'][1])


def test_clock_time_and_energy_enclosures_have_the_correct_angular_factor():
    result = units.energy_conversion((2, 2), generator_id='fixture-H', clock=clock((3, 3)))
    tlo, thi = fractions(result, 'seconds_per_source_time_interval')
    hlo, hhi = units.hbar_gev_seconds_interval()
    elo, ehi = fractions(result, 'energy_GeV_interval')
    assert elo == ehi == units.CESIUM_TRANSITION_ENERGY_GEV*F(2, 3)
    assert hlo*2/thi <= elo <= hhi*2/tlo
    assert not (hlo*2/(thi*6) <= elo <= hhi*2/(tlo*6))


def test_light_speed_normalization_does_not_set_a_time_unit():
    row = units.energy_conversion((1, 2), generator_id='specified-H')
    assert row['status'] == 'UNRESOLVED_ABSOLUTE_SCALE'
    for key in ('energy_GeV_interval', 'seconds_per_source_time_interval',
                'energy_equivalent_kg_interval', 'unique_physical_mass_GeV'):
        assert row[key] is None
    assert row['length_relation'] == 'length_unit_m=c_SI*T/v_model'
    assert row['matched_physical_signal_speed_m_s'] == '299792458'


def test_rescaling_preserves_c_and_changes_every_nonzero_energy():
    t, lam, v = F(3, 17), F(11, 5), F(7, 3)
    base = units.c_normalized_units((t, t), v)
    changed = units.c_normalized_units((lam*t, lam*t), v)
    length = F(base['length_unit_m_interval'][0])
    assert length*v/t == units.LIGHT_SPEED_M_S
    assert fractions(changed, 'length_unit_m_interval') == (lam*length, lam*length)
    assert fractions(changed, 'energy_unit_GeV_interval') == tuple(
        x/lam for x in fractions(base, 'energy_unit_GeV_interval'))
    assert changed['light_speed_matching_fixes_absolute_scale'] is False


def test_cs_time_rescaling_cannot_be_hidden_in_a_different_generator_id():
    with pytest.raises(ValueError, match='different generator'):
        units.energy_conversion((1, 2), generator_id='radial-H', clock=clock())
    with pytest.raises(ValueError, match='typed'):
        units.energy_conversion((1, 2), generator_id='fixture-H', clock={'gap': [1, 2]})


def test_zero_gap_conversion_has_exact_zero_energy_even_with_clock_uncertainty():
    row = units.energy_conversion((0, 0), generator_id='fixture-H', clock=clock())
    assert row['energy_GeV_interval'] == ['0', '0']
    assert row['physical_rest_mass_established'] is False


@pytest.mark.parametrize('bad', [(0, 1), (-1, 2), (2, 1), (1,), (0.1, 1), (True, 1)])
def test_invalid_clock_gap_is_rejected(bad):
    with pytest.raises(ValueError):
        clock(bad)


@pytest.mark.parametrize('name', ['source_gap_provenance', 'physical_process_provenance', 'same_time_provenance'])
def test_missing_identification_or_common_time_premise_is_rejected(name):
    with pytest.raises(ValueError):
        replace(clock(), **{name: ()})


def test_another_physical_transition_cannot_inherit_the_cs_frequency():
    with pytest.raises(ValueError, match='specified cesium'):
        replace(clock(), physical_process='unspecified oscillator')


@pytest.mark.parametrize('bad', [0, -1, F(-1, 2), 1.0, True])
def test_invalid_model_signal_speed_is_rejected(bad):
    with pytest.raises(ValueError):
        units.energy_conversion((1, 2), generator_id='H', model_signal_speed=bad)


def test_negative_energy_gap_cannot_be_presented_as_rest_energy():
    with pytest.raises(ValueError):
        units.energy_conversion((-1, 1), generator_id='H')


def test_radial_certificate_retains_an_unresolved_positive_coefficient():
    row = units.certificate()
    result = row['radial_charged_gap_scale_family']
    assert F(result['energy_GeV_times_seconds_per_source_time_interval'][0]) > 0
    assert result['energy_GeV_interval'] is None
    exterior = row['radial_exterior_one_body_rest_energy_scale_family']
    assert exterior['model_gap_interval'] == ['64', '64']
    assert exterior['energy_GeV_interval'] is None
    assert exterior['generator_id'] != result['generator_id']
    assert row['localized_charged_gap_is_relativistic_rest_energy'] is False
    assert row['accepted_same_generator_physical_clock_inputs'] == []
    assert row['actual_absolute_particle_mass_prediction'] is False
