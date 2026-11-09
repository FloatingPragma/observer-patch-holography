import copy
from fractions import Fraction as F
import json

import pytest

from . import source_mass_readout as mass


@pytest.fixture
def source():
    return json.loads((mass.flavor.OUTPUT).read_text())


def test_all_policies_original_coefficients_and_ratios(source):
    result = mass.summarize(source)
    for out, original in zip(result['policies'], source['current_gated_joint_spectrum']['policies']):
        for row, old in zip(out['channels'], original['channels']):
            assert row['ascending_dimensionless_action_mass_intervals'] == old['dimensionless_singular_values']
            boxes = [[F(x) for x in pair] for pair in old['dimensionless_singular_values']]
            if boxes[0][0] > 0:
                for pair, val in zip(row['ratios_to_smallest_interval'][1:], boxes[1:]):
                    assert F(pair[0]) == val[0]/boxes[0][1]
                    assert F(pair[1]) == val[1]/boxes[0][0]
            else:
                assert row['ratios_to_smallest_interval'] is None
    assert result['physical_pole_mass_table'] == []
    assert result['empirical_mass_inputs_consumed'] == []


def test_common_scale_cancels_in_ratios(source):
    original = mass.summarize(source)
    scaled = copy.deepcopy(source)
    for policy in scaled['current_gated_joint_spectrum']['policies']:
        for channel in policy['channels']:
            channel['dimensionless_singular_values'] = [[str(7*F(v)) for v in pair]
                                                      for pair in channel['dimensionless_singular_values']]
    changed = mass.summarize(scaled)
    for p, q in zip(original['policies'], changed['policies']):
        for a, b in zip(p['channels'], q['channels']):
            assert a['ratios_to_smallest_interval'] == b['ratios_to_smallest_interval']


def test_missing_policy_rejected(source):
    source['current_gated_joint_spectrum']['policies'].pop()
    with pytest.raises(ValueError, match='four'):
        mass.summarize(source)


@pytest.mark.parametrize('replacement', [True, 0, 'False'])
def test_unsupported_pole_promotion_rejected(source, replacement):
    source['scope']['physical_pole_masses_or_absolute_units_derived'] = replacement
    with pytest.raises(ValueError, match='scope'):
        mass.summarize(source)


def test_balanced_Koide_and_same_source_exclusion(source):
    import sympy as s
    # Positive square-root masses 1,1,4+sqrt(18) have exactly Q=2/3.
    z = 4+s.sqrt(18)
    assert s.simplify((2+z*z)/(2+z)**2) == s.Rational(2, 3)
    # Independent rational control (0,1,(2+sqrt(3))^2) is balanced;
    # outward rational enclosures must retain the balance point.
    r = mass.iv.sqrt_rational(F(3))
    upper_mass = (mass.iv.I(2)+r)**2
    bounded = mass.koide_interval([(F(),F()),(F(1),F(1)),(upper_mass.lo,upper_mass.hi)])
    assert not bounded['balanced_value_two_thirds_excluded']
    for policy in mass.summarize(source)['policies']:
        lepton = next(x for x in policy['channels'] if x['channel']=='charged_lepton')
        assert lepton['Koide_diagnostic']['balanced_value_two_thirds_excluded']


def test_duplicate_channel_rejected(source):
    channels = source['current_gated_joint_spectrum']['policies'][0]['channels']
    channels.append(copy.deepcopy(channels[0]))
    with pytest.raises(ValueError, match='channels'):
        mass.summarize(source)


def test_receipt_reuse_is_explicit_and_does_not_claim_cold_replay(monkeypatch):
    def forbidden():
        raise AssertionError('ordinary integration must identify receipt reuse explicitly')
    monkeypatch.setattr(mass.flavor, 'report', forbidden)
    value = mass.certificate()
    assert value['owning_source_receipt_authenticated'] is True
    assert value['owning_source_cold_scientific_replay_performed_here'] is False


def test_stale_parent_receipt_rejected(source, monkeypatch):
    source['source_bindings'] = {}
    monkeypatch.setattr(mass.strict_json, 'load', lambda path: source)
    monkeypatch.setattr(mass.flavor, 'source_bindings', lambda: {'expected': 'hash'})
    with pytest.raises(ValueError, match='stale'):
        mass.certificate()
