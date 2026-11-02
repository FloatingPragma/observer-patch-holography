"""Keep designs, finite gaps and physical predictions distinct."""
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import pytest

from tools import build_particle_mass_table as table


def test_scalar_and_spin_mass_addition_and_degeneracy():
    assert table.designed_masses((2, 0, 0, 1))['rest_energy_branches'] == ['1', '3']
    assert table.designed_masses((0, 0, 0, 1))['rest_energy_branches'] == ['1', '1']
    assert table.designed_masses((1, 0, 0, 1))['rest_energy_branches'] == ['0', '2']
    assert table.designed_masses((1, 0, 0, 1))['heavy_over_light_ratio'] is None
    for coefficients, expected in (
        ((2, 0, 0, 1), {'-3':4, '-1':4, '1':4, '3':4}),
        ((0, 0, 0, 1), {'-1':8, '1':8}),
        ((1, 0, 0, 1), {'-2':4, '0':8, '2':4}),
        ((0, 0, 0, 0), {'0':16}),
    ):
        levels = table.designed_masses(coefficients)['signed_zero_momentum_levels']
        assert {row['energy']:row['multiplicity'] for row in levels} == expected
        assert sum(row['multiplicity'] for row in levels) == 16


@pytest.mark.parametrize('ratio', [F(1), F(7, 2), F(206), F(3477)])
def test_arbitrary_ratios_are_designs_not_predictions(ratio):
    report = table.arbitrary_ratio_design(ratio)
    assert report['rest_energy_branches'] == ['1', str(ratio)]
    assert report['heavy_over_light_ratio'] == str(ratio)
    assert report['is_physical_prediction'] is False
    assert report['physical_species'] is None
    assert report['calibration_to_eV'] is None


def test_irrational_norm_is_not_rounded_into_an_exact_mass():
    report = table.designed_masses((1, 1, 1, 0))
    assert report['vector_norm_squared'] == '2'
    assert report['rest_energy_branches'] is None


def test_every_supplied_nonnegative_pair_has_an_exact_inverse_design():
    for a, b in ((F(0), F(0)), (F(3), F(0)), (F(7, 3), F(7, 3)), (F(13, 7), F(1, 9))):
        report = table.arbitrary_pair_design(a, b)
        assert report['rest_energy_branches'] == [str(b), str(a)]
        assert report['is_physical_prediction'] is False
    with pytest.raises(ValueError):
        table.arbitrary_pair_design(F(1), F(2))


@pytest.mark.parametrize('value', [True, 1., '1', None])
def test_nonexact_coefficient_inputs_rejected(value):
    with pytest.raises(ValueError):
        table.designed_masses((value, 0, 0, 1))


def test_malformed_fraction_cannot_smuggle_a_coefficient():
    bad = F(1, 2)
    bad._denominator = 0
    with pytest.raises(ValueError):
        table.designed_masses((bad, 0, 0, 1))


def test_outward_ratio_rounding_includes_exact_value():
    lo, hi = table.decimal_interval(F(1, 3), F(2, 3))
    assert F(lo) <= F(1, 3) < F(2, 3) <= F(hi)


def test_strict_json_rejects_duplicate_and_nan_claims():
    for text in ('{"mass":1,"mass":2}', '{"mass":NaN}'):
        with pytest.raises(ValueError):
            table.loads(text)


def test_source_manifest_detects_changed_bytes(tmp_path):
    pins = {}
    for name in table.INPUTS:
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(name.encode())
        pins[name] = hashlib.sha256(name.encode()).hexdigest()
    path = tmp_path/table.MANIFEST
    path.write_text(json.dumps({'schema': 'oph.particle_mass_table.sources.v1', 'reviewed_sources': pins}))
    (tmp_path/'tools/build_particle_mass_table.py').write_text('fixture')
    table.checked_snapshot(tmp_path)
    (tmp_path/table.INPUTS[0]).write_text('changed theorem')
    with pytest.raises(ValueError, match='input changed'):
        table.checked_snapshot(tmp_path)


def test_real_complete_table_scope_and_all_charged_lines():
    value = table.build()
    table.validate(value)
    assert value['physical_prediction_count'] == 0
    assert value['all_physical_particle_masses_derived'] is False
    assert len(value['centered_finite_CAR']['all_lines']) == 228
    assert value['ordered_finite_spin']['states'] == 256
    assert value['classes'][2]['rest_energy'] is None
    assert all(row['physical_prediction'] is False for row in value['classes'])


def test_finite_gap_promotion_and_zero_as_false_rejected():
    value = table.build()
    changed = copy.deepcopy(value)
    changed['centered_finite_CAR']['finite_gaps_are_relativistic_masses'] = True
    with pytest.raises(ValueError):
        table.validate(changed)
    changed = copy.deepcopy(value)
    changed['all_physical_particle_masses_derived'] = 0
    with pytest.raises(ValueError):
        table.validate(changed)


def test_supplied_coefficients_cannot_be_promoted_to_measured_species():
    value = table.build()
    value['coefficient_designs'][0]['physical_species'] = 'electron'
    value['physical_prediction_count'] = 1
    with pytest.raises(ValueError):
        table.validate(value)


def test_real_csv_disagreement_is_rejected():
    snapshot = table.checked_snapshot(table.ROOT)
    path = 'model/centered_particle_spectrum/lines.csv'
    snapshot[path] = snapshot[path].replace(b'0.36168271180597', b'0.36168271180596', 1)
    with pytest.raises(ValueError, match='CSV'):
        table.assemble(snapshot)


def test_in_flight_source_change_rejected(tmp_path, monkeypatch):
    path = tmp_path/'source'
    path.write_bytes(b'after')
    monkeypatch.setattr(table, 'checked_snapshot', lambda root: {'source': b'before'})
    monkeypatch.setattr(table, 'assemble', lambda data: {})
    with pytest.raises(ValueError, match='during build'):
        table.build(tmp_path)
