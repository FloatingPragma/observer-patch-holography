import pytest
import sympy as s
from . import dressing as d, symmetry as y


def test_full_joined_symmetry_and_its_baryon_domain():
    result = y.certificate()
    assert result['full_spatial_h_particle_hole_identity_exact']
    assert result['all_three_local_quantum_channels_mapped_exactly']
    assert result['global_native_singlet_and_B_zero_domain_preserved']
    assert not result['whole_global_native_singlet_domain_preserved_without_B_zero']


def test_complete_positive_electric_energy_geometry():
    rows = d.geometry()
    assert len(rows) == 13
    assert rows[1]['positive_electric_coefficient'] == 0
    assert all(row['positive_electric_coefficient'] > 0 for j, row in enumerate(rows) if j != 1)
    assert all(len(row['unit_longitudinal_field']) == 42 for row in rows)


def test_every_quark_line_translates_and_every_original_probability_survives():
    data = d.shifted_records()
    assert data['uc_minus_dc_corresponding_line_shift'] > 0
    parent = d.validate_parent()['exact_local_response']
    for family, record in data['channels'].items():
        assert record['pump_success'] == record['pump_failure'] == '1/2'
        assert record['outside_support_probability'] == '0'
        for old, new in zip(parent[family]['complete_positive_records'], record['complete_energy_records'], strict=True):
            assert new['original_probability_interval'] == old['probability_interval']
            assert new['energy_transfer_interval'][0]-s.Rational(old['source_excitation_interval'][0]) == record['positive_electric_energy']


def test_hypercharge_conservation_and_rejected_nonpositive_or_fitted_float_coupling():
    assert d.column_charge_controls()['all_thirty_elementary_charge_identities_exact']
    for coupling in (0, -1, 0.1):
        with pytest.raises(ValueError):
            d.shifted_records(coupling=coupling)
    with pytest.raises(ValueError):
        d.shifted_records(site=1)
