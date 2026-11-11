import numpy as np
from scipy import linalg as la
import sympy as s
import pytest

from . import bridge as b
from . import source as m


def test_full_native_polynomial_and_all_contact_components():
    result = b.local_operator_checks()
    assert result['all_native_global_current_commutators'] == 12
    assert result['all_fermion_pair_contacts_checked'] == 900
    assert result['nonzero_oriented_native_contacts'] > 8


def test_exact_lepton_control_and_failed_transfer_at_unit_quark_channels():
    result = b.lepton_control_embedding()
    assert result['complete_intertwined_source_columns'] == 320
    assert result['all_unit_channel_quark_creation_norm_squared_on_blank'] > 0
    assert not result['empty_quarks_reduce_all_three_unit_channel_member']


def test_matrix_free_columns_equal_full_local_polynomial_on_arbitrary_occupations():
    source = m.one_cell()
    polynomial = m.one_cell_polynomial()
    masks = (0, 1, 1 << 29, (1 << 30)-1, 0x15555555, 0x2aaaaaaa, 193739)
    for mask in masks:
        for scalar in range(5):
            expected = {key: complex(value) for key, value in m.car.column(polynomial, mask, scalar).items()}
            assert source.column(mask, scalar) == expected


def test_actual_multisite_full_CAR_columns_are_Hermitian_and_keep_scalar_transport():
    source = m.Source(((0, 0, -1), (0, 0, 0), (0, 0, 1)))
    assert source.modes == 90 and source.scalar_dimension == 13
    assert source.scalar_root_residual < 1e-9
    seeds = [(0, 0), ((1 << 90)-1, 0), (1 << 9, source.scalar_label(0, 1)),
             ((1 << 3) | (1 << 44), source.scalar_label(2, 4))]
    for key in seeds:
        for target, value in source.column(*key).items():
            assert abs(source.column(*target).get(key, 0)-np.conjugate(value)) < 1e-11
    scalar = source.scalar_label(0, 1)
    vacuum_transport = source.column(0, scalar)
    assert abs(vacuum_transport[(0, source.scalar_label(2, 1))]) > 0


def test_same_compact_scalar_and_fermion_sources_and_all_mode_contacts():
    result = b.spatial_contact()
    assert result['source']['CAR_modes'] == 10290
    assert result['source']['scalar_dimension'] == 1373
    assert len(result['complete_mode_records']) == 343
    assert result['raw_mode_kernel'].shape == (343, 4, 4)
    assert result['fermion_band_metric_error'] < 1e-10
    for j in (0, 1, 79, 342):
        field = result['scalar_modes'][:, j]
        value = m.compact_source()
        # Independent literal microscopic spin/site pair and full band frame.
        epsilon = np.array([[0., 1.], [-1., 0.]])
        local = la.block_diag(*(phase*f*epsilon for phase, f in zip(value.stagger, field, strict=True)))
        expected = result['W'].conj().T@local@result['W'].conj()
        assert la.norm(expected-result['raw_mode_kernel'][j]) < 1e-10
    assert result['complete_mode_records'][0]['complete_exterior_leakage'] > 0
    assert not result['complete_band_is_reducing_under_local_quantum_Higgs_interaction']


def test_one_state_reference_clock_contract_retains_all_branches_and_boundaries():
    result = b.reference_and_records()
    assert not result['local_gauge_field_or_local_Gauss_constraint_supplied']
    assert not result['charged_excitation_energies_presumed_positive']
    assert not result['one_cell_lepton_normalization_reused']
    assert not result['free_scalar_vacuum_normalization_reused_for_interacting_ground']


def test_exact_full_source_scalar_pulse_propagates_to_distant_pair():
    from .encounter import execute
    result = execute()
    assert result['source']['CAR_modes'] == 90
    assert result['all_three_unit_channels_executed']
    assert result['receiver_leading_power'] == 2
    assert result['maximum_fermion_occupation_in_each_exact_power'][-1] == 6
    expected = 8*(s.Rational(2, 3)-s.sqrt(2)/2)**2
    assert s.simplify(s.sympify(result['exact_receiver_leading_second_power_norm_squared'])-expected) == 0
    lower, upper = result['receiver_probability_bounds_numeric']
    assert 0 < lower < result['receiver_Taylor_probability'] < upper
    assert not result['preparation_is_interacting_ground']
    assert not result['full_spatial_ground_or_rest_mass_inferred']


def test_exact_pulse_refuses_float_coefficients_and_inexact_source_custody(monkeypatch):
    from .encounter import scaled_coefficient
    with pytest.raises(ArithmeticError):
        scaled_coefficient(s.Float(0.5))
    original = b.ROOT
    monkeypatch.setattr(b, 'ROOT', original/'nonexistent_source')
    with pytest.raises(FileNotFoundError, match='missing declared joined-source input'):
        b.bindings()


def test_complete_operator_packet_rejects_missing_and_truncated_scalar_modes():
    from .__main__ import check_archive, equal
    contact = b.spatial_contact()
    names = ('W', 'site_kernel', 'scalar_modes', 'raw_mode_kernel', 'canonical_mode_kernel')
    packet = {name: contact[name] for name in names}
    check_archive(packet, contact)
    with pytest.raises(ArithmeticError, match='inventory'):
        check_archive({name: value for name, value in packet.items() if name != 'scalar_modes'}, contact)
    with pytest.raises(ArithmeticError, match='shape'):
        check_archive({**packet, 'raw_mode_kernel': packet['raw_mode_kernel'][:-1]}, contact)
    assert not equal(True, 1)
    assert not equal({'receiver_probability': 3e-19}, {'receiver_probability': 0.0})


def test_bosonic_spectrum_is_complete_on_its_actual_free_reference():
    result = b.bosonic_source_content()
    assert len(result['all_scalar_suboperator_frequencies']) == 343
    assert result['complete_scalar_suboperator_rank'] == 1373
    assert result['scalar_complex_field_coordinates'] == 686
    assert result['dynamical_gauge_link_or_field_operator_count'] == 0
    assert not result['free_reference_is_joint_interacting_ground']
    assert not result['identical_clock_reference_embedding_of_quantum_link_source_proved']
    assert not result['W_Z_photon_gluon_mass_records_supplied']
