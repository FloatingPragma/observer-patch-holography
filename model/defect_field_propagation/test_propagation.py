from copy import deepcopy
import math

import numpy as np
import pytest

from . import propagation as p


def test_complete_native_geometry_and_exact_component_partition():
    assert p.geometry()['cycle'] == list(p.CYCLE)
    parts = p.components()
    assert sorted(i for row in parts for i in row) == list(range(90))
    assert sorted(map(len, parts)) == [1]*6+[2]*24+[6]*6
    owner = {i: k for k, row in enumerate(parts) for i in row}
    for i, j in zip(*np.nonzero(p.local_pairing()), strict=True):
        assert owner[i] == owner[j]


def test_ring_derivative_fourier_sign_and_both_doublers():
    D = p.derivative()
    for m in range(12):
        k = 2*math.pi*m/12
        v = np.exp(1j*k*np.arange(12))/math.sqrt(12)
        np.testing.assert_allclose(D@v, math.sin(k)*v, atol=2e-15)
    assert np.count_nonzero(abs(np.linalg.eigvalsh(D)) < 1e-12) == 2
    record, spectrum = p.dispersion()
    assert spectrum.shape == (12, 180)
    assert np.count_nonzero(abs(spectrum) < 1e-10) == 24
    assert record['all_doublers_retained']


def test_component_ring_has_same_complete_momentum_spectrum():
    indices = p.components()[0]
    h, delta = p.block(indices)
    actual = np.linalg.eigvalsh(p.nambu(h, delta))
    mass = abs(p.local_pairing()[indices[0], indices[1]])
    energies = [math.sqrt(math.sin(2*math.pi*m/12)**2+mass**2) for m in range(12)]
    expected = sorted(sign*energy for energy in energies for sign in (-1, 1) for _ in (0, 1))
    np.testing.assert_allclose(actual, expected, atol=2e-11, rtol=1e-12)


def test_source_symmetric_product_converges_to_same_mass_and_kinetic_sum():
    indices = min((r for r in p.components() if len(r) == 2),
                  key=lambda r: abs(p.local_pairing()[r[0], r[1]]))
    errors = []
    for steps in (1, 2, 4):
        source, record = p.execute_block(indices, 1, duration=.4, steps=steps)
        errors.append(record['sum_action_nambu_error'])
        assert record['effective_duration'] == .5
        assert source.resources['retained_work_fermion_modes'] == 24
        assert source.resources['routing_fswap_calls'] > 0
        assert all(item['diagnostics']['contracted_fock_columns'] == 4
                   for item in source.local_catalog.values())
    assert errors[1] < .3*errors[0] and errors[2] < .3*errors[1]


def test_nonzero_same_port_collision_cannot_be_discarded():
    # Independent four-mode CAR oracle: modes are (site0 a,b, site1 a,b).
    def annihilate(mask, mode):
        if not mask >> mode & 1:
            return None
        return mask ^ (1 << mode), (-1)**((mask & ((1 << mode)-1)).bit_count())
    initial = (1 << 0) | (1 << 3)
    mid, sign = annihilate(initial, 0)
    creation_sign = (-1)**((mid & ((1 << 2)-1)).bit_count())
    final = mid | (1 << 2)
    assert final == (1 << 2) | (1 << 3)
    assert sign*creation_sign == 1
    # A nonzero native-ring hopping matrix element makes this an actual
    # collision transition of the proposed action, not a hypothetical mode.
    h, _ = p.block(p.components()[0])
    assert h[2, 0] != 0


def test_local_catalog_rejects_scalar_phase_even_if_majorana_map_unchanged():
    result, _ = p.execute_block(p.components()[0])
    changed = deepcopy(result.local_catalog)
    row = next(iter(changed.values()))['source_fock_unitary']
    row['real'] = (-np.array(row['real'])).tolist()
    row['imag'] = (-np.array(row['imag'])).tolist()
    with pytest.raises(ArithmeticError):
        p.audit_local_catalog(changed)


def test_analytic_product_bounds_cover_independent_one_body_error():
    result, record = p.execute_block(p.components()[0], duration=.001, steps=2)
    assert record['source_product_full_Fock_operator_bound'] > 0
    assert 0 < record['ideal_product_to_static_sum_full_Fock_operator_bound'] <= 2
    # A full-Fock operator bound controls the Heisenberg CAR map with factor2.
    assert record['sum_action_nambu_error'] <= 2*(
        record['ideal_product_to_static_sum_full_Fock_operator_bound']
        +record['source_product_full_Fock_operator_bound'])


def test_kinetic_does_not_preserve_static_onsite_positive_energy_chart():
    indices = p.components()[0]
    h, delta = p.block(indices)
    mass = p.nambu(np.zeros_like(h), delta)
    kinetic = p.nambu(h, np.zeros_like(delta))
    assert np.linalg.norm(mass@kinetic-kinetic@mass) > 1e-5


@pytest.fixture(scope='module')
def wave():
    return p.wavepacket()


def test_source_wavepacket_and_current_keep_all_purifiers_and_pointer_outcomes(wave):
    record, arrays = wave
    assert record['retained_work_modes'] == record['retained_purifying_modes'] == 24
    assert abs(arrays['occupation_excess_profile'][2:]).max() > 1e-5
    for row in record['records']:
        assert abs(row['source_current_readout']-row['bond_current']) < 2e-8
        for pointer in row['pointer_records']:
            assert len(pointer['probabilities']) == 2
            assert abs(sum(pointer['probabilities'])-1) < 1e-12
            assert min(pointer['probabilities']) > 0
    assert arrays['reference_covariance_with_purifiers'].shape == (96, 96)


@pytest.mark.parametrize('bad', [True, 2, -2, '0'])
def test_invalid_native_clock_current_rejected(bad):
    with pytest.raises(ValueError):
        p.execute_block(p.components()[0], bad)
