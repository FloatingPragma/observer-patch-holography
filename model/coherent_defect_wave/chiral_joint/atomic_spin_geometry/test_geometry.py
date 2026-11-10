"""Independent spin assembly, full spectra and hostile source transfers."""
import numpy as np
import pytest
import sympy as s
from . import cube, momentum, transfer


def test_full_twisted_source_and_every_magnetization_sector():
    plain, twisted, unitary = cube.full_operators()
    assert unitary.H*twisted*unitary == plain
    for down in range(9):
        basis, block = cube.magnetization_block(down)
        assert plain.extract(basis, basis) == block
    for bad in (-1, 9, 2.0, True):
        with pytest.raises(ValueError):
            cube.magnetization_block(bad)


def test_complete_exact_spectrum_against_independent_dense_diagonalization():
    # Independent tensor-product assembly in a reversed spin-bit convention.
    pauli = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
             np.array([[1, 0], [0, -1]]))
    operators = []
    for site in range(8):
        row = []
        for sigma in pauli:
            matrix = np.array([[1]], dtype=complex)
            for target in range(8):
                matrix = np.kron(matrix, sigma/2 if target == site else np.eye(2))
            row.append(matrix)
        operators.append(row)
    h = np.zeros((256, 256), dtype=complex)
    for left, right, _ in cube.cube_edges():
        for axis in range(3):
            h += operators[left][axis]@operators[right][axis]
    actual = np.linalg.eigvalsh(h)
    receipt = cube.spectrum(); offset = 0
    for row in receipt['complete_spectrum']:
        low, high = map(lambda x: float(s.Rational(x)), row['energy_interval'])
        count = row['full_multiplicity']
        assert np.all(actual[offset:offset+count] >= low-2e-12)
        assert np.all(actual[offset:offset+count] <= high+2e-12)
        offset += count
    assert offset == 256 and receipt['ground_rank'] == 1 and receipt['first_excited_rank'] == 3


def test_actual_momentum_does_not_become_rotated_spin_momentum():
    result = momentum.certificate()
    rows = result['complete_first_triplet_translation_records']
    assert result['ground_translation_characters'] == [1, 1, 1]
    assert len(rows) == 8
    assert sorted(row['momenta_in_pi_units'] for row in rows if row['first_triplet_rank']) == [[0, 0, 1], [0, 1, 0], [1, 0, 0]]


def test_full_interacting_transfer_requires_positive_exchange_and_remainder(monkeypatch):
    result = transfer.certificate()
    assert s.Rational(result['reported_nonzero_hopping_coefficient']) == s.Rational(1, 2**29)
    assert s.Rational(result['full_physical_ground_gap_interval_at_reported_coefficient'][0]) > 0
    assert result['unique_full_physical_interacting_ground_proved']
    assert not result['lowest_gap_is_a_zero_momentum_particle_rest_mass']
    filter_law = result['formal_full_carrier_ground_filter']
    assert s.Rational(filter_law['success_probability']) + s.Rational(filter_law['failure_probability']) == 1
    original = transfer.exchange_certificate
    def invalid_exchange():
        row = dict(original()); row['exact_J_interval'] = ['0', '1']; return row
    monkeypatch.setattr(transfer, 'exchange_certificate', invalid_exchange)
    with pytest.raises(ValueError, match='positive exact rational second-order gap'):
        transfer.certificate()


def test_missing_source_or_dependency_cannot_be_pinned(monkeypatch):
    from . import certificate as custody
    monkeypatch.setattr(custody, 'DEPENDENCIES', (*custody.DEPENDENCIES, 'missing-spin-reference.json'))
    with pytest.raises(FileNotFoundError, match='missing declared spin geometry input'):
        custody.input_hashes()
