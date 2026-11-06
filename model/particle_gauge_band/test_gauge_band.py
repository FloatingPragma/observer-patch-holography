"""Complete native replay and mutation-sensitive algebra controls."""
from fractions import Fraction

import numpy as np
import sympy as s

from . import gauge_band as band
from . import quantum_link_completion as link
from . import native_link_projector_witness as native_car
from .verify import compare


def test_full_literal_gauge_band_and_polar_transport(tmp_path):
    receipt = band.execute(tmp_path)
    assert receipt['exact_bounds']['raw_record_annihilator_count'] == 180
    assert receipt['source_graph']['sites'] == 343
    assert receipt['source_graph']['edges'] == 882
    for block in receipt['curved_link_matter_blocks']:
        assert 0 < block['exterior_to_core_probability_ratio_max'] <= float(Fraction(block['exterior_to_core_probability_ratio_upper_exact']))
    assert receipt['lepton_BdG_execution']['raw_family_unordered_projector_pairs_executed'] == 91
    assert receipt['lepton_BdG_execution']['original_scalar_band_pairing_leakage_Frobenius'] > 1e-3
    assert receipt['gauge_covariance_execution']['full_BdG_operator_error'] < 1e-8


def test_exact_rotor_completion_all_blocks(tmp_path):
    receipt = link.execute(tmp_path)
    ranks = {row['block']: row['completed_combined_rank'] for row in receipt['completed_coordinate_blocks']}
    assert ranks == {'u^c': 18, 'L': 2, 'Q': 30, 'd^c': 3, 'e^c': 10}
    assert receipt['total_exact_link_covariance_checks'] == 1416
    assert receipt['link_Hamiltonian_commutant']['independent_energy_differences_after_common_offset'] == 5


def test_literal_native_link_noninteger_gauss_singlet(tmp_path):
    receipt = native_car.execute(tmp_path)
    assert receipt['noninteger_occupation_Gauss_singlet']['all_native_current_endpoint_checks'] == 72
    assert receipt['scope']['noninteger_naive_Fock_occupation_proved_exactly']


def test_changing_coordinate_charge_fails_covariance():
    masks, blocks, currents = link.exterior_data()
    tensors, _ = link.clebsch_gordan(masks, blocks, currents)
    indices, _, left, _ = link.link_representation(blocks, currents)
    _, operators = link.coordinate_operators(6, blocks, indices, tensors)
    U = operators[0, 0]
    mutant = s.SparseMatrix(U)
    mutant[0, indices[6, 0, 0]] = 1
    failures = 0
    for T, L in zip(currents, left, strict=True):
        value = T[blocks[6][0], blocks[6][0]]
        assert link.zero(L*U-U*L+value*U)
        failures += not link.zero(L*mutant-mutant*L+value*mutant)
    assert failures > 0


def test_polar_map_moves_band_and_preserves_particle_hole():
    plus = np.eye(8, dtype=complex)[:, :2]
    changed = plus.copy()
    changed[:, 0] = .8*plus[:, 0]+.6*np.eye(8)[:, 2]
    old = {'plus': plus, 'minus': band.ph(plus)}
    new = {'plus': changed, 'minus': band.ph(changed)}
    result = band.small_polar_transport(old, new)
    Z, U = result['Z'], result['U_small']
    full = np.eye(8)+Z@(U-np.eye(U.shape[0]))@Z.conj().T
    np.testing.assert_allclose(full@plus@plus.conj().T@full.conj().T, changed@changed.conj().T, atol=2e-14)
    np.testing.assert_allclose(full.conj().T@full, np.eye(8), atol=2e-14)
    assert result['particle_hole_transport_error'] < 2e-14


def test_exact_band_bounds_and_replay_rejects_mutation():
    bounds = band.exact_constants()
    assert Fraction(bounds['three_projector_alignment_minus_identity_norm_upper']) < 1
    assert Fraction(bounds['compressed_grading_square_lower']) > 0
    try:
        compare({'rank': 360}, {'rank': 358})
    except AssertionError:
        pass
    else:
        raise AssertionError('replay accepted a changed exact rank')


def test_replay_rejects_nonfinite_numeric_diagnostics():
    for left, right in [(0.0, float('nan')), (float('inf'), float('inf'))]:
        try:
            compare(left, right)
        except AssertionError:
            pass
        else:
            raise AssertionError('replay admitted a nonfinite diagnostic')
