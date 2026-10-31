"""Independent CAR, full-current and invariant-replay adversarial checks."""
import copy
import json

import numpy as np
import pytest

from ..common_source.chiral_candidate import action
from . import benchmark as report
from .quasiparticles import (MODE_COUNT, current_probe_data, particle_hole,
                            quasiparticles, transform_one_body_probe, validate_chart)


def test_full_action_and_charge_reconstructed_without_projecting_negative_or_zero_modes():
    chart = quasiparticles()
    basis = validate_chart(chart['positive_basis'], chart['energies'], chart['charge6'])
    assert basis.shape == (180, 180)
    assert np.count_nonzero(chart['energies'] == 0) == 6
    energies = np.r_[chart['energies'], -chart['energies']]
    assert np.linalg.norm((basis*energies)@basis.conj().T-chart['bdg']) < 2e-10
    charges = np.r_[chart['charge6'], -chart['charge6']]/6
    assert np.linalg.norm((basis*charges)@basis.conj().T-chart['original_nambu_charge']) < 2e-11
    assert np.array_equal(particle_hole(particle_hole(basis)), basis)
    assert np.linalg.norm(chart['positive_basis'].conj().T@particle_hole(chart['positive_basis'])) < 2e-11


def test_degenerate_chart_rotation_preserves_every_invariant_block():
    rng = np.random.default_rng(6041140)
    chart = quasiparticles()
    rotated = chart['positive_basis'].copy()
    for block in chart['projectors']:
        indices = block['indices']
        z = rng.normal(size=(len(indices), len(indices)))+1j*rng.normal(size=(len(indices), len(indices)))
        unitary, _ = np.linalg.qr(z)
        rotated[:, indices] = rotated[:, indices]@unitary
        assert np.linalg.norm(rotated[:, indices]@rotated[:, indices].conj().T-block['projector']) < 2e-11
    validate_chart(rotated, chart['energies'], chart['charge6'])
    assert np.linalg.norm(rotated-chart['positive_basis']) > 1


@pytest.mark.parametrize('damage', ['zero_mode_deleted', 'wrong_hole', 'changed_energy', 'changed_charge', 'NaN'])
def test_false_canonical_charts_rejected(damage):
    chart = quasiparticles()
    w, e, q = (chart[key].copy() for key in ('positive_basis', 'energies', 'charge6'))
    if damage == 'zero_mode_deleted':
        w[:, np.flatnonzero(e == 0)[0]] = 0
    elif damage == 'wrong_hole':
        w[:, 0] = particle_hole(w[:, 0])
    elif damage == 'changed_energy':
        e[0] += 1e-7
    elif damage == 'changed_charge':
        q[0] *= -1
    else:
        w[0, 0] = np.nan
    with pytest.raises(ValueError):
        validate_chart(w, e, q)


def test_independent_two_mode_jordan_wigner_pairing_cost_and_number_boundary():
    chart = quasiparticles()
    q = np.diag(chart['original_nambu_charge'])[:90].real
    # A charged-lepton original annihilator has one pairing partner. This
    # independently constructs its four-dimensional many-body Hamiltonian.
    index = next(i for i in np.flatnonzero(q == 1)
                 if np.count_nonzero(chart['pairing'][i]) == 1)
    partner = int(np.flatnonzero(chart['pairing'][index])[0])
    amplitude = chart['pairing'][index, partner]
    annihilators = []
    for mode in (0, 1):
        c = np.zeros((4, 4), complex)
        for mask in range(4):
            if mask >> mode & 1:
                c[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
        annihilators.append(c)
    a, b = annihilators
    h = amplitude*a.conj().T@b.conj().T + amplitude.conjugate()*b@a
    number = a.conj().T@a+b.conj().T@b
    parity = np.diag([1, -1, -1, 1])
    m = abs(amplitude)
    assert np.allclose(np.linalg.eigvalsh(h), [-m, 0, 0, m], atol=1e-13, rtol=0)
    assert np.linalg.norm(h@number-number@h) > m
    assert np.array_equal(h@parity, parity@h)
    assert np.linalg.norm(h@(q[index]*a.conj().T@a+q[partner]*b.conj().T@b)
                          -(q[index]*a.conj().T@a+q[partner]*b.conj().T@b)@h) < 1e-13


def test_all_twelve_probes_retained_and_exact_higgs_stabilizer_exhausts_preserving_span():
    data = current_probe_data()
    assert data['exact_higgs_orbit_rank'] == 3
    assert data['exact_stabilizer_dimension'] == 9
    assert len(data['probes']) == 12
    for probe in data['probes']:
        assert probe['full'].shape == (180, 180)
        assert np.linalg.norm(probe['full']-probe['full'].conj().T) < 2e-11
        assert probe['pair_block_norm'] > .9
    blocks = np.column_stack([p['pair_creation'].ravel() for p in data['probes']])
    assert np.linalg.norm(blocks@data['stabilizer_coefficients']) < 2e-11
    assert data['anomalous_singular_values'][2] > 2
    assert data['anomalous_singular_values'][3] < 2e-11


def test_electric_charge_and_spin_are_number_preserving_rest_probes():
    charge = transform_one_body_probe(action.spin_electric_charge())
    assert charge['pair_block_norm'] < 2e-11
    assert abs(charge['quasiparticle_vacuum_offset']) < 2e-11
    assert np.linalg.norm(charge['number_preserving']-quasiparticles()['electric_charge']) < 2e-11
    pauli = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))
    spin = [transform_one_body_probe(np.kron(s/2, np.eye(45))) for s in pauli]
    assert all(p['pair_block_norm'] < 2e-11 for p in spin)
    generators = [p['number_preserving'] for p in spin]
    assert np.linalg.norm(sum(s@s for s in generators)-.75*np.eye(MODE_COUNT)) < 2e-11
    assert np.linalg.norm(generators[0]@generators[1]-generators[1]@generators[0]-1j*generators[2]) < 2e-11


def test_transformed_bare_number_retains_pair_blocks_and_vacuum_offset():
    number = transform_one_body_probe(np.eye(MODE_COUNT))
    # 84 paired original modes have half occupation in this vacuum; the
    # explicitly retained six massless annihilators stay empty.
    assert abs(number['quasiparticle_vacuum_offset']-42) < 2e-11
    assert number['pair_block_norm'] > 9
    assert abs(np.trace(number['number_preserving']).real-6) < 2e-11


@pytest.mark.parametrize('bad', [np.eye(89), np.full((90, 90), np.nan), np.triu(np.ones((90, 90)))])
def test_invalid_probes_rejected(bad):
    with pytest.raises(ValueError):
        transform_one_body_probe(bad)


def test_actual_nine_charged_input_scales_compared_after_diagonalization():
    rows = report.charged_mass_comparison()
    assert len(rows) == 18
    assert sum(row['multiplicity'] for row in rows) == 84
    for row in rows:
        assert abs(row['diagonalized_energy_GeV']-row['independent_y_times_v_over_sqrt2_GeV']) < 2e-10
    lepton = [r['diagonalized_energy_GeV'] for r in rows if r['channel'] == 'charged_lepton' and r['electric_charge_sixths'] == -6]
    assert lepton[2]/lepton[0] > 3000


def test_same_mass_weyl_clifford_and_inertial_constraint():
    data = report.mass_shell_certificate()
    assert data['same_mass_kernel_and_all_180_directions']
    assert not data['native_port_transport_is_this_momentum']
    assert not data['rest_quasiparticle_basis_or_zero_half_continues_to_finite_p']


def test_receipt_replay_accepts_roundoff_but_rejects_spectrum_domain_and_hash_changes():
    actual = report.certificate()
    report.compare_receipt(actual, json.loads(report.RECEIPT.read_text()))
    noisy = copy.deepcopy(actual)
    noisy['charged_energy_blocks'][0]['diagonalized_energy_GeV'] += 1e-13
    report.compare_receipt(actual, noisy)
    for key in ('energy', 'domain', 'hash'):
        wrong = copy.deepcopy(actual)
        if key == 'energy':
            wrong['charged_energy_blocks'][0]['diagonalized_energy_GeV'] += 1e-7
        elif key == 'domain':
            wrong['domain']['all_Nambu_directions'] -= 12
        else:
            wrong['source_sha256'][next(iter(wrong['source_sha256']))] = '0'*64
        with pytest.raises(ValueError):
            report.compare_receipt(actual, wrong)
