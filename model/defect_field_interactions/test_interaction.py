"""Independent full-CAR and source-column tests of the driven extension."""
from decimal import Decimal
import math

import numpy as np
import pytest
from scipy.linalg import expm

from ..defect_field_benchmark.quasiparticles import quasiparticles
from . import interaction as model


def _annihilator(mode):
    value = np.zeros((4, 4), complex)
    for mask in range(4):
        if mask >> mode & 1:
            value[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
    return value


def test_actual_higgs_channels_have_opposite_creation_phase_conventions():
    rows = model.channel_pairings()
    assert {name: r['creation_phase_sign'] for name, r in rows.items()} == {'up': -1, 'down': 1, 'charged_lepton': 1}
    assert np.linalg.norm(model.higgs_pairing(0)-quasiparticles()['pairing']) < 2e-10
    # A uniform phase on all creation coefficients is an incorrect drive.
    phase = .7
    assert np.linalg.norm(model.higgs_pairing(phase)-np.exp(1j*phase)*model.higgs_pairing(0)) > 100


def test_complete_pair_factorization_preserves_charge_and_static_projectors():
    paired = model.paired_chart()
    base = quasiparticles()
    assert len(paired['pairs']) == 42
    occupied = [p[key] for p in paired['pairs'] for key in ('left', 'right')]
    assert len(set(occupied)) == 84
    assert set(occupied).isdisjoint(paired['zero_modes'])
    for block in base['projectors']:
        columns = paired['positive_basis'][:, block['indices']]
        assert np.linalg.norm(columns@columns.conj().T-block['projector']) < 2e-11
    q = base['original_nambu_charge']
    for phase in (0., .2, math.pi/2):
        h = model.nambu(model.higgs_pairing(phase))
        assert np.linalg.norm(h@q-q@h) < 2e-10


@pytest.mark.parametrize('energy,phase,duration', [(.102759757, 0., 1.), (.102759757, math.pi/2, 1.25), (2.862404953, -.7, .75)])
def test_source_gate_equals_independent_jordan_wigner_evolution(energy, phase, duration):
    a, b = _annihilator(0), _annihilator(1)
    number = a.conj().T@a+b.conj().T@b
    h = energy*math.cos(phase)*(number-np.eye(4)) + energy*math.sin(phase)*(a.conj().T@b.conj().T+b@a)
    assert np.linalg.norm(h-model.pair_hamiltonian(energy, phase)) < 1e-13
    actual, _ = model.pair_gate(energy, phase, duration)
    assert np.linalg.norm(actual-expm(-1j*duration*h), 2) < model.SOURCE_GUARD
    charge = a.conj().T@a-b.conj().T@b
    parity = np.diag([1, -1, -1, 1])
    assert np.linalg.norm(actual@charge-charge@actual) < model.SOURCE_GUARD
    assert np.linalg.norm(actual@parity-parity@actual) < model.SOURCE_GUARD
    if phase:
        assert np.linalg.norm(h@number-number@h) > .01


def test_routed_source_gate_preserves_all_intervening_fermionic_signs():
    annihilators = []
    for mode in range(4):
        c = np.zeros((16, 16), complex)
        for mask in range(16):
            if mask >> mode & 1:
                c[mask ^ (1 << mode), mask] = (-1)**(mask & ((1 << mode)-1)).bit_count()
        annihilators.append(c)
    a, b = annihilators[0], annihilators[3]
    energy, phase, duration = .102759757, .7, 1.25
    h = energy*math.cos(phase)*(a.conj().T@a+b.conj().T@b-np.eye(16))
    h += energy*math.sin(phase)*(a.conj().T@b.conj().T+b@a)
    word = model.routed_pair_word(energy, phase, duration, 0, 3)
    actual = model.compiler.execute(word, 4)
    assert np.linalg.norm(actual-expm(-1j*duration*h), 2) < model.SOURCE_GUARD
    # An ordinary qubit pair rotation is wrong on odd spectator occupation.
    naive = model.compiler.execute(model.compiler.remap(model.pair_word(energy, phase, duration), [0, 3]), 4)
    assert np.linalg.norm(naive-actual, 2) > .1


@pytest.mark.parametrize('indices', [(0, 66), (27, 78)])
def test_source_bare_chart_joins_both_static_mass_and_driven_phase(indices):
    chart = model.bare_component_chart(indices)
    rotation = chart['source_unitary']
    z = quasiparticles()['pairing'][indices]
    for phase in (0., -.8, 1.2):
        coefficient = z*np.exp(1j*chart['creation_phase_sign']*phase)
        wanted = np.zeros((4, 4), complex)
        wanted[3, 0], wanted[0, 3] = coefficient, coefficient.conjugate()
        actual = rotation@model.pair_hamiltonian(chart['energy'], phase)@rotation.conj().T
        assert np.linalg.norm(actual-wanted) < model.SOURCE_GUARD


@pytest.mark.parametrize('initial', [0, 1, 2, 3])
def test_complete_occupation_pointer_retains_every_initial_sector(initial):
    energy, phase, duration = .102759757, math.pi/2, 1.
    observed = model.observed_pair(energy, phase, duration, initial)
    p = math.sin(energy*duration)**2
    wanted = ([1-p, 0, 0, p], [0, 1, 0, 0], [0, 0, 1, 0], [p, 0, 0, 1-p])[initial]
    assert np.max(abs(observed['probabilities']-wanted)) < model.SOURCE_GUARD
    assert observed['amplitudes'].shape == (4, 4)


def test_driven_full_number_law_retains_underflow_tails_and_leaves_N2():
    law = model.extended_number_law(1., math.pi/2)
    values = list(map(Decimal, law['probabilities_decimal']))
    assert len(values) == 503
    assert all(value > 0 for value in values)
    assert min(values) < Decimal('1e-308')
    assert law['particle_numbers'][0] == 2
    assert law['particle_numbers'][-1] == 1006
    assert law['expected_particle_number'] > 100
    assert abs(law['expected_particle_number']-(2+2*sum(law['factor_probabilities']))) < 1e-10


def test_static_limit_does_not_create_particles():
    law = model.extended_number_law(1., 0.)
    assert law['expected_particle_number'] == 2
    assert Decimal(law['probability_staying_N2_decimal']) == 1


@pytest.mark.parametrize('current', [-1, 0, 1])
def test_one_clock_multiplier_applies_to_entire_driven_action(current):
    law = model.extended_number_law(.4, math.pi/2, current)
    assert law['effective_duration'] == .4*(1+current/4)
    first_energy = model.paired_chart()['pairs'][1]['energy']
    assert abs(law['factor_probabilities'][0]-math.sin(first_energy*law['effective_duration'])**2) < 1e-13


def test_actual_complete_repair_records_destroy_drive_inverse_echo():
    energy = .102759757
    rows = model.repair_echo(energy)
    expected = 1-2*math.sin(energy)**2*math.cos(energy)**2
    for row in rows.values():
        assert row['all_record_field_amplitudes'].shape == (32, 16)
        assert abs(row['observed_return_probability']-expected) < 8e-7
        assert row['observed_return_probability'] < .99


def test_full_projector_echo_composes_without_postmeasurement_Gaussian_assumption():
    from ..defect_field_bridge.full_domain import channel

    energy, duration = .102759757, 3.
    g, _ = model.pair_gate(energy, math.pi/2, duration)
    unitary = np.kron(g, g)
    psi = unitary[:, 0]
    rho = np.outer(psi, psi.conj())
    both_empty = math.cos(energy*duration)**4
    expected = both_empty**2+(1-both_empty)**2
    for law in ('C', 'E'):
        after = channel(rho, sites=2, flavors=2, law=law, seam=(0, 1))
        echo = unitary.conj().T@after@unitary
        assert abs(echo[0, 0].real-expected) < 8e-7
    full = model.full_repair_echo()
    assert .1 < full['two_endpoint_vacuum_probability'] < .9
    assert full['inverse_drive_return_probability'] < .9
    assert not full['postmeasurement_Gaussian_assumption']


@pytest.mark.parametrize('args', [(-1., 0., 1.), (1., float('nan'), 1.), (1., 0., float('inf'))])
def test_invalid_drive_rejected(args):
    with pytest.raises(ValueError):
        model.pair_gate(*args)
