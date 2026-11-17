"""Independent physical-fiber and native-reference checks of the clock join."""
from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from tools import strict_json
from . import acceptance_audit as audit
from ..bound_pair.fiber import PairFiber, clock_subspace
from ..bound_pair.protection import PeriodicPair


def test_exact_native_eigenvalues_and_independent_tensor_multiplicities():
    levels = audit.native_levels()
    assert levels == {0: 45, 1: 20, 2: 1}
    expanded = [level for level, multiplicity in levels.items() for _ in range(multiplicity)]
    for copies in (1, 2):
        measured = {}
        for row in product(expanded, repeat=copies):
            measured[sum(row)] = measured.get(sum(row), 0)+1
        assert audit.tensor_multiplicities(copies) == measured
    for copies in (3, 10, 100):
        result = audit.tensor_multiplicities(copies)
        assert sum(result.values()) == 66**copies
        assert set(result) == set(range(2*copies+1))


@pytest.mark.parametrize('side', [4, 6, 8])
def test_integer_certificate_against_actual_physical_pair_fiber(side):
    moments = audit.exact_pair_moments()
    f = PairFiber([0, 0, 0], side=side)
    state = f.d.copy()
    for expected in moments:
        assert abs(np.vdot(f.d, state)-float(expected)) < 2e-14
        assert f.exchange_residual(state) < 2e-14
        state = f.free_apply(state)
    pole = f.pole()['energy']
    assert F(13, 12) < pole < F(3, 2)
    assert np.linalg.norm(f.apply(f.eigenstate())-pole*f.eigenstate()) < 2e-13
    assert max(clock_subspace(side)['residuals'].values()) < 2e-13


def test_rayleigh_witness_in_independent_complete_position_car_action():
    space = PeriodicPair()
    f = PairFiber([0, 0, 0], side=4)
    contact = space.fiber_state([0, 0, 0], f.d)
    td = space.free_h@contact
    assert abs(np.vdot(td, td)-3/32) < 2e-14
    assert abs(np.vdot(contact, td)) < 2e-14
    # Actual compression of the complete8128-dimensional Hamiltonian.
    basis = np.column_stack((contact, td/np.linalg.norm(td)))
    compressed = basis.conj().T@(space.hamiltonian@basis)
    expected = np.array([[1, np.sqrt(3/32)], [np.sqrt(3/32), 0]])
    assert np.max(abs(compressed-expected)) < 2e-14
    assert np.linalg.eigvalsh(compressed)[-1] > 13/12
    certificate = audit.spectral_certificate()
    assert certificate['rayleigh_polynomial_at_lower'] == '-1/288'
    assert certificate['centered_generator_residual_over_U_strict_lower'] == '1/24'


def test_complete_reference_preserving_recurrence_countercontrol():
    energy = PairFiber([0, 0, 0], side=4).pole()['energy']
    # Two physical orthogonal field eigenstates: eta dark reference and K=0 pair.
    initial = np.array([1, 1j])/np.sqrt(2)
    evolved = np.exp(-2j*np.pi*np.array([0, energy]))*initial
    orthogonal_outcome = 1-abs(np.vdot(initial, evolved))**2
    distance = np.linalg.norm(np.outer(evolved, evolved.conj())-
                              np.outer(initial, initial.conj()), ord='nuc')/2
    assert orthogonal_outcome > 1/36
    assert distance > 1/6
    # The full native tensor generator has a scalar recurrence, for every copy count.
    for copies in (1, 3, 32):
        assert np.max(abs(np.exp(2j*np.pi*np.arange(2*copies+1))-1)) < 6e-14
    # Appending an untouched quantum reference preserves this witness distance.
    reference = np.array([[.6, .2j], [-.2j, .4]])
    delta = np.outer(evolved, evolved.conj())-np.outer(initial, initial.conj())
    assert abs(np.linalg.norm(np.kron(delta, reference), ord='nuc')/2-distance) < 2e-14


def test_compression_of_energy_expectations_is_not_clock_intertwining():
    energy = PairFiber([0, 0, 0], side=4).pole()['energy']
    native = np.diag([0, 1, 2])
    v = np.array([[1, 0], [0, np.sqrt(2-energy)], [0, np.sqrt(energy-1)]])
    field = np.diag([0, energy])
    assert np.linalg.norm(v.T@native@v-field) < 2e-14
    assert np.linalg.norm(native@v-v@field, 2) > 1/24


def test_supplied_nonproduct_gibbs_reference_is_an_explicit_escape():
    energy = PairFiber([0, 0, 0], side=4).pole()['energy']
    h = np.array([0, 1, energy])
    beta = math.log(25/8)
    reference = np.exp(-beta*h)
    reference /= reference.sum()
    modular = -np.log(reference)
    assert np.max(abs(modular-modular[0]-beta*h)) < 2e-15
    assert abs((modular[2]-modular[0])/beta-round(energy)) > 1/12
    # This thermal spectrum violates the excluded product-reference premise.


def test_exact_hilbert_space_capacity_is_necessary_but_not_an_encoder():
    for side in (4, 6, 8, 12, 16):
        result = audit.field_capacity(side)
        dimension = math.comb(2*side**3, 2)
        count = result['native_H66_factors_at_least']
        assert 66**(count-1) < dimension <= 66**count
        assert not result['isometry_constructed']
        assert not result['dimension_capacity_solves_modular_obstruction']
    assert audit.field_capacity(4)['complete_N2_dimension'] == 8128
    assert audit.field_capacity(4)['native_H66_factors_at_least'] == 3


@pytest.mark.parametrize('bad', [True, 0, -1, 1.5, 257])
def test_original_bounded_copy_count_required(bad):
    with pytest.raises(ValueError):
        audit.tensor_multiplicities(bad)


def test_source_custody_failure_is_not_silently_refreshed():
    audit.check_bindings()
    path = next(iter(audit.SOURCE_BINDINGS))
    with pytest.raises(ValueError, match='changed source'):
        audit.check_bindings(bindings={path: '0'*64})
    with pytest.raises(ValueError, match='missing or external'):
        audit.check_bindings(bindings={'../outside.py': '0'*64})


def test_acceptance_report_does_not_falsify_the_whole_program():
    result = audit.report()
    assert result['accepted_gate_changes'] == []
    assert result['canonical_premises_discharged'] == []
    assert not result['whole_program_falsified']
    assert result['physical_comparison_count_added'] == 0
    assert result['official_model'] == 'CCF-LT-2'
    assert not result['applies_to_official_model']
    assert 'independent native product-reference' in result['source_domain']
