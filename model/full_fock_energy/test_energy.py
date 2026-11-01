"""Independent CAR assembly and hostile complete-energy controls."""
from fractions import Fraction
import math
import numpy as np
import pytest
from scipy.linalg import expm

from ..full_fock_field.field import FieldGeometry, FullFock
from ..field_stress.energy import GaussianInteger, commutator
from ..observer.source_scalar import compiler
from ..field_thermal.thermal import ThermalField
from .energy import certify, charged_probe_certificate, grand_canonical, local_energies
from .instrument import execute, expectation, pointer_work


def ladder(modes, at):
    result = np.ones((1, 1), complex)
    for k in reversed(range(modes)):
        factor = (np.array([[0, 1], [0, 0]]) if k == at else
                  np.diag([1, -1]) if k < at else np.eye(2))
        result = np.kron(result, factor)
    return result


def independent_local(geometry):
    a = [ladder(geometry.modes, k) for k in range(geometry.modes)]
    sigma = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
             np.diag([1, -1]))
    locations = list(geometry.sites)
    result = [geometry.interaction*(a[2*x].conj().T@a[2*x])@
              (a[2*x+1].conj().T@a[2*x+1]) for x in range(geometry.volume)]
    for x, coordinate in enumerate(locations):
        for axis in range(3):
            neighbor = list(coordinate); neighbor[axis] += 1
            if tuple(neighbor) not in locations:
                continue
            y = locations.index(tuple(neighbor))
            bond = np.zeros_like(result[0])
            for i in range(2):
                for j in range(2):
                    hop = 1j*geometry.speed/2*sigma[axis][i, j]*(a[2*x+i].conj().T@a[2*y+j])
                    bond += hop+hop.conj().T
            result[x] += bond/2; result[y] += bond/2
    return result


@pytest.mark.parametrize('sides', [(2, 1, 1), (1, 2, 1), (1, 1, 2), (2, 2, 1), (3, 1, 1)])
def test_exact_complete_operators_against_independent_Jordan_Wigner(sides):
    geometry = FieldGeometry(sides)
    report, op = certify(geometry)
    expected = independent_local(geometry)
    assert sum(report['number_sector_dimensions']) == report['full_Fock_dimension']
    assert report['complete_operator_identities_exact']
    assert not report['number_sector_projection_applied']
    for actual, target in zip(op.local, expected, strict=True):
        np.testing.assert_array_equal(actual.matrix(op.scale).toarray(), target)
    for (x, y), current in op.currents.items():
        target = 1j*(expected[x]@expected[y]-expected[y]@expected[x])
        np.testing.assert_array_equal(current.matrix(op.scale**2).toarray(), target)


def test_range_two_current_is_necessary_for_a_regional_boundary():
    _, op = certify(FieldGeometry((3, 1, 1)))
    assert (0, 2) in op.currents
    regional, outward = op.region((0,))
    assert (commutator(op.total, regional).times_i()+outward).is_zero()
    truncated = outward-op.currents[0, 2]
    assert not (commutator(op.total, regional).times_i()+truncated).is_zero()


def test_interaction_is_retained_and_current_sign_matches_time_derivative():
    geometry = FieldGeometry((2, 1, 1))
    _, op = certify(geometry)
    _, free = certify(FieldGeometry((2, 1, 1), interaction=0))
    assert not (op.currents[0, 1]-free.currents.get((0, 1), GaussianInteger.zero(16))).is_zero()
    random = np.random.default_rng(94)
    state = random.normal(size=(16, 3))+1j*random.normal(size=(16, 3))
    state /= np.linalg.norm(state)
    h = sum(independent_local(geometry))
    dt = 1e-4
    before, after = expm(1j*dt*h)@state, expm(-1j*dt*h)@state
    local = op.local[0].matrix(op.scale).toarray()
    derivative = (expectation(after, local)-expectation(before, local))/(2*dt)
    outward = expectation(state, op.currents[0, 1].matrix(op.scale**2))
    assert abs(outward) > 1e-4
    assert abs(derivative+outward) < 1e-9


def test_every_occupation_sector_and_between_sector_coherences_are_present():
    geometry = FieldGeometry((2, 1, 1))
    _, op = certify(geometry)
    # Do not only test N=2: the complete matrix is compared to every block.
    number = np.array([mask.bit_count() for mask in range(16)])
    h = op.total.matrix(op.scale).toarray()
    for sector in range(5):
        indices = np.flatnonzero(number == sector)
        assert len(indices) == math.comb(4, sector)
        np.testing.assert_array_equal(h[np.ix_(indices, number != sector)], 0)
    assert h[15, 15] == 2  # filled state omitted by a two-particle projection
    assert h[0, 0] == 0
    psi = np.arange(1, 17)+1j*np.arange(16)
    psi /= np.linalg.norm(psi)
    evolved = expm(-.3j*h)@psi
    assert all(abs(evolved[k]) > 0 for k in range(16))
    assert abs(np.vdot(evolved, h@evolved)-np.vdot(psi, h@psi)) < 1e-14


@pytest.mark.parametrize('mode', range(4))
def test_complete_number_pointer_work_identity_and_Gibbs_backaction(mode):
    geometry = FieldGeometry((2, 1, 1))
    work, scale = pointer_work(geometry, mode)
    case = ThermalField(geometry)
    target = case.target()
    rho, h = target['reference'], target['hamiltonian']
    a = ladder(4, mode); n = a.conj().T@a; empty = np.eye(16)-n
    after = n@rho@n+empty@rho@empty
    direct = np.trace(h@(after-rho)).real
    assert direct > .0001
    np.testing.assert_allclose(direct, np.trace(rho@work.matrix(scale).toarray()).real, atol=1e-15)
    # The contact term commutes with the read, but hopping does not.
    atomic, _ = pointer_work(FieldGeometry((2, 1, 1), speed=0), mode)
    assert atomic.is_zero()


def test_source_continuation_and_both_pointer_energies_retain_reference():
    geometry = FieldGeometry((2, 1, 1))
    random = np.random.default_rng(7)
    initial = random.normal(size=(16, 3))+1j*random.normal(size=(16, 3))
    initial /= np.linalg.norm(initial)
    report, evolved, read = execute(geometry, initial)
    assert report['source_complete_evolution_error_diagnostic'] < float(Fraction(report['ideal_split_half_diamond_upper']))
    assert report['source_success_energy_drift'] < float(Fraction(report['ideal_unitary_energy_drift_upper']))
    assert all(p > .1 for p in report['pointer_branch_probabilities'])
    assert abs(report['executed_pointer_success_energy_change']) > .0001
    assert not report['success_normalization_applied']
    assert not report['abort_output_energy_assigned']
    assert not report['physical_total_energy_balance_closed']
    np.testing.assert_allclose(read.success.conj().T@read.success+read.failure_gram+
                               evolved.failure_gram, initial.conj().T@initial, atol=2e-11)
    assert report['evolution_resources']['native_pulses'] > 100


def test_leaking_pointer_exports_remain_unassigned_and_do_not_normalize(monkeypatch):
    import model.full_fock_energy.instrument as instrument
    original = instrument.execute_complete
    # Only modify the pointer execution. The intended instrument comparison
    # must reject its missing weight even though normalized successes agree.
    def leaking_pointer(word, bits, initial):
        if bits != 5:
            return original(word, bits, initial)
        native = compiler.native_pulse
        def leakage(kind, angle):
            u = native(kind, angle)
            generator = np.zeros((6, 6), complex)
            generator[0, 4] = generator[4, 0] = 1
            return expm(.02j*generator)@u
        with monkeypatch.context() as patch:
            patch.setattr(compiler, 'native_pulse', leakage)
            return original(word, bits, initial)
    monkeypatch.setattr(instrument, 'execute_complete', leaking_pointer)
    initial = np.ones((16, 1))/4
    with pytest.raises(ArithmeticError, match='complete source energy-accounting'):
        execute(FieldGeometry((2, 1, 1)), initial)


def test_changed_global_action_is_rejected(monkeypatch):
    original = FullFock.hamiltonian
    def changed(self, geometry):
        return original(self, geometry)+.03125*self.number(0)
    monkeypatch.setattr(FullFock, 'hamiltonian', changed)
    with pytest.raises(ArithmeticError, match='actual full-Fock action'):
        certify()


def test_domains_fail_closed_and_zero_hamiltonian_is_supported():
    report, op = certify(FieldGeometry((1, 1, 1), speed=0, interaction=0))
    assert report['nonzero_energy_current_pairs'] == 0 and op.total.is_zero()
    with pytest.raises(ValueError):
        certify(FieldGeometry((3, 2, 1)))
    with pytest.raises(ValueError, match='denominator'):
        certify(FieldGeometry((2, 1, 1), speed=.1))
    for sites in ((True,), (0, 0), (2,), (0.5,)):
        with pytest.raises(ValueError):
            local_energies().region(sites)
    for initial in (np.ones(16), np.ma.array(np.ones(16)/4, mask=True), np.full(16, np.nan)):
        with pytest.raises(ValueError):
            execute(FieldGeometry((2, 1, 1)), initial)
    with pytest.raises(ValueError):
        execute(FieldGeometry((2, 1, 1)), np.ones(16)/4, duration=.125)


def test_general_binary_rational_parameters_are_exact_not_rounded():
    geometry = FieldGeometry((2, 1, 1), speed=3/16, interaction=7/8)
    report, op = certify(geometry)
    assert report['local_energy_denominator'] == 64
    for local, expected in zip(op.local, independent_local(geometry), strict=True):
        np.testing.assert_array_equal(local.matrix(op.scale).toarray(), expected)


def test_source_prepared_same_reference_energy_report_keeps_all_stages():
    from .__main__ import build
    result = build()
    read = result['prepared_reference_energy_read']
    assert read['source_preparation_executed']
    assert read['preparation_resources']['native_pulses'] == 10249
    assert abs(read['complete_preparation_evolution_read_probability']-1) < 2e-10
    assert read['predicted_pointer_energy_change'] > .001
    assert read['preparation_export_abort_probability'] >= 0
    assert not read['abort_output_energy_assigned']
    assert not result['physical_mass_or_gravity_premises_discharged']


@pytest.mark.parametrize('sides', ((2, 1, 1), (1, 2, 1), (1, 1, 2), (3, 1, 1), (2, 2, 1)))
@pytest.mark.parametrize('mu', (Fraction(-1, 3), Fraction(0), Fraction(1, 2)))
def test_exact_number_grand_energy_currents_and_changing_probes(sides, mu):
    bundle = grand_canonical(FieldGeometry(sides), mu)
    assert charged_probe_certificate(bundle)['charge_changing_CAR_probes_checked'] == 2*bundle.field.geometry.modes
    h = sum(independent_local(bundle.field.geometry))
    a = [ladder(bundle.field.geometry.modes, i) for i in range(bundle.field.geometry.modes)]
    n = [a[2*x].conj().T@a[2*x]+a[2*x+1].conj().T@a[2*x+1]
         for x in range(bundle.field.geometry.volume)]
    np.testing.assert_allclose(bundle.grand.total.matrix(bundle.grand.scale).toarray(),
                               h-float(mu)*sum(n), atol=2e-15)
    for (x, y), current in bundle.number_currents.items():
        # Independent local continuity identifies the oriented bond current.
        direct = 1j*(h@n[y]-n[y]@h)
        if len(n) == 2:
            np.testing.assert_allclose(current.matrix(bundle.field.scale).toarray(), direct, atol=2e-15)
        distance = sum(abs(p-q) for p, q in zip(bundle.field.geometry.sites[x], bundle.field.geometry.sites[y]))
        assert distance == 1
    # Every region has both exact number and exact grand-energy balance.
    for mask in range(1 << len(n)):
        region = tuple(x for x in range(len(n)) if mask >> x & 1)
        k_local, j_out = bundle.grand.region(region)
        assert (commutator(bundle.grand.total, k_local).times_i()+j_out).is_zero()


def test_dimer_charge_gaps_are_not_bare_H_addition_and_removal_costs():
    # Independent JW matrices and number-block eigensolvers avoid attaching
    # an arbitrary basis vector to a degenerate charged doublet.
    geometry = FieldGeometry((2, 1, 1))
    h = sum(independent_local(geometry))
    number = np.diag([mask.bit_count() for mask in range(16)])
    mu = .5
    k = h-mu*number
    eigenvalues, eigenvectors = np.linalg.eigh(k)
    ground = eigenvectors[:, 0]
    reference_H = np.vdot(ground, h@ground).real
    reference_N = np.vdot(ground, number@ground).real
    gap = math.sqrt(17)/8-1/16
    for sector, charge in ((1, -1), (3, 1)):
        indices = np.array([mask for mask in range(16) if mask.bit_count() == sector])
        values, vectors = np.linalg.eigh(k[np.ix_(indices, indices)])
        low = np.zeros((16, 2), complex)
        low[indices] = vectors[:, :2]
        final_H = np.trace(low.conj().T@h@low).real/2
        delta_H = final_H-reference_H
        delta_K = values[:2].mean()-eigenvalues[0]
        delta_N = sector-reference_N
        assert abs(delta_N-charge) < 1e-14
        assert abs(delta_K-gap) < 1e-14
        assert abs(delta_K-(delta_H-mu*delta_N)) < 1e-14
        assert abs(delta_H-(gap+mu*charge)) < 1e-14
        assert abs(delta_H-delta_K) > .49
    assert gap-mu < 0 < gap  # A positive hole K-gap need not cost positive H.


def test_chemical_potential_domain_rejects_inexact_or_unbounded_inputs():
    for mu in (.5, True, Fraction(1, 65537)):
        with pytest.raises(ValueError):
            grand_canonical(chemical_potential=mu)
