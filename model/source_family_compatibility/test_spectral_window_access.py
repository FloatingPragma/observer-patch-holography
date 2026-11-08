"""Independent complete-source and error-budget controls for Lorentzian access."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import spectral_window_access as m


@pytest.fixture(scope='module')
def control():
    return m.finite_control()


@pytest.fixture(scope='module')
def numerical_source():
    h, rho, a, _ = m.dimer_source()
    h = np.asarray(h, dtype=complex)
    rho = np.asarray(rho, dtype=complex)
    a = np.asarray(a, dtype=complex)
    identity = np.eye(16)
    centered = h-identity/2
    clock = (identity-1j*centered/8)@np.linalg.inv(identity+1j*centered/8)
    power = identity.copy(); rows = []; maps = []
    weights = m.geometric_times(1, F(1, 4), 16)
    for j, weight in enumerate(weights['weights']):
        for label, op in (('addition', a.conj().T), ('removal', a)):
            evolved = power@op@power.conj().T
            for sign in (1, -1):
                k = np.sqrt(float(weight))*(op+sign*evolved)/2
                rows.append((j, label, sign)); maps.append(k)
        power = clock@power
    for label, op in (('addition', a.conj().T), ('removal', a)):
        rows.append((None, label, 0))
        maps.append(np.sqrt(float(weights['dummy_weight']))*op)
    return h, rho, a, rows, maps


def test_all_actual_source_records_are_complete_on_every_input(numerical_source, control):
    _, rho, _, rows, maps = numerical_source
    assert len(rows) == len(maps) == control['complete_record_count'] == 66
    effect = sum(k.conj().T@k for k in maps)
    assert np.linalg.norm(effect-np.eye(16), 2) < 2e-14
    assert np.linalg.norm(np.kron(effect, np.eye(3))-np.eye(48), 2) < 2e-14
    for row, meta, k in zip(control['all_records'], rows, maps, strict=True):
        assert (row['time_index'], row['charge'], row['sign']) == meta
        p = np.trace(rho@k.conj().T@k).real
        lo, hi = row['original_probability_interval']
        assert float(lo)-2e-14 <= p <= float(hi)+2e-14
    # Dropping the dummy or the other charge changes the original experiment.
    assert np.linalg.norm(sum(k.conj().T@k for k in maps[:-2])-np.eye(16), 2) > .02
    assert np.linalg.norm(sum(k.conj().T@k for row, k in zip(rows, maps, strict=True)
                              if row[1] == 'addition')-np.eye(16), 2) > .49


def test_actual_ground_score_compares_complete_spectral_measure(numerical_source, control):
    h, rho, a, rows, maps = numerical_source
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    ground = eigenvectors[:, 0]
    energies = eigenvalues-eigenvalues[0]
    weights = sum(np.abs(eigenvectors.conj().T@op@ground)**2 for op in (a, a.conj().T))
    exact = np.sum(weights/(1+energies**2))
    score = sum(row[2]*np.trace(rho@k.conj().T@k).real
                for row, k in zip(rows, maps, strict=True))
    lo, hi = control['exact_source_Lorentzian_interval']
    assert float(lo) <= exact <= float(hi)
    lo, hi = control['source_score_interval']
    assert float(lo) <= score <= float(hi)
    assert abs(exact-score) <= float(control['certified_score_error'])
    assert control['source_window']['original_total_charge_weight_lower'] > F(1, 3)
    assert weights[energies <= 1].sum() >= float(control['source_window']['original_total_charge_weight_lower'])


def test_excited_input_is_not_called_a_ground_energy_measure(numerical_source):
    h, _, a, rows, maps = numerical_source
    excited = np.eye(16)[:, 1]
    score = sum(row[2]*np.linalg.norm(k@excited)**2 for row, k in zip(rows, maps, strict=True))
    eig, vec = np.linalg.eigh(h)
    wrong = sum(np.sum(np.abs(vec.conj().T@op@excited)**2/(1+(eig-eig[0])**2))
                for op in (a, a.conj().T))
    assert abs(score-wrong) > .3


def test_geometric_kernel_error_is_global_and_retains_tail():
    for eta in (.125, 1., 7.):
        for dt in (.01, .25, 1.):
            r = 1/(1+eta*dt)
            for energy in np.r_[np.linspace(-20, 20, 401), 2*np.pi/dt]:
                q = ((1-r)/(1-r*np.exp(-1j*dt*energy))).real
                target = eta*eta/(eta*eta+energy*energy)
                assert abs(q-target) <= dt*abs(energy)/2+1e-13
                finite = sum((1-r)*r**j*np.cos(j*dt*energy) for j in range(16))
                assert abs(q-finite) <= r**16+1e-13
    t = m.geometric_times(1, F(1, 4), 16)
    assert sum(t['weights'], F())+t['dummy_weight'] == 1
    assert t['complete_records'] == 66


def test_complete_boundary_identity_includes_arbitrary_reference(numerical_source):
    h, _, a, _, _ = numerical_source
    eig, vec = np.linalg.eigh(h)
    u = (vec*np.exp(-.7j*eig))@vec.conj().T
    v = (vec*np.exp(-.4j*eig))@vec.conj().T
    left = np.zeros((16, 16), complex); right = left.copy()
    for op in (a, a.conj().T):
        difference = u@op@u.conj().T-v@op@v.conj().T
        right += difference.conj().T@difference/2
        for sign in (1, -1):
            delta = sign*difference/2
            left += delta.conj().T@delta
    assert np.linalg.norm(left-right) < 1e-15
    assert np.linalg.norm(np.kron(left-right, np.eye(2))) < 2e-15


def test_cayley_error_checked_against_full_source_clock(numerical_source):
    h = numerical_source[0]-np.eye(16)/2
    eig, vec = np.linalg.eigh(h)
    for n in (1, 2, 4):
        dt = F(1, 4)
        one = (np.eye(16)-1j*float(dt)*h/(2*n))@np.linalg.inv(np.eye(16)+1j*float(dt)*h/(2*n))
        actual = np.linalg.matrix_power(one, 15*n)
        ideal = (vec*np.exp(-1j*float(15*dt)*eig))@vec.conj().T
        bound = m.cayley_service_bound(15*dt, dt, F(5, 8), n)
        assert np.linalg.norm(actual-ideal, 2) < float(bound['unitary_operator_error'])


def test_one_pressure_observable_respects_magnetic_cell_density():
    # Atomic centered source: every local ground is singly occupied.
    # Its algebraic projection observable is one on that entire ground space.
    x = m.pressure_score_interval(8, 0, F(1, 4), (0, 0), (16, 16), (-16, -16), (0, 0, 0), stride=2)
    assert x['perturbation_density'] == F(1, 8)
    assert x['score_interval'] == (F(1), F(1))
    # A two-site-span diagnostic misses the last cell along that direction.
    y = m.pressure_score_interval(8, 0, F(1, 4), (0, 0), (12, 12), (-12, -12), (2, 0, 0), stride=2)
    assert y['omitted_translate_fraction'] == F(1, 4)
    assert y['score_interval'] == (F(1, 2), F(1))
    # An atomic spin sign distinguishes degenerate grounds, so it cannot certify
    # a universal positive expectation just by selecting a preferred one.
    z = m.pressure_score_interval(8, 0, F(1, 4), (0, 0), (-16, -16), (-16, -16), (0, 0, 0), stride=2)
    assert z['score_interval'] == (F(-1), F(1))


def test_pressure_to_original_window_is_one_joined_certificate():
    # The actual c=U=0 CAR source has a zero clock and its complete Ramsey
    # score is (1-dummy_weight)I. Its full ground space is retained.
    z = 1-m.geometric_times(1, F(1, 100), 1000)['dummy_weight']
    x = m.pressure_window_certificate(8, 0, F(1, 4), (0, 0), (16*z, 16*z),
        (-16*z, -16*z), (0, 0, 0), stride=2, eta=1, step=F(1, 100),
        count=1000, second_moment_upper=0, first_moment_upper=0,
        radius=2, support_radius=0, support_sites=1, clock_service_error=0,
        energy_upper=1)
    assert x['strict_positive_original_weight_certified']
    assert x['window']['original_total_charge_weight_lower'] > F(9, 10)
    assert x['score']['ground_score_certified_by_pressure']
    assert not x['score']['same_ground_marginal_required']
    assert not x['score']['a_finite_prepared_state_marginal_has_been_identified']
    assert not x['thermodynamic_ground_preparation_claimed']


def test_default_fixed_source_pressure_instance_is_explicit():
    x = m.default_terminating_instance()
    assert x['local_second_moment'] == F(35, 128)
    assert x['exact_Jensen_Lorentzian_lower'] == F(128, 163)
    assert x['analytic_Markov_original_window_lower'] == F(93, 128)
    assert not x['operational_floor_improves_analytic_moment_bound']
    assert x['time_error_upper'] < F(1, 256)
    assert x['uniform_score_lower'] > F(3, 4)
    assert x['published_strict_original_window_floor'] == F(15, 32)
    assert x['cut_radius'] == 124 and x['complete_records'] == 8194
    assert x['Cayley_subdivisions_per_time_step'] == 3047219812
    assert x['algebraic_complete_clock_error'] <= F(1, 512)
    assert x['pressure_search_termination_proved'] and not x['pressure_search_executed']
    assert not x['native_large_clock_tape_executed'] and not x['shrinking_window_uniform_weight_proved']


def test_window_conversion_preserves_total_mass_and_includes_kernel():
    assert m.window_from_score((F(3, 4), 1), 1, 1)['original_total_charge_weight_lower'] == F(1, 2)
    assert m.window_from_score((F(1, 4), F(1, 2)), 1, 1)['original_total_charge_weight_lower'] == 0
    x = m.window_from_score((1, 1), 1, 1)
    assert x['zero_energy_component_allowed'] and not x['separate_charge_equality_assumed']
    assert x['physical_mass_assigned'] is None


def test_complete_errors_enter_before_any_record_selection():
    base = dict(measured=(F(3, 4), F(4, 5)), eta=1, step=F(1, 100), count=1000,
                second_moment_upper=1, first_moment_upper=1, radius=256,
                support_radius=0, support_sites=1, coupling=F(1, 8))
    x = m.complete_score_interval(**base)
    y = m.complete_score_interval(**base, marginal_error=F(1, 100), preparation_error=F(1, 100),
                                   service_error=F(1, 100), repair_error=F(1, 100), packet_error=F(1, 100))
    assert y['original_score_error']-x['original_score_error'] == F(1, 10)
    assert y['original_Lorentzian_interval'][0] < x['original_Lorentzian_interval'][0]
    assert not x['conditional_normalization_applied']


def test_clock_homogeneity_preserves_exact_budgets():
    x = m.time_error(F(1, 4), F(1, 8), 40, F(25, 64), F(5, 8))
    y = m.time_error(2, F(1, 64), 40, 25, 5)
    assert x == y
    assert m.cayley_service_bound(5, F(1, 8), F(5, 8), 7) == m.cayley_service_bound(F(5, 8), F(1, 64), 5, 7)


def test_fraction_and_container_cache_detachment(control):
    before = repr(control)
    out = m.finite_control()
    out['source_c']._numerator = 999
    out['all_records'][0]['original_probability_interval'][0]._numerator = 999
    out['all_records'].clear()
    assert repr(m.finite_control()) == before


@pytest.mark.parametrize('args', [(0, 1, 4), (1, 0, 4), (1, 1, 0), (True, 1, 4), (1., 1, 4)])
def test_invalid_time_settings_rejected(args):
    with pytest.raises(ValueError):
        m.geometric_times(*args)


def test_invalid_moment_and_pressure_inputs_rejected():
    with pytest.raises(ValueError):
        m.time_error(1, 1, 4, 2, 1)
    with pytest.raises(ValueError):
        m.pressure_score_interval(6, 1, 1, (0, 0), (0, 0), (0, 0), (0, 0, 0), stride=4)
    with pytest.raises(ValueError):
        m.pressure_score_interval(8, 1, 1, (0, 0), (0, 0), (0, 0), (8, 0, 0))
    with pytest.raises(ArithmeticError):
        m.pressure_score_interval(8, 0, 1, (0, 0), (1000, 1000), (-1000, -1000), (0, 0, 0))
