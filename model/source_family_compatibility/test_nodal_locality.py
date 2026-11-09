"""Source, frame, probability and uniform-scaling controls for nodal locality."""
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest
import sympy as s

from . import nodal_locality as m
from .hole_dispersion import gamma_matrices
from .local_spectral import dimer_source
from .local_matter import annihilator


def test_integer_logarithm_is_outward_at_exact_thresholds():
    for exponent in range(-80, 81):
        power = F(2**exponent) if exponent >= 0 else F(1, 2**(-exponent))
        assert m.ceil_log2(power) == exponent
        assert m.ceil_log2(power*F(999, 1000)) == exponent
        assert m.ceil_log2(power*F(1001, 1000)) == exponent+1


@pytest.mark.parametrize('call', [
    lambda: m.ceil_log2(0), lambda: m.ceil_log2(.5),
    lambda: m.fourier_tail(True, 1), lambda: m.fourier_tail(0, 1),
    lambda: m.fourier_cutoff(1, 0), lambda: m.fourier_cutoff(1, 2),
    lambda: m.lorentzian_difference(-1, 1), lambda: m.lorentzian_difference(1, 0),
    lambda: m.lorentzian_difference(1, 1, accuracy=2),
    lambda: m.conditional_window(0, 1), lambda: m.conditional_window(F(10, 17), 1),
    lambda: m.conditional_window(F(1, 10), 0),
    lambda: m.conditional_window(F(9, 17), 1),
])
def test_rejects_unproved_domains_and_inexact_inputs(call):
    with pytest.raises(ValueError):
        call()


def test_actual_cut_bound_implies_matrix_tail_without_clustering():
    # Independent exact infinite shell sums, rather than a copy of the producer.
    n = s.symbols('n', integer=True, positive=True)
    for radius in (1, 2, 4, 9, 23):
        exact_tail = s.summation(n**4/s.Integer(4)**n, (n, radius+1, s.oo))
        assert exact_tail <= s.Rational(21, 2**radius)
    for n in range(2, 40):
        # 2B, 16 matrix entries per row and the actual shell count.
        factor = F(2*6*16*(24*n*n+2)*(4*n-3)**2, 2**(2*n-2))
        assert factor <= F(319488*n**4, 4**n)
    for c in (F(1, 8), F(1), F(3, 2)):
        for time in (F(0), F(1, 100), F(1, 4), F(2)):
            result = m.fourier_cutoff(time, F(1, 1000), c)
            assert result['original_response_tail_upper'] <= F(1, 1000)
            analytic = m.FOURIER_TAIL_CONSTANT*math.exp(49*float(c*time))*2.**(-result['cell_radius'])
            assert analytic <= .001*(1+1e-13)
    assert m.fourier_tail(1, 0) == 0


def test_matrix_bernstein_interpolation_is_an_operator_identity():
    for degree in range(1, 13):
        theta = (2*np.arange(1, 2*degree+1)-1)*np.pi/(2*degree)
        weights = (-1.)**np.arange(2*degree)/(4*degree*np.sin(theta/2)**2)
        assert abs(np.abs(weights).sum()-degree) < 2e-12
        for frequency in range(-degree, degree+1):
            derived = np.dot(weights, np.exp(1j*frequency*theta))
            assert abs(derived-1j*frequency) < 4e-12
        rng = np.random.default_rng(degree)
        coefficients = rng.normal(size=(2*degree+1, 3, 3))+1j*rng.normal(size=(2*degree+1, 3, 3))
        frequencies = np.arange(-degree, degree+1)
        values = np.einsum('jl,lmn->jmn', np.exp(1j*np.outer(theta, frequencies)), coefficients)
        exact_derivative = np.einsum('l,lmn->mn', 1j*frequencies, coefficients)
        assert np.linalg.norm(np.einsum('j,jmn->mn', weights, values)-exact_derivative) < 2e-11
        assert np.linalg.norm(exact_derivative, 2) <= degree*max(np.linalg.norm(x, 2) for x in values)+2e-11


def test_actual_magnetic_frame_is_glued_not_falsely_periodic():
    gamma = gamma_matrices()
    parities = np.array([(i//4, (i//2)%2, i%2) for i in range(8)])
    k = np.array([.21, -.37, .52])
    def physical(point):
        return sum(np.cos(point[j])*gamma[j] for j in range(3))
    def frame(point):
        return np.diag(np.exp(-1j*parities@point))
    def periodic(point):
        d = frame(point)
        return d.conj().T@physical(point)@d
    for axis in range(3):
        shifted = k.copy(); shifted[axis] += np.pi
        gluing = np.diag((-1)**parities[:, axis])
        assert np.linalg.norm(physical(shifted)-gluing@physical(k)@gluing) < 2e-15
        assert np.linalg.norm(periodic(shifted)-periodic(k)) < 5e-15
        # Treating the physical cos-Gamma matrix as periodic loses source transport.
        assert np.linalg.norm(physical(shifted)-physical(k), 2) > 1.5


@pytest.fixture(scope='module')
def finite():
    return m.finite_control()


def test_complete_interacting_source_has_contraction_and_nonunitary_memory(finite):
    h, rho, _, _ = dimer_source()
    h, rho = np.array(h, dtype=complex), np.array(rho, dtype=complex)
    modes = [np.array(annihilator(2, j), dtype=complex) for j in range(4)]
    eigenvalues, vectors = np.linalg.eigh(h)
    clock = (vectors*np.exp(1j*eigenvalues*.7))@vectors.conj().T
    response = np.empty((4, 4), dtype=complex)
    for i, ai in enumerate(modes):
        evolved = clock@ai@clock.conj().T
        for j, aj in enumerate(modes):
            response[i, j] = np.trace(rho@(evolved@aj.conj().T+aj.conj().T@evolved))
    assert np.linalg.norm(response, 2) <= 1+1e-13
    assert np.linalg.eigvalsh(np.eye(4)-response.conj().T@response).min() > .1
    assert finite['retarded_second_moment_minus_first_square'] == F(1, 4)
    assert finite['actual_Fock_dimension'] == 16


def test_original_measure_keeps_both_charges_and_uses_same_ground(finite):
    h, _, a, _ = dimer_source()
    h, a = np.array(h, dtype=complex), np.array(a, dtype=complex)
    e, v = np.linalg.eigh(h); ground = v[:, 0]
    energies = e-e[0]
    add = np.abs(v.conj().T@a.conj().T@ground)**2
    remove = np.abs(v.conj().T@a@ground)**2
    assert abs(add.sum()-.5) < 2e-14 and abs(remove.sum()-.5) < 2e-14
    assert abs((add+remove).sum()-1) < 2e-14
    q = np.sum((add+remove)/(1+energies**2))
    lower, upper = finite['original_Lorentzian_interval']
    assert float(lower) <= q <= float(upper)
    assert abs(np.sum(add/(1+energies**2))-q/2) < 2e-14
    assert not finite['finite_dimer_replaces_thermodynamic_reference']


def test_uniform_dyadic_envelope_handles_logarithm_jumps():
    for power in range(1, 50, 3):
        maximum = F(1, 2**power)
        for ratio in [F(1), F(999, 1000), F(1, 2), F(1, 3), F(1, 10000)]:
            d = maximum*ratio
            direct = d*(8+4*m.ceil_log2(F(m.FOURIER_TAIL_CONSTANT)/d))
            assert direct <= m.dyadic_log_remainder(power)
        # Probe just below every possible ceiling discontinuity in this range.
        for shift in range(1, 10):
            d = F(m.FOURIER_TAIL_CONSTANT, 2**(power+23+shift))*F(999999, 1000000)
            if d <= maximum:
                assert d*(8+4*m.ceil_log2(F(m.FOURIER_TAIL_CONSTANT)/d)) <= m.dyadic_log_remainder(power)


def test_conditional_plan_has_uniform_original_and_measured_margins():
    for z in (F(1, 1000), F(1, 100), F(1, 10), F(9, 20)):
        for k in (F(1, 3), F(2)):
            plan = m.conditional_window(z, k, F(3, 2))
            for ratio in (F(1), F(9, 10), F(1, 2), F(1, 3), F(1, 10000)):
                spacing = plan['spacing_upper']*ratio
                eta = spacing*plan['physical_Lorentzian_scale']/F(3, 2)
                actual_bound = m.lorentzian_difference(spacing*k, eta)
                assert actual_bound <= plan['joint_frequency_momentum_error_upper'] <= z/4
            assert plan['original_window_weight_lower'] >= z/2
            assert plan['inferred_original_window_weight_lower'] >= z/4
            assert not plan['actual_nodal_atom_lower_proved']
            assert plan['physical_energy_window'][0] == 0


def test_weak_continuity_alone_does_not_supply_the_locality_modulus():
    # A comparison measure delta_(sqrt(|k|)) is weakly continuous with a node
    # atom, but violates the source modulus. It is not a claimed source state.
    d = F(1, 10**16); eta = 10**6*d
    defect = 1-eta*eta/(eta*eta+d)
    assert defect > F(99, 100)
    assert m.lorentzian_difference(d, eta) < F(1, 1000)
    # A fixed gapped scalar measure has no nodal atom and loses every bounded
    # physical window when energy is multiplied by 1/a.
    for spacing in (F(1, 100), F(1, 1000)):
        assert 1/spacing > 10


def test_certificate_is_compact_detached_and_does_not_assert_atom_or_mass():
    from .packet import serialize
    first = m.certificate()
    expected = json.dumps(serialize(first), sort_keys=True)
    assert len(expected) < 12000
    first['actual_source_control']['original_Lorentzian_interval'][0]._numerator = -999
    assert json.dumps(serialize(m.certificate()), sort_keys=True) == expected
    assert first['actual_nodal_atom_positivity_proved'] is False
    assert first['physical_mass_assigned'] is None
    assert first['sharp_dispersion_proved'] is False
    assert first['growing_native_instrument_executed'] is False


def test_atom_input_is_checked_against_exact_paired_star_bound():
    scale = 2**40
    low = F(math.isqrt(577*scale*scale), scale)
    high = low+F(1, scale)
    assert low*low < 577 < high*high
    lower = (31+low)/(95+low)
    upper = (31+high)/(95+high)
    assert lower < upper < F(7, 15)
    assert m.conditional_window(lower, 1)['assumed_nodal_atom_lower'] == lower
    with pytest.raises(ValueError):
        m.conditional_window(upper, 1)
    # Unsigned squaring would incorrectly reject sufficiently small floors.
    assert m.conditional_window(F(1, 10000), 1)['actual_nodal_atom_lower_proved'] is False
