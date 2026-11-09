"""Independent source, Hamilton-equation and finite-lattice controls."""
from fractions import Fraction as F
import itertools

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import characteristics as c


def test_all_twelve_native_generators_and_electric_casimir():
    d = c.native_generators()
    assert len(d['generators']) == 12
    assert d['gram'] == s.diag(*([2]*11+[s.Rational(10, 3)]))
    assert s.Matrix(d['native_port_coefficients']).rank() == 12
    for t in d['generators']:
        assert t == t.H and s.trace(t) == 0
    assert sum(x*x for x in (6*d['generators'][-1]).diagonal()) == 120
    assert s.trace(d['casimir']) == 12


def test_all_gauge_directions_in_ordered_plaquette_hessian():
    d = c.native_generators()
    for j, t in enumerate(d['generators']):
        result = c.plaquette_hessian((t, s.zeros(15), s.zeros(15), s.zeros(15)))
        assert result['cost_second_derivative'] == d['gram'][j, j]/15
    # Noncommuting link variables: each pair contributes its mixed trace,
    # including the cancellation for a linearized pure gauge plaquette.
    x, y = d['generators'][0], d['generators'][1]
    assert x*y-y*x != s.zeros(15)
    assert c.plaquette_hessian((x, y, -x, -y))['quadratic_cost_coefficient'] == 0
    result = c.plaquette_hessian((x+y, 2*x-y, x, -y))
    assert result['quadratic_cost_coefficient'] == s.Rational(17, 15)


def test_ordered_exponential_remainder_with_noncommuting_native_links():
    gs = c.native_generators()['generators']
    xs = (gs[0]+gs[8]/3, gs[1]/2+gs[11], -gs[2]/5, -gs[9]/7)
    coefficient = float(c.plaquette_hessian(xs)['quadratic_cost_coefficient'])
    norm_upper = sum(max(sum(abs(x[i, j]) for j in range(15)) for i in range(15)) for x in xs)
    epsilon = F(1, 100)
    product = np.eye(15, dtype=complex)
    for x in xs:
        product = product@expm(1j*float(epsilon)*np.array(x, dtype=complex))
    actual = 1-np.trace(product).real/15
    error = abs(actual-float(epsilon)**2*coefficient)
    bound = c.weak_field_remainder(F(str(norm_upper))*epsilon)
    assert error <= float(bound)+2e-15
    assert float(bound) < 1e-4
    assert c.weak_field_remainder(0) == 0


def test_nonabelian_zero_momentum_directions_keep_quartic_interactions():
    gs = c.native_generators()['generators']
    result = c.constant_nonabelian_loop(gs[8], gs[9])
    assert result['linear_through_cubic_cost_coefficients'] == [0,0,0]
    assert result['quartic_cost_coefficient'] == s.Rational(1,15)
    assert result['nonzero_nonlinear_magnetic_interaction']
    assert not result['interacting_particle_scattering_derived']
    # Hypercharge commutes with every nonabelian source direction.
    assert c.constant_nonabelian_loop(gs[11],gs[8])['quartic_cost_coefficient'] == 0


def test_inverse_trace_metric_and_hamilton_equations_fix_factor_fifteen():
    a, e, k, index, q, p = s.symbols('a e k index q p', positive=True)
    # The mode has forward-difference squared norm2, e.g. z_x=i.
    electric = e*p*p/(2*a*index)
    magnetic = k*index*2*q*q/(30*a)
    dq = s.diff(electric, p)
    dp = -s.diff(magnetic, q)
    ddq = s.diff(dq, p)*dp
    assert s.simplify(ddq+2*e*k*q/(15*a*a)) == 0
    inherited = c.clock_family()
    assert set(inherited['all_twelve_speed_squared']) == {F(1, 105)}
    assert not inherited['same_clock_gauge_Dirac_characteristic_match']
    assert c.clock_family('1/100')['all_twelve_speed_squared'] == [F(1, 1500)]*12
    matched = c.clock_family(electric_scale='105')
    assert matched['all_twelve_speed_squared'] == [1]*12
    assert matched['same_clock_gauge_Dirac_characteristic_match']


def test_coordinate_mixing_keeps_all_trace_cross_terms():
    gram = c.native_generators()['gram']
    change = s.eye(12)
    change[0, 11] = s.Rational(2, 3)
    change[8, 3] = s.Rational(5, 7)
    changed = change.T*gram*change
    assert changed[0, 11] != 0 and changed[3, 8] != 0
    assert changed.inv()*changed == s.eye(12)
    assert changed.inv() == change.inv()*gram.inv()*change.inv().T
    # Treating the changed diagonal entries as the full metric is false.
    wrong_inverse = s.diag(*(1/changed[i, i] for i in range(12)))
    assert wrong_inverse*changed != s.eye(12)


def test_charge_coordinate_rescaling_preserves_canonical_pair_and_energy():
    data = c.charge_normalization()
    assert data['trace_index_rescaled'] == 120
    assert data['electric_coefficient_rescaled'] == F(1, 120)
    theta, pi, scale, index = s.symbols('theta pi scale index', positive=True)
    assert s.diff(theta/scale, theta)*s.diff(scale*pi, pi) == 1
    assert s.simplify((scale*pi)**2/(scale**2*index)-pi*pi/index) == 0
    assert s.simplify(scale**2*index*(theta/scale)**2-index*theta**2) == 0


@pytest.mark.parametrize('phases',[(s.I,-1,1),(-1,-1,-1),
                                  ((3+4*s.I)/5,(5-12*s.I)/13,1),(1,1,1)])
def test_lattice_curl_keeps_longitudinal_gauss_null_and_zero_modes(phases):
    result = c.lattice_symbol(phases)
    d, K = result['difference'], result['hessian']
    assert c.clean(K*d) == s.zeros(3, 1)
    assert c.clean(d.H*K) == s.zeros(1, 3)  # Gauss constraint is conserved.
    expected = 0 if phases == (1,1,1) else 2
    assert result['rank'] == expected
    assert result['zero_momentum_flat_coordinates'] == (3 if expected == 0 else 0)


def test_exact_real_space_curl_spectrum_on_periodic_three_cube():
    cells = list(itertools.product(range(3), repeat=3))
    index = {x:i for i,x in enumerate(cells)}
    derivatives = []
    for axis in range(3):
        D = s.zeros(27)
        for x, row in index.items():
            y = list(x); y[axis] = (y[axis]+1)%3
            D[row, index[tuple(y)]] = 1
            D[row, row] -= 1
        derivatives.append(D)
    Dx,Dy,Dz = derivatives
    zero = s.zeros(27)
    curl = s.BlockMatrix([[zero,-Dz,Dy],[Dz,zero,-Dx],[-Dy,Dx,zero]]).as_explicit()
    gradient = s.Matrix.vstack(*derivatives)
    K = curl.T*curl
    assert curl*gradient == s.zeros(81,27)
    assert K.to_DM().rank() == 2*(27-1)
    eigenvalues = np.linalg.eigvalsh(np.asarray(K,dtype=float))
    expected = [0]*29+[3]*12+[6]*24+[9]*16
    assert np.allclose(eigenvalues, expected, atol=2e-13, rtol=0)


def test_native_dirac_and_gauge_lattice_errors_against_direct_spectra():
    assert c.free_symbol((0,0,0))['energy_squared'] == 0
    assert c.free_symbol((s.pi/2,0,0))['energy_squared'] == 2
    for p in ((1,0,0),(0.2,-0.7,1.1),(2,1,-0.5)):
        p = np.asarray(p)
        norm2 = np.dot(p,p)
        for a in (0.5,0.125,0.03125):
            k = a*p
            gauge = np.sum(4*np.sin(k/2)**2)/a**2
            wilson = (np.sum(np.sin(k)**2)+np.sum(1-np.cos(k))**2)/a**2
            assert -2e-13 <= norm2-gauge <= a*a*norm2**2/12+2e-13
            assert abs(wilson-norm2) <= 7*a*a*norm2**2/12+2e-13
            group_velocity2 = np.sum(np.sin(k)**2)/np.sum(4*np.sin(k/2)**2)
            assert group_velocity2 <= 1+2e-14


def test_constant_coframe_has_same_gauge_and_dirac_principal_metric():
    k = s.Matrix([1,2,-1])
    frames = (s.eye(3),2*s.eye(3),
              s.Matrix([[1,s.Rational(1,3),0],[0,1,s.Rational(1,5)],[0,0,1]]))
    for frame in frames:
        data = c.coframe_principal(frame,k)
        K = data['gauge_principal_matrix_at_matched_clock']
        q = data['characteristic_squared']
        t = s.symbols('t')
        assert s.factor((t*s.eye(3)-K).det()) == s.factor(t*(t-q)**2)
        assert K*k == s.zeros(3,1)
        assert K.rank() == 2
    assert c.coframe_principal(2*s.eye(3),k)['characteristic_squared'] == s.Rational(3,2)


@pytest.mark.parametrize('kwargs',[{'electric_scale':0},{'magnetic_scale':-1},
                                  {'electric_scale':1.0},{'magnetic_scale':'0'}])
def test_invalid_clock_family_is_rejected(kwargs):
    with pytest.raises(ValueError): c.clock_family(**kwargs)


def test_certificate_does_not_promote_quadratic_comparison_to_particles():
    report = c.certificate()
    assert report['native_basis']['all_twelve_generators_used']
    assert report['scope']['same_native_magnetic_Hessian_and_electric_metric']
    assert not report['scope']['weak_field_localization_of_ground_state_proved']
    assert not report['scope']['nonperturbative_interacting_relativistic_limit_proved']
    assert not report['scope']['named_particle_mass_or_running_coupling_prediction']
