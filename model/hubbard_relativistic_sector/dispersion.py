"""Rational finite-rapidity bounds on the exact Hubbard Bessel dispersion."""
from fractions import Fraction as F
from functools import lru_cache

from ..common_source.mass_response import gap_certificate as arithmetic

I = arithmetic.RationalInterval


def _rational(value):
    if type(value) not in (int, F): raise ValueError('exact rational argument required')
    return F(value)


@lru_cache(maxsize=128)
def hyperbolic(value):
    value = _rational(value)
    negative = arithmetic.exp_negative_bounds(I.enclose(abs(value)))
    positive = 1/negative
    cosine, sine = (positive+negative)/2, (positive-negative)/2
    return cosine, sine if value >= 0 else -sine


def bessel_ratio(order, argument, terms=12):
    """Enclose K_order(z)/(sqrt(pi/(2z))*exp(-z)), order0 or1.

    Integrating Taylor's global one-sided remainder for (1+x)^(nu-1/2)
    against Gamma(nu+1/2) gives alternating adjacent truncation bounds.
    No floating Bessel evaluation or asymptotic equality is used.
    """
    z = _rational(argument)
    if type(order) is not int or order not in (0, 1) or z <= 0:
        raise ValueError('order0/1 and positive argument required')
    if type(terms) is not int or not 1 <= terms <= 32:
        raise ValueError('bounded positive Taylor order required')
    coefficient, total = F(1), F(1)
    for n in range(1, terms+1):
        coefficient *= F(4*order*order-(2*n-1)**2, 8*n)/z
        total += coefficient
    following = coefficient*F(4*order*order-(2*terms+1)**2, 8*(terms+1))/z
    lo, hi = min(total, total+following), max(total, total+following)
    # Intersect with independent global tangent bounds.
    if order == 0: lo, hi = max(lo, 1-F(1, 8)/z, F(0)), min(hi, F(1))
    else: lo, hi = max(lo, F(1)), min(hi, 1+F(3, 8)/z)
    if lo > hi: raise ArithmeticError('inconsistent exact Bessel bounds')
    return I(lo, hi)


def uniform_bounds(A, rapidity_radius):
    """All |theta|<=radius; E/M and vp/M errors from cosh and sinh."""
    A, radius = _rational(A), _rational(rapidity_radius)
    if A < 1 or radius < 0 or radius >= A:
        raise ValueError('A>=1 and0<=rapidity radius<A required')
    q = arithmetic.exp_negative_bounds(I.enclose(2*(A-radius)))
    energy_relative = (1+F(3, 8)/A)/(1-q)-1
    momentum_relative = I.enclose(max(F(1, 8)/A, (q/(1-q)).upper))
    cosine, sine = hyperbolic(radius)
    shell = cosine**2*(2*energy_relative+energy_relative**2)+sine**2*(2*momentum_relative+momentum_relative**2)
    return {'A': str(A), 'rapidity_radius': str(radius),
            'energy_error_over_M_cosh_radius': energy_relative.as_json(),
            'momentum_error_over_M_sinh_radius': momentum_relative.as_json(),
            'absolute_mass_shell_residual_over_M_squared': I(F(0), shell.upper).as_json()}


def point(A, theta, *, bessel_terms=12, odd_terms=8):
    """Exact enclosures of E/M and vp/M on the chosen running family."""
    A, theta = _rational(A), _rational(theta)
    if A < 1 or abs(theta) >= A or type(odd_terms) is not int or not 1 <= odd_terms <= 16:
        raise ValueError('A>=1,|theta|<A and bounded odd-series order required')
    radius = abs(theta); energy, momentum = I.enclose(0), I.enclose(0)
    for n in range(odd_terms):
        m = 2*n+1
        exponential = arithmetic.exp_negative_bounds(I.enclose((m-1)*A))
        weight = exponential/(m*arithmetic.sqrt_bounds(I.enclose(m)))
        cosine, sine = hyperbolic(m*radius)
        energy += weight*bessel_ratio(1, m*A, bessel_terms)*cosine
        momentum += weight*bessel_ratio(0, m*A, bessel_terms)*sine
    q = arithmetic.exp_negative_bounds(I.enclose(2*(A-radius)))
    cosine, sine = hyperbolic(radius)
    remainder = q**odd_terms/(1-q)
    energy = I(energy.lower, energy.upper+(cosine*(1+F(3, 8)/A)*remainder).upper)
    momentum = I(momentum.lower, momentum.upper+(sine*remainder).upper)
    if theta < 0: momentum = -momentum
    if theta == 0: momentum = I.enclose(0)
    # Interval squaring via multiplication is safe here, including negative p.
    square_p = I(min(momentum.lower**2, momentum.upper**2), max(momentum.lower**2, momentum.upper**2))
    if momentum.lower <= 0 <= momentum.upper: square_p = I(F(0), square_p.upper)
    shell = energy**2-square_p
    return {'A': str(A), 'theta': str(theta), 'E_over_M': energy.as_json(),
            'v_p_over_M': momentum.as_json(), 'mass_shell_over_M_squared': shell.as_json(),
            'Bessel_Taylor_terms': bessel_terms, 'odd_Bessel_terms': odd_terms}


def anchored_mass():
    """The stated flow through a=1,t=1/16,U=1; not a uniquely forced flow."""
    pi = arithmetic.pi_bounds()
    scale = arithmetic.exp_negative_bounds(pi/8)/pi
    return {'anchor': {'a': '1', 'v': '1/8', 't': '1/16', 'U': '1', 'A': 'pi/8'},
            'chosen_trajectory_rest_energy_exact': 'M0=exp(-pi/8)/pi',
            'chosen_trajectory_rest_energy': scale.as_json(),
            'chosen_trajectory_rest_mass_exact': 'm0=M0/v^2=64*exp(-pi/8)/pi',
            'chosen_trajectory_rest_mass': (64*scale).as_json(),
            'anchor_substitution': 'sqrt(2*pi*A0)=pi/2 =>a(A0)=1,t(A0)=1/16,U(A0)=1',
            'uniqueness_counterfamily': 'a_M(A)=[4v exp(-A)/(M sqrt(2pi A))]*[1+(M/M0-1)exp(-(A-A0))], A>=A0',
            'counterfamily_property': 'Every M>0 gives the same anchor and limit rest energy M; the bracket stays positive and tends to1',
            'trajectory_selected_by_unchanged_source': False}


def certificate():
    return {'exact_charge_dispersion': {
                'E': '(8t/pi) sum_(odd m>=1) K1(mA) cosh(m theta)/m',
                'relative_p': '(4/(pi*a)) sum_(odd m>=1) K0(mA) sinh(m theta)/m',
                'parameters': 'A=2pi*t/U,theta=A*sin(k),|k|<pi/2,v=2t*a',
                'centered_generator': 'K=H-(U/2)N; constant UL/4 cancels from excitation energies',
                'primary_source': 'https://arxiv.org/html/cond-mat/9410043'},
            'running_family': 'a=4v exp(-A)/(M sqrt(2pi A)),t=v/(2a),U=pi*v/(a*A)',
            'limit': {'E': 'M cosh(theta)', 'p': '(M/v) sinh(theta)',
                      'mass_shell': 'E^2-v^2 p^2=M^2', 'rest_mass': 'm=M/v^2'},
            'uniform_bounds': [uniform_bounds(A, 1) for A in (4, 8, 16, 32, 64, 128)],
            'finite_parameter_dispersion': [point(A, theta) for A in (4, 8, 16, 32, 64)
                                             for theta in (F(0), F(1, 2), F(1))],
            'anchored_trajectory': anchored_mass(),
            'scaling_trajectory_is_additional_input': True,
            'physical_SI_or_species_mass_calculated': False,
            'isolated_local_CAR_pole_or_finite_ring_single_quasiparticle_proved': False}
