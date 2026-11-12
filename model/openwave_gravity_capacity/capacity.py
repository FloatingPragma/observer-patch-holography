"""Amplitude-free consequences of the conditional capacity density law."""
from fractions import Fraction as F
import sympy as s


def capacity_from_expansion(a, expansion_squared, omega_m, omega_r, omega_dark):
    """Return N(a)/N(1) from a declared flat Friedmann comparison history."""
    a, e2, om, radiation, od = map(F, (a, expansion_squared, omega_m, omega_r, omega_dark))
    if a <= 0 or e2 <= 0 or min(om, radiation) < 0 or od <= 0:
        raise ValueError('positive expansion/capacity branch and nonnegative matter required')
    if om+radiation+od != 1:
        raise ValueError('flat reference fractions must sum to one')
    remainder = e2-om/a**3-radiation/a**4
    if remainder <= 0:
        raise ValueError('dark density residual must be strictly positive')
    return od/remainder


def expansion_from_capacity(a, capacity_ratio, omega_m, omega_r, omega_dark):
    a, ratio, om, radiation, od = map(F, (a, capacity_ratio, omega_m, omega_r, omega_dark))
    if a <= 0 or ratio <= 0 or min(om, radiation) < 0 or od <= 0:
        raise ValueError('positive capacity history required')
    if om+radiation+od != 1:
        raise ValueError('flat reference fractions must sum to one')
    return om/a**3+radiation/a**4+od/ratio


def two_epoch_density_ratio(n1, n2):
    n1, n2 = F(n1), F(n2)
    if min(n1, n2) <= 0:
        raise ValueError('positive capacities required')
    return n2/n1


def symbolic_certificate():
    a, n0, kappa = s.symbols('a N0 kappa', positive=True)
    exponent = s.symbols('epsilon', real=True)
    n = n0*a**exponent
    rho = kappa/n
    w = s.simplify(-1-a*s.diff(rho, a)/(3*rho))
    if s.simplify(w-(-1+exponent/3)) != 0:
        raise ArithmeticError('capacity drift sign changed')
    c = s.symbols('c', positive=True)
    scale_degeneracy = s.simplify((c*kappa)/(c*n)-rho) == 0
    return {
        'density_law': 'rho_DE(a)=kappa/N(a)',
        'continuity_assumption': 'd rho/d ln(a)=-3(1+w)rho',
        'drift_law': 'w=-1+(1/3)d ln(N)/d ln(a)',
        'power_family': 'N(a)=N0*a^epsilon',
        'power_family_w': str(w),
        'fixed_capacity_w': '-1',
        'fixed_capacity_CPL': ['-1', '0'],
        'monotone_capacity_implies_no_phantom': True,
        'counterexample_decreasing_capacity_exponent': '-1',
        'counterexample_w': '-4/3',
        'simultaneous_N_kappa_rescaling_invisible_to_density': scale_degeneracy,
        'capacity_amplitude_selected': False,
        'physical_horizon_or_Friedmann_dynamics_derived_here': False,
    }


def certificate():
    om, radiation, od = F(1, 4), F(1, 16), F(11, 16)
    rows = []
    # These rational profiles are algebraic controls; no measured data enter.
    for exponent in (-1, 0, 1, 2):
        for a in (F(1, 2), F(1), F(2)):
            nratio = a**exponent
            e2 = expansion_from_capacity(a, nratio, om, radiation, od)
            recovered = capacity_from_expansion(a, e2, om, radiation, od)
            if recovered != nratio:
                raise ArithmeticError('capacity/expansion inverse failed')
            rows.append({'a': str(a), 'capacity_exponent': exponent,
                         'capacity_ratio': str(nratio), 'E_squared': str(e2),
                         'recovered_capacity_ratio': str(recovered),
                         'w': str(-1+F(exponent, 3))})
    return {
        'symbolic': symbolic_certificate(),
        'epoch_invariant': 'N(a)*(E(a)^2-Omega_m/a^3-Omega_r/a^4)=Omega_DE*N(1)',
        'normalization_free_inputs': ['scale-factor ratio', 'matter and radiation fractions',
                                     'dimensionless expansion history'],
        'absolute_cosmic_capacity_or_Lambda_computed': False,
        'comparison_data_read': False,
        'algebraic_control_parameters': {'Omega_m': str(om), 'Omega_r': str(radiation),
                                         'Omega_DE': str(od)},
        'control_histories': rows,
        'de_sitter_conditional_relation': 'Lambda*ell_P^2=3*pi/N',
        'de_sitter_remaining_inputs': ['source-selected global record capacity',
                                      'same-carrier horizon/record identification',
                                      'independent physical Planck scale'],
        'fixed_N_does_not_identify_N_value': True,
        'raw_dimension_or_local_record_count_is_cosmic_capacity': False,
    }
