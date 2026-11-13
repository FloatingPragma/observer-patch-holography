"""Exact free symbols of the declared Pauli-root/Wilson prescription.

The Maxwell symbol is the quadratic small-field symbol of the same electric
and magnetic coefficients. This calculation does not replace compact-field
interactions by that quadratic approximation.
"""
from fractions import Fraction as F
from itertools import product
import numpy as np
from .operators import SPEED_SQUARED


def symbols(momentum, spacing, mass, regulator=F(1, 16)):
    p = np.asarray(momentum, dtype=float)
    a, m, r = float(spacing), float(mass), float(regulator)
    if p.shape != (3,) or a <= 0 or m <= 0 or r < 0:
        raise ValueError('three momenta, positive spacing/mass, nonnegative regulator')
    c = float(SPEED_SQUARED)**.5
    momentum_square = c*c*np.dot(np.sin(a*p), np.sin(a*p))/(a*a)
    e0 = np.sqrt(m*m+momentum_square)
    wilson = r*np.sum(1-np.cos(a*p))/a
    photon = 2*c*np.linalg.norm(np.sin(a*p/2))/a
    target = np.sqrt(m*m+c*c*np.dot(p,p))
    pnorm = np.linalg.norm(p)
    return {'positive_energy': e0+wilson, 'negative_comparator_energy': -e0-wilson,
            'photon_quadratic_energy': photon, 'continuum_Dirac_energy': target,
            'continuum_photon_energy': c*pnorm, 'Wilson_energy': wilson,
            'Dirac_absolute_error_bound': c*a*a*pnorm**3/6+r*a*pnorm**2/2,
            'Maxwell_absolute_error_bound': c*a*a*pnorm**3/24}


def corner_records(spacing, mass, regulator=F(1, 16)):
    """All eight corners; the exact lift is 2 r times the corner weight / a."""
    a, m, r = F(spacing), F(mass), F(regulator)
    if a <= 0 or m <= 0 or r < 0:
        raise ValueError('positive spacing/mass, nonnegative regulator')
    return [{'corner': list(k), 'mass_plus_lift_exact': str(m+2*r*sum(k)/a),
             'lift_exact': str(2*r*sum(k)/a)} for k in product((0,1), repeat=3)]


def refinement_records(masses):
    result = []
    for name, mass in masses:
        for a in (F(1), F(1,2), F(1,4), F(1,8)):
            for p in ((.1,.2,.3), (.4,0.,0.), (.25,.25,.25)):
                values = symbols(p,a,mass)
                assert abs(values['positive_energy']-values['continuum_Dirac_energy']) <= values['Dirac_absolute_error_bound']+1e-12
                assert abs(values['photon_quadratic_energy']-values['continuum_photon_energy']) <= values['Maxwell_absolute_error_bound']+1e-15
                result.append({'family': name, 'spacing_exact': str(a), 'momentum': list(p), **values})
    return {'common_speed_squared_exact': str(SPEED_SQUARED),
            'same_electric_magnetic_normalization': True,
            'records': result,
            'Wilson_corner_records': corner_records(1,1),
            'zero_regulator_corner_control': corner_records(1,1,0),
            'scope': 'free Pauli-root matter and quadratic Maxwell symbols; no interacting Lorentz or signed-Dirac transfer'}
