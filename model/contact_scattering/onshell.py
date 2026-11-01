"""Zero-total-momentum contact threshold of the unchanged Weyl/contact field.

The rational certificate encloses the *infinite lattice* free contact density
for 0 < |E| < 2c, with all eight Weyl nodes. It does not substitute a finite
broadened resolvent for an on-shell limit. See ONSHELL.md for the wave-operator
theorem, its hypotheses, and the distinction from physical QFT scattering.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from ..interacting_field.lattice import real


def _fraction(value, name):
    if type(value) not in (int, F):
        raise ValueError(f'exact rational {name} required')
    return F(value)


@lru_cache(maxsize=128)
def sphere_coefficient(degree):
    """Exact angular coefficient A_n; the proof gives 0 < A_n <= 1."""
    if type(degree) is not int or not 0 <= degree <= 64:
        raise ValueError('integer coefficient degree 0..64 required')
    answer = F(0)
    for i in range(degree+1):
        for j in range(degree-i+1):
            k = degree-i-j
            binomial = math.prod(math.comb(2*a, a) for a in (i, j, k))
            moment = F(math.factorial(degree)*math.prod(math.factorial(2*a) for a in (i, j, k)),
                       math.factorial(2*degree+1)*math.prod(math.factorial(a) for a in (i, j, k)))
            answer += F(binomial, 4**degree)*moment
    return answer


def _atan_interval(inverse, terms):
    x = F(1, inverse)
    partial = sum(((-1)**n*x**(2*n+1)/F(2*n+1) for n in range(terms)), F(0))
    other = partial+(-1)**terms*x**(2*terms+1)/F(2*terms+1)
    return min(partial, other), max(partial, other)


def pi_interval():
    """Machin identity and alternating series, rather than a rounded pi input."""
    a, b = _atan_interval(5, 32), _atan_interval(239, 10)
    return 16*a[0]-4*b[1], 16*a[1]-4*b[0]


def threshold_enclosure(energy, *, speed=F(1, 8), degree=24):
    e, c = _fraction(energy, 'energy'), _fraction(speed, 'speed')
    if c <= 0 or not 0 < abs(e) < 2*c:
        raise ValueError('c > 0 and 0 < |E| < 2c required; threshold and van Hove endpoints excluded')
    sphere_coefficient(degree)  # Validate the explicit arithmetic resource bound.
    r2 = (e/(2*c))**2
    partial = sum((sphere_coefficient(n)*r2**n for n in range(degree+1)), F(0))
    remainder = r2**(degree+1)/(1-r2)
    p_lo, p_hi = pi_interval()
    factor = e*e/(4*c**3)
    return {'density_lower':factor*partial/p_hi**2,
            'density_upper':factor*(partial+remainder)/p_lo**2,
            'angular_partial':partial, 'angular_tail_bound':remainder,
            'squared_radius':r2, 'degree':degree, 'node_count':8}


def angular_density(energy, *, speed=.125, order=48):
    """Independent smooth spherical quadrature, diagnostic floating arithmetic.

    This uses no sphere-coefficient series. The exact rational enclosure is
    supplied by threshold_enclosure, not by a comparison of two quadratures.
    """
    e, c = real(energy, 'energy'), real(speed, 'speed')
    if c <= 0 or not 0 < abs(e) < 2*c:
        raise ValueError('c > 0 and 0 < |E| < 2c required')
    if type(order) is not int or not 4 <= order <= 256:
        raise ValueError('quadrature order 4..256 required')
    z, weights = np.polynomial.legendre.leggauss(order)
    phi = 2*np.pi*(np.arange(2*order)+.5)/(2*order)
    radial = np.sqrt(1-z*z)[:, None]
    coordinates = (radial*np.cos(phi), radial*np.sin(phi),
                   np.broadcast_to(z[:, None], (order, 2*order)))
    jacobian = np.ones((order, 2*order))
    for coordinate in coordinates:
        jacobian /= np.sqrt(1-(e/(2*c)*coordinate)**2)
    angular_mean = float(weights@jacobian.mean(axis=1)/2)
    return e*e/(4*np.pi**2*c**3)*angular_mean


def scattering_eigenvalue(principal_green, contact_density, *, interaction=1.):
    """Rank-one bright-channel S from an admitted real boundary value and rho.

    This algebraic helper does not compute or certify that its inputs are
    the boundary values of this lattice. ONSHELL.md proves their existence
    and positive density on the open threshold interval. Dark ac channels
    have eigenvalue one; point spectrum is outside the scattering domain.
    """
    a = real(principal_green, 'principal-value Green function')
    rho = real(contact_density, 'positive contact density')
    u = real(interaction, 'repulsion', nonnegative=True)
    if rho <= 0:
        raise ValueError('positive absolutely continuous density required')
    real_part, imaginary = 1-u*a, math.pi*u*rho
    if not math.isfinite(real_part) or not math.isfinite(imaginary):
        raise ValueError('boundary products exceed the floating evaluation domain')
    scale = max(abs(real_part), abs(imaginary))
    if scale == 0:
        raise ValueError('boundary denominator is numerically unresolved')
    x, y = real_part/scale, imaginary/scale
    return complex((x-1j*y)/(x+1j*y))


def report():
    rows = []
    for e in (F(1, 64), F(1, 16), F(1, 8), F(3, 16)):
        exact = threshold_enclosure(e)
        numeric = angular_density(float(e))
        rows.append({'energy':str(e),
                     'exact_enclosure':{key:str(value) if isinstance(value, F) else value
                                        for key, value in exact.items()},
                     'density_interval_display':[float(exact['density_lower']), float(exact['density_upper'])],
                     'independent_angular_quadrature':numeric})
    root = Path(__file__).resolve().parents[2]
    sources = ('model/contact_scattering/onshell.py',
               'model/contact_scattering/test_onshell.py',
               'model/contact_scattering/ONSHELL.md',
               'model/contact_scattering/scattering.py',
               'model/interacting_field/lattice.py', 'model/bound_pair/fiber.py')
    return {'schema':'oph.lattice_contact_threshold.v1',
            'source_bindings':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                               for name in sources},
            'total_momentum':[0,0,0], 'speed':'1/8',
            'scope':'INFINITE_FIXED_SPACING_TWO_FERMION_CONTACT_RESPONSE',
            'density_certificate':'exact rational series and Machin pi enclosure',
            'rows':rows, 'eight_nodes_retained':True,
            'wave_operators':'complete on absolutely continuous subspaces by trace-class theorem',
            'interacting_density_threshold':'rho_U(E) / (E^2/(4*pi^2*c^3)) tends to 1 as E tends to 0',
            'dark_zero_energy_point_subspace_retained':True,
            'physical_massless_particle_identified':False,
            'many_body_or_chiral_continuum_scattering_established':False}


def main():
    import argparse
    from pathlib import Path
    from tools import strict_json
    from ..observer.scalar_response.experiment import compare_receipt
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    current = report()
    if args.check is not None:
        compare_receipt(strict_json.loads(args.check.read_text()), current)
        print('Exact lattice contact threshold receipt reproduced')
    else:
        print(json.dumps(current, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
