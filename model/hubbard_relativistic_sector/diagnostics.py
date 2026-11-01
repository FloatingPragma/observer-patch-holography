"""Independent floating checks; excluded from the exact dispersion bounds."""
import math

import numpy as np
from scipy.integrate import quad
from scipy.special import kv, jv


def fixed_source():
    t, U, a = 1/16, 1., 1.
    A = 2*math.pi*t/U
    odd = np.arange(1, 1002, 2)
    rows = []
    for theta in (0., .1, .2, .3):
        energy = float((8*t/math.pi)*sum(kv(1, odd*A)*np.cosh(odd*theta)/odd))
        vp = float((8*t/math.pi)*sum(kv(0, odd*A)*np.sinh(odd*theta)/odd))
        k = math.asin(theta/A); u = U/(4*t)
        # Independent original dressed-dispersion integrals in Melzer's
        # normalization. u=4 at this source makes the tails benign.
        def efun(x):
            if x == 0: return .5
            return jv(1, x)*math.cos(x*math.sin(k))/x*2*math.exp(-2*u*x)/(1+math.exp(-2*u*x))
        def pfun(x):
            if x == 0: return math.sin(k)
            return jv(0, x)*math.sin(x*math.sin(k))/x*2*math.exp(-2*u*x)/(1+math.exp(-2*u*x))
        alternative_e = 2*t*(u-math.cos(k)+quad(efun, 0, math.inf, epsabs=1e-13)[0])
        alternative_vp = 2*t*(k-quad(pfun, 0, math.inf, epsabs=1e-13)[0])
        residual = max(abs(energy-alternative_e), abs(vp-alternative_vp))
        if residual > 2e-12: raise ArithmeticError('independent dressed-dispersion diagnostic mismatch')
        rows.append({'theta': theta, 'E': energy, 'relative_p': vp/(2*t*a), 'v_p': vp,
                     'E_squared_minus_v_squared_p_squared': energy**2-vp**2,
                     'independent_integral_residual': residual})
    return {'t': t, 'U': U, 'a': a, 'v': 2*t*a, 'A': A, 'rows': rows,
            'floating_diagnostic_only': True, 'diagnostic_admission': 2e-12,
            'conclusion': 'the fixed-source lattice curve is not identified with a constant-v relativistic mass shell'}
