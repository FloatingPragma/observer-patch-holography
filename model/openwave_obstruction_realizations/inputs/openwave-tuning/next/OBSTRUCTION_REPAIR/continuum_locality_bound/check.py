"""Exact envelope arithmetic and independent matrix-form controls."""
from fractions import Fraction as F
from pathlib import Path
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    eta = F(33, 256)
    A = F(4913, 4096)
    assert (1+eta)**3-1 < F(1, 2)
    # sin(pi/12)^2=(2-sqrt(3))/4 > 1/16, since sqrt(3)<7/4.
    assert eta < F(1, 4) and F(7, 4)**2 > 3
    assert F(223, 256)**6 > F(1, 8)
    assert A < F(3, 2) and A*A < 2
    C3 = 320*16**6
    example = F(72*C3, 100000**3)
    assert example == F(386547056640, 10**15)
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]], complex)
    Z = np.diag([1, -1]).astype(complex)
    alpha = [np.kron(x, X) for x in (X, Y, Z)]
    beta = np.kron(I, Z)
    eye = np.eye(4)

    def algebraic(p):
        q = p@p
        omega = np.sqrt(1+q**3+0j)
        D = sum(x*a for x, a in zip(p, alpha))
        return np.sqrt((omega+1)/(2*omega))*eye+beta@D*q/np.sqrt(2*omega*(omega+1))

    largest_disk_norm, largest_real_error = 0., 0.
    count = 0
    for r in (0., .125, .5, 1., 2., 8., 64.):
        for c in (-1., -.5, 0., .5, 1.):
            p = np.array([r*math.sqrt(1-c*c), 0., r*c], complex)
            theta = math.atan(r**3)/2
            D = sum(x*a for x, a in zip(p, alpha))
            trig = math.cos(theta)*eye+(math.sin(theta)/r if r else 0)*beta@D
            largest_real_error = max(largest_real_error, float(np.linalg.norm(trig-algebraic(p), 2)))
            for k in range(64):
                z = p.copy()
                z[2] += max(1., r)/16*np.exp(2j*math.pi*k/64)
                norm = float(np.linalg.norm(algebraic(z), 2))
                assert norm < 3
                largest_disk_norm = max(largest_disk_norm, norm)
                count += 1
    assert largest_real_error < 1e-13
    result = {'schema': 'oph.continuum_locality.explicit_bound.v1',
        'exact_envelope_inequalities_pass': True,
        'full_three_dimensional_C3_rational_upper': C3,
        'rapidity_limit': 'log(2)', 'example_gap_in_R': 100000,
        'cross_observer_norm_upper_exact': str(example),
        'cross_observer_norm_upper': float(example),
        'matrix_disk_diagnostic_count': count,
        'matrix_disk_diagnostic_max_norm': largest_disk_norm,
        'real_formula_max_error': largest_real_error,
        'numerical_samples_are_not_the_proof': True,
        'physical_gravity_solved': False}
    (HERE/'RESULT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
