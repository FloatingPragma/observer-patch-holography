import itertools

import mpmath as mp
import numpy as np
import pytest

from . import diagnostic


def test_deformed_source_ray_equals_upper_half_plane_resolvent_boundary():
    a = 1/65536
    for omega, p in (([1., 0., 0.], [.7, -.3, .2]),
                     (np.array([1., 2., -3.])/np.sqrt(14), [-.6, .1, .9])):
        k = a*np.asarray(p)
        for s in (.8-.2j, 1.2-.15j):
            original = diagnostic.original_inner_ray(omega, k, a, s)
            deformed = diagnostic.deformed_inner_ray(omega, k, a, s)
            assert deformed == pytest.approx(original, abs=2e-11)
        boundary = diagnostic.deformed_inner_ray(omega, k, a, 1.)
        wrong = diagnostic.deformed_inner_ray(omega, k, a, 1., orientation=-1)
        assert boundary.imag < -a
        assert wrong.imag > a
        assert wrong == pytest.approx(boundary.conjugate(), abs=2e-11)


def test_actual_four_spin_source_resolvent_and_all_eight_edge_charts():
    with mp.workdps(70):
        a = mp.mpf(1)/65536
        c = mp.sqrt(3)/a**2
        k = [a*mp.mpf(x) for x in ('0.7', '-0.4', '0.2')]
        u = [a*mp.mpf(x) for x in ('0.2', '-0.6', '0.5')]
        matrices = [mp.matrix([[0, 1], [1, 0]]), mp.matrix([[0, -1j], [1j, 0]]), mp.matrix([[1, 0], [0, -1]])]
        def tensor(left, right):
            return mp.matrix([[left[i//2, j//2]*right[i % 2, j % 2] for j in range(4)] for i in range(4)])
        singlet = mp.matrix([0, 1/mp.sqrt(2), -1/mp.sqrt(2), 0])
        C, S = [mp.cos(x/2) for x in k], [mp.sin(x/2) for x in k]
        R = mp.sqrt(sum(x*x for x in C))
        z = mp.mpc(1, -.2)
        t = R-a*a*z/(2*mp.sqrt(3))
        energy = 2*c*t
        lam = R/mp.sqrt(3)*a*a*z-a**4*z*z/12
        shape = sum(mp.cos(k[j])*u[j]**2 for j in range(3))
        shape += sum(C[j]*S[j]*u[j]*mp.sqrt(1-u[j]**2) for j in range(3))**2/t**2
        numerator = 1-sum(S[j]**2*u[j]**2 for j in range(3))/t**2
        formula = t/(2*c)*numerator/(shape-lam)
        for signs in itertools.product((-1, 1), repeat=3):
            q = [signs[j]*mp.pi/2+mp.asin(u[j]) for j in range(3)]
            H0 = mp.zeros(4)
            for j, sigma in enumerate(matrices):
                H0 -= c*(mp.sin(k[j]/2+q[j])*tensor(sigma, mp.eye(2))+
                          mp.sin(k[j]/2-q[j])*tensor(mp.eye(2), sigma))
            actual = (singlet.H*mp.lu_solve(energy*mp.eye(4)-H0, singlet))[0]
            assert abs(actual-formula) < mp.mpf('1e-55')
