"""Check a finite CMI lower bound from native two-contact transport through a buffer."""

from math import log1p, pi, sinh, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh, expm

from .coefficients import scaled_native_contact


def spectral_function(matrix, function):
    values, vectors = eigh(matrix)
    return (vectors * function(values)) @ vectors.conj().T


def entropy(state):
    values = np.linalg.eigvalsh(state)
    return -float(values @ np.log(values))


def native_path(strength, beta, time):
    h = scaled_native_contact(strength) / 4
    path = np.kron(h, np.eye(2)) + np.kron(np.eye(2), h)
    state = expm(-beta * path)
    state /= np.trace(state)
    ab = state.reshape(4, 2, 4, 2).trace(axis1=1, axis2=3)
    bc = state.reshape(2, 4, 2, 4).trace(axis1=0, axis2=2)
    b = bc.reshape(2, 2, 2, 2).trace(axis1=1, axis2=3)
    reference = np.kron(ab, np.eye(2)/2)
    reduced_reference = np.kron(b, np.eye(2)/2)
    cmi = entropy(ab) + entropy(bc) - entropy(b) - entropy(state)

    # Basis bits are A,B,C. A single load travels from 100 to 001.
    transfer = float(abs(expm(-1j*time*path)[1, 4])**2)
    eigen = np.linalg.eigvalsh(state)
    reference_eigen = np.linalg.eigvalsh(reference)
    m = reference_eigen[0] / eigen[-1]
    big_m = reference_eigen[-1] / eigen[0]
    modular_time = abs(time / beta)
    q = transfer * sqrt(eigen[0]) * pi / (4 * sinh(pi * modular_time))
    lower = q*q / (4*log1p(4*big_m/(m*q)))
    return state, bc, reference, reduced_reference, cmi, transfer, lower


def integrated_deficit(state, bc, reference, reduced_reference):
    """Evaluate the positive resolvent-square identity independently of the entropies."""
    d_bc = len(bc)
    root = spectral_function(state, np.sqrt)
    inverse_root = spectral_function(bc, lambda e: 1 / np.sqrt(e))
    units = np.eye(d_bc**2).reshape(d_bc**2, d_bc, d_bc)
    vmap = np.column_stack([(np.kron(np.eye(2), u @ inverse_root) @ root).ravel()
                            for u in units])
    # Row-major vectorization of X -> reference @ X @ state^{-1}.
    big = np.kron(reference, np.linalg.inv(state).T)
    small = np.kron(reduced_reference, np.linalg.inv(bc).T)
    np.testing.assert_allclose(vmap.T @ vmap, np.eye(d_bc**2), atol=2e-13, rtol=0)
    np.testing.assert_allclose(vmap.T @ big @ vmap, small, atol=2e-13, rtol=0)
    eb, ub = eigh(big)
    es, us = eigh(small)
    v = spectral_function(bc, np.sqrt).ravel()
    full_vector = ub.T @ vmap @ v
    small_vector = us.T @ v
    transport = ub.T @ vmap @ us

    # Independently check the imaginary modular power against the matrix unitaries.
    unitary = spectral_function(reference, lambda e: e**1j) @ spectral_function(state, lambda e: e**(-1j))
    local = spectral_function(reduced_reference, lambda e: e**1j) @ spectral_function(bc, lambda e: e**(-1j))
    left = ub @ (eb**1j * full_vector - transport @ (es**1j * small_vector))
    right = ((unitary - np.kron(np.eye(2), local)) @ root).ravel()
    np.testing.assert_allclose(left, right, atol=2e-13, rtol=0)

    def integrand(s):
        diff = full_vector / (eb + s) - transport @ (small_vector / (es + s))
        return float(((eb + s) * diff) @ diff)

    value, error = quad(integrand, 0, np.inf, epsabs=1e-12, epsrel=1e-10)
    return value, error


if __name__ == "__main__":
    for contrast, beta in ((8, 1/16), (4, 1/16), (8, -1/16)):
        records = native_path(contrast, beta, abs(beta))
        cmi, transfer, lower = records[4:]
        integral, error = integrated_deficit(*records[:4])
        print(f"contrast {contrast}, beta {beta:g}, p {transfer:.9g}, "
              f"CMI lower bound {lower:.9g}, CMI {cmi:.9g}, "
              f"resolvent integral {integral:.9g} ± {error:.1g}")
