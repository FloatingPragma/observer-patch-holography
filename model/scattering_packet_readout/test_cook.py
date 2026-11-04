"""Independent continuous Friedrichs-model control of the two-tail echo."""
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh
from scipy.special import roots_legendre


def test_two_tail_cook_comparison_against_independent_continuous_contact_model():
    # Uniform free spectral measure on [-1,1]; its Green function and source
    # scattering phase are closed-form, independently of the lattice code.
    normalization = np.sqrt(3003/1024)

    def derivative_norms(energy):
        x = energy
        A = normalization*(1-4*x*x)**3/np.sqrt(2)
        Ap = -24*normalization*x*(1-4*x*x)**2/np.sqrt(2)
        App = normalization*(-24+576*x*x-1920*x**4)/np.sqrt(2)
        b, bp, bpp = 2-np.arctanh(x), -1/(1-x*x), -2*x/(1-x*x)**2
        d = np.pi/2
        phase = (b-1j*d)/(b+1j*d)
        first = 2j*d*bp/(b+1j*d)**2
        second = 2j*d*(bpp/(b+1j*d)**2-2*bp*bp/(b+1j*d)**3)
        return abs(App), abs(second*A+2*first*Ap+phase*App)

    constants = [quad(lambda x: derivative_norms(x)[i], -.5, .5,
                      points=[-1/np.sqrt(20), 1/np.sqrt(20)], epsabs=1e-10)[0]
                 for i in range(2)]
    assert constants[1] > constants[0] > 0
    observed = []
    for size in (128, 256):
        energy, weights = roots_legendre(size)
        contact = np.sqrt(weights/2)
        full = np.diag(energy)+.5*np.outer(contact, contact)
        values, vectors = eigh(full)
        packet = np.sqrt(weights)*normalization*np.maximum(0, 1-4*energy**2)**3
        packet /= np.linalg.norm(packet)
        green = .5*np.log((1+energy)/(1-energy))-.5j*np.pi
        denominator = 2-green
        phase = denominator.conjugate()/denominator
        errors = []
        for T in (4, 8, 16, 24):
            free = np.exp(1j*T*energy)
            actual = free*(vectors@(np.exp(-2j*T*values)*(vectors.T@(free*packet))))
            error = np.linalg.norm(actual-phase*packet)
            # Finite quadrature is a diagnostic of the continuous theorem.
            assert error < .5*sum(constants)/T
            errors.append(error)
            if T == 24:
                assert np.linalg.norm(actual-phase.conjugate()*packet) > 1.9
        assert all(a > b for a, b in zip(errors, errors[1:]))
        assert errors[-1] < .004
        observed.append(errors)
    np.testing.assert_allclose(observed[0], observed[1], atol=1e-6, rtol=0)
