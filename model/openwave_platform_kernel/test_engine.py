"""Checks that exercise native complement and the integrator contract."""
from fractions import Fraction as F
import numpy as np
from .export import source_packet


def test_complete_native_symbol_and_untuned_mass_intervals():
    p=source_packet()
    assert p['Nambu_rows_per_site']==180
    assert p['native_CAR_modes_per_site']==90
    assert p['massless_internal_coefficients_retained']==39
    for mass, interval in zip(p['native_mass_coefficients'],p['rational_mass_intervals']):
        assert abs(mass-float((F(interval[0])+F(interval[1]))/2))<1e-16
    assert len(p['native_H0_and_three_kinetic_sparse_matrices'][0])==24


def test_RK4_contractivity_and_operator_error():
    for steps in (10,20,40):
        x=np.linspace(-2/steps,2/steps,1001)
        r=1-1j*x-x*x/2+1j*x**3/6+x**4/24
        assert max(abs(r))<=1+2e-16
        assert max(abs(r**steps-np.exp(-1j*x*steps)))<=2**5/(120*steps**4)+1e-14
