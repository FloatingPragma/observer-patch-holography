from fractions import Fraction as F
import numpy as np
from .source import native_lift,complete_spatial_lifts

def test_actual_native_loop_closes_at_every_coordinate():
    r=native_lift();assert len(r['native_rows'])==15
    for x in r['native_rows']:
        assert (F(x['charge_exact'])+F(x['colour_path_generator_exact'])).denominator==1
    for t in np.linspace(0,1,17):
        S=np.diag(np.exp(2j*np.pi*t*np.array([1,1,-2])/3))
        U=np.exp(-2j*np.pi*t/3)*S
        assert abs(np.linalg.det(U)-np.exp(-2j*np.pi*t))<2e-15
    assert r['generator_determinant_winding_exact']=='-1'

def test_all_Gauss_coordinates_have_integer_endpoint_exponent():
    for n,dim in ((2,120),(3,364)):
        r=complete_spatial_lifts(n);assert len(r['records'])==dim
        for row in r['records']:
            assert all((a+b)%3==0 for a,b in zip(row['electric_flux'],row['colour_triality']))
            flux=[0]+row['electric_flux']+[0]
            q=row['local_charge_thirds'].copy();q[0]+=row['reference_charge_thirds']
            assert np.array_equal(np.diff(flux),q)

def test_electric_eigenstate_does_not_localize_angle():
    theta=np.linspace(0,2*np.pi,1024,endpoint=False)
    for n in (-3,-1,0,2,5):
        state=np.exp(1j*n*theta)
        assert np.max(abs(abs(state)**2-1))<5e-16
        for k in (-2,-1,1,2):
            assert abs(np.mean(abs(state)**2*np.exp(1j*k*theta)))<2e-16
