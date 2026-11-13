"""Exact compact Bloch operator and full native complement controls."""
from fractions import Fraction as F
from functools import lru_cache
import sys
from pathlib import Path
import numpy as np
import sympy as sy
from scipy.linalg import eigh
from scipy.optimize import brentq
from model.openwave_hadron_refinement_flow.positivity import inputs,declared_kappa
from .bloch import symbol_data,symbol,coupling_square


def constants(dU=F(0),dm=F(0)):
    p,_,unit=inputs();kappa=declared_kappa();p=dict(p);p['U']+=dU;p['mbar']+=dm
    common=3*p['mbar']+p['delta']/2+3*p['t']
    a=common+p['A']+F(3,4)*(p['Js']+p['Jw']+p['spin_gradient']+p['eta'])
    low=common+5*p['A']/9+8*p['U']/3+F(3,4)*(p['Js']+p['Jw']+p['spin_gradient']+p['eta'])
    high=low+2*(p['Js']+p['Jw'])+4*p['spin_gradient']
    return a,low,high,p,kappa,unit


def exact_rest(dU=F(0),dm=F(0)):
    a,b,c,p,k,u=constants(dU,dm);t=sy.Rational(p['t']);root=sy.sqrt(3)*t
    # Even inversion block, with a separate complete odd block.
    H=sy.Matrix([[a,-root,-root],[-root,c-t/2,-3*t/2],[-root,-3*t/2,b-t/2]])
    odd=sy.Matrix([[c+t/2,3*t/2],[3*t/2,b+t/2]])
    J=sy.Matrix([[sy.Rational(5,3),0,0],[0,1,sy.Rational(2,3)],[0,sy.Rational(2,3),1]])
    return H,odd,J


def certify_formula():
    d=symbol_data();a,b,c,p,k,u=constants();t=float(p['t'])
    expected=np.array([[float(c),0,-t/2,-3*t/2],[0,float(b),-3*t/2,-t/2],[-t/2,-3*t/2,float(c),0],[-3*t/2,-t/2,0,float(b)]])
    assert abs(d['a']-float(a))<2e-11
    assert np.max(abs(d['b']-expected))<2e-11
    assert np.max(abs(d['vl']-np.array([-np.sqrt(1.5)*t]*2+[0,0])))<2e-13
    T0=3*t*t*np.eye(4);T1=np.zeros((4,4));T1[2:,:2]=t*t/2*np.array([[1,3],[3,1]])
    assert np.max(abs(d['T0']-T0))<1e-14
    assert np.max(abs(d['T1']-T1))<1e-14
    # Singular values of the off-diagonal2x2 matrix are2t² andt²;
    # phase multiplication is unitary for every real k.
    assert sy.Matrix([[1,3],[3,1]]).eigenvals()=={4:1,-2:1}
    return {'full_compact_native_formula_verified':True,'exact_coupling_norm_squared':str(5*p['t']**2),
            'complete_coupling_square_spectrum':[str(j*p['t']**2)for j in(1,2,4,5)],
            'argument':'Finite support and translation covariance extend the checked complete four-cell CAR identity to every Bloch fibre.'}


def full_bounds(dU=F(0),dm=F(0)):
    a,b,c,p,kappa,unit=constants(dU,dm)
    common=3*p['mbar']+p['delta']/2-3*p['t']
    split=common+5*p['A']/9+8*p['U']/3+4*kappa/3+F(3,4)*(p['Js']-p['spin_gradient'])+F(11,4)*(p['Jw']-p['eta'])
    three=common+p['A']/3+4*p['U']+2*kappa/3+F(9,4)*(p['Js']+p['Jw']-p['spin_gradient']-p['eta'])
    C=min(split,three)
    # Gershgorin bound on every5x5 compact symbol, using sqrt(3/2)<5/4.
    top=max(a+5*p['t'],c+F(9,2)*p['t'],b+F(9,2)*p['t'])
    assert C>top
    b2=5*p['t']**2;epsilon=b2/(C-top)
    H,O,J=exact_rest(dU,dm)
    ev,V=eigh(np.array(H).astype(float));odd=eigh(np.array(O).astype(float),eigvals_only=True)
    v=V[:,0];e0=float(ev[0]);gap=min(ev[1],odd[0])-e0
    assert gap>float(epsilon)
    angle=float(epsilon)/(gap-float(epsilon))
    ga=float(v@np.array(J).astype(float)@v)
    variance=float(np.linalg.norm((np.array(J).astype(float)-ga*np.eye(3))@v))
    leakage=float(b2/(C-F(str(e0)))**2)
    # P J Q=0 because this actual current preserves the occupation span.
    gaerr=2*variance*angle+6*angle*angle+6*leakage
    return {'Q_floor':C,'all_P_upper':top,'uniform_energy_error':epsilon,
            'ground_P_energy':e0,'GeV_interval':[float(unit)*(e0-float(epsilon)),float(unit)*e0],
            'gA_P':ga,'gA_interval':[ga-gaerr,ga+gaerr],
            'P_eigenvector_sine_bound':angle,'full_Q_probability_upper':leakage,
            'source_parameters':p,'scope':'Infinite one-dimensional native three-quark colour-singlet Bloch fibre, k=0, S=I=1/2; no3D Lorentz or scattering identification.'}


def flow():
    target=float(F('1.2753'))
    def residual(U):
        H,O,J=exact_rest(F(str(U)))
        e,V=eigh(np.array(H).astype(float));v=V[:,0]
        return float(v@np.array(J).astype(float)@v)-target
    assert residual(-1)*residual(1)<0
    dU=F(round(brentq(residual,-1,1,xtol=1e-13)*2**40),2**40)
    p,kappa,unit=inputs()
    raw=full_bounds(dU)
    targetE=F('.93956542194')/unit
    # Calibrate the centre of the proved full-source mass enclosure.
    dm=F(round((float(targetE)-raw['ground_P_energy']+float(raw['uniform_energy_error'])/2)/3*2**40),2**40)
    result=full_bounds(dU,dm)
    return {'degree_counterterm_U':str(dU),'degree_counterterm_mass':str(dm),
            'calibration_inputs':['inherited neutron mass','inherited axial magnitude'],
            'two_cell_operator_change_exact':'0',
            'bulk_result':result,'boundary':'Other bulk hadrons are not inferred from preserving the two-cell calibration operator.'}


if __name__=='__main__':
    import json
    print(json.dumps({'formula':certify_formula(),'fixed_coefficients':full_bounds(),'calibration_flow':flow()},default=str),flush=True)
