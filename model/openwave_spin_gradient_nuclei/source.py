"""One positive native action with five disclosed empirical calibrations."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
import hashlib
import inspect
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import sympy as sy
from scipy.linalg import eigh
from scipy.optimize import brentq
from model.openwave_axial_state.source import components as base_components
from model.openwave_calibrated_hadrons.source import angular4,integral
from model.openwave_native_composites.spectrum import one_body
from model.openwave_baryon_weak_formfactors.native import formfactors
from model.openwave_hadron_spin_constraints.constraints import block
from model.openwave_nuclear_distance import nuclear

HERE=Path(__file__).resolve().parent
INPUTS=json.loads((HERE/'inputs.json').read_text())
DECLARED={'Js':F(1,4),'Jw':F(1,4),'t':F(1,8),'kappa_colour':F(4)}
GRADIENT=F(5,2)

@lru_cache(None)
def components(kind):
    d=base_components(kind);modes=d['modes'];basis=d['basis'];index={v:i for i,v in enumerate(modes)}
    global4=angular4(modes,basis,(0,1),True)
    G4=(2*d['parts']['S4']-global4).tocsr()
    raising=np.zeros((len(modes),len(modes)),dtype=np.int64)
    for species,x,s,f,c in modes:
        if s==0:raising[index[species,x,0,f,c],index[species,x,1,f,c]]=1 if x==0 else -1
    R=integral(one_body(raising,basis));X=R+R.T;Y=R-R.T
    z=np.array([sum((1 if x==0 else -1)*(1-2*s)for i,(_,x,s,_,_)in enumerate(modes)if mask>>i&1)for mask in basis],dtype=np.int64)
    assert (G4-(X@X-Y@Y+sp.diags(z*z,dtype=np.int64))).nnz==0
    assert (G4-G4.T).nnz==0
    parts={**d['parts'],'G4':G4}
    reduced={**d['integer_reduced_parts'],'G4':(d['Binteger'].T@G4@d['Binteger']).toarray()}
    hparts={**d['hparts'],'G4':reduced['G4']/np.sqrt(np.outer(d['Gram'],d['Gram']))}
    assert np.max(abs(G4@d['B']-d['B']@hparts['G4']))<3e-11
    return {**d,'parts':parts,'integer_reduced_parts':reduced,'hparts':hparts,'positive_gradient_identity_exact':True}

def weights(p):
    return {'N':p['mbar'],'D':p['delta']/2,'K':p['t'],'C48':p['U']/48,
            'S4':p['Js']/4,'I4':p['Jw']/4,'B9':p['A']/9,'Ig4':p['eta']/4,
            'Cg48':p['kappa_colour']/48,'G4':p['spin_gradient']/4}

def matrix(kind,p):
    d=components(kind)
    return sum(float(v)*d['hparts'][k]for k,v in weights(p).items())

def energy(kind,p,q=1,twice_spin=1):
    P=block(kind,q,twice_spin)
    return float(eigh(P.conj().T@matrix(kind,p)@P,eigvals_only=True)[0])

def provisional(D,A,eta,gradient):
    return {**DECLARED,'U':3*(D+4*A/9)/8,'A':A,'eta':eta,'spin_gradient':gradient,'mbar':0.,'delta':0.}

@lru_cache(None)
def fitted_split_coefficient(gradient):
    target=float(INPUTS['axial_calibration_target_magnitude'])
    def f(D):return formfactors(matrix('baryon',provisional(D,0.,0.,gradient)))['gA']-target
    assert f(-2)*f(2)<0
    return brentq(f,-2.,2.,xtol=2e-13)

@lru_cache(None)
def floating_parameters(gradient,eta):
    unit=float(INPUTS['electron_unit_GeV']);z={k:float(v)/unit for k,v in INPUTS['mass_calibration_targets_GeV'].items()}
    D=fitted_split_coefficient(gradient);target=(z['proton']+z['neutron'])/2-1.5*z['charged_pion']
    def p(A):return provisional(D,A,eta,gradient)
    def f(A):return energy('baryon',p(A))-1.5*energy('meson',p(A),1,0)-target
    assert f(100)*f(5000)<0
    A=brentq(f,100.,5000.,xtol=1e-10);result=p(A)
    result['mbar']=(z['charged_pion']-energy('meson',result,1,0))/2
    result['delta']=F(INPUTS['mass_calibration_targets_GeV']['neutron'])/F(INPUTS['electron_unit_GeV'])-F(INPUTS['mass_calibration_targets_GeV']['proton'])/F(INPUTS['electron_unit_GeV'])
    return result

def nuclear_matrix(p,S,I,parity):
    c=nuclear.channel(S,I,parity)
    if not c['rank']:return np.zeros((0,0))
    parts={**c['parts'],'G4':2*c['parts']['S4']-4*S*(S+1)*np.eye(c['rank'])}
    return sum(float(v)*parts[k]for k,v in weights(p).items())

def nuclear_energy(p,S=1,I=0,parity=1):
    H=nuclear_matrix(p,S,I,parity)
    return None if not len(H)else float(eigh(H,eigvals_only=True)[0])

@lru_cache(None)
def fit_parameters(gradient=GRADIENT):
    gradient=F(gradient);unit=float(INPUTS['electron_unit_GeV'])
    threshold=sum(float(INPUTS['mass_calibration_targets_GeV'][k])for k in('proton','neutron'))/unit
    binding=float(INPUTS['deuteron_binding_calibration_target_GeV'])/unit
    def f(eta):return threshold-nuclear_energy(floating_parameters(float(gradient),eta))-binding
    assert f(0)*f(10)<0
    eta=brentq(f,0.,10.,xtol=1e-10)
    p=floating_parameters(float(gradient),eta)
    p={k:(v if isinstance(v,F)else F(round(float(v)*2**44),2**44))for k,v in p.items()}
    assert p['mbar']>p['delta']/2>0 and all(p[k]>0 for k in('U','A','eta','Js','Jw','t','kappa_colour'))
    assert p['spin_gradient']>=0
    return p


CALIBRATION_RECIPE_FUNCTIONS=('components','weights','matrix','energy','provisional',
    'fitted_split_coefficient','floating_parameters','nuclear_matrix','nuclear_energy','fit_parameters')

def validate_declared_calibrations(record):
    """Authenticate the exact calibrated actions without rerunning a floating fit."""
    assert record['schema']=='oph.declared_five_input_rational_calibrations.v1'
    assert record['calibration_inputs_sha256']==hashlib.sha256((HERE/'inputs.json').read_bytes()).hexdigest()
    assert record['declared_inputs']==INPUTS
    assert record['binary_quantization_places']==44
    assert record['selected_gradient']==str(GRADIENT)
    assert set(record['recipe_functions_sha256'])==set(CALIBRATION_RECIPE_FUNCTIONS)
    for name in CALIBRATION_RECIPE_FUNCTIONS:
        assert record['recipe_functions_sha256'][name]==hashlib.sha256(inspect.getsource(globals()[name]).encode()).hexdigest()
    menu=INPUTS['architecture_selection']['diagnostic_gradient_menu']
    assert set(record['actions'])==set(menu)
    exact_delta=(F(INPUTS['mass_calibration_targets_GeV']['neutron'])-F(INPUTS['mass_calibration_targets_GeV']['proton']))/F(INPUTS['electron_unit_GeV'])
    result={}
    for gradient,raw in record['actions'].items():
        assert isinstance(raw,dict)
        assert set(raw)==set(DECLARED)|{'U','A','eta','spin_gradient','mbar','delta'}
        assert all(isinstance(v,str) and str(F(v))==v for v in raw.values())
        p={k:F(v)for k,v in raw.items()}
        assert all(p[k]==v for k,v in DECLARED.items())
        assert p['spin_gradient']==F(gradient) and p['delta']==exact_delta
        assert all((p[k]*2**44).denominator==1 for k in('U','A','eta','mbar'))
        assert p['mbar']>p['delta']/2>0 and all(p[k]>0 for k in('U','A','eta'))
        assert p['spin_gradient']>=0
        result[gradient]=p
    return result

@lru_cache(None)
def parameters(gradient=GRADIENT):
    """Return an input-bound exact declared action; fit_parameters is the recipe."""
    record=json.loads((HERE/'calibration_realizations.json').read_text())
    actions=validate_declared_calibrations(record)
    key=str(F(gradient))
    if key not in actions:raise ValueError('Gradient is outside the declared calibration menu; use fit_parameters for a new exploratory fit')
    return actions[key]


@lru_cache(None)
def sector(kind):
    p=parameters();d=components(kind);H=sum(float(v)*d['parts'][k]for k,v in weights(p).items());h=matrix(kind,p)
    assert np.linalg.norm(H@d['B']-d['B']@h)<3e-8
    return {**d,'H':H,'h':h,'parameters':p}

@lru_cache(None)
def exact_matrix(kind):
    d=components(kind);p=parameters();A=sy.zeros(len(d['Gram']))
    for k,v in weights(p).items():A+=sy.Rational(v.numerator,v.denominator)*sy.Matrix(d['integer_reduced_parts'][k].tolist())
    g=d['Gram']
    return sy.Matrix(len(g),len(g),lambda i,j:A[i,j]/sy.sqrt(int(g[i])*int(g[j])))

def baryon_states():
    return sector('baryon')
