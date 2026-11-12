"""Four empirical targets select one positive native split-site action."""
from functools import lru_cache
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
import sympy as sy
from model.openwave_calibrated_hadrons.source import components as old_components,weights as old_weights,integral
from model.openwave_native_composites.spectrum import one_body
from model.coherent_defect_wave.physical_attachment.chiral import representation
from itertools import product
import scipy.sparse as sp
from model.openwave_baryon_weak_formfactors.native import formfactors

HERE=Path(__file__).resolve().parent
INPUTS=json.loads((HERE/'inputs.json').read_text())
DECLARED={'Js':F(1,4),'Jw':F(1,4),'t':F(1,8),'eta':F(1),'kappa_colour':F(1)}

def colour_global48(modes,basis):
 subsets,_,g,_=representation();qids=[[subsets.index((c,w))for c in range(3)]for w in(3,4)]
 aids=[[subsets.index(tuple(c for c in range(3)if c!=missing))for missing in range(3)],
       [subsets.index(tuple(c for c in range(5)if c!=missing))for missing in range(3)]]
 index={x:i for i,x in enumerate(modes)};out=sp.csr_matrix((len(basis),len(basis)),dtype=np.int64)
 for axis,gen in enumerate([2*a for a in g[:7]]+[2*(g[6]+2*g[7])]):
  imaginary=any(sy.im(x)!=0 for x in gen);raw=gen/sy.I if imaginary else gen;A=np.zeros((len(modes),len(modes)),dtype=np.int64)
  for species,x,s,f,c in modes:
   ids=(qids if species==0 else aids)[f]
   for cc in range(3):
    value=raw[ids[c],ids[cc]];assert value.is_Integer
    A[index[species,x,s,f,c],index[species,x,s,f,cc]]=int(value)
  J=integral(one_body(A,basis));out+=(4 if axis==7 else 12)*(-1 if imaginary else 1)*(J@J)
 return out

@lru_cache(None)
def components(kind):
 d=old_components(kind);C=colour_global48(d['modes'],d['basis']);assert(C@d['Binteger']).nnz==0
 return {**d,'parts':{**d['parts'],'Cg48':C},'integer_reduced_parts':{**d['integer_reduced_parts'],'Cg48':np.zeros_like(d['integer_reduced_parts']['N'])},'hparts':{**d['hparts'],'Cg48':np.zeros_like(d['hparts']['N'])}}

def weights(p):return {**old_weights(p),'Cg48':p.get('kappa_colour',F(1))/48}

def matrix(kind,p):
 d=components(kind)
 return sum(float(v)*d['hparts'][k]for k,v in weights(p).items())

def fit_at_colour(U,eta=1,Js=.25,Jw=.25,t=.125):
 unit=F(INPUTS['electron_unit_GeV']);z={k:F(v)/unit for k,v in INPUTS['mass_calibration_targets_GeV'].items()}
 delta=z['neutron']-z['proton'];difference=float((z['neutron']+z['proton'])/2-F(3,2)*z['charged_pion'])
 p={'U':U,'eta':eta,'Js':Js,'Jw':Jw,'t':t,'mbar':0,'delta':0,'A':0,'kappa_colour':F(1)}
 def energy(kind,A):
  P=components(kind)['P_calibration'];return float(eigh(P.conj().T@matrix(kind,{**p,'A':A})@P,eigvals_only=True,subset_by_index=(0,0))[0])
 def objective(A):return energy('baryon',A)-1.5*energy('meson',A)-difference
 assert objective(0)<0<objective(20000)
 A=brentq(objective,0.,20000.,xtol=1e-10);mbar=(float(z['charged_pion'])-energy('meson',A))/2
 assert mbar>float(delta)/2
 return {**p,'A':A,'mbar':mbar,'delta':delta}

@lru_cache(None)
def parameters():
 target=float(INPUTS['axial_calibration_target_magnitude'])
 def objective(U):return formfactors(matrix('baryon',fit_at_colour(U)))['gA']-target
 assert objective(224)<0<objective(240)
 colour=brentq(objective,224.,240.,xtol=2e-11)
 p=fit_at_colour(F(round(colour*2**36),2**36))
 return {k:(F(round(v*2**36),2**36)if isinstance(v,(float,np.floating))else F(v))for k,v in p.items()}

@lru_cache(None)
def sector(kind):
 p=parameters();d=components(kind);H=sum(float(v)*d['parts'][k]for k,v in weights(p).items());h=matrix(kind,p)
 assert np.linalg.norm(H@d['B']-d['B']@h)<2e-8
 return {**d,'H':H,'h':h,'parameters':p}

def baryon_states():return sector('baryon')

def exact_matrix(kind):
 p=parameters();d=components(kind);g=d['Gram'];n=len(g);A=sy.zeros(n)
 for name,weight in weights(p).items():A+=sy.Rational(weight.numerator,weight.denominator)*sy.Matrix(d['integer_reduced_parts'][name].tolist())
 return sy.Matrix(n,n,lambda i,j:A[i,j]/sy.sqrt(int(g[i])*int(g[j])))
