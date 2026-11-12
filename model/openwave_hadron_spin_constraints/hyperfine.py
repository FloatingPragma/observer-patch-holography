"""Positive completion of the actual local colour-spin current square."""
from functools import lru_cache
from fractions import Fraction as F
from itertools import product
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from model.openwave_axial_state.source import components,matrix,parameters,INPUTS
from model.openwave_baryon_weak_formfactors.native import formfactors
from .constraints import block,UNIT,TARGETS

@lru_cache(None)
def hyperfine(kind):
 d=components(kind)
 if kind=='baryon':return 4*d['hparts']['B9']-6*d['hparts']['N']+2/3*d['hparts']['S4']
 P=np.diag([int(row[0][0]==row[1][0])for row in d['labels']]);return 4/3*P@d['hparts']['S4']@P

def H(kind,p,h):return matrix(kind,p)+h*hyperfine(kind)
def energy(kind,p,h,q,s):
 P=block(kind,q,s);return float(eigh(P.conj().T@H(kind,p,h)@P,eigvals_only=True)[0])

@lru_cache(None)
def fit(h):
 p0=parameters();targetbar=(float(INPUTS['mass_calibration_targets_GeV']['proton'])+float(INPUTS['mass_calibration_targets_GeV']['neutron']))/(2*UNIT);pion=float(INPUTS['mass_calibration_targets_GeV']['charged_pion'])/UNIT
 def at_U(U):
  p={**p0,'U':U,'A':0.,'mbar':0.,'delta':0.}
  def at_A(A):return {**p,'A':A}
  def obj(A):return energy('baryon',at_A(A),h,1,1)-1.5*energy('meson',at_A(A),h,1,0)-(targetbar-1.5*pion)
  if obj(0)>=0:return None
  A=brentq(obj,0.,10000.,xtol=1e-9);m=(pion-energy('meson',at_A(A),h,1,0))/2
  return {**p,'A':A,'mbar':m,'delta':p0['delta']}
 def obj(U):
  p=at_U(U)
  if p is None:return float('nan')
  return formfactors(H('baryon',p,h))['gA']-1.2753
 values=[]
 for U in(0.,64.,128.,256.,512.,1024.,2048.,4096.):
  z=obj(U);values.append((U,z))
  if len(values)>1 and np.isfinite(z)and np.isfinite(values[-2][1])and z*values[-2][1]<0:break
 brackets=[(a[0],b[0])for a,b in zip(values,values[1:])if np.isfinite(a[1])and np.isfinite(b[1])and a[1]*b[1]<0]
 if not brackets:return {'feasible':False,'scan':values,'coefficient':h}
 U=brentq(obj,*brackets[0],xtol=1e-9);p=at_U(U);return {'feasible':bool(p['mbar']>float(p['delta'])/2 and p['A']>0 and U>=0),'parameters':p,'coefficient':h,'scan':values}

def record(h):
 d=fit(h)
 if not d['feasible']:return d
 p=d['parameters'];vals={name:energy(kind,p,h,q,s)*UNIT for name,kind,q,s in(('proton','baryon',1,1),('neutron','baryon',0,1),('charged_pion','meson',1,0),('neutral_scalar','meson',0,0),('rho_charged','meson',1,2),('Delta_plus','baryon',1,3))};ff=formfactors(H('baryon',p,h));return {**d,'parameters':{k:float(v)for k,v in p.items()},'energies_GeV':vals,'gA':ff['gA'],'gV':ff['gV'],'all_baryon_levels':eigh(H('baryon',p,h),eigvals_only=True).tolist(),'all_meson_levels':eigh(H('meson',p,h),eigvals_only=True).tolist()}

@lru_cache(None)
def native_hyperfine_check(kind):
 from model.openwave_native_composites.spectrum import native,one_body
 import scipy.sparse as sp
 d=components(kind);modes=d['modes'];basis=d['basis'];index={v:i for i,v in enumerate(modes)};qg,ag,*_=native();pauli=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])];T=sp.csr_matrix((len(basis),len(basis)),dtype=complex)
 for x in(0,1):
  n=np.array([sum((b>>i)&1 for i,m in enumerate(modes)if m[1]==x)for b in basis]);local=sp.diags(4*n*n,dtype=complex)
  for a,sigma in product(range(8),pauli):
   A=np.zeros((len(modes),len(modes)),complex)
   for species,cell,s,f,c in modes:
    if cell!=x:continue
    for ss,cc in product(range(2),range(3)):A[index[species,x,s,f,c],index[species,x,ss,f,cc]]=(qg if species==0 else ag)[f][a][c,cc]*sigma[s,ss]
   J=one_body(A,basis);local-=J@J
  T+=local
 error=float(np.max(abs(T@d['B']-d['B']@hyperfine(kind))));assert error<2e-11
 return {'full_fixed_number_dimension':len(basis),'complete_colour_dimension':len(d['Gram']),'full_native_intertwiner_error':error,'positive_completion':'sum_x[4N_x²−sum_(a,k)(dGamma(T_a tensor sigma_k)_x)²]','full_Fock_positivity':'On each local N sector Cauchy gives sum||sum_jA_j psi||²<=N sum_j,a,k||A_j psi||²=4N²||psi||².'}

@lru_cache(None)
def hyperfine_certificate():
 rows=[record(h)for h in(0.,1.,10.,50.,100.,120.,150.)]
 for row in rows:
  if 'scan'in row:row['scan']=[[u,(z if np.isfinite(z)else None)]for u,z in row['scan']]
 return {'native_source_checks':{kind:native_hyperfine_check(kind)for kind in('baryon','meson')},'controls':rows,
  'target_rho_requires_h_at_least':3*((float(TARGETS['rho_vector'])-float(INPUTS['mass_calibration_targets_GeV']['charged_pion']))/UNIT-.5)/32,
  'exact_meson_spin_gap_bound':'0<=m_rho-m_pi<=(32h/3+2Js) in source units; min-max compares the same four orbital states.',
  'calibration_scope':'Each feasible h control refits the same four empirical inputs. The positive-A continuation supplies no calibrated rho solution in the recorded search; failed scans are not a theorem excluding every possible native hyperfine architecture.',
  'unperformed_claims':'No fitted rho or Delta success is inferred from an infeasible source or from a negative A coefficient.'}
