import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import numpy as np
import sympy as s
from scipy.linalg import eigh,solve_sylvester
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.openwave_local_gauss_quantum.gauss import finite_CAR_gauss_control,exact_charge_tree,charged_block
from model.openwave_local_gauss_quantum.__main__ import pins as gauss_pins
from model.openwave_infinite_current.current import norm_controls
from model.openwave_infinite_current.__main__ import pins as current_pins
from model.openwave_spatial_completion.spatial import analytic,shell_source
from model.openwave_spatial_completion.__main__ import pins as spatial_pins
from model.openwave_quantum_rotational.source import matter_exact
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def current_checks():
 m=matter_exact();B0=np.array(m['B'],complex);E=np.array(m['external'][2],complex);M=[np.array(x,complex)for x in m['derivative']];rows=[]
 for point in((0,0,0),(1,-2,3),(9,0,-4),(100,-50,20)):
  B=B0+np.sqrt(2)/8*sum(x*A for x,A in zip(point,M));v,U=eigh(B);K=(U*np.sqrt(256+v*v))@U.conj().T
  first=solve_sylvester(K,K,B@E+E@B);second=solve_sylvester(K,K,2*(E@E-first@first))
  assert np.linalg.norm(first,2)<np.sqrt(5)and np.linalg.norm(second,2)<141/400
  powers=[np.eye(20)]
  for _ in range(12):powers.append(powers[-1]@B)
  polynomial=sum(float(s.binomial(s.Rational(1,2),k))/16**(2*k-1)*sum(powers[j]@E@powers[2*k-1-j]for j in range(2*k))for k in range(1,7))
  coefficient=14*33/(2048*16**13);bound=coefficient*np.sqrt(5)*np.linalg.norm(B,2)**13;error=np.linalg.norm(first-polynomial,2)
  # Independent finite fibre diagnostic; primary interval proofs retain
  # their own arithmetic contract, rather than borrowing this tolerance.
  assert error<bound+2e-12
  rows.append({'field_coordinate':list(point),'derivative_norm':float(np.linalg.norm(first,2)),
   'second_derivative_norm':float(np.linalg.norm(second,2)),'ordered_polynomial_derivative_error':float(error),'global_derivative_remainder_bound':float(bound)})
 # Exact scalar Taylor coefficient and derivative signs.
 x=s.symbols('x',nonnegative=True);f=s.sqrt(256+x);derivative=s.diff(f,x,7)
 assert derivative.is_positive and s.simplify(derivative.subs(x,0)/factorial(7))==s.Rational(33,2048*16**13)
 return {'fibre_diagnostics':rows,'exact_norm_LDL_controls':norm_controls(),
  'proof_review':['uniform fibre Hilbert-Schmidt derivative bound extends to the unbounded field multiplier','second derivative Sylvester denominator is at least32','photon-coordinate moment26 is bounded by13 successive support-raising steps','exact Jz projection selects one infinite-source ground line','complete finite-time error includes all photon and spectator-reference outcomes']}

def spatial_checks():
 a=analytic();x,L=s.symbols('x L',positive=True)
 radial=s.simplify(4*s.pi*s.integrate(x**4*s.exp(-(1+2*L**2/5)*x**2),(x,0,s.oo))/s.pi**s.Rational(3,2))
 assert s.simplify(radial-s.Rational(3,2)/(1+2*L**2/5)**s.Rational(5,2))==0
 energy=16+s.Rational(3,64)/L**2+s.Rational(3,8192)*L**2/(1+s.Rational(2,5)*L**2)**s.Rational(5,2)+s.Rational(3,64)/L**4*(L**2+L**-2-2)
 assert s.limit(energy,L,s.oo)==16
 rows=[shell_source(k,1,1.)for k in(2,3,4)];assert all(r['complete_ground_Casimir_residual']<1e-9 for r in rows)
 assert a['untrapped_threshold_is_normalizable_eigenstate']is False and a['translation_invariant_physical_electron_dispersion_derived']is False
 return {'exact_Gaussian_integral_checked':True,'exact_release_limit':'16','complete_shell_controls':rows,
  'proof_review':['finite Hermite cores and uniform second moments justify square-root form recovery','positive forms give the Mosco liminf; common oscillator coercivity supplies compact sublevels forOmega>0','mass16 extension on the shell complement moves it to high oscillator energy','minmax and equivariant projection establish exactly two confined spin-half ground states','threshold nonattainment uses photon-vacuum projection and Fourier support atp0','no neutral-floor, Gauss-tree or physical translation-invariant dispersion transfer']}

def certificate():
 required={};receipts={}
 for name,fn in(('openwave_local_gauss_quantum',gauss_pins),('openwave_infinite_current',current_pins),('openwave_spatial_completion',spatial_pins)):
  path=ROOT/'model'/name/'results.json';r=json.loads(path.read_text());assert r['source_sha256']==fn();required.update(r['source_sha256']);receipts[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
 assert set(SOURCE_FILES)<=set(gauss_pins())
 for p in import_closure({str(p.relative_to(ROOT))for p in HERE.glob('*.py')}):required[p]=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
 control=finite_CAR_gauss_control();star=charged_block();assert star['star_exact_electric_energy_increment']=='1/16';assert star['star_field_energy_is_intrinsic_isolated_electron_rest_mass']is False
 return {'verdict':'PASS_INDEPENDENT_QUANTUM_DOMAIN_AND_GAUSS_REVIEW','source_sha256':{**required,**receipts},
  'local_Gauss':{'all16_words_and_complete_reference_records':control,'native_charge_domain':exact_charge_tree(),'star_conditional_GeV':star['conditional_GeV_using_unchanged_electron_unit'],
   'reviewed_tree_is_orbital_graph_not_physical_spatial_locality':True,'root_compensator_and_positive_electric_energy_retained':True},
  'infinite_current':current_checks(),'spatial_completion':spatial_checks(),
  'direct_scientific_data_binding_repair':'new Gauss producer directly binds actual action.SOURCE_FILES consumed by native_data; older sealed native-charge source unchanged',
  'physical_electron_Gauss_Coulomb_or_spatial_dispersion_identity_inferred':False}

def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True);g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');a=p.parse_args();r=certificate();path=HERE/'results.json'
 if a.write:path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
 else:
  from model.openwave_larmor_dynamics.certificate import equivalent
  assert equivalent(json.loads(path.read_text()),r)
 print(r['verdict'])
if __name__=='__main__':main()
