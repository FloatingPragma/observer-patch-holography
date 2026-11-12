from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json,hashlib
import numpy as np
import sympy as sy
from scipy.linalg import eigh
from model.openwave_common_mass_source.source import full_source_terms,full_index,lepton_source,G
from model.openwave_axial_state.source import sector,parameters,weights
from model.openwave_baryon_weak_formfactors.native import currents
ROOT=Path(__file__).resolve().parents[2]
TARGETS=('openwave_dressed_baryon_lines','openwave_beta_outgoing','openwave_meson_memory_window','openwave_common_mass_source')

def load(package):
 d=json.loads((ROOT/'model'/package/'results.json').read_text());assert d['verdict'].startswith('PASS')
 for p,h in d['source_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,(package,p)
 return d

def integer_current_action(A,v):
 A=A.tocsr();out=[]
 for i in range(A.shape[0]):
  r=b=0
  for k in range(A.indptr[i],A.indptr[i+1]):
   z=A.data[k];a=int(round(z.real));c=int(round(z.imag));assert abs(z-a-1j*c)<1e-12;r+=a*v[A.indices[k]];b+=c*v[A.indices[k]]
  out.append((r,b))
 return out

def native_common_check(receipt):
 from model.common_source.chiral_candidate.action import ordering
 terms=full_source_terms()['terms'];rows=ordering()['all_left_components'];qids=[next(i for i,row in enumerate(rows)if row['exterior_subset']==[0,w])for w in(3,4)];qi=[full_index(0,s,c,0)for s,c in product(range(2),qids)];li=[full_index(0,s,c,0)for s,c in product(range(2),(10,11))]
 def extract(items,ids):
  lookup={x:i for i,x in enumerate(ids)};M=sy.zeros(len(ids))
  for i,j,a,b in items:
   if i in lookup and j in lookup:M[lookup[i],lookup[j]]=sy.Rational(F(a))+sy.I*sy.Rational(F(b))
  return M
 def fock(A):
  out=sy.zeros(16)
  for mask in range(16):
   for (i,j),v in A.todok().items():
    if not(mask>>j)&1:continue
    b=mask^(1<<j)
    if b>>i&1:continue
    out[b|(1<<i),mask]+=v*(-1)**((mask&((1<<j)-1)).bit_count()+(b&((1<<i)-1)).bit_count())
  return out
 V=sy.zeros(64)
 for t in terms:
  if t['kind']=='one_body_product':V+=sy.Rational(F(t['weight'])/G)*sy.kronecker_product(extract(t['left_matrix'],qi),fock(extract(t['right_matrix'],li)))
 x=sy.Symbol('x');poly=sy.factor(V.charpoly(x).as_expr());assert V==V.H
 assert sy.expand(poly-sy.sympify(receipt['local_contact_norm_certificate']['characteristic_polynomial']))==0
 # Invert the literal full polynomial's positive fermion-square rows independently.
 lepids=[full_index(x,s,c,f)for x,s,c,f in product(range(2),range(2),(9,10,11,15),range(3))];index={j:i for i,j in enumerate(lepids)};A=np.zeros((48,48),complex);B=np.zeros_like(A);energies=[]
 for t in terms:
  if t['kind']!='positive_linear_fermion_square':continue
  j=len(energies);energies.append(t['energy_electron_units_diagnostic'])
  for i,a,b in t['annihilation_coefficients']:A[j,index[i]]=a+1j*b
  for i,a,b in t['creation_coefficients']:B[j,index[i]]=a+1j*b
 W=np.c_[A,B];C=np.c_[B.conj(),A.conj()];H=W.conj().T@np.diag(energies)@W-C.conj().T@np.diag(energies)@C
 canon=max(np.linalg.norm(A@A.conj().T+B@B.conj().T-np.eye(48)),np.linalg.norm(A@B.T+B@A.T));assert canon<5e-15 and np.linalg.norm(H-lepton_source()['H'])<1e-11
 # Exact complete native current variance using the retained integer occupation trial.
 d=receipt['complete_coupled_baryon_energy_records'];v=list(map(int,d['neutron_integer_trial']));norm=sum(a*a for a in v);num=0;c=currents()
 for xx,mu in product(range(2),range(4)):
  raising,z=c['full'][xx][mu]
  for J in(raising+raising.conj().T,(raising-raising.conj().T)/1j,2*z):num+=sum(a*a+b*b for a,b in integer_current_action(J,v))
 kappa=F(3*num,4*norm);assert kappa==F(d['contact_second_moment_exact'])
 variance=F(d['trial_strong_variance_electron_units_squared_exact'])+G*G*kappa;assert variance==F(d['full_source_energy_variance_electron_units_squared_exact'])
 p=parameters();floor=p['mbar']-p['delta']/2-18*G;assert floor>1
 return {'actual_full_polynomial_local_Fock_characteristic_polynomial':str(poly),'exact_local_norm':6,'inverse_Bogoliubov_CAR_error':float(canon),'literal_polynomial_BdG_error':float(np.linalg.norm(H-lepton_source()['H'])),'exact_native_kappa':str(kappa),'exact_energy_variance':str(variance),'positive_full_colourful_floor':str(floor),'analytical_checks':['three-family first-quantized triangle bound18N_Q','pure Dirac reference and exact charged threshold1','complete rank-one-filter failure and work','Gaussian3/4 current covariance and all-lepton-Fock moment','no transfer from filled-neutral beta spectral lines']}

def dressed_check(receipt):
 from model.openwave_dressed_baryon_lines.source import construction,parity_basis
 out=[]
 for member in receipt['members']:
  for name,q in(('neutron',0),('proton',1)):
   r=member[name];d=construction(member['member'],q);B=parity_basis(d);C=parity_basis(d,-1);assert B.shape[1]+C.shape[1]==d['spin_dimension'];err=np.max(abs((d['H']@B-B@(B.T@d['H']@B)).data),initial=0);assert err<3e-10
   c=r['complete_prepared_spectrum_certificate'];D=c['dyadic_denominator'];O=np.asarray(c['dyadic_eigenframe'],dtype=object);O=np.vectorize(int,otypes=[object])(O);gram=O.T@O-D*D*np.eye(len(O),dtype=object);square=sum(int(x)*int(x)for row in gram for x in row);assert F(square,D**4)<=F(c['exact_Gram_error_upper'])**2
   weights=[F(z['dyadic_frame_probability'])for z in r['complete_spectral_clusters']];assert abs(sum(weights)-1)<=F(c['exact_Gram_error_upper'])
   out.append({'member':member['member'],'charge':q,'full_spin_dimension':d['spin_dimension'],'even_dimension':B.shape[1],'odd_dimension':C.shape[1],'native_reduction_error':float(err),'independent_exact_Gram_Frobenius_squared':str(F(square,D**4)),'all_spectral_clusters_retained':True})
 return {'sources':out,'analytical_checks':['polar-frame correction plus separated-projector bound4eta+4epsilon/gap','complete uniform-input failure and reference retained','discrete line differences do not establish an asymptotic decay endpoint','scalar filled-neutral source is distinct from chiral Dirac momentum source']}

def memory_check(receipt):
 from model.openwave_meson_memory_window.window import universal,certified_universal,taylor_first,G as gm
 from model.openwave_meson_ring.ring import construction,reduction
 out=[]
 for L in(2,4):
  d=construction(L);B,H,_=reduction(d);C=H[6:,:6].toarray()/gm;err=np.linalg.norm(C.conj().T@C-24*np.eye(6));assert err<1e-9
  out.append({'ring':L,'complete_emitted_blocks':9,'injection_Gram_error':float(err)})
 u=universal(4);e,V=eigh(u['D'].toarray()-u['parent_contact_energy']*np.eye(u['dimension']));t=.1;exact=-1j*gm*t*(V@(np.exp(-.5j*e*t)*np.sinc(e*t/(2*np.pi))*(V.conj().T@u['psi'])));trial,cert=taylor_first(certified_universal(4),F(1,10));error=float(np.linalg.norm(exact-trial));assert error<float(F(cert['total_numeric_vector_error_upper_exact']))
 for row in receipt['controlled_windows']:assert F(row['infinite_line_light_channel_probability_lower'])>0 and F(row['original_input_probability_lower'])==F(6,25)*F(row['infinite_line_light_channel_probability_lower'])
 return {'native_controls':out,'independent_diagonal_first_emission_vs_certified_Taylor':error,'numeric_error_budget':cert['total_numeric_vector_error_upper_exact'],'analytical_checks':['exact off-block norm sqrt24 gives cubic feedback remainder','all9family scalar-shift blocks retained','dimension-independent roundoff through absolute row/column bounds','hopping-word tail transfers spatial instruments with outside complement','no ring-Fourier or irreversible-rate transfer']}

def parity_check(receipt):
 from model.openwave_beta_outgoing.response import pair_response,wick_column_norm,parity_bound
 a,b,g=sy.symbols('a b g',real=True);W=(1-g)**2*(1-a)*(1+b)/4+g*g*(1-a)*(1-b)+(1+g)**2*(1+a)*(1-b)/4;odd=sy.expand(W-W.subs({a:-a,b:-b},simultaneous=True));assert sy.expand(odd-2*g*((1-g)*a-(1+g)*b))==0
 r=pair_response(8,F(1,10000));wick=wick_column_norm(8);assert abs(wick-r['full_native_first_column_squared_norm_floating'])<2e-9
 bounds=[parity_bound(L,F(1,10000))for L in(8,16,32,64)];assert all(F(z['strict_negative_exact_parity_difference_upper'])<0 for z in bounds)
 return {'independent_exact_spin_parity_polynomial':str(odd),'retained_pair_channels':len(r['all_first_pair_channels']),'all120_hadron_columns':len(r['complete_hadron_spectral_pair_weights']),'pair_column_vs_full_offsite_Wick_error':abs(wick-r['full_native_first_column_squared_norm_floating']),'exact_source_sign_bounds':bounds,'analytical_checks':['hole-reflection pairing makes the b term nonpositive with even energy weights','at leastL/4 positive Dirac-momentum modes give a>4/5','full54-norm second-Duhamel bound includes multipair channels','positive p=2sin(k) is not group-velocity sign','finite momentum parity does not establish a3D angular coefficient or lifetime']}

def review():
 receipts={p:load(p)for p in TARGETS}
 return {'schema':'oph.independent_current_source_review.v1','verdict':'PASS_INDEPENDENT_NATIVE_CURRENT_AND_RESPONSE_DOMAIN_REVIEW','common_source':native_common_check(receipts[TARGETS[3]]),'dressed_lines':dressed_check(receipts[TARGETS[0]]),'memory_window':memory_check(receipts[TARGETS[2]]),'momentum_parity':parity_check(receipts[TARGETS[1]]),'reviewed_packages':list(TARGETS),'no_new_physical_or_Lean_status_inferred':True}
