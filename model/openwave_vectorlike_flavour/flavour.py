"""Native hypercharge mixing and complete vectorlike free mass branches."""
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations_with_replacement,combinations
from decimal import Decimal,localcontext
from pathlib import Path
import json
import mpmath as mp
import sympy as s
import numpy as np
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.lepton_mass_compiler.compiler import electron_calibration
from model.lepton_mass_compiler.analyzer import rational_effect_intervals

M=1;G=1

@lru_cache(None)
def ward_certificate():
 subsets,_,currents,vertices=representation();colour=sum((t*t for t in currents[:6]),s.zeros(15))+currents[6]**2+(currents[6]+2*currents[7])**2/3
 weak=sum((t*t for t in currents[8:11]),s.zeros(15));ids=[i for i in range(15) if colour[i,i]!=0 and weak[i,i]==0]
 assert len(ids)==6
 R=[t.extract(ids,ids) for t in currents];Y=R[-1]
 ext=[s.diag(t,-r.T,r) for t,r in zip(currents,R)]
 mixing=s.zeros(27);heavy=s.zeros(27)
 for k,i in enumerate(ids):
  mixing[i,15+k]=mixing[15+k,i]=Y[k,k]
  heavy[15+k,21+k]=heavy[21+k,15+k]=1
 for t in ext:
  assert t.T*mixing+mixing*t==s.zeros(27)
  assert t.T*heavy+heavy*t==s.zeros(27)
 checks=0
 for i,j,k in combinations_with_replacement(range(12),3):
  anti=[-r.T for r in R]
  z=s.trace(R[i]*(R[j]*R[k]+R[k]*R[j]))+s.trace(anti[i]*(anti[j]*anti[k]+anti[k]*anti[j]))
  assert s.simplify(z)==0;checks+=1
 Q=ext[10]+ext[11]
 assert Q.T*mixing+mixing*Q==s.zeros(27) and Q.T*heavy+heavy*Q==s.zeros(27)
 return {'native_colourful_weak_singlet_indices':ids,'native_right_hypercharges':[str(Y[k,k]) for k in range(6)],
  'all_current_mass_Ward_checks':24,'vectorlike_symmetric_cubic_anomaly_cancellations':checks,
  'added_all_left_internal_coordinates_per_family':12,'full_internal_coordinates_per_family':27,
  'three_family_two_spin_CAR_per_site':162,'charged_Dirac_channels_per_family':13,'neutral_internal_modes_per_family':1,
  'complete_heavy_symmetric_tensor_rank':heavy.rank(),'mixing_tensor_rank':mixing.rank(),
  'architecture':'A full conjugate/same-representation pair of native colourful weak singlets is added; hypercharge mixing G and vectorlike mass M are both1. No experimental quark masses enter.'}

def branches(m,y,g=G,heavy=M):
 with mp.workdps(65):
  m,y,g,heavy=map(mp.mpf,(m,y,g,heavy));a=m*m+g*g*y*y+heavy*heavy;disc=mp.sqrt(a*a-4*m*m*heavy*heavy)
  low=mp.sqrt(2*m*m*heavy*heavy/(a+disc));high=mp.sqrt((a+disc)/2)
  return low,high

def ordering_proof():
 a,z=s.symbols('a z',positive=True);root=s.sqrt(a*a-4*z)
 # In the domain a>2sqrt(z), the small singular value squared decreases with a.
 derivative=s.diff(2*z/(a+root),a)
 assert s.simplify(derivative+(a/root-1)/2)==0
 return {'mass_matrix':'[[m,0],[G*Y_right,M]]','A':'m^2+G^2*Y_right^2+M^2',
  'low_squared':'2*m^2*M^2/(A+sqrt(A^2-4*m^2*M^2))','high_squared':'(A+sqrt(A^2-4*m^2*M^2))/2',
  'determinant_product':'m_light*m_heavy=m*M','strict_ordering':'0<m_up_light<m_down_light<m for G>0,M>0,m>0',
  'bare_family_limit_ratio':'sqrt((M^2+G^2/9)/(M^2+4*G^2/9))',
  'physical_fit_boundary':'Unequal branches follow from a declared vectorlike action, not a source-forced Standard Model Yukawa or a measured quark-mass fit.'}

def complete_channel(m,y):
 mass=np.array([[m,0.],[float(y),1.]])
 # Independent full chiral/spin symbol, retaining both singular branches.
 p=np.array([.13,-.07,.11]);sig=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])];sp=sum(x*t for x,t in zip(p,sig))
 h=np.block([[-np.kron(np.eye(2),sp),np.kron(mass,np.eye(2))],[np.kron(mass.T,np.eye(2)),np.kron(np.eye(2),sp)]])
 w=np.linalg.eigvalsh(h);sv=np.linalg.svd(mass,compute_uv=False);expected=np.sort(np.repeat(np.r_[-np.sqrt(sv*sv+sum(p*p)),np.sqrt(sv*sv+sum(p*p))],2))
 assert np.max(abs(w-expected))<1e-13
 # Original light-left preparation is measured in every full rest spectral output.
 rest=np.block([[np.zeros((2,2)),mass],[mass.T,np.zeros((2,2))]])
 e,v=np.linalg.eigh(rest);probs=abs(v[0])**2
 assert abs(sum(probs)-1)<1e-13
 return {'full_spin_chiral_dimension':8,'momentum':p.tolist(),'dispersion_max_error':float(np.max(abs(w-expected))),
  'rest_original_leg_spectral_energies':[round(float(x),14) for x in e],
  'rest_complete_original_leg_probabilities':[round(float(x),14) for x in probs],
  'spin_reference_and_other_native_occupations':'retained as exact free complements',
  'full_CAR_addition_records':addition_records(mass)}

def addition_records(mass):
 rest=np.block([[np.zeros((2,2)),mass],[mass.T,np.zeros((2,2))]])
 h=np.kron(rest,np.eye(2));e,v=np.linalg.eigh(h);negative=v[:,e<0];positive=v[:,e>0]
 seed=np.eye(8)[:,0];added=positive@(positive.T@seed);success=float(np.vdot(added,added).real)
 assert abs(success-.5)<1e-12
 time=1.;added=v@(np.exp(-1j*time*e)*(v.T@added))/np.sqrt(success)
 W=np.column_stack((negative,added));assert np.linalg.norm(W.conj().T@W-np.eye(5))<1e-12
 probabilities=np.zeros(256)
 for ids in combinations(range(8),5):
  dirac_mask=sum(1<<i for i in ids);native_mask=dirac_mask&15
  for r in range(2):
   for ss in range(2):
    if not dirac_mask>>(4+2*r+ss)&1:native_mask|=1<<(4+2*r+1-ss)
  probabilities[native_mask]=abs(np.linalg.det(W[list(ids),:]))**2
 assert abs(sum(probabilities)-1)<1e-12
 for mask,p in enumerate(probabilities):
  if p>1e-18:assert (mask&15).bit_count()-(mask>>4).bit_count()==1
 return {'full_active_Fock_dimension':256,'complete_unconditional_native_occupation_probabilities':[round(float(success*p),14) for p in probabilities],
  'creation_failure_probability':.5,'exact_creation_success_probability':'1/2','time':time,
  'reference':'Filled negative free Dirac bands, equivalently the zero-charge native pairing ground; idle apparatus carries opposite added charge',
  'native_charge_on_success':'q, on every nonzero original native occupation outcome'}

def floor_certificate():
 boxes=rational_effect_intervals()['effect_intervals_k0_k1_k2']
 boxes=[tuple(map(F,b)) for b in boxes];ehi=boxes[1][1];out=[]
 for k in (1,2,0):
  for q in (F(2,3),F(1,3)):
   a=ehi*ehi*q*q;c=q*q;mlo,mhi=boxes[k]
   diagonal=mlo*mlo+c-a
   determinant=mlo*mlo*(1-a)-a*(1+c-a)
   assert 0<a<1 and diagonal>0 and determinant>0
   assert 2*a<mlo*mlo+1+c
   out.append({'character':k,'charge_magnitude':str(q),'first_principal_minor_lower':str(diagonal),'determinant_lower':str(determinant)})
 return {'method':'Exact rational Sylvester bounds on M^T M-(e_light |q|)^2 I; all charged branches retained','rows':out,'analyzer_effect_intervals':[[str(x) for x in b] for b in boxes]}

def certificate():
 data=electron_calibration();rows=data['rows'];out=[]
 with mp.workdps(60):
  unit=mp.mpf(data['source_energy_unit_GeV']);light=mp.mpf(rows[0]['compiled_source_rest_energy'])
  for family,row in enumerate(rows,1):
   m=mp.mpf(row['compiled_source_rest_energy'])
   for channel,y,q in [('up','-2/3',mp.mpf(2)/3),('down','1/3',-mp.mpf(1)/3)]:
    yy=mp.mpf(-2)/3 if channel=='up' else mp.mpf(1)/3
    lo,hi=branches(m,yy)
    assert 0<lo<m<hi and lo>=light*abs(q)
    out.append({'family':family,'channel':channel,'native_right_hypercharge':y,'native_electric_charge':mp.nstr(q,20),
     'bare_analyzer_mass':mp.nstr(m,35),'light_source_gap':mp.nstr(lo,35),'heavy_source_gap':mp.nstr(hi,35),
     'light_gap_in_electron_units':mp.nstr(lo/light,35),'heavy_gap_in_electron_units':mp.nstr(hi/light,35),
     'conditional_light_GeV':mp.nstr(lo*unit,35),'conditional_heavy_GeV':mp.nstr(hi*unit,35),
     'complete_channel_records':complete_channel(float(m),float(yy))})
  dn=mp.mpf(out[1]['light_gap_in_electron_units'])-mp.mpf(out[0]['light_gap_in_electron_units'])
  return {'verdict':'PASS_NATIVE_VECTORLIKE_UP_DOWN_SPLITTING','native_current_certificate':ward_certificate(),'exact_charged_floor_certificate':floor_certificate(),'exact_mass_ordering':ordering_proof(),
   'source_parameters':{'heavy_mass':M,'hypercharge_mixing':G,'neutral_Higgs_background_norm':1},'calibration':data,'quark_branches':out,
   'light_down_minus_up_electron_units':mp.nstr(dn,35),'light_down_minus_up_conditional_GeV':mp.nstr(dn*mp.mpf(rows[0]['electron_anchored_GeV']),35),
   'extra_branches_retained':6,'free_charged_Fock_floor':'H-Eground >= m_e*abs(Q) at the stated parameters; all charged branches checked, neutral zeros retained',
   'free_charge_floor_method':'For each 2x2 block M^T M-m_e^2*q^2 I, positivity follows from its positive first diagonal and determinant; rational analyzer intervals certify both.',
   'scope':'The added vectorlike module and hypercharge coupling are declared architecture. All heavy states, neutral spectators and lattice valleys remain. No interacting quark pole, MSbar mass, physical heavy-particle identity, measured Yukawa fit or uniqueness theorem is inferred.'}
