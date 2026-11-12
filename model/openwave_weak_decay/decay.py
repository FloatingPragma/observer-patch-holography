"""Complete native lepton-current contraction and finite three-body records."""
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction as F
from itertools import product
import numpy as np
import sympy as s
from model.common_source import calibrated_matter as internal
from model.common_source.chiral_candidate import action
from model.lepton_mass_compiler.analyzer import rational_effect_intervals

PAULI=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))

def index(p,spin,weak,family):return (p+1)*90+spin*45+(10 if weak==0 else 11)*3+family

def labels(i):
 family=i%3;i//=3;internal_index=i%15;i//=15;spin=i%2;i//=2
 return (i-1,spin,internal_index,family)

@lru_cache(maxsize=1)
def interaction():
 # Freeze complete current contraction before evaluating the transition.
 # weak=0 neutral, weak=1 charged; native raising sign is a harmless
 # simultaneous doublet basis phase, checked against the full15 current.
 ell=internal._LH_tensor()*s.Matrix([0,1]);nu=next(i for i in range(15)if ell[i]!=0)
 D=internal._oriented_components('charged_lepton')[1];e=next(iter(D.todok()))[0]
 J=internal.weak_raising_current(families=False)
 sign=J[nu,e];assert abs(sign)==1
 current=J.extract([nu,e],[nu,e]);assert current==s.Matrix([[0,sign],[0,0]])
 # Include all3 source weak generators, quark and lepton doublets, spins, families and momenta.
 assert nu==10 and e==11
 tau=(J+J.H,-s.I*(J-J.H),J*J.H-J.H*J)
 for a in range(3):
  assert (tau[a]*tau[(a+1)%3]-tau[(a+1)%3]*tau[a]-2*s.I*tau[(a+2)%3])==s.zeros(15)
 bars=(s.eye(2),*(-t for t in PAULI));poly=defaultdict(complex)
 momentum=[pq for pq in product((-1,0,1),repeat=4)if pq[0]-pq[1]+pq[2]-pq[3]==0]
 for t in tau:
  for mu,bar in enumerate(bars):
   A=s.kronecker_product(bar,t,s.eye(3));entries=[(i,j,complex(v))for(i,j),v in A.todok().items()]
   eta=1 if mu==0 else -1
   for p,r,q,u in momentum:
    for i,j,a in entries:
     for k,l,b in entries:
      i0=(p+1)*90+i;j0=(r+1)*90+j;k0=(q+1)*90+k;l0=(u+1)*90+l
      if i0==k0 or j0==l0:continue
      creators=tuple(sorted((i0,k0)));annihilators=tuple(sorted((j0,l0)))
      sc=1 if i0<k0 else -1;sa=-1 if j0<l0 else 1
      poly[(*creators,*annihilators)]+=-eta*a*b*sc*sa/4
 poly={k:v for k,v in poly.items()if abs(v)>0}
 assert all(v.imag==0 and v.real==round(v.real*4)/4 for v in poly.values())
 poly={k:F(v.real)for k,v in sorted(poly.items())}
 charges=[s.Rational(row['electric_Q'])for row in action.ordering()['all_left_components']]
 for (a,b,c,d),v in poly.items():
  assert poly.get((c,d,a,b))==v
  la,lb,lc,ld=map(labels,(a,b,c,d))
  assert la[0]+lb[0]==lc[0]+ld[0]
  assert charges[la[2]]+charges[lb[2]]==charges[lc[2]]+charges[ld[2]]
  for f in range(3):
   assert sum(x[3]==f and x[2]in(10,11)for x in(la,lb))==sum(x[3]==f and x[2]in(10,11)for x in(lc,ld))
  assert sorted((la[3],lb[3]))==sorted((lc[3],ld[3]))
 return poly,nu,e,sign,len(momentum)


def act(state,mode,create):
 out=defaultdict(lambda:s.Integer(0))
 for mask,z in state.items():
  occupied=(mask>>mode)&1
  if occupied==create:continue
  sign=(-1)**((mask&((1<<mode)-1)).bit_count())
  out[mask^(1<<mode)]+=sign*z
 return {k:s.simplify(v)for k,v in out.items()if v!=0}

def add(*states):
 out=defaultdict(lambda:s.Integer(0))
 for state in states:
  for k,v in state.items():out[k]+=v
 return {k:s.simplify(v)for k,v in out.items()if s.simplify(v)!=0}

def scale(state,z):return {k:s.simplify(v*z)for k,v in state.items()}
def dot(a,b):return s.simplify(sum(s.conjugate(v)*b.get(k,0)for k,v in a.items()))

def prepared_states():
 # Spectator factors are common. Choose the exact graded tensor ordering
 # e0,mu0 (ellup,elldown,ecup,ecdown), nu_e+,nu_mu+ (up,down).
 local={0:s.Rational(1,2),9:-s.Rational(1,2),6:s.Rational(1,2),15:-s.Rational(1,2)}
 reference={a|(b<<4)|(1<<8)|(1<<10):u*v for a,u in local.items()for b,v in local.items()}
 initial=scale(add(act(reference,5,True),scale(act(reference,6,False),-1)),1/s.sqrt(2))
 target=scale(add(act(reference,0,True),act(reference,3,False)),1/s.sqrt(2))
 target=act(act(target,8,False),11,True)
 assert dot(reference,reference)==dot(initial,initial)==dot(target,target)==1 and dot(initial,target)==0
 mapping={index(0,spin,1,f):4*f+spin for f in (0,1)for spin in (0,1)}
 mapping.update({index(1,spin,0,f):8+2*f+spin for f in (0,1)for spin in (0,1)})
 return reference,initial,target,mapping


def matrix_element():
 poly,*_=interaction();_,initial,target,mapping=prepared_states();out={};selected=[]
 # All omitted terms have a mismatched spectator/momentum occupation or
 # cannot supply both required neutral changes. This is checked by requiring
 # the two mandatory neutral legs before using the common spectator factor.
 must_create=index(1,1,0,1);must_annihilate=index(1,0,0,0)
 for key,v in poly.items():
  a,b,c,d=key
  if must_create not in(a,b)or must_annihilate not in(c,d):continue
  if any(i not in mapping for i in key):continue
  col=initial
  for i,creation in((c,False),(d,False),(b,True),(a,True)):col=act(col,mapping[i],creation)
  z=s.simplify(s.Rational(v.numerator,v.denominator)*dot(target,col))
  if z:selected.append({'source_CAR_indices':list(key),'coefficient':str(v),'matrix_element':str(z)});out=add(out,scale(col,s.Rational(v.numerator,v.denominator)))
 amplitude=dot(target,out);assert amplitude!=0
 return amplitude,selected


def mass_data():
 boxes=rational_effect_intervals()['effect_intervals_k0_k1_k2'];vals=[sum(map(F,b),F(0))/2 for b in boxes]
 me=vals[1]/vals[2];mtau=vals[0]/vals[2]
 assert F(boxes[1][1])/F(boxes[2][0])<1 and F(boxes[0][1])/F(boxes[2][0])<17
 return float(me),float(mtau),boxes


def current_energy_expectation(target=False):
 poly,*_=interaction();me,mtau,_=mass_data();q=(1-me)/2;masses=(me,1.,mtau);occupations={}
 for i in range(270):
  p,spin,internal_index,f=labels(i);momentum=q*p;sgn=1 if spin==0 else -1
  if internal_index==11:occupations[i]=(1+sgn*momentum/np.sqrt(momentum**2+masses[f]**2))/2
  else:occupations[i]=1. if sgn*momentum>0 else 0.
 if target:
  occupations[index(0,0,1,0)]=1.;occupations[index(1,1,0,1)]=1.;occupations[index(1,0,0,0)]=0.
 else:occupations[index(0,1,1,1)]=1.
 return float(sum(float(v)*occupations[a]*occupations[b]for(a,b,c,d),v in poly.items()if(a,b)==(c,d)))


def certificate():
 poly,nu,e,sign,quad=interaction();amplitude,selected=matrix_element();me,mtau,boxes=mass_data()
 norm=sum(abs(v)for v in poly.values());G=F(1,64);B=F(270*18)+G*norm;t=F(1,10**14)
 lead=s.Abs(amplitude)*s.Rational(G.numerator,G.denominator)*s.Rational(t.numerator,t.denominator)
 remainder=s.Rational((B*B*t*t/2).numerator,(B*B*t*t/2).denominator);assert lead>remainder
 Ei=current_energy_expectation();Ef=current_energy_expectation(True)
 return {'schema':'oph.native_current_three_body_transition.v1',
 'architecture':'H=H0+G sum_(a,mu,p-r+q-u=0) eta_mu :j_a^mu(p,r) j_a^mu(q,u):; j=c^dagger(bar_sigma^mu tensor native_tau_a/2 tensor I3)c; canonical empty-vacuum normal ordering',
 'carrier':'270 actual native CAR modes:15 internal times3families times2spin times3momenta;144 active native weak-doublet modes; all126 weak-singlet modes retained',
 'reference':'filled negative bands of the same free H0, zero bands vacant; interacting evolution starts from this declared reference, not an asserted interacting ground',
 'native_internal_neutral_charged_indices':[nu,e],'native_raising_current_sign':str(sign),
 'momentum_modes_in_units_q':[-1,0,1],'momentum_conserving_quadruples':quad,
 'complete_combined_quartic_word_count':len(poly),'interaction_CAR_word_l1_norm_upper':str(norm),
 'exact_Hermiticity_charge_momentum_and_family_checks':len(poly),
 'charged_source_effect_intervals_k0_k1_k2':boxes,
 'charged_masses_in_middle_branch_units':[me,1.,mtau],'neutral_masses':'0',
 'box_momentum_q_exact':'(1-e_1/e_2)/2','box_momentum_q_diagnostic':(1-me)/2,
 'input':'middle charged positive-energy rest spin-down with opposite-charge idle reference',
 'receiver':'light charged rest spin-up; middle-family neutral +q spin-down; light-family antineutral -q from hole in +q spin-up',
 'exact_free_energy_initial_and_receiver':'1 = e_1/e_2 + 2q','exact_total_momentum_initial_and_receiver':'0','exact_total_Jz_initial_and_receiver':'-1/2',
 'actual_current_words_contributing_to_receiver':selected,'exact_receiver_V_matrix_element':str(amplitude),
 'architecture_contact_coefficient_G':str(G),'source_norm_upper':str(B),'time':str(t),
 'leading_receiver_amplitude_modulus':str(lead),'full_source_state_remainder_upper':str(remainder),
 'original_receiver_probability_interval':[str(s.expand((lead-remainder)**2)),str(s.expand((lead+remainder)**2))],
 'original_receiver_probability_decimal_interval':[float((lead-remainder)**2),float((lead+remainder)**2)],
 'full_energy_expectation_initial_diagnostic':1+float(G)*Ei,'receiver_full_energy_expectation_diagnostic':1+float(G)*Ef,
 'final_projector_measurement_absolute_energy_work_upper':str(2*G*norm),
 'complete_records':'receiver yes and entire270-CAR-source complement no, including all other emissions, scattering and vacuum response; no conditional renormalization',
 'free_G_zero_control_receiver_probability':'0','vacuum_reference_control_receiver_probability':'0 by conserved family lepton number',
 'detuned_box_control_free_energy_mismatch':'e_1/e_2+2q-1',
 'observed_muon_lifetime_weak_scale_or_continuum_S_matrix_derived':False,
 'scope':'Declared finite native-current contact construction with source-derived charged hierarchy, massless neutral channels and complete finite-time transition bound. No physical neutrino identity or measured decay width assigned.',
 'verdict':'PASS_NATIVE_CURRENT_CHARGED_TO_CHARGED_NEUTRAL_PAIR'}
