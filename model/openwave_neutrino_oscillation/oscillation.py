"""Native family geometry, charged-current frames and complete oscillation laws."""
from functools import lru_cache
import numpy as np
import sympy as s
from model.native_flavour_transport.transport import source,word_unitary,numeric,record,zero
from model.common_source import calibrated_matter as internal
from model.common_source.chiral_candidate import action

@lru_cache(maxsize=1)
def construction():
 y=source()['y'];V=word_unitary((2,));M=s.expand(1)*V.conjugate()*y*V.H
 M=M.applyfunc(s.expand);assert zero(M-M.T)
 labels=action.ordering()['all_left_components'];ell=internal._LH_tensor()*s.Matrix([0,1])
 neutral=next(i for i in range(15)if ell[i]!=0);assert abs(ell[neutral])==1
 D=internal._oriented_components('charged_lepton')[1]
 charged,singlet=next(iter(D.todok()));sign=D[charged,singlet]
 assert labels[neutral]['electric_Q']=='0' and labels[charged]['electric_Q']=='-1'
 k=s.kronecker_product(D,y)+s.kronecker_product(D,y).T+s.kronecker_product(ell*ell.T,M)
 Q=action.electric_charge();assert zero(Q.T*k+k*Q)
 eps=s.Matrix([[0,1],[-1,0]]);pair=s.kronecker_product(eps,k.conjugate());assert zero(pair+pair.T)
 # The inherited LH invariant supplies the complete internal Ward equation.
 ward=0;LH=internal._LH_tensor()
 for generator,fundamental in zip(action.generators(families=False),action._representation()[1]):
  assert zero(generator.T*LH+LH*fundamental[3:,3:]);ward+=1
 masses,F=np.linalg.eigh(numeric(y));F=np.real_if_close(F).real
 W=np.eye(45,dtype=complex)
 for i,block in [(charged,F),(singlet,float(sign)*F),(neutral,numeric(V)@F)]:W[3*i:3*i+3,3*i:3*i+3]=block
 transformed=W.T@numeric(k)@W
 assert np.max(abs(transformed[3*neutral:3*neutral+3,3*neutral:3*neutral+3]-np.diag(masses)))<2e-13
 assert np.max(abs(transformed[3*charged:3*charged+3,3*singlet:3*singlet+3]-np.diag(masses)))<2e-13
 current=W.conj().T@numeric(internal.weak_raising_current())@W
 B=current[3*neutral:3*neutral+3,3*charged:3*charged+3]
 # Retain the actual internal sign of the current; its common phase cancels.
 qsign=complex(internal.weak_raising_current(families=False)[neutral,charged]);assert abs(qsign)==1
 C=B.conj().T/qsign.conjugate()
 expected=F.T@numeric(V)@F;assert np.max(abs(C-expected))<2e-13
 assert np.max(abs(C.conj().T@C-np.eye(3)))<2e-13
 return dict(y=y,V=V,M=M,k=k,masses=masses,C=C,ward=ward,neutral=neutral,charged=charged,
             current_sign=qsign.real,frame_residual=float(np.max(abs(C-expected))))


def records(C,masses,momentum,time,initial):
 if momentum<=0 or time<0 or initial not in range(3):raise ValueError('positive momentum, nonnegative time and one native family required')
 C=np.asarray(C,complex);m=np.asarray(masses,float)
 if C.shape!=(3,3)or m.shape!=(3,)or np.min(m)<0 or np.max(abs(C.conj().T@C-np.eye(3)))>1e-12:raise ValueError('complete unitary family frame and nonnegative masses required')
 E=np.sqrt(momentum**2+m*m);w2=(E+momentum)/(2*E)
 K=np.diag(np.sqrt(w2))@C.conj().T
 evolution=np.diag(np.exp(-1j*(E-E[0])*time));T=K.conj().T@evolution@K
 prep=K[:,initial];out=T[:,initial];p=list(map(float,abs(out)**2));success=float(np.vdot(prep,prep).real)
 failure0=1-success;failure1=success-sum(p)
 probabilities=p+[failure0,failure1]
 assert min(probabilities)>-3e-14 and abs(sum(probabilities)-1)<2e-13
 ideal=C@evolution@C.conj().T;ideal_p=list(map(float,abs(ideal[:,initial])**2))+[0.,0.]
 variation=sum(abs(a-b)for a,b in zip(probabilities,ideal_p))/2
 r=float(max(1-w2));assert variation<=2*r+2e-13
 return {'momentum':momentum,'source_time':time,'initial_weak_family':initial,'positive_energies':E.tolist(),
         'positive_band_left_chirality_weights':w2.tolist(),
         'complete_original_probabilities':probabilities,'labels':['family0','family1','family2','preparation_failure','weak_return_complement'],
         'ideal_unitary_oscillation_probabilities':ideal_p,'complete_law_TV_from_ideal':variation,
         'proved_uniform_TV_upper':2*r}


def certificate():
 d=construction();C=d['C'];m=d['masses'];rows=[];contrast=0.
 for p in(1.,3.,10.):
  for t in(0.,1.,10.,100.):
   for a in range(3):
    r=records(C,m,p,t,a);r['antineutral_complete_probabilities']=records(C.conjugate(),m,p,t,a)['complete_original_probabilities']
    contrast=max(contrast,max(abs(np.array(r['complete_original_probabilities'])-r['antineutral_complete_probabilities'])));rows.append(r)
 identity=[records(np.eye(3),m,3.,100.,a)for a in range(3)]
 assert all(all(abs(r['complete_original_probabilities'][b])<1e-26 for b in range(3)if b!=a)for a,r in enumerate(identity))
 phases=np.diag(np.exp(1j*np.array([.2,.7,1.1])))
 phase_error=max(np.max(abs(np.array(records(C@phases,m,3.,100.,a)['complete_original_probabilities'])-records(C,m,3.,100.,a)['complete_original_probabilities']))for a in range(3))
 assert phase_error<2e-13 and contrast>1e-9
 B=abs(C)**2;s13=B[0,2];J=float(np.imag(C[0,0]*C[1,1]*np.conj(C[0,1]*C[1,0])))
 return {'schema':'oph.native_current_neutrino_oscillation.v1',
  'source_family_Y':record(d['y']),'neutral_symmetric_coefficient':record(d['M']),
  'declared_port_word':[2],'declared_pulse_phase':'-i','neutral_higgs_component_norm':'1',
  'source_neutral_masses_in_native_units':m.tolist(),'actual_charged_current_mixing_real':C.real.tolist(),'actual_charged_current_mixing_imag':C.imag.tolist(),
  'model_mixing_parameters':{'sin2_theta13':float(s13),'sin2_theta12':float(B[0,1]/(1-s13)),'sin2_theta23':float(B[1,2]/(1-s13)),'J':J},
  'full_native_LH_Ward_identities':d['ward'],'full45_mass_kernel_and90_CAR_pairing_retained':True,
  'actual_current_frame_residual':d['frame_residual'],'complete_records':rows,
  'maximum_executed_neutral_antineutral_probability_contrast':contrast,
  'zero_port_phase_control_has_no_family_appearance':True,'Majorana_column_phase_control_error':float(phase_error),
  'complete_instrument':'K=diag(sqrt(w_i)) C^dagger; F0=sqrt(I-K^dagger K); successful positive-band evolution followed by K^dagger and F1=sqrt(I-K K^dagger). All outputs retained.',
  'exact_spectral_definition':'lambda_i are the ordered roots of det(xI-Y); F is a real orthonormal ordered eigenframe; C=F^T[I+(-i-1)P2]F; Mnu=V*Y V^dagger.',
  'imported_PMNS_angles_mass_splittings_or_neutrino_GeV_scale':False,
  'measured_PMNS_fit_or_observed_neutrino_identification':False,
  'scope':'Declared free same-Higgs neutral architecture from native family Y and one coherent port phase; complete chirality-weighted weak-current records, with no empirical neutrino parameters.',
  'verdict':'PASS_NATIVE_FAMILY_NEUTRINO_OSCILLATION_RECORDS'}
