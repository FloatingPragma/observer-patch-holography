"""Exact reducing state witnesses and complete-source current enclosures."""
from functools import lru_cache
from fractions import Fraction as F
from math import isqrt
import numpy as np
import sympy as sy
import scipy.sparse as sp
from scipy.linalg import eigh
from model.openwave_calibrated_hadrons.source import angular4,integral
from model.openwave_six_quark.certificate import spectral_enclosure
from model.openwave_baryon_weak_formfactors.native import currents,formfactors
from .fragments import fragment_controls
from .source import parameters,sector,components,exact_matrix,INPUTS,fit_at_colour,matrix,weights

def sqrt_upper(x,bits=60):
 den=2**bits;n=(x.numerator*den*den+x.denominator-1)//x.denominator;q=isqrt(n)
 if q*q<n:q+=1
 return F(q,den)

@lru_cache(None)
def reduced_line():
 d=sector('baryon');c=currents();n=len(d['basis']);eye=sp.eye(n,dtype=np.int64,format='csr')
 S4=angular4(d['modes'],d['basis'],(0,1),True);I4=angular4(d['modes'],d['basis'],(0,1),False)
 ids=[i for i,row in enumerate(d['labels'])if sum(w==0 for x,s,w in row)==1 and sum(s==0 for x,s,w in row)==2]
 B=(15*eye-S4)@(15*eye-I4)@(eye+integral(c['full_parity']))@d['Binteger'][:,ids]
 G=sy.Matrix((B.T@B).toarray().tolist());_,pivots=G.rref();B=B[:,list(pivots)].tocsr();Gram=sy.Matrix((B.T@B).toarray().tolist());L,D=Gram.LDLdecomposition(hermitian=True);assert all(x>0 for x in D.diagonal())
 assert (S4@B-3*B).nnz==0 and(I4@B-3*B).nnz==0 and(integral(c['full_parity'])@B-B).nnz==0
 A=sy.zeros(len(pivots));parts={}
 for k,coefficient in weights(parameters()).items():
  parts[k]=(B.T@d['parts'][k]@B).toarray();A+=sy.Rational(coefficient.numerator,coefficient.denominator)*sy.Matrix(parts[k].tolist())
 Li=L.inv();T=Li*A*Li.T;g=list(D.diagonal());H=sy.Matrix(len(g),len(g),lambda i,j:T[i,j]/sy.sqrt(g[i]*g[j]))
 X=integral(sum(c['full'][x][0][0]for x in range(2))).T@integral(sum(c['full'][x][3][0]for x in range(2)))
 assert (B.T@(X-X.T)@B).nnz==0
 transform=np.asarray(Li.T.evalf(18),float)@np.diag([1/np.sqrt(float(x))for x in g]);return dict(d=d,B=B,Gram=Gram,H=H,transform=transform,X=X,parts=parts)

def exact_action(A,v):
 A=A.tocsr();return [sum(int(A.data[k])*int(v[A.indices[k]])for k in range(A.indptr[i],A.indptr[i+1]))for i in range(A.shape[0])]

@lru_cache(None)
def ground_certificate():
 r=reduced_line();d=r['d'];H=r['H'];e,U=eigh(np.asarray(H.evalf(18),float));coef=r['transform']@U[:,0];integers=[round(float(x)*2**45)for x in coef];v=exact_action(r['B'],integers);norm=sum(x*x for x in v)
 Hv=[F(0)for _ in v]
 for k,c in weights(parameters()).items():
  for i,a in enumerate(exact_action(d['parts'][k],v)):Hv[i]+=c*a
 energy=sum(x*y for x,y in zip(v,Hv))/norm;residual=sum((x-energy*y)**2 for x,y in zip(Hv,v))/norm
 spec=spectral_enclosure(H);gap=F(spec['all_eigenvalue_intervals'][1][0])-energy;assert gap>0
 current=F(sum(a*b for a,b in zip(v,exact_action(r['X'],v))),norm);sine=sqrt_upper(residual)/gap;error=12*sine
 target=F(INPUTS['axial_calibration_target_magnitude']);unc=F(INPUTS['axial_quoted_uncertainty'])
 assert abs(current-target)+error<unc/100
 coloc=[]
 for mask in d['basis']:
  n0=sum((mask>>i)&1 for i,m in enumerate(d['modes'])if m[1]==0);coloc.append(n0 in(0,3))
 pcol=F(sum(x*x for x,yes in zip(v,coloc)if yes),norm)
 probs=[F(x*x,norm)for x in v];N=len(v)
 return {'exact_selected_sector_dimension':H.rows,'quantum_numbers':'native colour0,S=I=1/2, charge0, Sz=+1/2, even cell parity','complete_selected_spectrum':spec,
  'integer_native_occupation_trial':list(map(str,v)),'integer_trial_norm_squared':str(norm),'trial_energy_rational':str(energy),'trial_residual_squared_rational':str(residual),
  'gap_above_trial_lower_rational':str(gap),'ground_vector_sine_error_upper':str(sine),'trial_native_axial_expectation_rational':str(current),
  'native_ground_axial_interval_rational':[str(current-error),str(current+error)],'native_ground_axial_interval':[float(current-error),float(current+error)],
  'colocated_probability_rational':str(pcol),'colocated_ground_probability_interval':[float(max(F(0),pcol-2*sine)),float(min(F(1),pcol+2*sine))],
  'complete_trial_original_records':{'uniform_full_N3_input_dimension':N,'filter_success_probability':str(F(1,N)),
    'success_native_occupation_probabilities':[str(p/N)for p in probs],'failure_native_occupation_probabilities':[str((1-p)/N)for p in probs],
    'all_other_Fock_occupations':'retained as orthogonal zero-probability outcomes for this input','complete_nonselective_filter_and_occupation_mean_work':'0 on the uniform full occupation input',
    'exact_ground_complete_record_transfer_l1_bound':str(4*sine/N),'complete_filter_reference_stable_channel_bound':str(4*sine),'arbitrary_external_reference_transfer':True,
    'scope':'Exact integer state preparation at source time0; spectral bound transfers its full occupation records to the actual stationary ground line. All failures are retained.'}}

@lru_cache(None)
def spectra():
 out={}
 for kind in('baryon','meson'):
  d=sector(kind);H=exact_matrix(kind);charges=np.rint(np.diag(d['q']).real).astype(int);rows=[]
  for q in sorted(set(charges)):
   ids=np.flatnonzero(charges==q).tolist();others=np.flatnonzero(charges!=q).tolist();assert all(H[i,j]==0 for i in ids for j in others)
   rows.append({'charge':int(q),'source_indices':ids,'certificate':spectral_enclosure(H.extract(ids,ids))})
  out[kind]=rows
 return out

def certificate():
 p=parameters();unit=F(INPUTS['electron_unit_GeV']);g=ground_certificate();s=spectra();fits=[]
 for name,kind,q in(('proton','baryon',1),('neutron','baryon',0),('charged_pion','meson',1)):
  row=next(x for x in s[kind]if x['charge']==q);interval=[unit*F(v)for v in row['certificate']['all_eigenvalue_intervals'][0]];target=F(INPUTS['mass_calibration_targets_GeV'][name]);assert interval[0]<=target<=interval[1]
  fits.append({'comparison':name,'target_GeV':str(target),'same_source_energy_GeV_interval':list(map(str,interval)),'role':'calibration input'})
 f=formfactors(sector('baryon')['h']);controls=[]
 for U in(0,224,240,2048):
  par=fit_at_colour(U);ff=formfactors(matrix('baryon',par));controls.append({'colour_coefficient':U,'other_coefficients':{k:float(v)for k,v in par.items()},'gV':ff['gV'],'gA':ff['gA'],'scope':'Three masses are refitted on each diagnostic source member; this is not one fixed-parameter trajectory.'})
 return {'schema':'oph.native_split_site_axial_calibration.v1','verdict':'PASS_COMPLETE_NATIVE_STATE_AXIAL_REPAIR_WITH_FOUR_DISCLOSED_INPUTS',
  'full_native_CAR_modes':60,'input_ancestry':INPUTS,'parameters_rational':{k:str(v)for k,v in p.items()},'complete_colour_spectra':s,'three_mass_fits':fits,
  'current_ground_certificate':g,'floating_current_comparison':{k:f[k]for k in('gV','gA','current_leakage_operator_norm')},'colour_regime_controls':controls,'coloured_fragment_controls':fragment_controls(),
  'positive_full_Fock_floor':'H>=min(m_up,m_down)N_colourful+N_charged_lepton; neutral spectators retained. Every source term is positive.',
  'native_current_changed':False,'effective_axial_multiplier_installed':False,'fourth_empirical_input':'The measured magnitude1.2753 selects U and the resulting split-site state; this agreement is a calibration, not a prediction.',
  'local_scalar_obstruction':'The exact colocated nucleon coefficient stays5/3. The selected source uses split-site orbital/spin/flavour correlations outside that multiplicity-one local irrep.',
  'scope':'Complete finite positive native action, source current and calibrated spin-half rest-band records. No physical continuum, renormalized QCD current, outgoing lifetime or full hadron spectrum is identified.'}
