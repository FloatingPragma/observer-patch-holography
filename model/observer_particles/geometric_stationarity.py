#!/usr/bin/env python3
"""Actual native entropy envelope, area dictionary, and radial backreaction.

No empirical inputs. The geometry is the existing static radial-product
class; proposed comparison variations and bulk-entropy dictionaries are
explicit. Curvature and matter stress are derived independently.
"""
from __future__ import annotations
from evidence_contract import public_receipt
import sys
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins
sys.dont_write_bytecode=True
import hashlib,json,math
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];PUBLIC=ROOT
sys.path.insert(0,str(PUBLIC));sys.path.insert(0,str(PUBLIC/'code'))
from model.native_equilibration import common_variation as cv,cut_completion as cut
from model.observer_time_fields import continuum_dynamics as field
from model.relative_observer_geometry import intrinsic_geometry as geo
from quantum_information import entropy_response

OUT=HERE/'results/geometric_stationarity.json';ART=HERE/'results/geometric_stationarity_arrays.npz'
SOURCES=(
 'code/geometry/ENTROPY_FIRST_LAW_AUDIT.md',
 'Lean/Thermodynamics/EinsteinPremiseLink.lean','Lean/Thermodynamics/EquationOfState.lean',
 'Lean/Geometry/EinsteinBranchInhabitant.lean','Lean/Geometry/CrossingReadAreaLaw.lean',
 'Lean/Geometry/SourceOrderEinsteinComposition.lean',
 'Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean',
 'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
 'Lean/ObserverPatchHolography/EinsteinBranch/FiniteCapGeneratorSplit.lean',
 'Lean/ObserverPatchHolography/EinsteinBranch/CollarPremiseDerivation.lean',
 'Lean/ObserverPatchHolography/EinsteinBranch/FixedCapacityWLaw.lean',
 'code/particles/calibration/derive_d10_gravity_shared_edge_entropy_bridge.py',
 'code/causal_refinement/crossing_read_area_law.py',
 'paper/tex_fragments/CROSSING_READ_AREA_LAW.tex','paper/tex_fragments/PAIRED_SOURCE_GRAVITY.tex',
 'model/relative_observer_geometry/README.md','Lean/docs/LIBRARY_GUIDE.md','Lean/docs/PROOF_INDEX.md')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean(x):
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(tuple,list)):return [clean(v) for v in x]
 if isinstance(x,(np.ndarray,s.MatrixBase)):return clean(x.tolist())
 if isinstance(x,(np.integer,np.floating)):return x.item()
 if isinstance(x,F):return str(x)
 if isinstance(x,s.Basic):return str(x)
 if isinstance(x,complex):return [x.real,x.imag]
 return x

def logm(a):
 v,u=np.linalg.eigh((a+a.conj().T)/2);assert min(v)>0
 return (u*np.log(v))@u.conj().T

def entropy(a):return float(-np.trace(a@logm(a)).real)
def logledger(x):return sum(float(v)*math.log(int(k)) for k,v in x.items())

def native_envelope():
 # Execute the strongest public finite source envelope, not a supplied MaxEnt ansatz.
 actual=cv.run()
 w=cv.common.reference();rows=[]
 for extra in combinations(range(2,12),2):
  cap=tuple(sorted((*cv.CORE,*extra)))
  cert=cut.cut_certificate(cap)
  boundary=sum((a in cap)!=(b in cap) for a,b in cut.quantum_mixing.seams())
  row=cv.region_entropy(cap,w)
  rows.append((boundary,cap,row,cert['internal_native_seams']))
 lo=min(rows,key=lambda r:(r[0],r[1]));hi=max(rows,key=lambda r:(r[0],r[1]))
 assert lo[0]<hi[0] and lo[2]['counted_edge']==hi[2]['counted_edge']
 assert all(row[2]==lo[2] for row in rows)
 # Same selected reference, same core moments, same cap cardinality, distinct boundary.
 return {'source_envelope':actual,
  'all_connected_complement_four_port_core_caps':len(rows),
  'same_cap_cardinality_sector_counts_and_entropy_ledger_exact':True,
  'minimum_crossing_seams':{'cap':lo[1],'count':lo[0],'complement_internal_seams':lo[3]},
  'maximum_crossing_seams':{'cap':hi[1],'count':hi[0],'complement_internal_seams':hi[3]},
  'shared_edge_entropy_exact_prime_logs':lo[2]['counted_edge'],
  'shared_edge_entropy_nats':logledger(lo[2]['counted_edge']),
  'boundary_area_proxy_no_go':'No nonzero universal coefficient equates this counted edge entropy to crossing-seam count across these actual selected-reference caps.',
  'scope':'This is a finite 2D native boundary proxy test, not rejection of whole-screen 3D boundary area or every physical area dictionary.'}

def metric_area_volume():
 mu,aa,bb,dd,L,z=s.symbols('mu aa bb dd L z',positive=True,real=True)
 screen=s.Matrix([[aa,bb],[bb,(bb*bb+dd)/aa]])
 spatial=s.diag(1,screen);C=mu*s.diag(-1,spatial)
 N2=mu*s.sqrt(dd);g=s.simplify(s.sqrt(-C.det())*C.inv())
 area=s.simplify(s.sqrt(g[2:4,2:4].det()))
 volume=s.simplify(s.sqrt(g[1:4,1:4].det()))
 assert area==mu and s.simplify(volume-mu*s.sqrt(N2))==0
 # A geometric comparison family within the actual static product action class.
 # Domain coordinates/radial boundary held fixed; coefficients really vary.
 muz=mu*s.exp(-2*z/3);screenz=s.exp(2*z)*screen
 Cz=muz*s.diag(-1,1,screenz)
 gz=s.simplify(s.sqrt(-Cz.det())*Cz.inv())
 Az=s.simplify(s.sqrt(gz[2:4,2:4].det()))
 Vz=s.simplify(L*s.sqrt(gz[1:4,1:4].det()))
 assert s.simplify(Vz-L*volume)==0
 assert s.simplify(s.diff(Az,z).subs(z,0)+2*mu/3)==0
 # Read actual native point coefficient; do not invent its spatial derivatives.
 read=geo.native_chart_increment(level=1,attempts=2)
 actual=s.Matrix([[s.Rational(v) for v in row] for row in read['averaged_screen_tensor']])
 actual_metric=geo.wave_density_metric(s.diag(1,actual).tolist())
 av=s.simplify(s.sqrt(actual_metric['metric'][2:4,2:4].det()))
 vv=s.simplify(s.sqrt(actual_metric['metric'][1:4,1:4].det()))
 assert av==1
 # At fixed coordinate beta each angularly excited oscillator's entropy falls
 # when its angular stiffness rises. This sign holds termwise in full Fock space.
 api=field.product_mode_source_control(6,radial_mode=1,radial_points=3,length=F(2),inverse_temperature=F(1),steps=1)
 a=float(s.N(s.sympify(api['angular_frequency_expression']),60));c=float(s.N(s.sympify(api['radial_frequency_expression']),60))
 omega=math.hypot(a,c);r=math.tanh(omega/2);p=(1-r)/2
 spin_derivative=-2*a*a*p*(1-p)
 bose_n=1/math.expm1(omega);bose_derivative=-a*a*bose_n*(bose_n+1)
 assert spin_derivative<0 and bose_derivative<0
 return {'metric_identity':str(g),'lapse_squared':str(N2),
  'whole_screen_area_density_exact':str(area),'spatial_volume_density_exact':str(volume),
  'screen_area_is_independent_of_radial_coordinate':True,
  'fixed_counting_measure_implies_fixed_physical_volume':False,
  'native_point_readout':read,'comparison_screen_tensor':clean(actual),
  'comparison_metric':clean(actual_metric['metric']),
  'comparison_mu_one_area_density':str(av),'comparison_mu_one_volume_density':str(vv),
  'fixed_volume_comparison':{
   'family':'mu(z)=mu0 exp(-2z/3), a_screen(z)=exp(2z) a0, radial coordinate interval unchanged',
   'proper_volume_exactly_constant':True,'area_log_derivative':'-2/3',
   'lapse_log_derivative':'2/3','field_operator':'B(z)=-d_r^2+exp(2z) A',
   'beta_convention':'same model-coordinate inverse temperature beta=1; local proper beta changes with lapse',
   'same_native_Hilbert_dimension_selects_this_variation':False,
   'two_spin_code_entropy_derivative':spin_derivative,
   'same_ideal_boson_mode_entropy_derivative':bose_derivative,
   'full_Fock_thermal_entropy_derivative':'-beta^2 sum_modes angular_eigenvalue*n(omega)*(n(omega)+1) < 0',
   'trial_generalized_entropy':'S_thermal(cavity)+boundary_area/(4 Newton), with fixed positive Newton',
   'trial_generalized_entropy_derivative_strictly_negative':True,
   'reason':'Every thermal-mode term is nonpositive, at least one is negative, and each radial-end area derivative is negative.',
   'scope':'A genuine fixed-proper-volume cavity comparison inside the inherited metric/action class. It is not a demonstrated source-selected variation or a small-ball vacuum-entanglement entropy identification.'}}

def radial_curvature():
 t,r,x,y=s.symbols('t r x y',real=True);coords=(t,r,x,y)
 f=s.Function('f')(x,y);u=s.Function('u')(x,y);v=s.Function('v')(x,y);w=s.Function('w')(x,y)
 g=f*s.diag(-1,1,s.Matrix([[u,v],[v,w]]));inv=g.inv()
 @lru_cache(None)
 def Gamma(i,j,k):
  return sum(inv[i,l]*(s.diff(g[l,k],coords[j])+s.diff(g[l,j],coords[k])-s.diff(g[j,k],coords[l]))/2 for l in range(4))
 def Ricci(i,j):
  return sum(s.diff(Gamma(k,i,j),coords[k])-s.diff(Gamma(k,i,k),coords[j])+
    sum(Gamma(k,k,l)*Gamma(l,i,j)-Gamma(k,j,l)*Gamma(l,i,k) for l in range(4)) for k in range(4))
 rtt,rrr,rtr=Ricci(0,0),Ricci(1,1),Ricci(0,1)
 assert s.simplify(rtt+rrr)==0 and s.simplify(rtr)==0
 A,L=s.symbols('amplitude length',positive=True,real=True)
 phi=A*s.sin(s.pi*r/L)*s.cos(s.pi*t/L)
 assert s.simplify(s.diff(phi,t,2)-s.diff(phi,r,2))==0
 plus=s.simplify((s.diff(phi,t)+s.diff(phi,r))**2)
 minus=s.simplify((s.diff(phi,t)-s.diff(phi,r))**2)
 point={t:0,r:L/4};stress=s.simplify(plus.subs(point));assert stress==A*A*s.pi*s.pi/(2*L*L)
 assert s.simplify(minus.subs(point)-stress)==0
 return {'metric_class':'g=f(x,y)*diag(-1,1,h_screen(x,y)), f>0, h_screen positive; arbitrary smooth screen dependence',
  'Ricci_tt_plus_Ricci_rr_exact':0,'Ricci_tr_exact':0,
  'both_radial_null_Einstein_contractions_exact':0,
  'cosmological_metric_term_changes_radial_null_balance':False,
  'actual_restored_angular_constant_mode':'A 1=0; phi=amplitude sin(pi r/L) cos(pi t/L)',
  'same_product_wave_equation_and_Dirichlet_boundary_exact':True,
  'independent_coherent_stress_difference':'delta T_mu_nu=d_mu phi d_nu phi - g_mu_nu*(d phi)^2/2',
  'radial_null_stress_at_t0_rLover4':str(stress),
  'same_metric_reference_and_excited_state_can_both_satisfy_Einstein_with_positive_coupling':False,
  'positive_coupling_residual_at_length2_amplitude1':'-pi^3 Newton',
  'first_order_parameter_warning':'Single coherent displacement changes stress at order amplitude^2; this alone is not a failure of a first-order entropy law.',
  'first_order_state_family':'Equal classical mixture of +/-sqrt(epsilon) displacement gives delta rho=-[F,[F,sigma]]/2 per epsilon and delta T_kk linear in epsilon.',
  'required_new_response':'For this fixed-reference scalar increment with other sources fixed, a nonzero geometric response outside the static radial-product ansatz is necessary: delta G_kk=8 pi Newton delta T_kk. Independently changing compensating stress is not excluded.',
  'scope':'This is exact in each smooth description in the inherited class, not a pointwise curvature assigned to a nonsmooth native limit or a no-go for dynamic OPH geometries.'}

def source_mixture():
 # Actual restored constant angular mode, public field source at beta1.
 api=field.product_mode_source_control(0,radial_mode=1,radial_points=3,length=F(2),inverse_temperature=F(1),steps=1)
 omega=float(s.N(s.sympify(api['joint_frequency_expression']),60));pol=api['state_clock_parameters']['polarization'];rf=float(pol)
 ident=np.eye(2);X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex)
 rotation=[field.source_compiler.pulse('Y',-float(field._source_dyadic_scalar(s.pi/4)),[0])]
 rotations=sum((field.source_compiler.remap(rotation,[j]) for j in range(2)),[])
 prep=field._thermal_code_preparation(pol,2)+rotations
 blank=np.zeros(16,complex);blank[0]=1;purified=field.source_compiler.execute(prep,4,blank)
 amp=purified.reshape(4,4);sigma=amp.T@amp.conj()
 H=np.kron(omega*(ident-X)/2,ident)+np.kron(ident,omega*(ident-X)/2)
 Fop=math.sqrt(omega/rf)*(np.kron(Y,ident)+np.kron(ident,Y))/2
 comm=lambda a,b:a@b-b@a
 tangent=-comm(Fop,comm(Fop,sigma))/2
 ds=entropy_response.entropy_tangent(sigma,tangent);de=float(np.trace(H@tangent).real)
 assert ds>0 and abs(ds-de)<2e-10 and abs(np.trace(tangent))<1e-12
 rows=[];arrays={'mixture_reference_sigma':sigma,'mixture_H':H,'mixture_generator_F':Fop,'mixture_tangent':tangent}
 for eps in (1e-4,1e-5,1e-6):
  states=[];marginals=[];errors=[]
  for sign in (-1,1):
   word=field.spin_weyl_word([[0,sign*math.sqrt(eps*omega/rf)]],copies=2)
   out=field.source_compiler.execute(prep+word,4,blank)
   states.append(out/math.sqrt(2));a=out.reshape(4,4);marginals.append(a.T@a.conj())
   errors.append(abs(np.linalg.norm(out)**2-1))
  rho=sum(marginals)/2
  ds_eps=(entropy(rho)-entropy(sigma))/eps
  cq=math.log(2)+(entropy(marginals[0])+entropy(marginals[1]))/2
  assert abs(np.trace(rho)-1)<1e-12 and max(errors)<1e-12
  assert abs(cq-(math.log(2)+entropy(sigma)))<2e-12
  rows.append({'epsilon':eps,'all_coin_probabilities':[float(np.vdot(v,v).real) for v in states],
    'unread_work_entropy_secant':ds_eps,'first_order_entropy_prediction':ds,
    'secant_derivative_error':ds_eps-ds,
    'full_coin_plus_work_CQ_entropy_change':cq-math.log(2)-entropy(sigma),
    'energy_secant':float(np.trace(H@(rho-sigma)).real)/eps,
    'word_normalization_error':max(errors)})
  arrays['mixture_states_eps_'+str(eps)]=np.array(states)
 assert abs(rows[-1]['secant_derivative_error'])<abs(rows[0]['secant_derivative_error'])
 return {'domain':'actual public two-copy constant-angular/radial-mode code, plus two purification references and retained classical displacement sign',
  'fixed_ideal_H_and_reference_beta':1,'frequency':omega,'state_variation':'rho_epsilon=1/2(U_+ sigma U_+*+U_- sigma U_-*), U_+/-=exp(+/- i sqrt(epsilon) F)',
  'source_words_executed_for_both_signs':True,'all_outcomes_retained':True,
  'analytic_first_order_entropy_derivative':ds,'analytic_first_order_energy_derivative':de,
  'trace_zero_tangent_error':float(abs(np.trace(tangent))),
  'positive_unread_work_derivative_at_fixed_background_geometry':True,
  'complete_coin_work_CQ_entropy_derivative':0,
  'finite_source_checks':rows,
  'scope':'Finite spin-code entropy, not continuum local-ball entanglement. Choice of which records belong to bulk matters; no sign-conditioned success filter or memory erasure is used.'},arrays

def warped_focusing():
 # A minimal dynamic extension, not a new selected source action.
 u,v,theta,azimuth=s.symbols('u v theta azimuth',real=True)
 coords=(u,v,theta,azimuth);R=s.Function('R')(u,v)
 g=s.Matrix([[0,-1,0,0],[-1,0,0,0],[0,0,R**2,0],[0,0,0,R**2*s.sin(theta)**2]])
 inv=g.inv()
 @lru_cache(None)
 def Gamma(i,j,k):
  return s.simplify(sum(inv[i,l]*(s.diff(g[l,k],coords[j])+s.diff(g[l,j],coords[k])-s.diff(g[j,k],coords[l]))/2 for l in range(4)))
 def Ricci(i,j):
  return s.simplify(sum(s.diff(Gamma(k,i,j),coords[k])-s.diff(Gamma(k,i,k),coords[j])+sum(Gamma(k,k,l)*Gamma(l,i,j)-Gamma(k,j,l)*Gamma(l,i,k) for l in range(4)) for k in range(4)))
 ruu=Ricci(0,0);rvv=Ricci(1,1)
 assert s.simplify(ruu+2*s.diff(R,u,2)/R)==0
 assert s.simplify(rvv+2*s.diff(R,v,2)/R)==0
 assert all(Gamma(i,0,0)==0 for i in range(4))
 expansion=2*s.diff(R,u)/R
 assert s.simplify(s.diff(expansion,u)+expansion**2/2-2*s.diff(R,u,2)/R)==0
 return {'declared_extension':'ds^2=-2 du dv+R(u,v)^2 dOmega^2; u at fixed v is affine null coordinate',
  'independent_null_curvature_uu':str(ruu),'independent_null_curvature_vv':str(rvv),
  'necessary_null_Einstein_response':'d^2R/du^2=-4 pi Newton R T_uu',
  'general_base_covariant_form':'nabla_k nabla_k R=-4 pi Newton R T_kk for affinely parametrized base-null k in a smooth spherical warped product',
  'area_expansion_equation':'d theta/du=-theta^2/2-8 pi Newton T_uu, theta=(1/Area)dArea/du',
  'metric_chosen_to_define_T':False,'full_Einstein_system_or_source_selector_solved':False,
  'scope':'Necessary backreaction target obtained from an independently specified metric and matter tensor. It supplies a local response equation, not an OPH derivation of its selection or of a full solution.'}


def source_hashes():
 paths={PUBLIC/p for p in SOURCES}
 assert all(p.is_file() for p in paths)
 for module in tuple(sys.modules.values()):
  path=getattr(module,'__file__',None)
  if path:
   p=Path(path).resolve()
   if p.is_relative_to(PUBLIC) and p.is_file():paths.add(p)
 return bound_pins({p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths)})

def run():
 native=native_envelope();print('Native common envelope and boundary control: PASS',flush=True)
 metric=metric_area_volume();curvature=radial_curvature();focus=warped_focusing();print('Area, volume and independent radial curvature: PASS',flush=True)
 mixture,arrays=source_mixture();np.savez_compressed(ART,**arrays)
 return {'status':'PASS_SCOPED_GEOMETRIC_ENTROPY_AND_BACKREACTION_TESTS','producer_sha256':sha(__file__),
  'source_sha256':bound_pins(source_hashes()),'artifacts':{ART.name:sha(ART)},'empirical_inputs':[],
  'native_common_entropy':native,'inherited_area_volume_dictionary':metric,
  'independent_radial_null_curvature_and_stress':curvature,'necessary_dynamic_area_response':focus,
  'actual_source_first_order_mixture':mixture,
  'additional_source_routes':{
    'common_entropy_envelope':'Executed actual C/E source; two conserved moments select a smooth family, shared by all1024core-containing caps. These moments are not physicalvolume.',
    'crossing_read_area_law':'Positive actual golden-source continuum crossing-area result under its population/read-radius/layer inputs; one nat/read and Planck calibration are declared. Not the240-port product-field source.',
    'Einstein_register_inhabitant':'Joint finite consistency and first-law type attachment; geometry and small-ball contractions are explicitly defined from stress. Not an independentcurvature solution.',
    'D10_shared_edge_entropy':'log(dimR3 dimR2)=logdimR3+logdimR2, conditional P/4/pixel cancellation. No physical area deformation/clock/scale source.',
    'collar_and_capacity':'Additional defect independence/cut spacing/capacity and physical dictionaries are premises, not equations selecting the current geometry response.'},
  'closure_assessment':{'native_finite_common_envelope_closed':True,'whole_screen_area_density_mapped_from_same_action':True,
    'fixed_native_count_or_core_moments_identified_with_fixed_physical_volume':False,
    'physical_Sgen_stationarity_derived':False,'static_product_metric_generic_matter_backreaction_compatible':False,
    'all_OPH_dynamic_geometries_excluded':False}}

if __name__=='__main__':
 result=clean(run());OUT.write_text(json.dumps(public_receipt(result),indent=2,sort_keys=True)+'\n', encoding='utf-8', newline='\n')
 print(result['status'])
