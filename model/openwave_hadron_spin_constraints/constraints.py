"""Same-current positive extensions, complete spectra and withheld comparisons."""
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from model.openwave_axial_state.source import parameters,matrix,components,INPUTS
from model.openwave_baryon_weak_formfactors.native import formfactors

UNIT=float(INPUTS['electron_unit_GeV']);TARGETS=INPUTS['inherited_mass_input_ancestry']['held_out_comparisons_GeV']

@lru_cache(None)
def block(kind,q,twice_spin):
 d=components(kind);spin=twice_spin*(twice_spin+2)/4
 Q=d['q']-q*np.eye(len(d['q']));S=d['spin']-spin*np.eye(len(d['spin']));e,v=eigh(Q@Q+S@S);return v[:,e<1e-8]

def extended_matrix(kind,p):
 d=components(kind);return matrix(kind,p)+float(p.get('global_spin',0))*d['spin']+float(p.get('charge_square',0))*d['q']@d['q']

def energy(kind,p,q,twice_spin):
 P=block(kind,q,twice_spin);return float(eigh(P.conj().T@extended_matrix(kind,p)@P,eigvals_only=True)[0])

@lru_cache(None)
def extension(gs=0.,chi=0.):
 p0=parameters();D=float(8*p0['U']/3-4*p0['A']/9);constant=float(3*p0['mbar']+p0['A']);target=float(INPUTS['mass_calibration_targets_GeV']['charged_pion'])/UNIT
 def par(A):return {**p0,'A':A,'U':3*(D+4*A/9)/8,'mbar':(constant-A-.75*gs-.5*chi)/3,'delta':float(p0['delta'])+chi,'global_spin':gs,'charge_square':chi}
 A=brentq(lambda a:energy('meson',par(a),1,0)-target,200.,1800.,xtol=1e-10);p=par(A);assert p['U']>0 and p['mbar']>p['delta']/2>0
 return p

def selected_energies(p):
 return {name:energy(kind,p,q,s)*UNIT for name,kind,q,s in(('proton','baryon',1,1),('neutron','baryon',0,1),('charged_pion','meson',1,0),('neutral_pion_candidate','meson',0,0),('rho_charged_candidate','meson',1,2),('Delta_plus_candidate','baryon',1,3))}

@lru_cache(None)
def certificate():
 p0=parameters();base=selected_energies(p0);delta0=float(p0['delta']);chi=((float(INPUTS['mass_calibration_targets_GeV']['charged_pion'])-float(TARGETS['neutral_pion']))/UNIT-delta0)/2
 cases=[]
 for charge in(0.,chi):
  for role in('no_spin_fit','Delta_fit','rho_fit'):
   gs=0.
   if role!='no_spin_fit':
    key='Delta_plus_candidate'if role=='Delta_fit'else'rho_charged_candidate';target=float(TARGETS['Delta_baryon'if role=='Delta_fit'else'rho_vector']);gs=brentq(lambda s:selected_energies(extension(s,charge))[key]-target,0.,900.,xtol=1e-10)
   p=extension(gs,charge);values=selected_energies(p);full={}
   for kind in('baryon','meson'):
    H=extended_matrix(kind,p);full[kind]=eigh(H,eigvals_only=True).tolist()
   f=formfactors(extended_matrix('baryon',p));assert abs(f['gA']-1.2753)<2e-8 and abs(f['gV']-1)<2e-10
   assert abs(values['proton']-float(INPUTS['mass_calibration_targets_GeV']['proton']))<2e-10 and abs(values['neutron']-float(INPUTS['mass_calibration_targets_GeV']['neutron']))<2e-10
   cases.append({'control':role+('_and_neutral_pion_fit'if charge else''),'parameters':{k:float(v)for k,v in p.items()},'rest_energies_GeV':values,'all_colour_singlet_levels':full,'gV':f['gV'],'gA':f['gA'],'calibration_ancestry':['proton','neutron','charged_pion','axial_magnitude']+(['neutral_pion']if charge else[])+(['Delta_baryon']if role=='Delta_fit'else['rho_vector']if role=='rho_fit'else[]),'physical_complete_hadron_assignment':False})
 # Exact operator inequalities give incompatible spin-coefficient intervals.
 p0delta=energy('baryon',p0,1,3)-energy('baryon',p0,1,1);needed_delta=(float(TARGETS['Delta_baryon'])-base['proton'])/UNIT
 delta_gs=(needed_delta-p0delta)/3;needed_rho=(float(TARGETS['rho_vector'])-base['charged_pion'])/UNIT
 rho_interval=[(needed_rho-.5)/2,needed_rho/2];assert delta_gs<rho_interval[0]
 from .hyperfine import hyperfine_certificate
 return {'native_colour_spin_completion':hyperfine_certificate(),'schema':'oph.native_hadron_spin_constraints.v1','verdict':'PASS_POSITIVE_NATIVE_CONTROLS_AND_COMMON_SPIN_OBSTRUCTION','base_source':'openwave_axial_state fixed four-input action',
  'baseline_rest_energies_GeV':base,'controls':cases,'exact_baryon_identity':'H_colour0=3*mbar+delta*(3-2Nu)/2+A+(8U/3-4A/9)*P_split+K+Js*sumSx²+Jw*sumIx²+eta*I²; added global_spin*S² and charge_square*Q² act as scalars at fixed S,Q.',
  'fixed_baryon_wavefunction_path':{'D':str(8*p0['U']/3-4*p0['A']/9),'delta_U_equals_delta_A_over_six':True,'native_axial_state_unchanged':True},
  'common_spin_obstruction':{'rho_minus_pion_exact_interval':'[2*global_spin,2*global_spin+2*Js]','Delta_minus_proton_exact_change':'3*global_spin','Delta_fit_global_spin':delta_gs,'rho_required_global_spin_interval':rho_interval,'incompatible':True},
  'charge_square_control':{'coefficient_from_fifth_empirical_input':chi,'exact_relations':'m_n-m_p=delta-chi; m_pi_charged-m_neutral_up_antiquark=delta+chi','leptons_affected':False,'physical_neutral_pion_wavefunction_identified':False},
  'arithmetic_scope':'All spectral numbers are deterministic floating diagnostics. The incompatibility follows exact operator inequalities with a gap above400 in the required coefficient intervals; no narrow numerical margin is used.',
  'full_Fock_source':'Each extension is a positive square of the actual global colourful spin or electric current, never a selected-composite projector. All full native Fock complements remain present.',
  'scope':'The common positive spin square cannot jointly repair rho and Delta. A separate charge square can fit the neutral scalar energy with a fifth empirical input, without identifying its observed isospin wavefunction. These controls do not replace the fixed four-input source.'}
