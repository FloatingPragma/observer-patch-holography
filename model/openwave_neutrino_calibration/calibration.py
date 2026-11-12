"""Explicit neutral family tensor in the native charged-family frame."""
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
import json,math
import numpy as np
from model.lepton_mass_compiler import analyzer

HERE=Path(__file__).resolve().parent

def inputs():return json.loads((HERE/'oscillation_inputs.json').read_text())
def sqrt_bounds(x,bits=160):
 x=F(x);n=math.isqrt((x.numerator*2**(2*bits))//x.denominator)
 return F(n,2**bits),F(n+1,2**bits)
def encoded(bounds):return list(map(str,bounds))

@lru_cache(None)
def calibration():
 d=inputs();boxes=analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2'];mid=[sum(map(F,b))/2 for b in boxes]
 charged_exact=[mid[1]/mid[2],F(1),mid[0]/mid[2]]
 electron=F(json.loads((HERE.parents[0]/'lepton_cycle/experimental_inputs.json').read_text())['electron_MeV']['value'])/1000
 unit=electron/charged_exact[0]
 s12,s23,s13=[F(d[k])for k in('sin2_theta12','sin2_theta23','sin2_theta13')]
 s,c=[math.sqrt(float(x))for x in(s12,s23,s13)],[math.sqrt(float(1-x))for x in(s12,s23,s13)]
 z=np.exp(1j*math.pi*float(F(d['delta_CP_degrees']))/180)
 R12=np.array([[c[0],s[0],0],[-s[0],c[0],0],[0,0,1]],complex)
 R13=np.array([[c[2],0,s[2]/z],[0,1,0],[-s[2]*z,0,c[2]]],complex)
 R23=np.array([[1,0,0],[0,c[1],s[1]],[0,-s[1],c[1]]],complex);R=R23@R13@R12
 m0=F(d['declared_lightest_mass_eV']);squares=[m0*m0,m0*m0+F(d['delta_m21_squared_eV2']),m0*m0+F(d['delta_m31_squared_eV2'])]
 source_squares=[x/(10**18*unit**2)for x in squares];delta=1-charged_exact[0]
 recoil_q_squared=(delta**4-2*delta**2*(source_squares[0]+source_squares[1])+(source_squares[0]-source_squares[1])**2)/(4*delta**2)
 assert max(source_squares)<F(5,10**10)**2 and recoil_q_squared>F(49,100)**2
 assert sum(sqrt_bounds(x)[1]for x in source_squares[:2])<delta
 bounds=[sqrt_bounds(x)for x in squares];physical=np.array([float(sum(b)/2)for b in bounds]);neutral=physical*1e-9/float(unit)
 assert np.linalg.norm(R.conj().T@R-np.eye(3))<1e-14
 assert np.max(neutral)<5e-10 and float(unit)<.106
 return {'R':R,'neutral':neutral,'charged':np.array(list(map(float,charged_exact))),
  'charged_exact':charged_exact,'middle_charged_unit_GeV_exact':unit,'mass_squared_eV2_exact':squares,
  'source_mass_squared_exact':source_squares,'source_recoil_q_squared_exact':recoil_q_squared,
  'physical_masses_eV':physical,'mass_eV_intervals':bounds,'empirical_electron_anchor_GeV':electron,'oscillation_inputs':d}

def neutral_family_tensor():
 """R diag(mu), creation-pair convention, middle charged mass equal to1."""
 d=calibration();return d['R']@np.diag(d['neutral'])

def rigid_obstruction():
 boxes=analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2'];tau,e,mu=[tuple(map(F,b))for b in boxes]
 ratio=((mu[0]**2-e[1]**2)/(tau[1]**2-e[0]**2),(mu[1]**2-e[0]**2)/(tau[0]**2-e[1]**2))
 d=inputs();sol=tuple(map(F,d['solar_3sigma_eV2']));atm=tuple(map(F,d['normal_atmospheric_3sigma_eV2']));inv=tuple(map(F,d['inverted_atmospheric_absolute_3sigma_eV2']))
 observed=(sol[0]/atm[1],sol[1]/atm[0]);inverted_observed=(sol[0]/inv[1],sol[1]/inv[0]);iratio=(1-ratio[1],1-ratio[0])
 assert ratio[1]<observed[0] and iratio[0]>inverted_observed[1]
 return {'normal_scale_independent_ratio_interval':encoded(ratio),'normal_empirical_box_ratio_interval':encoded(observed),
  'inverted_scale_independent_ratio_interval':encoded(iratio),'inverted_empirical_box_ratio_interval':encoded(inverted_observed),
  'one_scale_fits_both_mass_splittings':False,'normal_discrepancy_factor_lower':str(observed[0]/ratio[1]),
  'reason':'A common positive neutral scale cancels from both squared-mass ratios. Inverted ordering assigns the isolated light eigenvalue to m3 and the two heavier levels to m1,m2.'}

def catalogue():
 d=calibration();rows=[]
 for j,b in enumerate(d['mass_eV_intervals']):
  rows.append({'mass_eigenstate':j+1,'conditional_mass_eV':float(d['physical_masses_eV'][j]),
   'conditional_mass_GeV':float(d['physical_masses_eV'][j])*1e-9,'mass_eV_interval':encoded(b),
   'mass_GeV_interval':encoded((b[0]/10**9,b[1]/10**9)),
   'source_mass_squared_exact':str(d['mass_squared_eV2_exact'][j]/(10**18*d['middle_charged_unit_GeV_exact']**2)),
   'neutrino_and_antineutrino_same_mass':True,'active_and_sterile_coordinates_share_mass':True})
 s12=F(d['oscillation_inputs']['sin2_theta12']);s13=F(d['oscillation_inputs']['sin2_theta13'])
 beta_squared=sum(w*m for w,m in zip(((1-s13)*(1-s12),(1-s13)*s12,s13),d['mass_squared_eV2_exact']))
 return {'calibration_inputs':d['oscillation_inputs'],'middle_charged_unit_GeV':float(d['middle_charged_unit_GeV_exact']),
  'middle_charged_unit_GeV_exact':str(d['middle_charged_unit_GeV_exact']),
  'electron_anchor_GeV_exact':str(d['empirical_electron_anchor_GeV']),
  'source_recoil_q_squared_exact':str(d['source_recoil_q_squared_exact']),
  'exact_maximum_neutral_mass_bound':'1/2000000000','exact_recoil_q_lower_bound':'49/100',
  'charged_source_ratios_exact':list(map(str,d['charged_exact'])),'neutral_family_tensor_real':neutral_family_tensor().real.tolist(),
  'neutral_family_tensor_imag':neutral_family_tensor().imag.tolist(),'PMNS_real':d['R'].real.tolist(),'PMNS_imag':d['R'].imag.tolist(),
  'rows':rows,'sum_masses_eV_interval':encoded((sum(b[0]for b in d['mass_eV_intervals']),sum(b[1]for b in d['mass_eV_intervals']))),
  'effective_beta_mass_eV_diagnostic':math.sqrt(float(beta_squared)),'effective_beta_mass_eV_interval':encoded(sqrt_bounds(beta_squared)),
  'source_choices':{'mass_ordering':'normal','lightest_mass_eV':'0.001','Majorana_mass':'0 by conserved Dirac lepton number'},
  'oscillation_fit_is_calibration':True,'independent_neutrino_mass_or_mixing_prediction':False}
