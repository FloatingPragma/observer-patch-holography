"""One electron anchor for the declared native-current composite source."""
from decimal import Decimal,localcontext
from pathlib import Path
import json
from model.lepton_mass_compiler.compiler import electron_calibration
from model.openwave_native_composites.spectrum import selected,U,J,T,MASS
HERE=Path(__file__).resolve().parent

def certificate():
 with localcontext() as ctx:
  ctx.prec=40
  inherited=electron_calibration();old=inherited['rows'];anchor=Decimal(old[0]['electron_anchored_GeV'])
  light=Decimal(old[0]['compiled_source_rest_energy'])
  ratios=[Decimal(r['compiled_source_rest_energy'])/light for r in old]
  calibration={'energy_unit_GeV':str(anchor),'electron_empirical_anchor_GeV':str(anchor),'family_number_energies':[str(r) for r in ratios],
   'phase_ancestry':'theta=2/9 and family labels inherited from retrospective lepton-cycle architecture',
   'mass_architecture':'number energies reuse native analyzer effect ratios; this positive-Laplacian action is distinct from the Weyl/Higgs pairing action',
   'positive_interaction_parameters':{'U':U,'J':J,'t':T},'parameters_fitted_to_hadrons':False}
  rows=[]
  for family,(name,r) in enumerate(zip(('light','second','third'),ratios),1):
   q=r+Decimal(4)/3*Decimal(str(U))+Decimal(3)/4*Decimal(str(J))
   rows.append({'channel':name+' charged lepton','family':family,'source_energy':str(r),'conditional_GeV':str(anchor*r),'meaning':'constant-mode number energy; inherited electron-calibrated family value'})
   rows.append({'channel':name+' up/down quark or conjugate','family':family,'source_energy':str(q),'conditional_GeV':str(anchor*q),'meaning':'finite stationary excitation energy; not an MSbar quark mass or relativistic pole'})
  for kind,charge,spin,label in [('meson',1,0,'charged spin-zero meson'),('meson',0,0,'neutral spin-zero meson'),('meson',1,1,'charged spin-one meson'),('baryon',1,.5,'charge-one spin-half baryon'),('baryon',0,.5,'charge-zero spin-half baryon'),('baryon',2,1.5,'charge-two spin-three-halves baryon')]:
   r=selected(kind,charge,spin);e=Decimal(str(r['ground_energy']))
   rows.append({'channel':label,'source_energy':str(e),'conditional_GeV':str(anchor*e),'meaning':'bound finite-source energy; physical hadron identification unproved'})
  meson=Decimal(str(selected('meson',1,0)['ground_energy']));delta=ratios[1]-ratios[0]
  rows.append({'channel':'light/second-family spin-zero meson','source_energy':str(meson+delta),'conditional_GeV':str(anchor*(meson+delta)),'meaning':'binding eigenvectors unchanged by constant family number-energy shift; kaon-like assignment is conditional'})
  observed=json.loads((HERE/'comparison_inputs.json').read_text())
  obs={k:Decimal(v) for k,v in observed['GeV'].items()}
  predicted={r['channel']:Decimal(r['conditional_GeV']) for r in rows}
  comparisons=[]
  for target,label in [('pion_charged','charged spin-zero meson'),('proton','charge-one spin-half baryon'),('kaon_charged','light/second-family spin-zero meson')]:
   p=predicted[label];o=obs[target]
   comparisons.append({'identification_hypothesis':label+' = '+target,'predicted_GeV':str(p),'comparison_GeV':str(o),'ratio_to_comparison':str(p/o),'relative_error':str(p/o-1),'verdict':'REJECTED physical mass assignment for this fixed architecture'})
  exact={'proton_neutron_energy_difference':'0 for all positive U,J,t at equal light family masses',
   'light_second_meson_energy_difference_GeV':str(anchor*delta),
   'observed_charged_kaon_minus_pion_GeV':str(obs['kaon_charged']-obs['pion_charged']),
   'neutral_meson_singlet_removed':False,'independent_up_down_splitting':'absent by exact weak-flavour symmetry',
   'family_offset_cancels_binding':'H_(a,b)=H_(1,1)+(m_a+m_b-2m_1)I in the same meson sector'}
  assert all(Decimal(c['ratio_to_comparison'])<Decimal('.5') for c in comparisons)
  return {'verdict':'PASS: complete common-scale comparison rejects the stated physical hadron assignments','calibration':calibration,'energy_rows':rows,'comparison_only_inputs':observed,'physical_assignment_tests':comparisons,'exact_obstructions':exact,
   'bound_energy_upper_GeV':{'spin_zero_light_meson':str(anchor*Decimal('2.25')),'spin_half_light_baryon':str(anchor*Decimal('3.5625'))},
   'physical_particle_rest_mass_claim':False,'continuum_or_Lorentz_dispersion_proved':False}
