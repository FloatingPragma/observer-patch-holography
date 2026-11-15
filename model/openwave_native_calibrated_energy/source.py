"""Interaction from the actual corrected source and one exposed energy anchor."""
from pathlib import Path
from fractions import Fraction as F
import json
import numpy as np
from model.openwave_native_coherent_clock.source import rc,sel,evolve,pair_vectors,strings,infer_gap
from model.characteristic_constraints.bounds import trig_interval
ROOT=Path(__file__).resolve().parents[2]

def partial_transpose_binary(rho):
 # Matching order is spectator,binary,syndrome. Only the binary tensor
 # factor is transposed; the entire spectator/syndrome complement remains.
 x=rho.reshape(3,2,6,3,2,6);return x.transpose(0,4,2,3,1,5).reshape(36,36)

def run():
 h=tuple(F(x)for x in rc.logical_hamiltonian().diagonal());gaps=[h[2*j+1]-h[2*j]for j in range(3)]
 contrast=gaps[1]-gaps[0];assert contrast==-F(5,1024)
 time=F(256);steps=65536;rows=[]
 # One supplied joint source preparation, not a tensor product imposed
 # after interaction. Its code coordinates are binary+ x spectator(0+1).
 psi=(pair_vectors(0,'X')[0]+pair_vectors(1,'X')[0])/np.sqrt(2);rho=np.outer(psi,psi.conj())
 sin=trig_interval(abs(contrast)*time/2);ideal_neg=(sin[0]/2,sin[1]/2)
 for law in ('C','E'):
  out=evolve(rho,law,time,steps);err=F(rc.error_bound(law,time,steps)['trace_norm_upper_exact']);lower=max(F(0),ideal_neg[0]-err)
  values=np.linalg.eigvalsh(partial_transpose_binary(out));neg=float(sum(max(0.,-x)for x in values));assert lower>0 and neg>=float(lower)-2e-10
  zero=evolve(rho,law,F(0),steps);zero_neg=float(sum(max(0.,-x)for x in np.linalg.eigvalsh(partial_transpose_binary(zero))));assert zero_neg<2e-10
  rows.append({'law':law,'coherent_parameter_exact':str(time),'steps':steps,'complete_binary_vs_spectator_syndrome_negativity_lower_exact':str(lower),'ideal_negativity_interval_exact':strings(ideal_neg),'finite_source_trace_norm_error_exact':str(err),'actual_complete_negativity_diagnostic':neg,'phase_off_negativity_diagnostic':zero_neg})
 anchor=F(json.loads((ROOT/'model/openwave_single_valley_leptons/results.json').read_text())['charged_lepton_GeV'][0]['GeV_exact']);scale=anchor/min(map(abs,gaps));assert scale>0
 calibrated=[]
 for j,g in sorted(enumerate(gaps),key=lambda x:abs(x[1])):
  epsilon=F(rc.error_bound('E',time,steps)['trace_norm_upper_exact']);sine=trig_interval(g*time)
  interval=infer_gap((max(F(-1),sine[0]-epsilon),min(F(1),sine[1]+epsilon)),time)
  absolute=tuple(sorted(map(abs,interval)));assert interval[0]*interval[1]>0
  calibrated.append({'spectator':j,'signed_source_gap_exact':str(g),'absolute_source_gap_exact':str(abs(g)),'anchor_scaled_gap_GeV_exact':str(scale*abs(g)),'record_inferred_anchor_scaled_gap_GeV_interval_exact':strings(tuple(scale*x for x in absolute))})
 ratios=[F(x['absolute_source_gap_exact'])/min(map(abs,gaps))for x in calibrated];assert ratios==[F(1),F(3),F(12)]
 return {'schema':'oph.same_native_interaction_and_calibrated_gap.v1','verdict':'PASS_NATIVE_SPECTATOR_INTERACTION_AND_EXPOSED_COMMON_SCALE_TEST','shared_charge_clock_policy':'Use native_coherent_clock compensated inverse-feedback policy; uncorrected charged channel is not inherited.','same_source_energies_exact':strings(h),'signed_gaps_exact':strings(gaps),'interaction':{'exact_binary_spectator_coupling_difference':str(contrast),'decomposition':'h = I_binary tensor E0(s) + n_binary tensor Delta(s); Delta(s) is not scalar','generator_from_actual_complete_repair_compression':True,'target_interaction_Hamiltonian_installed':False,'two_hop_phase_architecture_supplied':True,'all_reference_and_syndrome_outputs_retained':True,'finite_entanglement_controls':rows},'energy_attachment':{'electron_anchor_GeV_exact':str(anchor),'anchor_role':'empirical electron energy applied to the smallest absolute source transition; species assignment unproved','one_common_GeV_per_source_unit_exact':str(scale),'rows':calibrated,'absolute_gap_ratios_exact':strings(ratios),'OpenWave_row11_comparison_ratios_exact':['1','207','3477'],'literal_three_lepton_ratio_match':False,'physical_rest_pole_or_particle_identification':False,'target_exposure':'retrospective; electron and charged-lepton hierarchy were exposed before this experiment','new_independent_GeV_prediction':False},'scope':{'one_binary_mode_with_three_spectator_values':True,'three_independent_lepton_species':False,'physical_seconds_or_work_scale_derived':False,'same_program_as_coherent_spatial_packet':False,'new_Nature_validation':False}}
