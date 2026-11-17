"""Frozen-data comparison, separate from source generation."""
from pathlib import Path
import json
ROOT=Path(__file__).parent
r=json.loads((ROOT/'RESULT.json').read_text())
# These are the local ROW_10_symmetry comparison packet's values. They
# are comparison data and are not imported by run.py or review.py.
reference={'sin2_12_13_23':[.308,.02215,.470],
 'asymmetric_one_sigma_minus':[.011,.00058,.013],
 'asymmetric_one_sigma_plus':[.012,.00056,.017],
 'splitting_ratio':.029805,'splittings_eV2':[7.49e-5,2.513e-3],
 'source':'Frozen ROW_10_symmetry/PROOF.md, NuFIT6.0 normal ordering IC24+SK comparison',
 'independent_current_dataset_verification':False}
angles=list(map(float,r['polar_orientation_sin2_12_13_23']))
pulls=[(a-b)/(lo if a<b else hi)for a,b,lo,hi in zip(angles,reference['sin2_12_13_23'],
 reference['asymmetric_one_sigma_minus'],reference['asymmetric_one_sigma_plus'])]
ratio=float(r['splitting_ratio'])
old=json.loads((ROOT/'inputs/chain_result.json').read_text())
real=next(q for q in old['candidates']if q['candidate']=='real_control')
out={'reference_comparison_inputs':reference,'source_chosen_with_problem_known':True,
 'source_PLAN_fixed_before_this_candidate_run':True,'parameters_retuned_after_comparison':False,
 'polar_angle_diagnostic':angles,'angle_pulls_against_frozen_packet':pulls,
 'all_three_within_three_sigma':all(abs(x)<=3 for x in pulls),
 'physical_light_current_exactly_unitary':False,
 'actual_heavy_current_weights':r['all_heavy_current_weights'],
 'full_pole_splitting_ratio':ratio,'ratio_over_reference':ratio/reference['splitting_ratio'],
 'overall_mass_rescaling_can_repair_splitting_ratio':False,
 'conditional_mass_eV':r['light_eV_conditional'],'conditional_splittings_eV2':r['splittings_eV2_conditional'],
 'zero_phase_control':{'complete_basis_equivalence_verified':True,'bright_Weyl_modes':33,
  'dark_Weyl_modes':30,'dark_Dirac_pairs':15,'dark_masses_exact':1,'dark_current_exact':0,
  'real_chain_orientation':real['unitary_orientation_sin2_12_13_23']},
 'source_lightness_mechanism_verified':True,'measured_neutrino_spectrum_and_mixing_solved':False,
 'physical_interacting_source_or_native_compiler_transfer':False}
(ROOT/'COMPARISON.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
