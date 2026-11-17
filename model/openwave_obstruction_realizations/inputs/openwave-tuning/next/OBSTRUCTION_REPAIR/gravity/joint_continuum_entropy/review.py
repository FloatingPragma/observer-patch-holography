import json
from pathlib import Path
from mixture_response import bound
outdir=Path('L40_h0p5')
result=json.loads((outdir/'RESULT.json').read_text())
mixtures=json.loads((outdir/'MIXTURE_BOUNDS.json').read_text())
allrows=[r for ball in mixtures['records'] for r in ball['records']]
if not all(r['low_rank']['numerical_gate_pass'] for r in allrows):raise ValueError('Do not form common entropy bounds from a failed determinant.')
epsilon=min(r['half_percent_controlled_interval']['epsilon'] for r in allrows)
common=[]
for ball in mixtures['records']:
    for row in ball['records']:
        interval=bound.finite_entropy_interval(row['actual_mixture_entropy_tangent'],epsilon,row['low_rank']['log_Z'])
        interval['relative_entropy_error_bound']=interval['relative_entropy_upper_bound']/interval['linear_entropy']
        common.append(dict(radius_cells=ball['radius_cells'],centre=row['centre'],**interval))
common_record=dict(status='ONE_COMMON_EXACTLY_NEUTRAL_MIXTURE_SOURCE',epsilon=epsilon,actual_total_energy=epsilon*mixtures['records'][0]['actual_branch_energy'],branch_records=[dict(label='vacuum',probability=1-epsilon,relative_charge=0),dict(label='four_particle_four_hole_Slater',probability=epsilon,relative_charge=0)],regional_observer_excludes_external_flag=True,records=common,analytic_entropy_bound=True,floating_point_interval_certified=False)
(outdir/'COMMON_MIXTURE_SOURCE.json').write_text(json.dumps(common_record,indent=2))
summaries=[]
for solution in result['solutions']:
    held=sorted(solution['held_out'],key=lambda x:abs(x['relative_residual']),reverse=True)
    middle=[r for r in held if r['radius']==7]
    axis=[r for r in held if r['centre'][1:]==[0,0]]
    summaries.append(dict(width=solution['width'],rank=solution['rank'],condition_number=solution['condition_number'],max_training_absolute_residual=max(abs(r['total_residual']) for r in solution['training']),max_held_relative=max(abs(r['relative_residual']) for r in held),max_held_axis_relative=max(abs(r['relative_residual']) for r in axis),unused_central_radius=[r for r in middle if r['centre']==[0,0,0]],max_unused_radius_relative=max(abs(r['relative_residual']) for r in middle),worst_held=held[:3],sigma_amplitudes=solution['sigma_amplitudes'],lapse_amplitudes=solution['lapse_amplitudes'],sigma_centre=solution['sigma_centre'],lapse_centre=solution['lapse_centre']))
review=dict(status='COMPLETED_REST_SLICE_PREPARATION_DIAGNOSTIC_WITH_HELD_OUT_FAILURE',solutions=summaries,uniform_gate=[dict(radius=b['physical_radius'],worst_uniform_contamination=max(r['uniform_residual_relative_to_this_metric'] for r in b['responses'])) for b in result['balls']],common_mixture_epsilon=epsilon,common_mixture_energy=common_record['actual_total_energy'],max_actual_matter_entropy_bound_fraction=max(r['relative_entropy_error_bound'] for r in common),periodic_mean_rows=[dict(radius=b['physical_radius'],values=[m['translation_mean_first_derivative'] for m in b['matter']]) for b in result['balls']],scope='Three separable axis-cosine modes, rest-slice finite cell observers and finite torus only. Failure does not exclude fuller geometry bases or gravity. No independent all-timelike-observer lapse equation or finite-geometry entropy pass.')
(outdir/'REVIEW.json').write_text(json.dumps(review,indent=2))
print(json.dumps(review,indent=2),flush=True)
