"""Resume the user-interrupted production calculation from real checkpoints."""
import json,time
from pathlib import Path
import numpy as np
from kernels import cell,response_kernels
from regions import build,aggregates,trace_response
from matter import prepare_matter
from mixture_response import response_with_bound

outdir=Path('L40_h0p5');start=time.process_time()
metadata=json.loads((outdir/'METADATA.json').read_text());cfg=metadata['configuration']
L=cfg['L'];h=cfg['h'];aliases=cfg['aliases']
records=[json.loads((outdir/f'BALL_r{r}.json').read_text()) for r in cfg['radii']]
source=cell.prepare(L,h,1.,0.,aliases,.5)
matter=prepare_matter(source,6.)
missing=sorted({j for data in records for j in [1,2,3] if not any(r['shift']==[j,0,0] for r in data['responses'])})
kernels={}
for j in missing:
    kernels[j]=response_kernels(L,h,1.,0.,aliases,.5,(j,0,0))
    print(json.dumps(dict(prepared_missing_shift=j,cpu_seconds=time.process_time()-start)),flush=True)

mixture_records=[]
for data in records:
    radius=data['radius_cells'];region=build(source,radius)
    required=[j for j in missing if not any(r['shift']==[j,0,0] for r in data['responses'])]
    if required:
        d,acc=aggregates(region,source,[(j,0,0) for j in required])
        for j,a in zip(required,acc):
            km,kn=kernels[j];raw=trace_response(a,km,source,d);lapse=trace_response(a,kn,source,d)
            average=float(np.sinc(j/L)*np.mean(np.cos(2*np.pi*region['cells'][:,0]*j/L)))
            fixed=raw.real+average*data['inward_shape_response']
            data['responses'].append(dict(shift=[j,0,0],q=2*np.pi*j/(L*h),metric_fixed_cells=raw.real,shape_volume_average=average,metric_fixed_proper_volume=fixed,lapse=lapse.real,max_trace_imaginary=max(abs(raw.imag),abs(lapse.imag)),uniform_residual_relative_to_this_metric=abs(data['uniform_fixed_volume_residual'])/max(abs(fixed),1e-15)))
        data['responses'].sort(key=lambda x:x['shift'])
        (outdir/f'BALL_r{radius}.json').write_text(json.dumps(data,indent=2))
    values=[response_with_bound(region,source,matter,c) for c in [(0,0,0),(4,0,0),(4,4,4)]]
    mixture_records.append(dict(radius_cells=radius,width=6.,actual_branch_energy=matter['actual_energy'],records=values))
    (outdir/'MIXTURE_BOUNDS.json').write_text(json.dumps(dict(status='ACTUAL_NEUTRAL_MIXTURE_ENTROPY_BRACKETS',external_flag_excluded_from_regional_algebra=True,source_branch_charges_equal=True,records=mixture_records),indent=2))
    print(json.dumps(dict(completed_radius=radius,cpu_seconds=time.process_time()-start,mixture_log_Z=[v['low_rank']['log_Z'] for v in values])),flush=True)
    if time.process_time()-start>260:raise RuntimeError('Completed radius checkpoint; resource ceiling requires another bounded resume.')

centres=[(0,0,0),(2,0,0),(4,0,0),(6,0,0),(8,0,0),(4,4,4)]
solutions=[]
for wi,W in enumerate(cfg['widths']):
    train=[];held=[]
    for data in records:
        for ci,c in enumerate(centres):
            factor=np.array([np.mean(np.cos(2*np.pi*np.array(c)*j/L)) for j in [1,2,3]])
            a=np.array([r['metric_fixed_proper_volume'] for r in data['responses']])*factor
            b=np.array([r['lapse'] for r in data['responses']])*factor
            row=dict(radius=data['radius_cells'],centre=c,A=np.r_[a,b],J=data['matter'][wi]['actual_centre_records'][ci]['first_derivative'])
            (train if data['radius_cells'] in [min(cfg['radii']),max(cfg['radii'])] and ci<3 else held).append(row)
    A=np.array([r['A'] for r in train]);J=np.array([r['J'] for r in train])
    solution,res,rank,sv=np.linalg.lstsq(A,-J,rcond=1e-12)
    def report(rows):
        return [dict(radius=r['radius'],centre=r['centre'],actual_matter_entropy_tangent=r['J'],predicted_geometry_entropy_tangent=float(r['A']@solution),total_residual=float(r['A']@solution+r['J']),relative_residual=float((r['A']@solution+r['J'])/max(abs(r['J']),1e-14))) for r in rows]
    solutions.append(dict(width=W,rank=int(rank),singular_values=sv.tolist(),condition_number=float(sv[0]/sv[-1]),sigma_amplitudes=solution[:3].tolist(),lapse_amplitudes=solution[3:].tolist(),sigma_centre=float(solution[:3].sum()),lapse_centre=float(solution[3:].sum()),sigma_supremum_bound=float(np.sum(abs(solution[:3]))),lapse_supremum_bound=float(np.sum(abs(solution[3:]))),training=report(train),held_out=report(held),global_stationarity=False,finite_geometry_entropy_test=False,physical_lapse_equation=False))
out=dict(status='ACTUAL_MULTI_BALL_COFRAME_LAPSE_TANGENT_SELECTION_COMPONENT',metadata=metadata,balls=records,solutions=solutions,resume_cpu_seconds=time.process_time()-start,interrupted_prior_run=True,Poisson_or_Einstein_update_inserted=False,lapse_relation_inserted=False,rest_slice_only=True,physical_G=False)
(outdir/'RESULT.json').write_text(json.dumps(out,indent=2))
print(json.dumps(dict(status=out['status'],cpu_seconds=out['resume_cpu_seconds'],solutions=[dict(width=x['width'],rank=x['rank'],condition=x['condition_number'],sigma_centre=x['sigma_centre'],lapse_centre=x['lapse_centre'],max_held_relative=max(abs(r['relative_residual']) for r in x['held_out'])) for x in solutions]),indent=2),flush=True)
