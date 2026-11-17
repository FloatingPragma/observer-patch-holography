import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.special import xlogy
from kernels import cell,s,response_kernels
from regions import build,aggregates,trace_coefficients,trace_response
from matter import prepare_matter,response

ap=argparse.ArgumentParser()
ap.add_argument('--L',type=int,default=40);ap.add_argument('--h',type=float,default=.5)
ap.add_argument('--radii',type=int,nargs='+',default=[6,7,8]);ap.add_argument('--widths',type=float,nargs='+',default=[4.,6.])
ap.add_argument('--aliases',type=int,default=1);ap.add_argument('--output',default='production')
args=ap.parse_args();outdir=Path(args.output);outdir.mkdir(exist_ok=True,parents=True)
start=time.process_time();source=cell.prepare(args.L,args.h,1.,0.,args.aliases,.5)
centres=[(0,0,0),(2,0,0),(4,0,0),(6,0,0),(8,0,0),(4,4,4)]
shifts=[(0,0,0),(1,0,0),(2,0,0),(3,0,0)]
prepared=[prepare_matter(source,W) for W in args.widths]
coeff=np.fft.fftn(source['cov'],axes=(0,1,2))
Ccell=cell.PB+np.einsum('...a,aij->...ij',coeff,cell.GAMMA)
ee=np.clip(np.linalg.eigvalsh(Ccell),0,1)
cell_bulk_density=float(np.mean(np.sum(-xlogy(ee,ee)-xlogy(1-ee,1-ee),axis=-1))/args.h**3)
del coeff,Ccell,ee
metadata=dict(configuration=vars(args),continuous_source_R=1.,mass=0.,twist=.5,observer_bulk_entropy_density=cell_bulk_density,observer_bulk_scope='The finite cell algebra traces subcell modes, so its global state need not be pure. This measured bulk entropy is not subtracted.',source_preparations=[{k:v for k,v in m.items() if not isinstance(v,np.ndarray)} for m in prepared],periodic_positive_energy_mean_obligation='Not solved by a nonzero-q fit; mean entropy responses are measured explicitly.')
(outdir/'METADATA.json').write_text(json.dumps(metadata,indent=2))
print(json.dumps(metadata),flush=True)
regions=[];records=[]
for radius in args.radii:
    region=build(source,radius);d,acc=aggregates(region,source,shifts)
    Km0=trace_coefficients(acc[0],source,'metric',d)
    Kshape=trace_coefficients(acc[0],source,'shape',d)
    data=dict(radius_cells=radius,physical_radius=radius*args.h,ell_over_R=radius*args.h,cell_count=len(region['cells']),full_observer_modes=4*len(region['cells']),sector_multiplicities=region['sizes'],entropy=region['entropy'],eigenranges=region['eigenranges'],floor_counts=region['floors'],bulk_entropy_estimate=cell_bulk_density*len(region['cells'])*args.h**3,uniform_metric_response=Km0.real,inward_shape_response=Kshape.real,uniform_fixed_volume_residual=(Km0+Kshape).real,matter=[],responses=[])
    for matter in prepared:
        values=[response(region,source,matter,c) for c in centres]
        mean=trace_response(acc[0],matter['average_kernel'],source,d)
        data['matter'].append(dict(width=matter['width'],actual_total_energy=matter['actual_energy'],actual_centre_records=values,translation_mean_first_derivative=mean.real))
    # Keep only traces/aggregates after the actual covariance entropy step.
    region['C'].clear();region['K'].clear()
    regions.append((d,acc));records.append(data)
    (outdir/f'BALL_r{radius}.json').write_text(json.dumps(data,indent=2))
    print(json.dumps(dict(radius=radius,entropy=data['entropy'],uniform_residual=data['uniform_fixed_volume_residual'],cpu_seconds=time.process_time()-start)),flush=True)
    if time.process_time()-start>250:raise RuntimeError('Resource ceiling approached after a completed ball checkpoint.')

for j,shift in enumerate(shifts[1:],1):
    km,kn=response_kernels(args.L,args.h,1.,0.,args.aliases,.5,shift)
    for data,(d,acc) in zip(records,regions):
        radius=data['radius_cells'];cells=np.array([c for c in __import__('itertools').product(range(-radius,radius+1),repeat=3) if sum(x*x for x in c)<=radius*radius])
        raw=trace_response(acc[j],km,source,d)
        lapse=trace_response(acc[j],kn,source,d)
        average=float(np.sinc(shift[0]/args.L)*np.mean(np.cos(2*np.pi*cells[:,0]*shift[0]/args.L)))
        fixed=raw.real+average*data['inward_shape_response']
        data['responses'].append(dict(shift=list(shift),q=2*np.pi*shift[0]/(args.L*args.h),metric_fixed_cells=raw.real,shape_volume_average=average,metric_fixed_proper_volume=fixed,lapse=lapse.real,max_trace_imaginary=max(abs(raw.imag),abs(lapse.imag)),uniform_residual_relative_to_this_metric=abs(data['uniform_fixed_volume_residual'])/max(abs(fixed),1e-15)))
        (outdir/f'BALL_r{radius}.json').write_text(json.dumps(data,indent=2))
    print(json.dumps(dict(completed_shift=shift,cpu_seconds=time.process_time()-start)),flush=True)
    if time.process_time()-start>280:raise RuntimeError('Resource ceiling approached after a completed Fourier response checkpoint.')

solutions=[]
for wi,W in enumerate(args.widths):
    train=[];held=[]
    for data in records:
        for ci,c in enumerate(centres):
            factor=np.array([np.mean(np.cos(2*np.pi*np.array(c)*j/args.L)) for j in [1,2,3]])
            a=np.array([r['metric_fixed_proper_volume'] for r in data['responses']])*factor
            b=np.array([r['lapse'] for r in data['responses']])*factor
            row=dict(radius=data['radius_cells'],centre=c,A=np.r_[a,b],J=data['matter'][wi]['actual_centre_records'][ci]['first_derivative'])
            (train if data['radius_cells'] in [min(args.radii),max(args.radii)] and ci<3 else held).append(row)
    A=np.array([r['A'] for r in train]);J=np.array([r['J'] for r in train])
    solution,res,rank,sv=np.linalg.lstsq(A,-J,rcond=1e-12)
    def results(rows):
        return [dict(radius=r['radius'],centre=r['centre'],actual_matter_entropy_tangent=r['J'],predicted_geometry_entropy_tangent=float(r['A']@solution),total_residual=float(r['A']@solution+r['J']),relative_residual=float((r['A']@solution+r['J'])/max(abs(r['J']),1e-14))) for r in rows]
    solutions.append(dict(width=W,rank=int(rank),singular_values=sv.tolist(),condition_number=float(sv[0]/sv[-1]),sigma_amplitudes=solution[:3].tolist(),lapse_amplitudes=solution[3:].tolist(),sigma_supremum_bound=float(np.sum(abs(solution[:3]))),lapse_supremum_bound=float(np.sum(abs(solution[3:]))),training=results(train),held_out=results(held),global_stationarity=False,finite_geometry_entropy_test=False))

out=dict(status='ACTUAL_MULTI_BALL_COFRAME_LAPSE_TANGENT_SELECTION_COMPONENT',metadata=metadata,balls=records,solutions=solutions,cpu_seconds=time.process_time()-start,Poisson_or_Einstein_update_inserted=False,lapse_relation_inserted=False,physical_G=False)
(outdir/'RESULT.json').write_text(json.dumps(out,indent=2))
print(json.dumps(dict(status=out['status'],cpu_seconds=out['cpu_seconds'],solutions=[dict(width=x['width'],rank=x['rank'],condition=x['condition_number'],sigma_bound=x['sigma_supremum_bound'],lapse_bound=x['lapse_supremum_bound'],max_held_relative=max(abs(r['relative_residual']) for r in x['held_out'])) for x in solutions]),indent=2),flush=True)
