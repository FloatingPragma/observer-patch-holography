"""Small physical-domain and input probes; no large wavepacket rerun."""
import hashlib
import importlib.util
import json
from pathlib import Path
import warnings
import zipfile
import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('stacked_physical_review_source',HERE/'stacked_wavepacket.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def declared_input_contract(n,mass,sigma,q0,times):
    if type(n) is not int or n<4:raise ValueError('integer grid count at least4 required')
    if not np.isfinite(mass) or mass<=0:raise ValueError('this positive-mode formula requires finite positive mass')
    if not np.isfinite(sigma) or sigma<=0:raise ValueError('finite positive momentum width required')
    q=np.asarray(q0,float);t=np.asarray(times,float)
    if q.shape!=(3,) or not np.all(np.isfinite(q)):raise ValueError('three finite momentum components required')
    if t.ndim!=1 or not len(t) or not np.all(np.isfinite(t)):raise ValueError('nonempty finite one-dimensional time samples required')
    return {'amplitude_components':n**3*111,'note':'Grid/width/time choices and resource budget must be frozen; negative time is only a declared unitary reversal, not elapsed proper time.'}


def zero_links_control():
    q=np.array([.75,.45,.1]);other=q.copy();other[2]=.31
    b=m.bloch_incidence(q,False);second=m.bloch_incidence(other,False)
    if not np.array_equal(b,second):raise ValueError('removed-link model still depends on q3')
    e0,_=m.positive_mode(q,.23,False);e1,_=m.positive_mode(other,.23,False)
    active0,_=m.positive_mode(q,.23,True);active1,_=m.positive_mode(other,.23,True)
    n=4;grid=m.band_grid(n,.23,False);wavevectors=grid['q']
    qq=np.stack(np.meshgrid(wavevectors,wavevectors,wavevectors,indexing='ij'),axis=-1)
    delta=(qq-q+np.pi)%(2*np.pi)-np.pi
    amplitude=np.exp(-np.sum(delta**2,axis=-1)/(4*.2**2))*np.exp(-1j*np.sum(qq*(n//2),axis=-1))
    amplitude/=np.linalg.norm(amplitude.ravel());base=amplitude[...,None]*grid['spinor']
    marginals=[]
    for time in (0.,.37):
        psi=np.fft.ifftn(base*np.exp(-1j*time*grid['energy'])[...,None],axes=(0,1,2),norm='ortho')
        marginals.append(np.sum(abs(psi)**2,axis=(0,1,3)))
    error=float(np.max(abs(marginals[0]-marginals[1])))
    if error>1e-12 or e0!=e1 or abs(active0-active1)<.01:raise ValueError('zero-links rejecting control failed')
    return {'grid':4,'mass':.23,'times':[0.,.37],
            'zero_links_q3_energy_difference':abs(e0-e1),
            'active_links_q3_energy_difference':abs(active0-active1),
            'zero_links_z_marginal_change':error,
            'scope':'Executed finite no-vertical-coupling control; a three-coordinate plot alone creates no propagation.'}


def input_probes():
    probes=[]
    calls=[('zero_mass_formula',lambda:m.positive_mode([0,0,0],0.)),
           ('negative_mass_at_zero_momentum',lambda:m.positive_mode([0,0,0],-.23)),
           ('zero_momentum_width',lambda:m.simulate(n=4,sigma=0.,times=[0.])),
           ('empty_time_samples',lambda:m.simulate(n=4,times=[]))]
    for name,call in calls:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                output=call()
        except ValueError as exc:
            probes.append({'case':name,'status':'input_rejected','exception':type(exc).__name__,'detail':str(exc)})
            continue
        except Exception as exc:
            probes.append({'case':name,'status':'late_untyped_input_failure','exception':type(exc).__name__,'detail':str(exc)})
            continue
        try:
            if name in ('zero_mass_formula','negative_mass_at_zero_momentum'):
                finite=bool(np.isfinite(output[0]) and np.all(np.isfinite(output[1])))
            else:
                finite=bool(json.dumps(output[0],allow_nan=False))
        except ValueError:
            finite=False
        probes.append({'case':name,'status':'finite_return_without_explicit_contract' if finite else 'nonfinite_return_accepted'})
    return probes


def run():
    published=HERE/'STACKED_WAVEPACKET.json';snapshot=json.loads(published.read_text()) if published.is_file() else None
    result={'schema':'oph.stacked_wavepacket.physical_assumption_review.v1',
            'zero_links':zero_links_control(),'source_input_probes':input_probes(),
            'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'stacked_wavepacket.py')},
            'large_grid_rerun':False,'source_edited':False,'official_rows_promoted':0}
    archive=HERE/'position_distributions.npz'
    if archive.is_file():
        with zipfile.ZipFile(archive) as packed:
            fields=[name for name in packed.namelist() if name.endswith('.npy')]
            spatial_shape=None
            if 'spatial_probability.npy' in fields:
                with packed.open('spatial_probability.npy') as stream:
                    version=np.lib.format.read_magic(stream)
                    reader=np.lib.format.read_array_header_1_0 if version==(1,0) else np.lib.format.read_array_header_2_0
                    spatial_shape=reader(stream)[0]
        result['spatial_archive']={'fields':fields,'spatial_probability_shape':spatial_shape,
                                   'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
                                   'inspection':'header only; no fullgrid probability tensor replay'}
    if snapshot:
        observed=snapshot['observables'];samples=observed['samples']
        result['published_snapshot']={'sha256':hashlib.sha256(published.read_bytes()).hexdigest(),
            'centroid_ray_error':observed['max_centroid_ray_distance'],
            'maximum_edge_shell_probability':observed['maximum_edge_shell_probability'],
            'final_position_std':samples[-1]['position_std'],
            'total_centroid_displacement':(np.array(samples[-1]['centroid'])-samples[0]['centroid']).tolist(),
            'central_vs_weighted_velocity_distance':float(np.linalg.norm(np.array(observed['central_group_velocity'])-observed['weighted_group_velocity'])),
            'proper_parameter_is_reference_interpretation':True}
    return result

if __name__=='__main__':
    result=run();(HERE/'PHYSICAL_ASSUMPTION_REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))
