"""Larger-volume coupling, correlated changes and autocorrelation diagnostics."""
from pathlib import Path
import numpy as np
from model.openwave_scale_coupling.coupling import jackknife,tree_coefficients,gauge_control
HERE=Path(__file__).resolve().parent

def influence(series,side):
    mean=series.mean(axis=0);k=tree_coefficients(side)
    D=np.zeros((3,5));D[0,0]=-1/mean[0]
    D[1,[0,1,2]]=[-1/mean[0],2/mean[1],-1/mean[2]]
    D[2,[2,3,4]]=[-1/mean[2],2/mean[3],-1/mean[4]]
    D/=k[:,None]
    return (series-mean)@D.T

def autocorrelation(x,maxlag=100):
    x=np.asarray(x,float);x=x-x.mean();n=len(x)
    if np.dot(x,x)==0:return {'tau_integrated':.5,'effective_records':n,'lags_used':0}
    gamma=np.array([np.dot(x[:n-lag],x[lag:])/n for lag in range(min(maxlag,n-1)+1)])
    pairs=[]
    for j in range(0,len(gamma)-1,2):
        value=(gamma[j]+gamma[j+1])/gamma[0]
        if value<=0:break
        pairs.append(min(value,pairs[-1])if pairs else value)
    tau=max(.5,-.5+sum(pairs))
    return {'tau_integrated':float(tau),'effective_records':float(n/(2*tau)),
        'lags_used':2*len(pairs),'method':'initial-positive monotone covariance-pair diagnostic; finite histories do not certify mixing'}

def run():
    files=sorted((HERE.parent/'openwave_scale_coupling/data').glob('*.npz'))+sorted((HERE/'data').glob('*.npz'))
    if len(files)!=6:raise ValueError('all six original and larger-box chains required')
    chains=[];by={}
    for path in files:
        with np.load(path,allow_pickle=False)as z:
            side=int(z['side']);x=z['loops'];assert x.shape==(400,5)
            assert float(z['beta'])==6 and np.isfinite(x).all()
            assert np.max(abs(x[:,0]-z['plaquette']))<2e-13
            err=gauge_control(z['final_links'],side);assert err<3e-12
            values=influence(x,side)
            chains.append({'path':str(path.relative_to(HERE.parents[1])),'side':side,'seed':int(z['seed']),'start':str(z['start']),
                'warmup':int(z['warmup']),'stride':int(z['stride']),
                'blocks':[jackknife(x,side,b)for b in (5,10,20,40,80)],
                'influence_autocorrelation':[autocorrelation(values[:,j])for j in range(3)],
                'first_half':jackknife(x[:200],side,20),'second_half':jackknife(x[200:],side,20),'final_gauge_error':err})
            by.setdefault(side,[]).append(x.copy())
    assert set(by)=={6,8,10}and all(len(v)==2 for v in by.values())
    combined=[];contrasts=[]
    for side,series in sorted(by.items()):
        # Fixed equal chain weights; block means never cross a chain seam.
        a,b=[jackknife(x,side,20)for x in series];mean=jackknife(np.concatenate(series),side,20)
        diff=np.array(a['g_squared'])-b['g_squared'];var=np.array(a['covariance'])+b['covariance']
        combined.append({'side':side,'mean':mean,'hot_cold_difference':diff.tolist(),'hot_cold_standard_errors':np.sqrt(np.diag(var)).tolist(),
            'scope':'independent starts and common-block jackknife; neither prove equilibration nor bound systematic errors'})
    for left,right in zip(combined,combined[1:]):
        a=left['mean'];b=right['mean'];diff=np.array(b['g_squared'])-a['g_squared'];cov=np.array(a['covariance'])+b['covariance']
        contrasts.append({'sides':[left['side'],right['side']],'larger_minus_smaller_g_squared':diff.tolist(),
            'independent_box_standard_errors':np.sqrt(np.diag(cov)).tolist()})
    return {'schema':'oph.native_colour_larger_volume_controls.v1','bare_beta':6,'bare_g_squared':1,
        'same_native_seam_action_and_loop_scheme':True,'new_raw_records':800,'total_retained_records':2400,
        'scales_lattice_lengths':[.5,1.5,2.5],'all_chains':chains,'same_bare_coupling_by_volume':combined,
        'adjacent_volume_contrasts':contrasts,'model_has_identified_particle_core':False,
        'physical_GeV_scale_or_continuum_beta_function':False,
        'interpretation':'Finite-scale coupling with explicit larger-volume, hot/cold, temporal-correlation and split-history diagnostics; no rigorous Monte Carlo mixing or continuum extrapolation.',
        'verdict':'PASS_LARGER_VOLUME_NATIVE_COUPLING_CONTROLS'}
