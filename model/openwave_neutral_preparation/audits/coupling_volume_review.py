"""Independent raw-loop, finite-Gaussian and correlated-estimator audit.

This reviewer does not import the producer's coupling or analysis routines.
It reconstructs closed paths, the Fourier current norm, common deletions and
influence covariances directly from the retained native configurations.
"""
from functools import lru_cache
from pathlib import Path
import hashlib,json
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
PKG=ROOT/'model/openwave_coupling_volume'
RECT=((1,1),(2,1),(2,2),(3,2),(3,3))

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

@lru_cache(None)
def geometry(side):
    coordinates=np.array(np.unravel_index(np.arange(side**4),(side,)*4)).T
    angle=2*np.pi*coordinates/side
    denominator=4*np.sum(np.sin(angle/2)**2,axis=1)
    values={}
    for r in range(4):
        for t in range(4):
            if not r or not t:values[r,t]=0.;continue
            # Oriented currents of a rectangle, with the zero mode omitted.
            ex=np.exp(1j*angle[:,0]);ey=np.exp(1j*angle[:,1])
            x_current=np.sum(np.stack([ex**j for j in range(r)]),axis=0)*(1-ey**t)
            y_current=np.sum(np.stack([ey**j for j in range(t)]),axis=0)*(ex**r-1)
            values[r,t]=2/3*np.sum((abs(x_current[1:])**2+abs(y_current[1:])**2)/denominator[1:])/side**4
    return np.array([values[r,r]+values[r-1,r-1]-2*values[r,r-1]for r in range(1,4)])

def transform(means,side):
    means=np.asarray(means)
    if np.any(means<=0)or not np.isfinite(means).all():raise ValueError('nonpositive or invalid loop means')
    log=np.log(means)
    return np.stack([-log[...,0],-log[...,2]-log[...,0]+2*log[...,1],-log[...,4]-log[...,2]+2*log[...,3]],axis=-1)/geometry(side)

def independent_blocks(data,side,block):
    n=len(data);count=n//block;assert n%block==0 and count>=4
    total=data.sum(axis=0)
    omitted=np.array([(total-data[i:i+block].sum(axis=0))/(n-block)for i in range(0,n,block)])
    estimates=transform(omitted,side);center=estimates.mean(axis=0)
    cov=sum(np.outer(v-center,v-center)for v in estimates)*(count-1)/count
    return transform(total/n,side),cov

def check_estimator(row,data,side):
    g,cov=independent_blocks(data,side,row['block_samples'])
    difference=np.max(abs(g-np.array(row['g_squared'])))
    error=np.max(abs(cov-np.array(row['covariance'])))
    assert difference<2e-12 and error<2e-12,(difference,error)
    assert np.max(abs(np.sqrt(np.diag(cov))-row['standard_errors']))<2e-11
    assert np.max(abs(geometry(side)-row['tree_geometry_coefficients']))<2e-14
    return max(float(difference),float(error))

def tau_from_full_correlation(x):
    n=len(x);x=x-x.mean();positive=np.correlate(x,x,mode='full')[n-1:n+100]/n
    if positive[0]==0:return .5,0
    previous=float('inf');summed=0.;pairs=0
    for i in range(0,len(positive)-1,2):
        pair=(positive[i]+positive[i+1])/positive[0]
        if pair<=0:break
        previous=min(previous,pair);summed+=previous;pairs+=1
    return max(.5,summed-.5),2*pairs

def paths(U,side):
    indices=np.arange(side**4);powers=side**np.arange(4)
    xyz=(indices[:,None]//powers)%side
    def advance(position,axis,direction):
        digit=(position//powers[axis])%side
        return position+(((digit+direction)%side)-digit)*powers[axis]
    out=[]
    for r,t in RECT:
        trace=0.
        for a in range(4):
            for b in range(4):
                if a==b:continue
                position=indices.copy();transport=np.broadcast_to(np.eye(3,dtype=complex),(side**4,3,3)).copy()
                for axis,sign,count in ((a,1,r),(b,1,t),(a,-1,r),(b,-1,t)):
                    for _ in range(count):
                        if sign==1:
                            transport=transport@U[position,axis];position=advance(position,axis,1)
                        else:
                            position=advance(position,axis,-1);transport=transport@U[position,axis].conj().swapaxes(1,2)
                assert np.array_equal(position,indices)
                trace+=np.trace(transport,axis1=1,axis2=2).real.mean()/3
        out.append(trace/12)
    return np.array(out)

def gauge_check(U,side):
    # Site-dependent diagonal SU(3) gauges, independently generated. This
    # supplements the producer's full random-SU(3) gauge check.
    n=side**4;site=np.arange(n);powers=side**np.arange(4)
    rng=np.random.default_rng(6161+side);angles=rng.uniform(-np.pi,np.pi,(n,2))
    diagonal=np.exp(1j*np.column_stack((angles,-angles.sum(axis=1))))
    changed=np.empty_like(U)
    for axis in range(4):
        digit=(site//powers[axis])%side
        end=site+(((digit+1)%side)-digit)*powers[axis]
        changed[:,axis]=diagonal[:,:,None]*U[:,axis]*diagonal[end,None,:].conj()
    baseline=paths(U,side);after=paths(changed,side)
    return baseline,float(np.max(abs(after-baseline)))

def main():
    result=json.loads((PKG/'results.json').read_text());pins=result['source_sha256']
    assert all(sha(ROOT/p)==value for p,value in pins.items())
    assert result['total_retained_records']==2400 and result['new_raw_records']==800
    assert result['bare_beta']==6 and result['bare_g_squared']==1
    assert not result['physical_GeV_scale_or_continuum_beta_function'] and not result['model_has_identified_particle_core']
    expected={(6,'cold'):(3101,200),(6,'hot'):(3102,300),(8,'cold'):(3103,300),(8,'hot'):(3104,400),(10,'cold'):(4101,400),(10,'hot'):(4102,500)}
    rows=[];seen=set();by={};maximum=0.
    for row in result['all_chains']:
        p=ROOT/row['path']
        with np.load(p,allow_pickle=False)as z:
            side=int(z['side']);start=str(z['start']);key=(side,start);assert key not in seen;seen.add(key)
            assert (int(z['seed']),int(z['warmup']))==expected[key]
            assert int(z['samples'])==400 and int(z['stride'])==2 and int(z['overrelaxations'])==2 and float(z['beta'])==6
            data=z['loops'];assert data.shape==(400,5)and np.isfinite(data).all()
            assert np.max(abs(z['plaquette']-data[:,0]))<2e-13
            assert row['side']==side and row['start']==start and row['seed']==int(z['seed'])and row['warmup']==int(z['warmup'])
            for b in row['blocks']:maximum=max(maximum,check_estimator(b,data,side))
            maximum=max(maximum,check_estimator(row['first_half'],data[:200],side),check_estimator(row['second_half'],data[200:],side))
            # Complex-step differentiation avoids the producer's analytic
            # influence formula and finite-difference cancellation.
            mean=data.mean(axis=0);h=1e-25;gradient=np.empty((3,5))
            for i in range(5):
                probe=mean.astype(complex);probe[i]+=1j*h
                gradient[:,i]=transform(probe,side).imag/h
            influence=(data-mean)@gradient.T
            for j,item in enumerate(row['influence_autocorrelation']):
                tau,lags=tau_from_full_correlation(influence[:,j]);assert abs(tau-item['tau_integrated'])<2e-11
                assert lags==item['lags_used']and abs(400/(2*tau)-item['effective_records'])<2e-9
            U=z['final_links'];assert U.shape==(side**4,4,3,3)
            unitary_error=float(np.max(abs(U.conj().swapaxes(-1,-2)@U-np.eye(3))))
            determinant_error=float(np.max(abs(np.linalg.det(U)-1)))
            assert unitary_error<2e-12 and determinant_error<2e-12
            baseline,gauge_error=gauge_check(U,side)
            assert gauge_error<2e-12 and np.max(abs(baseline-data[-1]))<2e-12
            rows.append({'side':side,'start':start,'records':len(data),'blocks_checked':[b['block_samples']for b in row['blocks']],
                'final_path_record_error':float(np.max(abs(baseline-data[-1]))),'independent_diagonal_gauge_error':gauge_error,
                'unitarity_error':unitary_error,'determinant_error':determinant_error})
            by.setdefault(side,[]).append(data.copy())
    assert seen==set(expected)and len(result['all_chains'])==6
    combined=[]
    for row in result['same_bare_coupling_by_volume']:
        side=row['side'];a,b=by[side];maximum=max(maximum,check_estimator(row['mean'],np.concatenate((a,b)),side))
        ga,ca=independent_blocks(a,side,20);gb,cb=independent_blocks(b,side,20)
        assert np.max(abs(ga-gb-row['hot_cold_difference']))<2e-12
        assert np.max(abs(np.sqrt(np.diag(ca+cb))-row['hot_cold_standard_errors']))<2e-11
        combined.append(independent_blocks(np.concatenate((a,b)),side,20))
    for row,(a,ca),(b,cb)in zip(result['adjacent_volume_contrasts'],combined,combined[1:]):
        assert np.max(abs(b-a-row['larger_minus_smaller_g_squared']))<2e-12
        assert np.max(abs(np.sqrt(np.diag(ca+cb))-row['independent_box_standard_errors']))<2e-11
    refs={**pins,'model/openwave_coupling_volume/results.json':sha(PKG/'results.json'),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}
    review={'schema':'oph.independent_larger_volume_coupling_review.v1','verdict':'PASS','input_sha256':refs,
        'producer_transitive_pins_checked':len(pins),'chains':rows,'retained_records':2400,
        'maximum_estimator_difference':maximum,'independent_methods':['closed oriented paths from links','site-dependent SU3 diagonal gauge','finite Fourier current norm','common-block leave-out estimates and covariance','complex-step influence gradients','full discrete autocorrelation','hot/cold and adjacent-volume covariance sums'],
        'scientific_boundary':'Same bare action and three finite volumes. Statistical errors are finite-history diagnostics, not mixing certificates; no continuum running, particle-core matching or physical GeV scale.'}
    assert all(sha(ROOT/p)==value for p,value in refs.items())
    Path(__file__).with_suffix('.json').write_text(json.dumps(review,indent=2,sort_keys=True)+'\n')
    print('PASS independent larger-volume coupling review')
if __name__=='__main__':main()
