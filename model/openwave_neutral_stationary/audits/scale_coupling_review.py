"""Independent finite-current, loop and blocked-statistics audit."""
from pathlib import Path
import hashlib
import json
import itertools
import numpy as np
import scipy.linalg as la
import sympy as s
from model.openwave_scale_coupling import coupling as producer

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TARGET=ROOT/'model/openwave_scale_coupling'


def independent_green(side):
    points=list(itertools.product(range(side),repeat=4));idx={p:i for i,p in enumerate(points)}
    edges=[]
    for p in points:
        for mu in range(4):
            q=list(p);q[mu]=(q[mu]+1)%side
            edges.append((idx[p],idx[tuple(q)]))
    incidence=np.zeros((len(points),len(edges)))
    for j,(a,b) in enumerate(edges):incidence[a,j]=-1;incidence[b,j]=1
    lap=incidence@incidence.T;p0=np.ones_like(lap)/len(points)
    green=la.solve(lap+p0,np.eye(len(points)),assume_a='pos')-p0
    assert la.norm(lap@green-(np.eye(len(points))-p0),np.inf)<1e-11
    return points,idx,incidence,green


def independent_current(side,r,t,points,idx):
    current=np.zeros((len(points),4));point=[0]*4
    for mu,length,sign in ((0,r,1),(1,t,1),(0,r,-1),(1,t,-1)):
        for _ in range(length):
            if sign<0:point[mu]=(point[mu]-1)%side
            current[idx[tuple(point)],mu]+=sign
            if sign>0:point[mu]=(point[mu]+1)%side
    assert point==[0]*4
    return current


def independent_loops(links,side):
    points=list(itertools.product(range(side),repeat=4))
    idx={p:sum(p[axis]*side**axis for axis in range(4)) for p in points}
    means=[]
    for r,t in producer.RECTANGLES:
        traces=[]
        for origin in points:
            for mu in range(4):
                for nu in range(4):
                    if mu==nu:continue
                    point=list(origin);matrix=np.eye(3,dtype=complex)
                    for axis,length,sign in ((mu,r,1),(nu,t,1),(mu,r,-1),(nu,t,-1)):
                        for _ in range(length):
                            if sign<0:point[axis]=(point[axis]-1)%side
                            link=links[idx[tuple(point)],axis]
                            matrix=matrix@(link if sign>0 else link.conj().T)
                            if sign>0:point[axis]=(point[axis]+1)%side
                    assert tuple(point)==origin
                    traces.append(np.trace(matrix).real/3)
        means.append(float(np.mean(traces)))
    return np.array(means)


def independent_jackknife(samples,k,block):
    def observable(x):
        w=x.mean(axis=0)
        return -np.array([np.log(w[0]),np.log(w[2])+np.log(w[0])-2*np.log(w[1]),
                          np.log(w[4])+np.log(w[2])-2*np.log(w[3])])/k
    nb=len(samples)//block;values=[]
    for i in range(nb):
        values.append(observable(np.concatenate((samples[:i*block],samples[(i+1)*block:]))))
    values=np.array(values);center=values.mean(axis=0)
    covariance=sum(np.outer(v-center,v-center) for v in values)*(nb-1)/nb
    return observable(samples),covariance


def run():
    generators=[]
    for a,b in ((0,1),(0,2),(1,2)):
        x=s.zeros(3);x[a,b]=x[b,a]=s.Rational(1,2);generators.append(x)
        y=s.zeros(3);y[a,b]=-s.I/2;y[b,a]=s.I/2;generators.append(y)
    generators.extend((s.diag(1,-1,0)/2,s.diag(1,1,-2)/(2*s.sqrt(3))))
    assert s.Matrix(8,8,lambda i,j:s.trace(generators[i]*generators[j]))==s.eye(8)/2
    assert sum((x*x for x in generators),s.zeros(3))==s.Rational(4,3)*s.eye(3)
    assert producer.fundamental_casimir()==s.Rational(4,3)
    tree=[]
    for side in (3,4):
        points,idx,incidence,green=independent_green(side)
        for r,t in producer.RECTANGLES:
            current=independent_current(side,r,t,points,idx)
            assert np.max(abs(incidence@current.reshape(-1)))==0
            value=(2/3)*sum(current[:,mu]@green@current[:,mu] for mu in range(4))
            expected=producer.tree_loop(side,r,t)
            assert abs(value-expected)<1e-12
            tree.append({'side':side,'rectangle':[r,t],'independent_value':float(value),
                         'absolute_difference':float(abs(value-expected))})
        assert abs(producer.tree_loop(side,1,1)-(1-side**-4)/3)<1e-14
    rng=np.random.default_rng(45177);side=3
    # Fresh noncommuting SU(3) links from QR; no production ensemble is reused.
    links=[]
    for _ in range(side**4*4):
        a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        q,r=la.qr(a);q=q@np.diag(np.diag(r)/abs(np.diag(r)))
        q[:,0]/=la.det(q);links.append(q)
    links=np.array(links).reshape(side**4,4,3,3)
    seam_error=max(abs((1-(np.trace(u).real+2)/5)-(3/5)*(1-np.trace(u).real/3))
                   for u in links.reshape(-1,3,3))
    assert seam_error<3e-16
    up,_=producer.K.neighbours(*(side,)*4)
    independent=independent_loops(links,side);direct=producer.loops(links,up)
    assert np.max(abs(independent-direct))<1e-14
    rows=[];combined={};betas=[]
    for path in sorted((TARGET/'data').glob('*.npz')):
        with np.load(path,allow_pickle=False) as z:
            samples=z['loops'];n=int(z['side']);betas.append(float(z['beta']))
            assert samples.shape==(400,5) and np.all(np.isfinite(samples))
            assert np.max(abs(samples[:,0]-z['plaquette']))<2e-13
            combined.setdefault(n,[]).append(samples.copy())
            k=producer.tree_coefficients(n);maximum=0.
            for block in (5,10,20,40,80):
                value,cov=independent_jackknife(samples,k,block)
                actual=producer.jackknife(samples,n,block)
                maximum=max(maximum,float(np.max(abs(value-actual['g_squared']))),
                            float(np.max(abs(cov-actual['covariance']))))
                assert maximum<1e-12
            rows.append({'raw_file':path.name,'side':n,'start':str(z['start']),
                         'independent_covariance_maximum_difference':maximum})
    assert betas==[6.]*4 and set(combined)=={6,8}
    volumes=[]
    for n,chains in sorted(combined.items()):
        value,cov=independent_jackknife(np.concatenate(chains),producer.tree_coefficients(n),20)
        delta=np.diff(value);transform=np.array([[-1,1,0],[0,-1,1]])
        err=np.sqrt(np.diag(transform@cov@transform.T))
        assert np.all(delta>5*err)
        volumes.append({'side':n,'independent_g_squared':value.tolist(),
                        'scale_differences':delta.tolist(),'blocked_errors':err.tolist()})
    rejected=[]
    for bad in ([.6,-.2,.1,-.05,.03],[.6,.2,0.,.05,.03],[.6,.2,.1,float('nan'),.03]):
        try:producer.creutz(bad)
        except ValueError:rejected.append(True)
        else:rejected.append(False)
    assert all(rejected),'nonpositive/nonfinite individual means must fail closed'
    paths=sorted(set(TARGET.glob('*.py'))|set(TARGET.glob('*.md'))|{TARGET/'results.json'}|set((TARGET/'data').glob('*.npz')))
    return {'schema':'oph.scale_coupling.independent_review.v1',
        'independent_fundamental_generator_trace_checks':64,
        'independent_fundamental_Casimir':'4/3','independent_seam_cost_maximum_difference':float(seam_error),
        'tree_current_green_controls':tree,'fresh_noncommuting_loop_maximum_difference':float(np.max(abs(independent-direct))),
        'fresh_closed_loop_traces':3**4*12*5,'raw_records_retained':1600,
        'independent_block_statistics':rows,'combined_volumes':volumes,
        'all_four_bare_betas_equal_six':True,'nonpositive_or_nonfinite_means_rejected':rejected,
        'finite_volume_long_loop_shift':float(volumes[0]['independent_g_squared'][2]-volumes[1]['independent_g_squared'][2]),
        'scope':'Finite declared pure-SU3 ensemble scale trend with Gaussian finite-volume normalization. No rigorous mixing, continuum beta function, physical GeV scale or PW119 particle-source equivalence follows.',
        'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'verdict':'PASS'}


if __name__=='__main__':
    result=run();(HERE/'scale_coupling_review.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['verdict'])
