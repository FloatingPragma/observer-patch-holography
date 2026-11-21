"""Propagate finite shot uncertainty of the direct trimer and local coupling instruments."""
import argparse
import gzip
import json

import numpy as np

from .experiment import HERE, SHOTS, canonical, same_receipt


B = 1500


def sample_rates(rng,counts,click_index):
    # sufficient statistics of the archived joint-word distribution, including all quiets.
    table=np.asarray(counts)
    arms,T,hand,matter=table.shape
    outcomes=np.array([(h//2,h%2,int(m==click_index)) for h in range(hand) for m in range(matter)])
    ps=table.reshape(arms,T,-1)/SHOTS
    mean=ps @ outcomes
    cov=np.einsum('atw,wi,wj->atij',ps,outcomes,outcomes)-mean[:,:,:,None]*mean[:,:,None,:]
    L=np.linalg.cholesky(cov/SHOTS+np.eye(3)*1e-22)
    normal=rng.normal(size=(B,arms,T,3))
    return mean[None]+np.einsum('atij,batj->bati',L,normal)


def times(rates,offset,a,b):
    xy=rates[:,:,:,:2].mean(axis=1)-offset[:,None,None]
    x,y=xy[:,:,0],xy[:,:,1]
    t=np.unwrap(np.arctan2(b[:,None]*x-a[:,None]*y,a[:,None]*x+b[:,None]*y),axis=1)
    return t-t[:,:1]


def slope(t,y):
    tx=t-t.mean(axis=1,keepdims=True)
    return np.sum(tx*y,axis=1)/np.sum(tx*tx,axis=1)


def ci(x):
    return np.quantile(x,[.025,.975]).tolist()


def run():
    rng=np.random.default_rng(202611213044)
    with gzip.open(HERE/'midpoint.json.gz','rt') as f:
        receipt=json.load(f)
    cal=rng.binomial(SHOTS,np.array(receipt['clock_calibration'])/SHOTS,size=(B,4))/SHOTS
    offset=cal.mean(axis=1)
    ac,bc=((cal[:,:2]-cal[:,2:])/2).T
    rates=sample_rates(rng,receipt['calibration']['joint_counts'],1)
    t=times(rates,offset,ac,bc)
    a=np.full(B,receipt['calibration']['hopping'])
    for _ in range(3):
        derivative=t*np.sin(2*a[:,None]*t)
        error=rates[:,0,:,2]-np.sin(a[:,None]*t)**2
        a+=np.sum(error*derivative,axis=1)/np.sum(derivative**2,axis=1)
    trace=(rates[:,1,:,2]-.5)+1j*(rates[:,2,:,2]-.5)
    d=abs(slope(t,np.unwrap(np.angle(trace),axis=1)))
    c=np.sqrt(a**3/(4*np.hypot(a,d/2)))
    rows=[]
    for line in receipt['lines']:
        if line['band']!=1:
            continue
        rates=sample_rates(rng,line['joint_counts'],2)
        t=times(rates,offset,ac,bc)
        y=rates[:,:,:,2]
        trace=(y[:,0]-y[:,1])+1j*(y[:,2]-y[:,3])
        omega=abs(slope(t,np.unwrap(np.angle(trace),axis=1)))
        if line['target']==0:
            m=omega
            continue
        if line['target']==.30:
            continue
        prediction=np.hypot(m,c*line['q'])
        rows.append(dict(q=line['q'], frequency_95=ci(omega), prediction_95=ci(prediction),
                         signed_fraction_error_95=ci(prediction/omega-1)))
    return dict(replicates=B, method='common binomial clock calibration, joint normal multinomial rates; fixed discrete q',
                hopping_95=ci(a),degree_95=ci(d),speed_95=ci(c),rest_95=ci(m),rows=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    r=canonical(run())
    output=HERE/'uncertainty.json'
    if args.write:
        with output.open('x') as f:
            json.dump(r,f,indent=2)
            f.write('\n')
    else:
        assert same_receipt(r,json.loads(output.read_text()))
    print(r)
