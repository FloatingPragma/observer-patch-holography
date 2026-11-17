from pathlib import Path
import hashlib,json
import numpy as np
from numpy.polynomial.legendre import leggauss

HERE=Path(__file__).resolve().parent
P=HERE.parents[1]/'work/20261115_codex_native_boost'
r=json.loads((P/'RECEIPT.json').read_text())

def cartesian(eps,v,n):
    q,w=leggauss(n)
    x,y,z=np.meshgrid(q,q,q,indexing='ij')
    wt=w[:,None,None]*w[None,:,None]*w[None,None,:]
    rs=x*x+y*y+z*z
    inside=rs<1
    f=np.maximum(1-rs,0)**6
    grad=np.array([-12*x*np.maximum(1-rs,0)**5,
                   -12*y*np.maximum(1-rs,0)**5,
                   -12*z*np.maximum(1-rs,0)**5])
    e=np.sqrt(1+eps*eps*rs)
    d=1+v*eps*x/e
    dg=np.array([v*eps*(1/e-eps*eps*x*x/e**3),
                 -v*eps**3*x*y/e**3,-v*eps**3*x*z/e**3])
    grad_q=grad/d-f*dg/d**2
    grad_eta=np.empty_like(grad_q)
    grad_eta[0]=grad_q[0]/d
    grad_eta[1]=grad_q[1]-v*eps*y/e*grad_eta[0]
    grad_eta[2]=grad_q[2]-v*eps*z/e*grad_eta[0]
    rest=np.sum(wt*grad[0]**2)/np.sum(wt*f*f)
    norm=np.sum(wt*f*f/d)
    widths=np.sqrt(np.array([np.sum(wt*d*g*g) for g in grad_eta])/norm/rest)
    widths[0]*=np.sqrt(1-v*v)
    return widths

comparisons=[]
for eps in (.5,.0625):
    a=cartesian(eps,.8,48);b=cartesian(eps,.8,64)
    target=next(x['width_ratios'] for x in r['rows']
                if x['velocity']==.8 and x['sigma_over_mass']==eps)
    err=max(abs(b[0]-target[0]),abs(b[1]-target[1]))
    assert err<1e-7
    assert abs(b[1]-b[2])<1e-12
    comparisons.append(dict(epsilon=eps,cartesian=b.tolist(),
        producer=target,max_error=err,order_change=float(np.max(abs(a-b)))))
out=dict(kind='same-group mathematical and independent Cartesian quadrature review',
    no_new_native_calls=True,comparisons=comparisons,
    quadrature_is_interval_certificate=False,
    proof_sha256=hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest(),
    receipt_sha256=hashlib.sha256((P/'RECEIPT.json').read_bytes()).hexdigest(),
    blocking_findings=[],G6='open')
(HERE/'REVIEW.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
