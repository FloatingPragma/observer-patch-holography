"""Exact signed integer Gram products using bounded10-bit limbs."""
import numpy as np

def exact_product(V,W):
    V=np.asarray(V,dtype=np.int64);W=np.asarray(W,dtype=np.int64);base=1<<10
    def limbs(A):
        sign=np.where(A<0,-1,1);B=np.abs(A);out=[]
        while np.any(B):out.append((B%base)*sign);B=B//base
        return out
    vl,wl=limbs(V),limbs(W);out=np.zeros((V.shape[1],W.shape[1]),dtype=object)
    assert V.shape[0]*base*base<(1<<62)
    for i,A in enumerate(vl):
        for j,B in enumerate(wl):out+=(A.T@B).astype(object)*(base**(i+j))
    return out

