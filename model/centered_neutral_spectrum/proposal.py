"""Untrusted floating proposals, retained only after exact admission."""
import numpy as np
from scipy import linalg
from . import source
BITS=60
ENERGY_BITS=56

def key(s,e,k=0,kind='global'):
    return f'{kind}_{s}_{e}_{k}'

def propose():
    result={}
    for s,e,k,kind in [(*x,'global') for x in source.labels()]+[(1,0,0,'probe'),(0,1,0,'probe')]:
        b=source.block(s,e,k,kind);name=key(s,e,k,kind)
        result[name+'_columns']=b['columns']
        if b['dimension']:
            en,v=linalg.eigh(b['numerator'].astype(float)/16,b['gram'].astype(float))
        else:en=np.empty(0);v=np.empty((0,0))
        if not np.isfinite(en).all() or not np.isfinite(v).all():raise ArithmeticError('nonfinite eigensolver proposal')
        result[name+'_energies']=np.array([round(float(x)*(1<<ENERGY_BITS)) for x in en],dtype=np.int64)
        result[name+'_vectors']=np.array([[round(float(x)*(1<<BITS)) for x in row] for row in v],dtype=np.int64).reshape(v.shape)
    from sympy import Matrix, ilcm
    from math import gcd
    from functools import reduce
    from ..centered_particle_spectrum import packet as old_packet
    old=old_packet.load_witness()
    g=source.parent.reduced('ground',tuple(int(x) for x in old['ground_columns']))
    for spin,eta in ((1,0),(0,1)):
        b=source.block(spin,eta,kind='probe');name=key(spin,eta,0,'probe')
        null=Matrix(b['numerator']+32*b['gram']).nullspace()
        if len(null)!=1:raise ArithmeticError('expected one exact dark line')
        z=null[0];z*=ilcm(*[x.q for x in z]);z/=reduce(gcd,[int(x) for x in z])
        t=source.transition(spin,g,b)
        pulled=Matrix(g['gram']).inv()*Matrix(t).T*z
        denominator=int(ilcm(*[x.q for x in pulled]));numerator=pulled*denominator
        result[name+'_dark']=np.array(list(z),dtype=np.int64)
        result[name+'_pullback']=np.array(list(numerator),dtype=np.int64)
        result[name+'_pullback_denominator']=np.array([denominator],dtype=np.int64)
    return result
