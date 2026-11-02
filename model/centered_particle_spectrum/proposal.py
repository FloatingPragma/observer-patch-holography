"""Untrusted floating proposal; only certificate.verify admits a result."""
import numpy as np
from scipy import linalg
from . import source
from .certificate import VECTOR_BITS,ENERGY_BITS


def propose():
    result={}
    for kind in ('ground','local'):
        block=source.reduced(kind)
        e,v=linalg.eigh(block['numerator'].astype(float)/16,block['gram'].astype(float))
        if not np.isfinite(e).all() or not np.isfinite(v).all():
            raise ArithmeticError('nonfinite eigensolver proposal')
        result[kind+'_columns']=np.array(block['columns'],dtype=np.int64)
        values=[int(round(x*(1<<ENERGY_BITS))) for x in e]
        vectors=[[int(round(x*(1<<VECTOR_BITS))) for x in row] for row in v]
        if max(abs(x) for x in values)>np.iinfo(np.int64).max:
            raise OverflowError('energy proposal outside signed integer witness range')
        if max(abs(x) for row in vectors for x in row)>np.iinfo(np.int64).max:
            raise OverflowError('vector proposal outside signed integer witness range')
        result[kind+'_energies']=np.array(values,dtype=np.int64)
        result[kind+'_vectors']=np.array(vectors,dtype=np.int64)
    return result
