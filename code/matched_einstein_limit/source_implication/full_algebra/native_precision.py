"""Fixed-point kernels for both parity sectors of the previously certified radar."""
from fractions import Fraction
from math import isqrt

import numpy as np
from scipy import sparse

from ..isometric_collar import certify as prior
from .station import iso


def generator(board, route, charge, odd):
    h = prior.integer_flight(board, route, charge)
    if odd:
        delta = 2*(iso.ranging.parity_blocks(charge,8,-1)-iso.ranging.parity_blocks(charge,8,1))
        np.testing.assert_array_equal(delta,delta.astype('int64'))
        delta = delta.astype('int64')
        ends = sorted(route[:2])
        for i in range(2):
            for j in range(2):
                co = np.kron(delta[i+1,j+1]-int(i==j)*delta[0,0],np.eye(2,dtype='int64'))
                dyad = sparse.csr_matrix(([1],([ends[i]],[ends[j]])),shape=(board.n,board.n))
                h += sparse.kron(dyad,co,format='csr')
    return h


def kernel_panels(h, stations, co):
    d = h.shape[0]
    row = int(abs(h).sum(axis=1).max())
    assert row < 64
    de = prior.DEGREE
    epsilon = Fraction(3*(1+isqrt(d))*de*(de+1)*(de+2),6*prior.UNIT)
    epsilon += Fraction(3*(de+1),prior.FRAC) + Fraction(1,10**59)
    assert 4*epsilon < Fraction(1,10**19)
    size = len(stations)
    if size % 4:
        raise ValueError("each enrolled station requires four coin columns")
    ports = size // 4
    shape = (d,size)
    high = np.zeros(shape, dtype=np.int64)
    low = high.copy()
    high[stations] = (prior.UNIT // prior.BASE)*np.eye(size,dtype=np.int64)
    out = np.zeros((2,12,size,size),dtype=object)
    prev_hi = prev_lo = None
    for k in range(de+1):
        assert row*(int(np.max(np.abs(high)))+1) < 2**62 and row*prior.BASE < 2**62
        val = high[stations].astype(object)*prior.BASE + low[stations]
        out[k%2] += (1 if k%4 in (0,3) else -1)*co[:,k,None,None]*val[None,:,:]
        if k < de:
            nh,nl = prior.divided_product(h,high,low,6 if k==0 else 5)
            if k:
                nl -= prev_lo
                nh += (nl >> prior.LIMB)-prev_hi
                nl &= prior.MASK
            prev_hi,prev_lo,high,low = high,low,nh,nl
    return (out//prior.UNIT).reshape(2,12,ports,4,ports,4).transpose(0,1,2,4,3,5)


def exact_effects():
    board = iso.pairing.formation.federation(1)
    route = iso.pairing.formation.region(board)
    stations = np.array([4*p+i for p in route for i in range(4)])
    coeff = prior.weights()
    pairs = []
    for odd in (False,True):
        total = 0
        for sign in (-1,1):
            u = kernel_panels(generator(board,route,sign,odd),stations,coeff)
            total = total + np.array(prior.direction_effects(u,12),dtype=object)
        pairs.append(total)
    error = prior.error_effect()
    assert error < Fraction(1,10**14)
    return np.stack(pairs,axis=1),8*prior.FRAC,error


def formed_state():
    """Normalized approximation to the even freely formed two-coin source."""
    board = iso.pairing.formation.federation(1)
    route = iso.pairing.formation.region(board)
    _, vec, error = prior.certified_kernels(board,route,1,prior.weights())
    assert error < 1e-20
    a,b = (prior.sp.Matrix(k.reshape(-1,4)) for k in vec)
    rho = a.T*a + b.T*b + prior.sp.I*(b.T*a-a.T*b)
    rho /= rho.trace()
    return rho
