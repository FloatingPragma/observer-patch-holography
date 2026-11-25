"""Reduced channels of the four observed loads and both complete head memories."""
from fractions import Fraction
from math import isqrt

import numpy as np

from .. import measure

prior = measure.bridge.native.prior
Q = 2**52
LIMB = 2**26


def fixed_product(x, scalar):
    """Exact componentwise floor of x*scalar/Q on two 26-bit limbs."""
    high, low = x//LIMB, x%LIMB
    a, b = divmod(scalar, LIMB)
    return high*a + (high*b + low*a + (low*b)//LIMB)//LIMB


def whole_columns(h, stations, launches, coefficients):
    dimension = h.shape[0]
    row = int(abs(h).sum(axis=1).max())
    degree = prior.DEGREE
    assert prior.UNIT//prior.BASE == Q and row<64 and dimension == 3840
    base_error = Fraction(3*(1+isqrt(dimension))*degree*(degree+1)*(degree+2), 6*prior.UNIT)
    base_error += Fraction(3*(degree+1), prior.FRAC) + Fraction(1, 10**59)
    assert 4*base_error < Fraction(1, 10**19)
    # At each Chebyshev term rounding T to Q, its coefficient to Q and
    # their product to Q costs <= (3+2+1)/Q per real component.
    column_error = base_error + Fraction(6*(degree+1)*(1+isqrt(dimension)), Q)
    assert 2*column_error < Fraction(12, 10**11)
    high = np.zeros((dimension, len(launches)), np.int64)
    low = high.copy()
    high[launches] = Q*np.eye(len(launches), dtype=np.int64)
    selected = np.zeros((2, len(stations), len(launches)), dtype=object)
    whole = np.zeros((2, dimension, len(launches)), dtype=np.int64)
    previous_high = previous_low = None
    for k in range(degree+1):
        assert np.max(abs(high)) <= 2*Q and row*(int(np.max(abs(high)))+1)<2**62 and row*prior.BASE<2**62
        sign = 1 if k%4 in (0,3) else -1
        coef = int(coefficients[k])
        assert abs(coef)<=3*prior.FRAC
        selected[k%2] += sign*coef*(high[stations].astype(object)*prior.BASE+low[stations])
        whole[k%2] += sign*fixed_product(high, coef*Q//prior.FRAC)
        # The check bounds the magnitude *after* the addition. With each
        # addend <=6Q and previous accumulation <=2^60, int64 never wraps.
        assert np.max(abs(whole[k%2]))<2**60
        if k < degree:
            nh, nl = prior.divided_product(h, high, low, 6 if k==0 else 5)
            if k:
                nl -= previous_low
                nh += (nl >> prior.LIMB)-previous_high
                nl &= prior.MASK
            previous_high,previous_low,high,low = high,low,nh,nl
    return selected//prior.UNIT,whole


def gram(re, im):
    """Return the exact coefficients of X -> sum V X V*, before division."""
    n = re.shape[1]
    r,z = (x.reshape(-1,4*n).astype(object) for x in (re,im))
    real = r.T@r+z.T@z
    skew = z.T@r
    return np.stack([x.reshape(n,4,n,4).transpose(0,2,1,3) for x in (real,skew-skew.T)])


def channels(board, route, rings, odd, coefficients):
    selected_indices = measure.bridge.beam.indices(route)
    launch_indices = measure.bridge.beam.indices(sum(rings, []))
    inside = []
    outside = 0
    for charge in (-1,1):
        h=measure.bridge.native.generator(board, route, charge, odd)
        selected,whole=whole_columns(h,selected_indices,launch_indices,coefficients)
        # The source launches a uniform classical mixture on ten orthogonal ports.
        ins=selected.reshape(2,16,10,4).transpose(0,2,1,3)
        outside_ports = np.setdiff1d(np.arange(board.n),route)
        out=whole.reshape(2,board.n,4,10,4)[:,outside_ports].transpose(0,1,3,2,4).reshape(2,-1,4,4)
        inside.append((gram(*ins),20*prior.FRAC**2))
        outside += gram(*out)
    return [*inside,(outside,20*Q**2)]
