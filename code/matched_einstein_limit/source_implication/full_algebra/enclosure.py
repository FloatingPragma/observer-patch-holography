"""Outward-rounded second derivatives of the four-station area and volume."""

from itertools import combinations

import mpmath as mp
import numpy as np

from .station import ORDER

iv = mp.iv


class Jet:
    def __init__(self, value, grad, hess):
        self.v, self.g, self.h = value, grad, hess

    def __add__(self, other):
        if isinstance(other, Jet):
            return Jet(self.v+other.v, self.g+other.g, self.h+other.h)
        return Jet(self.v+other, self.g, self.h)

    __radd__ = __add__

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, Jet):
            cross = np.outer(self.g, other.g)
            return Jet(self.v*other.v, self.g*other.v+other.g*self.v,
                       self.h*other.v+other.h*self.v+cross+cross.T)
        return Jet(self.v*other, self.g*other, self.h*other)

    __rmul__ = __mul__

    def unary(self, result, first, second):
        return Jet(result, self.g*first, self.h*first+np.outer(self.g,self.g)*second)

    def log(self):
        return self.unary(iv.ln(self.v), 1/self.v, -1/(self.v*self.v))

    def exp(self):
        e = iv.exp(self.v)
        return self.unary(e,e,e)

    def sqrt(self):
        if self.v.a <= 0:
            raise ValueError("interval contains a degenerate simplex")
        v = iv.sqrt(self.v)
        return self.unary(v,1/(2*v),-1/(4*self.v*v))


def geometry_jet(ranges):
    dims = len(ranges)
    zero = iv.mpf(0)
    g0 = np.array([zero] * dims, dtype=object)
    h0 = np.full((dims,dims), zero, dtype=object)
    logs = np.zeros((4,4), dtype=object)
    for i, (x,a,b) in enumerate(zip(ranges,*ORDER)):
        g = g0.copy(); g[i] = iv.mpf(1)
        logs[a,b] = Jet(x,g,h0).log()
    clocks = [sum(logs[a,b]-logs[b,a] for b in range(4))*.25 for a in range(4)]
    d = np.zeros((4,4),dtype=object)
    for a,b in combinations(range(4),2):
        d[a,b] = d[b,a] = (logs[a,b]+logs[b,a]-clocks[a]-clocks[b]).exp()
    def gram(a,b,c):
        return (d[a,b]+d[a,c]-d[b,c])*.5
    area = sum((d[a,b]*d[a,c]-gram(a,b,c)*gram(a,b,c)).sqrt()*.5
               for a,b,c in combinations(range(4),3))
    a,b,c = (d[0,i] for i in (1,2,3))
    p,q,r = (gram(0,i,j) for i,j in ((1,2),(1,3),(2,3)))
    det = a*b*c + 2*p*q*r - a*r*r - b*q*q - c*p*p
    volume = det.sqrt()*(iv.mpf(1)/6)
    return area,volume
