"""Actual integer routes and certified coordinates, with no dropped mode."""
from fractions import Fraction as F
from functools import lru_cache
import json
import numpy as np
import sympy as s
from providers import ROOT, source, native

D = 190
A, B = F(1, 2), F(3)


def interval_sqrt(x, bits):
    return native.sqrt_interval(F(x), bits)


def reciprocal(x):
    if x[0] <= 0:
        raise ValueError('positive interval required')
    return 1/x[1], 1/x[0]


def midpoint(x):
    return (x[0]+x[1])/2


@lru_cache(None)
def data():
    base = ROOT/'openwave-tuning/common_dynamics/quadratic_field'
    frame = json.loads((base/'SPECTRAL_FRAME.json').read_text())
    spectral = json.loads((base/'FOCK.json').read_text())['spectral_certificate']
    Q = [[F(x, frame['denominator']) for x in row] for row in frame['frame_integers']]
    lam = list(map(F, spectral['eigenvalues_exact']))
    _, incidence, faces, _, _, _ = source.field.cell_complex.build(2)
    den=frame['denominator']; limb=1<<22
    integer_frame=np.array(frame['frame_integers'],dtype=np.int64)
    curvature=np.array(faces*faces.T,dtype=np.int64)
    integer_lambda=np.array([int(x*den) for x in lam],dtype=np.int64)
    assert all(x*den==int(x*den) for x in lam)
    assert int(np.max(np.sum(abs(curvature),axis=1)))*int(np.max(abs(integer_frame)))<1<<63
    high=integer_frame//limb; low=integer_frame-high*limb
    assert 288*max(int(np.max(abs(high))),int(np.max(abs(low))))**2<1<<63
    gram=(high.T@high).astype(object)*limb**2
    gram+=((high.T@low).astype(object)+(low.T@high).astype(object))*limb
    gram+=(low.T@low).astype(object)-den**2*np.eye(D,dtype=object)
    residual=(curvature@integer_frame).astype(object)*den-integer_frame.astype(object)*integer_lambda[None,:].astype(object)
    gram2=F(sum(int(x)**2 for x in gram.flat),den**4)
    residual2=F(sum(int(x)**2 for x in residual.flat),den**4)
    eps=F(spectral['frame_Gram_error_upper']); residual_bound=F(spectral['frame_residual_upper'])
    assert gram2<=eps**2 and residual2<=residual_bound**2
    shift=residual_bound/min(lam)
    polar=F(spectral['polar_Gram_error_upper'])
    assert polar>=eps+2*(1+eps)*shift+shift*shift
    assert F(spectral['actual_curvature_operator_error_upper'])>=residual_bound+(F(1504,17)+max(lam))*(shift+polar)
    # Exact rational null vectors and exact orthogonalization. The isotropic
    # harmonic packet and p^2 Hamiltonian are invariant under frame rotation.
    raw = incidence.col_join(faces.T).nullspace()
    assert len(raw) == 3
    orthogonal = []
    for vector in raw:
        v = vector.copy()
        for u in orthogonal:
            v -= u*(u.dot(v)/u.dot(u))
        orthogonal.append(v)
    harmonics = [[F(x) for x in v] for v in orthogonal]
    norms = [sum(x*x for x in v) for v in harmonics]
    assert all(incidence*v == s.zeros(96,1) and faces.T*v == s.zeros(208,1) for v in orthogonal)
    # The reduced Laplacian fixes the gauge and yields exact minimum-energy
    # particulars for all eight one-charge root sectors.
    lap = incidence*incidence.T
    inverse = lap[1:,1:].inv(method='DM')
    sectors = []
    for cell in range(8):
        q = s.zeros(96,1); q[0] += 1; q[12*cell] -= 1
        potential = s.zeros(96,1); potential[1:,:] = inverse*q[1:,:]
        p = incidence.T*potential
        assert incidence*p == q and faces.T*p == s.zeros(208,1)
        sectors.append({'particular':[F(x) for x in p], 'coulomb':A*F(p.dot(p))/2})
    return {'Q':Q,'lambda':lam,'spectral':spectral, 'harmonics':harmonics,
            'harmonic_norms':norms, 'sectors':sectors,
            'integer_spectral_recheck':{'Gram_squared_exact':str(gram2),'residual_squared_exact':str(residual2),
                'frozen_upper_bounds_verified':True,'curvature_norm_majorant':'1504/17, frozen DYN16 geometric certificate'}}


@lru_cache(None)
def frequency(j, bits=80):
    return interval_sqrt(A*B*data()['lambda'][j], bits)


@lru_cache(None)
def qscale(j, bits=80):
    w = frequency(j,bits)
    return interval_sqrt(A/(2*w[1]), bits)[0], interval_sqrt(A/(2*w[0]), bits)[1]


@lru_cache(None)
def pscale(j, bits=80):
    w = frequency(j,bits)
    return interval_sqrt(w[0]/(2*A), bits)[0], interval_sqrt(w[1]/(2*A), bits)[1]


@lru_cache(None)
def current(cell, axis, sign):
    destination, path, J = source.matter.route(cell,axis,sign)
    j = tuple(map(int,J))
    d = data()
    support=[e for e,value in enumerate(j) if value]
    coexact = [sum(d['Q'][e][m]*j[e] for e in support) for m in range(D)]
    harmonic = [sum(v[e]*j[e] for e in support) for v in d['harmonics']]
    return {'destination':int(destination),'path':list(map(int,path)), 'integer_current':j,
            'coexact_rational_numerators':coexact,'harmonic_rational_numerators':harmonic}


@lru_cache(None)
def coefficient(cell, axis, sign, mode, bits=80):
    j = current(cell,axis,sign)
    if mode < D:
        return native.scale(qscale(mode,bits), j['coexact_rational_numerators'][mode], bits)
    k = mode-D
    root = interval_sqrt(2*data()['harmonic_norms'][k],bits)
    return native.scale(reciprocal(root), j['harmonic_rational_numerators'][k],bits)


def edge_coefficient(edge, mode, bits=80):
    if not 0 <= edge < 288 or not 0 <= mode < 193:
        raise ValueError('actual edge and retained mode required')
    if mode < D:
        return native.scale(pscale(mode,bits),data()['Q'][edge][mode],bits)
    k = mode-D
    root = interval_sqrt(2*data()['harmonic_norms'][k],bits)
    return native.scale(reciprocal(root),data()['harmonics'][k][edge],bits)
