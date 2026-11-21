"""Native 99-state closure, actual path addresses, and source compilation."""
from itertools import combinations

import numpy as np
from scipy.linalg import eigh_tridiagonal

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30_arrows.experiment import local_check
from model.executable_model import federation


K = candidate.KAPPA


def preparation():
    car = federation(2)
    gateways = {i: {} for i in range(car.cells)}
    for a, b in car.external:
        gateways[a // 12][b // 12] = a
        gateways[b // 12][a // 12] = b
    parent = {0: None}
    longest = []

    def visit(c, route):
        nonlocal longest
        if len(route) > len(longest):
            longest = route.copy()
        for other in sorted(gateways[c]):
            if other not in parent:
                parent[other] = c
                visit(other, route + [other])

    visit(0, [0])
    sites = []
    for i, c in enumerate(longest):
        entry = [gateways[c][longest[i-1]]] if i else sorted(gateways[c].values())
        exit_ = [gateways[c][longest[i+1]]] if i+1 < len(longest) else sorted(gateways[c].values())
        blocks = [(a, b, d, e) for a in entry for e in exit_ if a != e
                  for b in car.neighbours[a] if car.degree[b] == 5
                  for d in car.neighbours[e] if car.degree[d] == 5 and d in car.neighbours[b]]
        sites.extend(min(blocks, key=lambda x: (x[0], x[3], x[1], x[2])))
    sites = np.array(sites)
    assert len(set(sites)) == len(sites)
    assert all(b in car.neighbours[a] for a, b in zip(sites[:-1], sites[1:]))
    assert np.all(car.degree[sites].reshape(-1, 4) == [6, 5, 5, 6])
    clock_edges = [e for e in car.seams if e[0] // 12 == e[1] // 12 == longest[0]
                   and not set(e).intersection(sites)]
    clock = next((a,b) for a,b in combinations(clock_edges,2) if not set(a).intersection(b))
    # All terminal bright-coin projectors commute. Even at a gateway the incident
    # head count minus its calibrated launch value gives that port's occupancy change.
    heads = np.zeros(len(sites), dtype=int)
    heads[:-1] = 1
    launch_heads = heads.copy()
    for address in range(len(sites)):
        reads = heads-launch_heads
        reads[0] += 1
        assert reads.sum() == 1 and np.count_nonzero(reads) == 1 and reads[address] == 1
        if address+1<len(sites):
            heads[address] -= 1
            heads[address+1] += 1
    return car, longest, sites, clock


def column_checks():
    law = candidate.native_law
    err, readerr = local_check()
    # Unoccupied collars are prepared by W|0,0>. Check all possible touching load columns.
    zero = np.zeros((law.DIM,3))
    for c in (0,1):
        idx = (1 + law.LABELS.index((c,0))) * law.N_SRC
        for j, loads in enumerate(((0,0),(1,0),(0,1))):
            zero[idx + law.SRC.index(loads),j] = 1/np.sqrt(2)
    h0,h1 = law.local_term_parts()
    dark = float(np.max(np.abs((h0 + K*h1) @ zero)))
    blank = np.eye(law.DIM)[:,[law.SRC.index((0,0)),law.SRC.index((1,0))]]
    ref = float(np.max(abs(h1 @ blank - blank @ np.diag([0,1]))))
    assert max(err,readerr,dark,ref)<1e-11
    return dict(active=err, occupation_read=readerr, dark=dark, degree_reference=ref)


def matrices(car, sites, degree=True):
    diag = np.full(len(sites), -K)
    diag[[0,-1]] = -K/2
    if degree:
        diag += K*(6-car.degree[sites])
    off = np.full(len(sites)-1, -K/2)
    return diag,off


def apply(diag,off,v):
    out = diag*v
    out[:-1] += off*v[1:]
    out[1:] += off*v[:-1]
    return out


def bloch(q, degree=True):
    block = np.diag(np.full(4,-float(K))).astype(complex)
    if degree:
        block[1,1] += K
        block[2,2] += K
    for j in range(3):
        block[j,j+1] = block[j+1,j] = K/2
    block[0,3] = K/2*np.exp(-1j*q)
    block[3,0] = block[0,3].conjugate()
    return np.linalg.eigh(block)


def standing(M, q, band, diag, off):
    E,u = bloch(q)
    n = np.arange(M)[:,None]
    vector = np.exp(1j*q*n)*u[:,band] + np.exp(-1j*q*(1+n))*u[::-1,band]
    vector = vector.ravel()
    vector *= np.exp(-1j*np.angle(vector[np.argmax(abs(vector))]))
    assert np.linalg.norm(vector.imag) < 1e-8 and np.linalg.norm(vector.real)>1e-4
    vector = vector.real * (-1)**np.arange(4*M)
    vector /= np.linalg.norm(vector)
    residual = float(np.linalg.norm(apply(diag,off,vector) - E[band]*vector))
    assert residual < 1e-9, residual
    return vector,residual


def packet(M, q, width=0.045):
    n = np.arange(M)[:,None]
    wave = np.zeros((M,4), complex)
    points = np.arange(-8*M,8*M)*np.pi/(8*M)
    for k in points[np.abs(points-q)<8*width]:
        _,u = bloch(k)
        polarization = u[:,1]*np.exp(-1j*np.angle(u[0,1]))
        wave += np.exp(-((k-q)/(2*width))**2 + 1j*k*(n-(M-1)/2))*polarization
    wave = wave.ravel() * (-1)**np.arange(4*M)
    return wave/np.linalg.norm(wave)


def compile_source(target):
    """Apply serial two-port repair mixers, then local degree-phase pulses; release."""
    psi = np.zeros(len(target),complex)
    psi[0] = 1
    tails = np.cumsum(abs(target[::-1])**2)[::-1]
    for j in range(len(psi)-1):
        if tails[j] < 1e-29:
            break
        theta = 2*np.arctan2(np.sqrt(tails[j+1]),abs(target[j]))
        shared = np.expm1(1j*theta)*(psi[j]+psi[j+1])/2
        psi[j:j+2] += shared
    phases = np.angle(target*np.conj(psi))
    swaps = 0
    for j,phase in enumerate(phases):
        if j%4 in (0,3):
            neigh = j+1 if j%4 == 0 else j-1
            psi[[neigh,j]] = -psi[[j,neigh]]
            psi[neigh] *= np.exp(1j*phase)
            psi[[neigh,j]] = -psi[[j,neigh]]
            swaps += 2
        else:
            psi[j] *= np.exp(1j*phase)
    error = float(np.linalg.norm(psi-target))
    assert error<1e-9, error
    return psi, error


class FreeChannel:
    def __init__(self, diag, off):
        self.e,self.v = eigh_tridiagonal(diag,off)
        self.max_source_error = 0.

    def launch(self, target):
        psi,err = compile_source(target)
        self.max_source_error = max(self.max_source_error,err)
        return psi

    def clicks(self,psi,address,times):
        weight = self.v[address]*(self.v.T @ psi)
        return abs(weight @ np.exp(-1j*self.e[:,None]*times))**2

    def profile(self,psi,time):
        return abs(self.v @ ((self.v.T @ psi)*np.exp(-1j*self.e*time)))**2
