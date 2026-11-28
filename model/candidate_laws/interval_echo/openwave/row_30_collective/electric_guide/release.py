"""Local physical sign and dimer-number read after removing the four-dimer guards."""

from itertools import combinations

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1 import law


TIMES = (0, 1/16, 1/4, 1, 4, 16, 64)


def apparatus(board, path):
    words = [sum(1 << i for i in bits) for bits in combinations(range(8), 4)]
    ix = {w:i for i,w in enumerate(words)}
    h = np.zeros((len(words), len(words)))
    y = np.zeros(h.shape, dtype=complex)
    native = law.seam_term(8)
    for col, word in enumerate(words):
        def load(j):
            return 2*((word >> j) & 1)
        h[col,col] += sum(8*(6-board.degree[path[1+i]])*load(i)**2 for i in range(8))
        ports = [1, *(load(i) for i in range(8)), 1]
        for i, (a, b) in enumerate(zip(ports, ports[1:])):
            p = law.SRC.index((a,b))
            h[col,col] += native[p,p]
            if 0 < i < 8 and a != b:
                h[ix[word ^ (3 << (i-1))], col] += native[law.SRC.index((b,a)),p]
        if load(0) != load(1):
            y[ix[word ^ 3],col] = 1j if load(0)==2 else -1j
    good = np.array([all(((w >> 2*i) & 3) in (1,2) for i in range(4)) for w in words])
    return h, y, good


def probabilities(board, path):
    full, endpoint, enrolled = apparatus(board, path)
    panels = []
    for keep in (enrolled, np.ones(len(enrolled), bool)):
        good = enrolled[keep]
        y = endpoint[np.ix_(keep, keep)]
        e, q = np.linalg.eigh(full[np.ix_(keep, keep)])
        initial = q.T[:, good]
        signed = y[np.ix_(good, good)]
        rows = []
        for t in TIMES:
            u = q @ (np.exp(-1j*e*t)[:, None] * initial)
            retention = np.linalg.norm(u[good])**2/good.sum()
            valid = np.linalg.norm(y @ u)**2/good.sum()
            correlation = np.trace(u.conj().T @ y @ u @ signed).real/good.sum()
            rows.append([retention, (valid-correlation)/2, 1-valid, (valid+correlation)/2])
        panels.append(rows)
    return np.array(panels)


def counts(board, path, rng, shots):
    p = probabilities(board,path)
    bank = np.empty_like(p, dtype=np.int32)
    for i, panel in enumerate(p):
        for j, row in enumerate(panel):
            rounded = np.clip(row, 0, 1).round(12)
            bank[i,j,0] = rng.binomial(shots, rounded[0])
            vals = rounded[1:]
            bank[i,j,1:] = rng.multinomial(shots, vals/vals.sum())
    return bank


def read(bank, shots):
    if (bank.shape!=(2,len(TIMES),4) or bank.dtype.kind not in "iu" or np.min(bank)<0
            or np.max(bank)>shots or np.any(bank[:,:,1:].sum(axis=2)!=shots)):
        raise ValueError("Incomplete guard release")
    return [dict(t=t, guarded_retention=float(bank[0,i,0]/shots),
                 free_retention=float(bank[1,i,0]/shots),
                 guarded_response=float((bank[0,i,3]-bank[0,i,1])/shots),
                 free_response=float((bank[1,i,3]-bank[1,i,1])/shots)) for i,t in enumerate(TIMES)]
