"""Twisted (recharted) port Laplacian on the flat federation: lowest bands near q = 0, per sector.

Same lattice as check_bands.py (control gluing face_axis, glue ports P[k], template-identity frame rule: z = 1 on
glue seams in the common template gauge). Inside a carrier the hop a -> b carries the recharting phase
z_ab = phase of <chi_b|chi_a> for the section type (spin j, helicity m): flux 2m per face.
L(q) = D - A_twisted(q), a 24 x 24 Hermitian matrix; lambda = eigenvalues (relaxation rates in units of 1/(2M) per
attempt). Readings of a band lambda(q) = lambda_0 + c q^2 (both are declared clocks, see the theory session's note):
  first-order clock   H = L (or -log T): E = lambda
  recurrence clock    cos(omega) = 1 - lambda/(2M) per attempt, omega^2 ~ lambda/M: rest frequency^2 ~ lambda_0, speed^2 ~ c
"""
import json, sys
import numpy as np
sys.path.insert(0, '.')
from check_bands import N, NG, E, P, DEL, sections

def twisted(q, W):
    L = np.zeros((24, 24), complex)
    for blk in (0, 12):
        for a, b in E:
            z = np.vdot(W[b].conj(), W[a].conj())        # <chi_b|chi_a> with rows W[p] = <chi_p|
            z = z / abs(z)
            L[blk + a, blk + a] += 1; L[blk + b, blk + b] += 1
            L[blk + b, blk + a] -= z; L[blk + a, blk + b] -= np.conj(z)
    for k, d in enumerate(DEL):
        u, dn = P[k], 12 + P[k]
        ph = np.exp(-1j * (q @ d))
        L[u, u] += 1; L[dn, dn] += 1
        L[u, dn] -= ph; L[dn, u] -= np.conj(ph)
    return L

out = {}
for name, W in (('scalar', np.ones((12, 1), complex)), ('spinor convex (flux 1)', sections(0.5, 0.5, N)),
                ('spinor great (flux 7)', sections(0.5, 0.5, NG)), ('vector convex helicity 1 (flux 2)', sections(1, 1, N)),
                ('vector great helicity 1', sections(1, 1, NG))):
    d = W.shape[1]
    w0 = np.linalg.eigvalsh(twisted(np.zeros(2), W))
    print(f'{name}: lowest levels at q=0: {np.round(w0[:d + 3], 6)}')
    rec = {'q0': w0[:8].tolist(), 'dirs': []}
    for ang in (0, 30, 45, 90):
        e = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
        rows = []
        for qq in (0.04, 0.02, 0.01):
            w = np.linalg.eigvalsh(twisted(qq * e, W))[:d]
            rows.append(((w - w0[:d]) / qq ** 2).tolist())
        lin = [(np.linalg.eigvalsh(twisted(qq * e, W))[:d] - w0[:d]) / qq for qq in (0.02, 0.01)]
        rec['dirs'].append({'angle': ang, 'curvature': rows, 'linear': [l.tolist() for l in lin]})
        print(f'   direction {ang:2d} deg: (lambda(q) - lambda(0))/q^2 at q = 0.02, 0.01: {np.round(rows[1], 6)} {np.round(rows[2], 6)};  /q at 0.01: {np.round(lin[1], 6)}')
    out[name] = rec
json.dump(out, open('twisted_bands.json', 'w'), indent=1)
