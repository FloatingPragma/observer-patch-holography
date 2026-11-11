"""Recharted port Laplacian on the flat federation for two gluings and two frame rules.

Gluings (unit cell: one up and one down carrier, 24 ports; neighbour k of the up carrier is the down carrier in the
cell displaced by delta_k in (0, A1, A2)):
  bridge   one glue seam per neighbour: (up, P[k]) - (down, P[k]), P[k] the port across the control-face edge k
           (the control gluing face_axis of the simulation package)
  seam     the two carriers share the control-face edge that faces the neighbour: two glue seams per neighbour,
           (up, T[k+1]) - (down, T[k+2]) and (up, T[k+2]) - (down, T[k+1]), T the control face (crossed endpoints,
           as the geometry of an up and a down carrier implies; six control-face ports meet at every mesh vertex)
Frame rules for the recharting phase z on a glue seam:
  template  z = 1 (the same template sections in both carriers)
  ambient   every port has an axis in one common frame: the up carrier has its control axis along +z and its
            control-face vertex k pointing away from neighbour k; the down carrier is the up carrier turned by
            half a turn about z. z = phase of <chi(axis_b)|chi(axis_a)> on every seam.
Clock-free output: levels of L(q) = D - A_twisted(q) (kernel units) near q = 0.
"""
import json, sys
import numpy as np
sys.path.insert(0, '.')
from check_bands import N, E, F, P, TOP, spin_matrices

A1 = np.array([1.0, 0.0]); A2 = np.array([0.5, 3 ** 0.5 / 2]); DEL = (np.zeros(2), A1, A2)
tau = -(A1 + A2) / 3
dirs = [(tau + d) / np.linalg.norm(tau + d) for d in DEL]            # direction from the up carrier to neighbour k

def frame_up():
    n = N[list(TOP)].sum(0); n /= np.linalg.norm(n)
    # rotation taking the control axis to +z
    v = np.cross(n, [0, 0, 1.0]); s = np.linalg.norm(v); c = n[2]
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    R = np.eye(3) + K + K @ K * ((1 - c) / s ** 2)
    X = N @ R.T
    # turn about z so that control vertex 0 points away from neighbour 0
    a0 = np.arctan2(X[TOP[0]][1], X[TOP[0]][0]); t0 = np.arctan2(-dirs[0][1], -dirs[0][0])
    def rz(t): return np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])
    X = X @ rz(t0 - a0).T
    ok = all(np.allclose(X[TOP[k]][:2] / np.linalg.norm(X[TOP[k]][:2]), -dirs[k], atol=1e-9) for k in range(3))
    return X, ok, rz

def coherent(j, m, n):
    Jx, Jy, Jz = spin_matrices(j)
    w, U = np.linalg.eigh(n[0] * Jx + n[1] * Jy + n[2] * Jz)
    return U[:, int(np.argmin(abs(w - m)))]

def build(q, gluing, rule, j, m, AX):
    """AX: 24 x 3 ambient axes. Returns L(q)."""
    chi = [coherent(j, m, AX[p]) if j > 0 else np.ones(1) for p in range(24)]
    def z(a, b):
        if j == 0: return 1.0
        v = np.vdot(chi[b], chi[a]); return v / abs(v)
    L = np.zeros((24, 24), complex)
    def add(a, b, ph, zz):
        L[a, a] += 1; L[b, b] += 1; L[b, a] -= zz * ph; L[a, b] -= np.conj(zz * ph)
    for blk in (0, 12):
        for a, b in E: add(blk + a, blk + b, 1.0, z(blk + a, blk + b))
    for k, d in enumerate(DEL):
        ph = np.exp(1j * (q @ d))                                     # hop from the up cell to the cell displaced by d
        pairs = [(P[k], 12 + P[k])] if gluing == 'bridge' else [(TOP[(k + 1) % 3], 12 + TOP[(k + 2) % 3]), (TOP[(k + 2) % 3], 12 + TOP[(k + 1) % 3])]
        for a, b in pairs: add(a, b, ph, 1.0 if rule == 'template' else z(a, b))
    return L

X, ok, rz = frame_up()
if not ok:                                                            # template orientation is the mirror image: swap two neighbours
    DEL = (DEL[0], DEL[2], DEL[1]); dirs = [dirs[0], dirs[2], dirs[1]]
    X, ok, rz = frame_up()
AX = np.vstack([X, X @ rz(np.pi).T])
print('orientation consistent:', ok)
for g in ('bridge', 'seam'):
    k = 0
    pairs = [(P[k], 12 + P[k])] if g == 'bridge' else [(TOP[1], 12 + TOP[2])]
    a, b = pairs[0]
    print(f'gluing {g}: angle between the axes of two glued ports in the common frame: {np.degrees(np.arccos(np.clip(AX[a] @ AX[b], -1, 1))):.3f} degrees')
out = {}
for g in ('bridge', 'seam'):
    for rule in ('template', 'ambient'):
        for name, j, m in (('scalar', 0, 0), ('spinor', 0.5, 0.5), ('vector', 1, 1)):
            if name == 'scalar' and rule == 'ambient': continue
            d = int(2 * j + 1)
            w0 = np.linalg.eigvalsh(build(np.zeros(2), g, rule, j, m, AX))
            nb = d if d > 1 else 1
            rec = {'q0': w0[:8].tolist(), 'dirs': []}
            lines = []
            for ang in (0, 30, 60, 90):
                e = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
                wq = {qq: np.linalg.eigvalsh(build(qq * e, g, rule, j, m, AX))[:nb + 1] for qq in (0.04, 0.02, 0.01)}
                lin = (wq[0.01] - w0[:nb + 1]) / 0.01
                quad = (wq[0.02] - w0[:nb + 1]) / 0.02 ** 2
                rec['dirs'].append({'angle': ang, 'linear_at_0.01': lin.tolist(), 'quadratic_at_0.02': quad.tolist()})
                lines.append(f'      {ang:2d} deg: d lambda / q {np.round(lin, 6)}   d lambda / q^2 {np.round(quad, 5)}')
            out[f'{g}|{rule}|{name}'] = rec
            print(f'{g:6s} {rule:8s} {name:6s}: lowest levels at q=0 {np.round(w0[:nb + 3], 6)}')
            print('\n'.join(lines))
json.dump(out, open('gluing_bands.json', 'w'), indent=1)
