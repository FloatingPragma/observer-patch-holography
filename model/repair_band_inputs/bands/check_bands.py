"""Bands of the consensus (check) operators on the flat federation, by sector.

Flat honeycomb federation in the form of oph_defects.bloch.bloch_honeycomb (control gluing face_axis): unit cell =
one up and one down carrier (24 ports); glue seam k joins port P[k] of the up carrier to port P[k] of the down
carrier in the cell displaced by delta_k in (0, A1, A2). Frame rule on glue seams: template identity (the spin
direction relative to the carrier's own ports is kept), so a section has the same template amplitudes in both
carriers and the glue check is s(u, P_k) = s(d, P_k).

H_j(q) = sum over local checks inside carriers |c><c| + w_glue * sum over glue seams |g_k(q)><g_k(q)|,
g_k(q) = (|u,P_k> - exp(i q.delta_k) |d,P_k>)/sqrt2. Local checks: seams for scalars (then H = L/2, L the port
Laplacian), faces for spinors, rhombi for vectors. The zero modes at q = 0 are the global sections.
Reported: number of zero modes at q = 0, and the low branches E(q)/q^2 for small q in several directions.
"""
import itertools, json, sys
import numpy as np

PHI = (1 + 5 ** 0.5) / 2
PV = np.array([(-1, PHI, 0), (1, PHI, 0), (-1, -PHI, 0), (1, -PHI, 0), (0, -1, PHI), (0, 1, PHI),
               (0, -1, -PHI), (0, 1, -PHI), (PHI, 0, -1), (PHI, 0, 1), (-PHI, 0, -1), (-PHI, 0, 1)], float)
N = PV / np.linalg.norm(PV, axis=1)[:, None]
G = np.where(np.isclose(abs(PV), PHI), np.sign(PV) * (-1 / PHI), PV); NG = G / np.linalg.norm(G, axis=1)[:, None]
E = [(a, b) for a in range(12) for b in range(a + 1, 12) if abs(np.linalg.norm(PV[a] - PV[b]) - 2) < 1e-9]
ES = set(E)
adj = lambda a, b: (min(a, b), max(a, b)) in ES
F = [f for f in itertools.combinations(range(12), 3) if adj(f[0], f[1]) and adj(f[1], f[2]) and adj(f[0], f[2])]
RH = [[a, b] + [c for c in range(12) if adj(a, c) and adj(b, c)] for a, b in E]
TOP = (0, 11, 5)
assert tuple(sorted(TOP)) in [tuple(sorted(f)) for f in F]
P = []
for k in range(3):
    e = {TOP[(k + 1) % 3], TOP[(k + 2) % 3]}
    third = [v for f in F if e <= set(f) and set(f) != set(TOP) for v in f if v not in e]
    assert len(third) == 1; P.append(third[0])
A1 = np.array([1.0, 0.0]); A2 = np.array([0.5, 3 ** 0.5 / 2])
DEL = (np.zeros(2), A1, A2)

def spin_matrices(j):
    d = int(round(2 * j + 1)); m = np.array([j - i for i in range(d)])
    Jz = np.diag(m).astype(complex); Jp = np.zeros((d, d), complex)
    for i in range(1, d):
        Jp[i - 1, i] = np.sqrt(j * (j + 1) - m[i] * (m[i] + 1))
    return (Jp + Jp.conj().T) / 2, (Jp - Jp.conj().T) / (2j), Jz

def sections(j, m, axes):
    Jx, Jy, Jz = spin_matrices(j); rows = []
    for n in axes:
        w, U = np.linalg.eigh(n[0] * Jx + n[1] * Jy + n[2] * Jz)
        rows.append(U[:, int(np.argmin(abs(w - m)))].conj())
    return np.array(rows)

def local_checks(W, clusters):
    rows = []
    for S in clusters:
        u, s, vh = np.linalg.svd(W[S, :], full_matrices=True)
        r = int((s > 1e-9).sum())
        for c in u[:, r:].T:
            row = np.zeros(12, complex); row[S] = c.conj(); rows.append(row)
    return np.array(rows)

def H_of(q, C, wg=1.0):
    H = np.zeros((24, 24), complex)
    Hc = C.conj().T @ C
    H[:12, :12] = Hc; H[12:, 12:] = Hc
    for k, d in enumerate(DEL):
        g = np.zeros(24, complex); g[P[k]] = 1 / np.sqrt(2); g[12 + P[k]] = -np.exp(1j * (q @ d)) / np.sqrt(2)
        H += wg * np.outer(g, g.conj())
    return H

SECT = {
    'scalar (seam checks)': (np.ones((12, 1), complex), [list(e) for e in E]),
    'spinor convex (face checks)': (sections(0.5, 0.5, N), [list(f) for f in F]),
    'spinor great (face checks)': (sections(0.5, 0.5, NG), [list(f) for f in F]),
    'vector helicity 1 convex (rhombus checks)': (sections(1, 1, N), RH),
    'vector helicity 0 convex (rhombus checks)': (sections(1, 0, N), RH),
    'vector helicity 1 great (rhombus checks)': (sections(1, 1, NG), RH),
}
out = {'glue_ports': P}
print('glue ports', P, ' polar angle of glue ports from the control-face axis (degrees):',
      np.round(np.degrees(np.arccos(N[P] @ (N[list(TOP)].sum(0) / np.linalg.norm(N[list(TOP)].sum(0))))), 3))
for name, (W, cl) in SECT.items():
    C = local_checks(W, cl)
    d = W.shape[1]
    w0 = np.linalg.eigvalsh(H_of(np.zeros(2), C))
    nz = int((w0 < 1e-10).sum())
    rec = {'checks_per_carrier': int(len(C)), 'zero_modes_q0': nz, 'gap_above_q0': float(w0[nz]), 'branches': []}
    line = f'{name}: {len(C)} checks per carrier, zero modes at q=0: {nz}, next level {w0[nz]:.6f}'
    print(line)
    for ang in (0, 15, 30, 45, 90):
        e = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
        vals = []
        for qq in (0.02, 0.01, 0.005):
            w = np.linalg.eigvalsh(H_of(qq * e, C))[:nz]
            vals.append(w / qq ** 2)
        rec['branches'].append({'angle': ang, 'E_over_q2_at_q_0.02_0.01_0.005': [v.tolist() for v in vals]})
        print(f'   direction {ang:3d} deg: E/q^2 of the {nz} lowest branches at q=0.01: {np.round(vals[1], 6)}   at q=0.005: {np.round(vals[2], 6)}')
    out[name] = rec
json.dump(out, open('check_bands.json', 'w'), indent=1)
