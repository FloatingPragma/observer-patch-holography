"""Spectrum of the family walk in the walk clock, from the Bloch kernel (steps 5 and 6 at kernel level).

Family walk (declared): the address register of a quantum is that of the carrier that holds it: its sixty oriented
internal seams and both orientations of its three (bridge) or six (shared seam) glue seams, NA addresses. The step is
W = Q R with the public reversible seam operation Q and the reflection R over that register. Its blank-record
discriminant is 1 - L_G / NA with L_G the recharted port Laplacian of the federation (derivation in STEPS.md and in
the docstring of gluing_bands.py; the uniform normalization holds because every carrier has NA addresses). So every
band lambda_n(q) of L_G gives walk phases +-omega_n(q), cos omega = 1 - lambda / NA, and a blank-record excitation
A|x> has weight |<phi_n(q)|x>|^2 / 2 on each sign.

Reported per sector and gluing, in the one clock "walk steps": the massless speed v0 (loads), and for each of the
lowest branches at q = 0: rest phase omega_0, v^2 = d(omega^2)/d(q^2), v / v0, reduced Compton length v / omega_0
in mesh edges, the q up to which omega^2 = omega_0^2 + v^2 q^2 holds within one percent of the kinetic term, and
the weight of a port-localized excitation (averaged over the 24 ports of the cell) on the branch at q -> 0.
"""
import json, sys
import numpy as np
sys.path.insert(0, '.')
import gluing_bands as GB

def omega(lam, NA): return np.arccos(np.clip(1 - lam / NA, -1, 1))

def analyse(g, rule, name, j, m, NA, nlow):
    d = int(2 * j + 1)
    L0 = GB.build(np.zeros(2), g, rule, j, m, GB.AX)
    w0, U0 = np.linalg.eigh(L0)
    rows = []
    e = np.array([1.0, 0.0])
    qs = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8])
    lam = np.array([np.linalg.eigvalsh(GB.build(q * e, g, rule, j, m, GB.AX)) for q in qs])
    lam90 = np.array([np.linalg.eigvalsh(GB.build(q * np.array([0.0, 1.0]), g, rule, j, m, GB.AX)) for q in qs])
    for n in range(nlow):
        om0 = omega(w0[n], NA)
        om = omega(lam[:, n], NA); om90 = omega(lam90[:, n], NA)
        v2 = (om[0] ** 2 - om0 ** 2) / qs[0] ** 2
        flat = abs(lam[3, n] - w0[n]) < 1e-9
        qmax = None
        if not flat and v2 > 0:
            dev = abs((om ** 2 - om0 ** 2) / (v2 * qs ** 2) - 1)
            ok = qs[dev < 0.01]; qmax = float(ok.max()) if len(ok) else 0.0
        aniso = float(abs(om[4] ** 2 - om90[4] ** 2) / max(1e-15, abs(om[4] ** 2 - om0 ** 2))) if not flat else 0.0
        weight = float((abs(U0[:, n]) ** 2).mean())            # equal for every normalized eigenvector: 1/24; kept for the record
        rows.append({'level': float(w0[n]), 'omega0': float(om0), 'v2': float(v2), 'flat': bool(flat), 'qmax_1pct': qmax, 'anisotropy_at_q0.2': aniso})
    return rows

out = {}
for g, NA in (('bridge', 66), ('seam', 72)):
    sc = analyse(g, 'template', 'scalar', 0, 0, NA, 1)[0]
    v0 = sc['v2'] ** 0.5
    out[g] = {'addresses': NA, 'massless_speed_edges_per_step': v0, 'sectors': {}}
    print(f'gluing {g}: {NA} addresses per carrier; loads: massless, speed {v0:.6f} edges per walk step (v^2 = {sc["v2"]:.6e}), one percent up to q = {sc["qmax_1pct"]}')
    for rule in ('template', 'ambient'):
        for name, j, m, nlow in (('spinor', 0.5, 0.5, 3), ('vector', 1, 1, 4)):
            if g == 'seam' and rule == 'template': continue        # no global section for crossed endpoints with z = 1
            if g == 'bridge' and rule == 'ambient': continue       # identical to the template rule (one glue seam per neighbour)
            rows = analyse(g, rule, name, j, m, NA, nlow)
            out[g]['sectors'][f'{name}|{rule}'] = rows
            print(f'  {name}, frame rule {rule}:')
            for r in rows:
                if r['flat']: print(f'    level {r["level"]:.6f}: rest phase {r["omega0"]:.5f} rad per step, FLAT (does not move)')
                else:
                    v = np.sign(r['v2']) * abs(r['v2']) ** 0.5
                    print(f'    level {r["level"]:.6f}: rest phase {r["omega0"]:.5f}, v^2 = {r["v2"]:.4e} (v / v0 = {v / v0:+.4f}), reduced Compton length {abs(v) / r["omega0"]:.4f} edges, shell within 1 percent up to q = {r["qmax_1pct"]}, anisotropy {r["anisotropy_at_q0.2"]:.1e}')
json.dump(out, open('walk_spectrum.json', 'w'), indent=1)
