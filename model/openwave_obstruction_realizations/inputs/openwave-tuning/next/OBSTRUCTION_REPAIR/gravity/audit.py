"""Independent scope/response audit; never writes into producer directories.

Run with the theory repository's Python environment.  The finite-q checks use
continuous Bloch quadrature, not a periodic finite-volume mode when q L/2pi is
noninteger.  No thermodynamic error bound or Einstein realization is inferred.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PAULI = (
    np.array([[0, 1], [1, 0]], complex),
    np.array([[0, -1j], [1j, 0]], complex),
    np.diag([1., -1.]).astype(complex),
)
I2 = np.eye(2, dtype=complex)
I4 = np.eye(4, dtype=complex)
ALPHA = tuple(np.kron(p, PAULI[0]) for p in PAULI)
BETA = np.kron(I2, PAULI[2])


def dagger(a):
    return np.swapaxes(a.conj(), -1, -2)


def symbol(k, mu, axis=None):
    """G = i(U-U*)/2 and exact first momentum derivative."""
    k = np.atleast_2d(k)
    mass = math.cos(mu) * I4 - 1j * math.sin(mu) * BETA
    u = np.broadcast_to(mass, (len(k), 4, 4)).copy()
    du = np.zeros_like(u)
    for a in range(3):
        f = np.cos(k[:, a, None, None]) * I4 - 1j * np.sin(k[:, a, None, None]) * ALPHA[a]
        df = -np.sin(k[:, a, None, None]) * I4 - 1j * np.cos(k[:, a, None, None]) * ALPHA[a]
        du = f @ du + (df @ u if a == axis else 0)
        u = f @ u
    return .5j * (u - dagger(u)), .5j * (du - dagger(du))


def grid(L):
    n = np.arange(L)
    return 2 * np.pi * (np.stack(np.meshgrid(n, n, n, indexing="ij"), -1).reshape(-1, 3) + .5) / L


def quadratic_coefficient(k, mu, axis=0):
    h, dh = symbol(k, mu, axis)
    e, u = np.linalg.eigh(h)
    d = dagger(u) @ dh @ u
    eb, ea = e[:, :, None], e[:, None, :]
    gap = eb - ea
    mask = (eb > 0) & (ea < 0)
    first = np.zeros_like(d)
    first[mask] = ((eb + ea) * d)[mask] / (2 * gap[mask])
    chi2 = 2 * np.sum(np.abs(first[mask]) ** 2 / gap[mask]) / len(k)
    return float(chi2), float(np.min(np.abs(e))), h, e, u


def susceptibility(k, q, mu, h=None, e=None, u=None):
    if h is None:
        h, _ = symbol(k, mu)
        e, u = np.linalg.eigh(h)
    hq, _ = symbol(k + np.array([q, 0, 0]), mu)
    eq, uq = np.linalg.eigh(hq)
    vertex = .5 * (h + hq)
    w = dagger(uq) @ vertex @ u
    df = (eq < 0).astype(float)[:, :, None] - (e < 0).astype(float)[:, None, :]
    de = eq[:, :, None] - e[:, None, :]
    mask = df != 0
    return float(-np.sum(np.abs(w[mask]) ** 2 * df[mask] / de[mask]) / len(k))


def fixed_volume_geometry():
    """Metric exp(2 a r^2) delta_ij, R3(0)=-24a, same coordinate topology."""
    records = []
    for a in (-.04, .04):
        for ell in (.6, .3, .15):
            volume = lambda r: 4 * np.pi * quad(lambda x: x*x*math.exp(3*a*x*x), 0, r, epsabs=1e-13)[0]
            target = 4 * np.pi * ell**3 / 3
            r = brentq(lambda r: volume(r) - target, .5*ell, 1.5*ell, xtol=1e-14)
            actual = 4*np.pi*r*r*math.exp(2*a*r*r) - 4*np.pi*ell*ell
            expected = 16*np.pi*a*ell**4/5
            fixed_coordinate = 4*np.pi*ell*ell*(math.exp(2*a*ell*ell)-1)
            records.append(dict(a=a, ell=ell, radius=r, volume_residual=volume(r)-target,
                                area_contrast=actual, leading_area_contrast=expected,
                                ratio=actual/expected,
                                wrong_fixed_coordinate_ratio=fixed_coordinate/expected))
    return records


def fixed_graph_entropy(v=1., scale=1., L=32):
    """Diagnostic only: a declared Wilson chain, unchanged vertices/edges.

    Uniform overall scaling leaves its projector unchanged. Changing kinetic
    speed relative to the fixed Wilson regulator changes its actual entropy.
    This is not a source-derived geometric prescription or an Einstein model.
    """
    k = 2*np.pi*(np.arange(L)+.5)/L
    h = scale*(v*np.sin(k)[:,None,None]*PAULI[0] + (.1+1-np.cos(k))[:,None,None]*PAULI[2])
    e,u = np.linalg.eigh(h)
    p = (u*(e<0)[:,None,:]) @ dagger(u)
    xs = np.arange(L//2)
    dif = xs[:,None]-xs[None,:]
    corr = np.einsum('ijk,kab->ijab', np.exp(1j*dif[:,:,None]*k)/L, p)
    corr = corr.transpose(0,2,1,3).reshape(L, L)
    lam = np.clip(np.linalg.eigvalsh(corr), 1e-15, 1-1e-15)
    return float(-np.sum(lam*np.log(lam)+(1-lam)*np.log(1-lam)))


def main():
    out = {"scope": "Independent finite Bloch, spectral and geometric diagnostics; no gravity row solved."}
    out["inputs"] = {}
    paths = [
        "openwave-push/blockers/X1f_entanglement_area_response/stage2_lib.py",
        "openwave-push/blockers/X1f_entanglement_area_response/stage2_SOLUTION.md",
        "openwave-push/blockers/X1f_entanglement_area_response/stage2_VERDICT.md",
        "openwave-push/blockers/X1f_entanglement_area_response/stage2_part_chi_mu0p05_L32.json",
        "reverse-engineering-reality/model/matched_einstein_limit/PROOF.md",
    ]
    for path in paths:
        data = (ROOT/path).read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        snap = HERE/"inputs"/digest
        snap.parent.mkdir(exist_ok=True)
        snap.write_bytes(data)
        out["inputs"][path] = digest
    coeff = []
    for mu in (0., .05):
        for L in (12,20,28):
            c,gap,*_ = quadratic_coefficient(grid(L),mu)
            record = dict(mu=mu,L=L,chi_q2=c,gap_on_quadrature=gap)
            coeff.append(record)
            print(record,flush=True)
    out["analytic_quadratic_coefficients"] = coeff
    k = grid(20)
    c,gap,h,e,u = quadratic_coefficient(k,.05)
    out["finite_q_checks"] = [dict(q=q,chi_over_q2=susceptibility(k,q,.05,h,e,u)/q**2,
                                   analytic_chi_q2=c) for q in (.02,.01,.005,.002)]
    k32 = grid(32)
    q = 2*np.pi/32
    chi = susceptibility(k32,q,.05)
    ref = json.loads((ROOT/paths[3]).read_text())["chi"]["(1, 0, 0)"]["NN"]
    out["producer_reproduction"] = dict(L=32,mu=.05,q=q,ours=chi,producer=ref,absolute_error=abs(chi-ref))
    out["fixed_proper_volume_geometry"] = fixed_volume_geometry()
    s0 = fixed_graph_entropy()
    sm,sp = fixed_graph_entropy(.99),fixed_graph_entropy(1.01)
    out["fixed_graph_diagnostic"] = dict(reference_entropy=s0,uniform_scale_entropy=fixed_graph_entropy(scale=1.1),
                                        speed_minus_entropy=sm,speed_plus_entropy=sp,
                                        dS_dspeed=(sp-sm)/.02,scope="Declared one-dimensional regulator control only")
    checks = {
        "massless_q2_zero": max(r["chi_q2"] for r in coeff if r["mu"]==0)<1e-25,
        "massive_q2_nonzero": min(r["chi_q2"] for r in coeff if r["mu"]>.0)>1e-6,
        "massive_q2_not_roundoff": min(r["chi_q2"] for r in coeff if r["mu"]>.0)>1e10*max(r["chi_q2"] for r in coeff if r["mu"]==0),
        "finite_q_approaches_derivative": abs(out["finite_q_checks"][-1]["chi_over_q2"]/c-1)<.02,
        "producer_value_reproduced": abs(chi-ref)<1e-13,
        "fixed_volume_geometric_limit": max(abs(r["ratio"]-1) for r in out["fixed_proper_volume_geometry"] if r["ell"]==.15)<.01,
        "wrong_volume_rejected": min(abs(r["wrong_fixed_coordinate_ratio"]-1) for r in out["fixed_proper_volume_geometry"])>1,
        "overall_scaling_null": abs(out["fixed_graph_diagnostic"]["uniform_scale_entropy"]-s0)<1e-10,
        "relative_regulator_changes_entropy": abs(sp-sm)>.001,
    }
    out["checks"] = checks
    out["all_checks_pass"] = all(checks.values())
    (HERE/"RESULT.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2),flush=True)
    if not out["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
