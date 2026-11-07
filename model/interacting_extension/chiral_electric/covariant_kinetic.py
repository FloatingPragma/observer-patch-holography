#!/usr/bin/env python3
"""A finite gauge-configuration kinetic join for a source overlap projector.

This is a one-coordinate local patch of the source Euclidean overlap bundle,
not a many-body chiral gauge Hamiltonian or a 3+1D continuum construction.
The complete complementary subspace is retained in the finite-penalty parent.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import block_diag, eigh, polar

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
from ..dependencies.chiral_patch import (link_control, overlap_data, overlap_variation,
    projector_basis, canonical_transport, canonical_transport_variation)
from .._source import verify_sources, PIN


def adj(a):
    return a.conj().T


def small(x, name, limit=2e-9):
    value = float(x)
    if not np.isfinite(value) or value > limit:
        raise ArithmeticError(f"{name}: {value}")
    return value


def source_patch(repo=ROOT, charge=6):
    repo = Path(repo).resolve()
    if type(charge) is not int or charge not in (-4,-3,1,2,6):
        raise ValueError("charge must be a native integer hypercharge")
    names=("model/common_source/overlap_regulator.py",)
    verified=verify_sources(repo,names)
    if repo != ROOT:
        verify_sources(ROOT,names)
    from model.common_source import overlap_regulator as source
    gamma, _ = source.gamma_matrices()
    original = link_control(charge)
    # theta is the total holonomy in direction0 on the side-two lattice.
    def data(theta):
        links = original.copy()
        links[:, 0] *= np.exp(1j * charge * theta / 2)
        d = overlap_data(source, charge, links_override=links)
        derivative = np.zeros_like(links)
        derivative[:, 0] = 1j * charge / 2 * links[:, 0]
        _, dp = overlap_variation(d, derivative, gamma)
        d["dP"] = dp
        return d
    d0 = data(0.)
    p0 = d0["P"]
    v0 = projector_basis(p0)
    def frame(theta):
        d = data(theta)
        q, aux = canonical_transport(d["P"], p0)
        dq = canonical_transport_variation(p0, d["dP"], aux)
        d["V"] = q @ v0
        d["dV"] = dq @ v0
        return d
    return frame, {"commit": PIN, "charge_q6Y": charge,
        "source_sha256": verified,
        "configuration": "flat holonomy in direction0 of source side-two Euclidean stencil",
        "temporal_spin_boundary": "antiperiodic",
        "original_toron_per_link": [.17, -.11, .07, .03],
        "source_selected_gauge_inertia_or_patch_boundary": False}


def differential_checks(frame, theta=.037):
    d = frame(theta); p, dp, v, dv = (d[k] for k in ("P", "dP", "V", "dV"))
    n, r = v.shape
    ident = np.eye(n)
    connection = p @ dp - dp @ p
    berry = adj(v) @ dv
    phi = adj(dv) @ (ident-p) @ dv
    # dP + [connection,P] = 0. The opposite sign is the parallel-transport generator.
    checks = {
        "connection_antihermitian": small(np.linalg.norm(connection+adj(connection), "fro"), "Kato connection anti-Hermitian"),
        "covariant_projector_derivative": small(np.linalg.norm(dp+connection@p-p@connection, "fro"), "covariant derivative of P"),
        "frame_connection_antihermitian": small(np.linalg.norm(berry+adj(berry), "fro"), "Berry connection"),
        "frame_P_derivative": small(np.linalg.norm(dp@v+p@dv-dv, "fro"), "d(PV)=dV"),
        "Born_Huang_hermitian": small(np.linalg.norm(phi-adj(phi), "fro"), "Born-Huang adjoint"),
        "Born_Huang_identity": small(np.linalg.norm(phi-(adj(dv)@dv+berry@berry), "fro"), "Born-Huang sign"),
    }
    phi_eigenvalues = eigh((phi+adj(phi))/2, eigvals_only=True)
    if phi_eigenvalues[0] < -2e-12:
        raise ArithmeticError("Born-Huang positivity")
    # Independent centered derivative: checks the sign-kernel and Kato Frechet derivative.
    derivative_rows = []
    for step in (1e-3, 5e-4, 2.5e-4):
        plus, minus = frame(theta+step), frame(theta-step)
        fd = (plus["V"]-minus["V"])/(2*step)
        derivative_rows.append({"step": step,
            "frame_first_derivative_error": float(np.linalg.norm(fd-dv, "fro")),
            "projector_first_derivative_error": float(np.linalg.norm(
                (plus["P"]-minus["P"])/(2*step)-dp, "fro"))})
    if derivative_rows[-1]["frame_first_derivative_error"] > 1e-5:
        raise ArithmeticError("Independent chiral-frame derivative")
    return {"theta": theta, "ambient_dimension": n, "chiral_rank": r,
        "Wilson_gap": float(d["gap"]), "bare_electric_leakage_norm": float(np.linalg.norm((ident-p)@dv, 2)),
        "Born_Huang_minimum_eigenvalue": float(phi_eigenvalues[0]),
        "Born_Huang_maximum_eigenvalue": float(phi_eigenvalues[-1]),
        "checks": checks, "centered_derivative_check": derivative_rows}, d


def grid_model(frame, points, halfwidth=.16, inertia=.01):
    if type(points) is not int or points<1:
        raise ValueError("points must be a positive integer")
    if any(isinstance(x,bool) or not np.isscalar(x) or not np.isfinite(x) or x<=0
           for x in (halfwidth,inertia)):
        raise ValueError("halfwidth and kinetic coefficient must be finite and positive")
    step = 2*halfwidth/(points+1)
    theta = -halfwidth + step*np.arange(1,points+1)
    data = [frame(x) for x in theta]
    frames = [d["V"] for d in data]
    vbig = block_diag(*frames)
    n, r = frames[0].shape
    lap = (2*np.eye(points)-np.eye(points,k=1)-np.eye(points,k=-1))*inertia/step**2
    parent = np.kron(lap,np.eye(n))
    compressed = adj(vbig) @ parent @ vbig
    pbig = vbig @ adj(vbig)
    overlaps = [adj(frames[j])@frames[j+1] for j in range(points-1)]
    parallel = np.eye(points*r)*2*inertia/step**2
    parallel = parallel.astype(complex)
    for j,c in enumerate(overlaps):
        unitary,_ = polar(c)
        parallel[j*r:(j+1)*r,(j+1)*r:(j+2)*r] = -inertia/step**2*unitary
        parallel[(j+1)*r:(j+2)*r,j*r:(j+1)*r] = -inertia/step**2*adj(unitary)
    # Independently assemble every block from neighboring overlaps.
    assembly = np.eye(points*r,dtype=complex)*2*inertia/step**2
    for j,c in enumerate(overlaps):
        assembly[j*r:(j+1)*r,(j+1)*r:(j+2)*r] = -inertia/step**2*c
        assembly[(j+1)*r:(j+2)*r,j*r:(j+1)*r] = -inertia/step**2*adj(c)
    checks = {"complete_compression_assembly_error": small(np.linalg.norm(compressed-assembly,"fro"),"Complete kinetic compression"),
        "positive_projected_kinetic_minimum": float(eigh(compressed,eigvals_only=True,subset_by_index=(0,0))[0]),
        "projector_error": small(np.linalg.norm(pbig@pbig-pbig,"fro"),"Grid projector")}
    return dict(points=points, step=step, theta=theta, V=vbig, P=pbig, parent=parent,
        compressed=compressed, parallel=parallel, data=data, checks=checks,
        inertia=inertia,halfwidth=halfwidth)


def frame_covariance(grid, seed=62062):
    nsite=grid["points"];r=grid["data"][0]["V"].shape[1]
    rng=np.random.default_rng(seed)
    qs=[]
    for _ in range(nsite):
        z=rng.normal(size=(r,r))+1j*rng.normal(size=(r,r))
        q,_=np.linalg.qr(z);qs.append(q)
    qbig=block_diag(*qs)
    changed=grid["V"]@qbig
    actual=adj(changed)@grid["parent"]@changed
    expected=adj(qbig)@grid["compressed"]@qbig
    cov=small(np.linalg.norm(actual-expected,"fro"),"Local frame covariance",1e-8)
    bare=np.kron((2*np.eye(nsite)-np.eye(nsite,k=1)-np.eye(nsite,k=-1)),
                 np.eye(r))*grid["inertia"]/grid["step"]**2
    error=float(np.linalg.norm(adj(qbig)@bare@qbig-bare,2))
    return {"random_seed":seed,"exact_frame_covariance_error":cov,
            "bare_identical_link_non_covariance_norm":error,
            "frame_redundancy_is_not_a_full_gauge_Ward_identity":True}


def complete_penalty(grid, count=8):
    parent,p,v=grid["parent"],grid["P"],grid["V"]
    q=np.eye(len(p))-p
    a=eigh(grid["compressed"],eigvals_only=True,subset_by_index=(0,count-1))
    cross=q@parent@v
    bnorm=float(np.linalg.norm(cross,2))
    parentnorm=float(eigh(parent,eigvals_only=True,subset_by_index=(len(parent)-1,len(parent)-1))[0])
    rows=[];packet={}
    for factor in (10.,40.,160.):
        penalty=factor*parentnorm
        h=parent+penalty*q
        values,vectors=eigh(h,subset_by_index=(0,count-1))
        bound=bnorm**2/(penalty-a[-1])
        errors=a-values
        if errors.min() < -1e-8 or errors.max()>bound+1e-8:
            raise ArithmeticError("Complete-complement Schur energy bracket")
        residual=float(np.linalg.norm(h@vectors-vectors*values,2))
        rows.append({"penalty_factor_times_parent_norm":factor,
            "penalty":penalty,"full_dimension":len(parent),
            "compressed_dimension":grid["compressed"].shape[0],
            "first_eight_parent_energies":values.tolist(),
            "maximum_difference_from_compression":float(errors.max()),
            "analytic_Schur_bound_float64":bound,
            "maximum_complement_amplitude":float(np.linalg.norm(q@vectors,2)),
            "eigensystem_residual":residual})
        if factor==160.:
            packet={"finite_penalty_H":h,"finite_penalty_energies":values,
                    "finite_penalty_vectors":vectors}
    return {"cross_block_norm":bnorm,"parent_kinetic_norm":parentnorm,
            "compressed_first_eight":a.tolist(),"runs":rows,
            "full_complement_retained":True,
            "finite_penalty_is_declared_not_source_selected":True},packet


def main(repo=ROOT,out=None):
    out=Path(out) if out is not None else None
    if out is not None:
        out.mkdir(parents=True,exist_ok=True)
    frame,source=source_patch(repo)
    differential,d=differential_checks(frame)
    grids=[]
    for n in (5,9,17):
        grid=grid_model(frame,n)
        e=eigh(grid["compressed"],eigvals_only=True,subset_by_index=(0,7))
        ep=eigh(grid["parallel"],eigvals_only=True,subset_by_index=(0,7))
        grids.append({"points":n,"step":grid["step"],
            "first_eight_compressed":e.tolist(),
            "first_eight_connection_only":ep.tolist(),
            "ground_difference_from_connection_only":float(e[0]-ep[0]),
            "minimum_Wilson_gap":min(float(x["gap"]) for x in grid["data"]),
            "maximum_plaquette_defect":max(float(x["plaquette"]) for x in grid["data"]),
            "checks":grid["checks"]})
        if n==9:central=grid
    covariance=frame_covariance(central)
    penalty,arrays=complete_penalty(central)
    result={"scope":"single native charge q=6 on a one-holonomy Euclidean overlap configuration patch",
        "source":source,"differential":differential,"grids":grids,
        "frame_covariance":covariance,"complete_penalty":penalty,
        "declared_controls":{"halfwidth":.16,"gauge_kinetic_coefficient":.01,
            "configuration_boundary":"Dirichlet; every interior component retained",
            "finite_penalty_limit_is_declared_completion":True},
        "floating_point_values_are_not_interval_certificates":True,
        "no_particle_mass_or_full_interacting_continuum_claim":True}
    if out is not None:
        np.savez_compressed(out/"chiral_electric_packet.npz",
            theta=central["theta"],V=central["V"],P=central["P"],
            parent_kinetic=central["parent"],compressed_kinetic=central["compressed"],
            connection_only_kinetic=central["parallel"],
            local_P=d["P"],local_dP=d["dP"],local_V=d["V"],local_dV=d["dV"],**arrays)
        (out/"diagnostics.json").write_text(json.dumps(result,indent=2)+"\n")
        print(json.dumps({"differential":differential,
            "ground_grid_values":[(r["points"],r["first_eight_compressed"][0]) for r in grids],
            "covariance":covariance,"complete_penalty":penalty},indent=2))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path,default=ROOT)
    p.add_argument("--out",type=Path,default=HERE)
    a=p.parse_args();main(a.repo,a.out)
