"""Small same-cone manufactured-field test of the dressed Maxwell/scalar action.

The analytic uniform estimate lives in WHITNEY_SPATIAL_CONSISTENCY.tex. This
deterministic probe checks three refinements, all scheduled results, analytic
first variations, and the actual potential-dependent interpolation. It does
not simulate an observer history or prove trajectory/quantum convergence.
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path
import re

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "runtime/whitney_spatial_consistency_receipt.json"
SCHEMA = "oph.whitney_spatial_consistency.v1"
SCOPE = "MANUFACTURED_SMOOTH_ACTION_VARIATION_REFINEMENT__NO_TRAJECTORY_OR_PHYSICAL_IDENTIFICATION"
PARAMETERS = {"charge": 0.7, "mass_squared": 0.4, "quartic": 0.2,
              "time_interval": [0.0, 0.5], "levels": [0, 1, 2]}
PIN_PATHS = ["Lean/Screen/SeamCurrentEdge30Moment.lean",
             "Lean/ObserverPatchHolography/CoreAxioms.lean",
             "Lean/Screen/WhitneySpatialConsistency.lean",
             "code/electromagnetism/whitney_spatial_consistency.py",
             "code/electromagnetism/verify_whitney_spatial_consistency.py",
             "code/electromagnetism/test_whitney_spatial_consistency.py"]
# Every row is coefficient, followed by powers of t,x,y,z. These are fixed
# manufactured data, not fitted to the measured errors.
POLYNOMIALS = {
    "Ax": [(0.2,0,0,0,0),(0.1,1,0,0,0),(0.15,0,0,1,0),(0.07,0,1,0,1)],
    "Ay": [(-0.1,0,0,0,0),(0.09,1,1,0,0),(0.13,0,0,0,1),(0.03,0,2,0,0)],
    "Az": [(0.05,0,0,0,0),(0.11,0,1,0,0),(0.04,1,0,1,0),(0.06,0,0,1,1)],
    "phi": [(0.1,0,0,0,0),(0.2,0,1,0,0),(-0.1,0,0,1,0),(0.05,0,0,0,2),
            (0.03,1,0,0,0),(0.03,1,1,0,0)],
    "u": [(1.,0,0,0,0),(0.1,0,1,0,0),(0.08,0,0,2,0),(0.03,1,0,0,1)],
    "v": [(0.2,0,0,0,0),(0.12,0,0,0,1),(0.05,0,1,1,0),(0.04,1,1,0,0)],
    "Bx": [(0.1,0,0,0,0),(0.09,0,0,1,0),(0.04,1,0,0,1)],
    "By": [(-0.07,0,0,0,0),(0.05,0,0,0,2),(0.03,1,0,0,0)],
    "Bz": [(0.08,0,0,0,0),(0.06,0,1,1,0)],
    "eta": [(0.03,0,0,0,0),(0.04,0,0,1,0),(0.02,1,1,0,0)],
    "zu": [(0.05,0,0,0,0),(0.06,0,0,1,0),(0.02,1,0,0,0)],
    "zv": [(0.04,0,1,0,0),(0.03,0,0,0,2),(-0.01,1,0,0,1)],
}
EDGES = tuple(combinations(range(4), 2))


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n").encode("ascii")


def polynomial(name, time, points, derivative=None):
    coordinates = [np.broadcast_to(time, points.shape[:-1]), *np.moveaxis(points, -1, 0)]
    result = np.zeros(points.shape[:-1])
    for coefficient, *powers in POLYNOMIALS[name]:
        if derivative is not None:
            coefficient *= powers[derivative]
            if not coefficient:
                continue
            powers[derivative] -= 1
        term = np.full(points.shape[:-1], coefficient)
        for coordinate, power in zip(coordinates, powers, strict=True):
            if power:
                term *= coordinate**power
        result += term
    return result


def vector(prefix, time, points, derivative=None):
    return np.stack([polynomial(prefix+axis, time, points, derivative) for axis in "xyz"], axis=-1)


def scalar(real, imaginary, time, points, derivative=None):
    return polynomial(real, time, points, derivative)+1j*polynomial(imaginary, time, points, derivative)


def original_mesh():
    source = (ROOT/"Lean/Screen/SeamCurrentEdge30Moment.lean").read_text(encoding="utf-8")
    source = source.split("def portVector",1)[1].split("theorem portVector_positivePort",1)[0]
    golden = (1+np.sqrt(5))/2
    literals = {"0":0., "1":1., "-1":-1., "φ":golden, "-φ":-golden}
    vertices = np.array([[0.,0.,0.]]+[[literals[v.strip()] for v in row.split(",")]
                        for row in re.findall(r"!\[([^\[\]]+)\]", source)])
    source = (ROOT/"Lean/ObserverPatchHolography/CoreAxioms.lean").read_text(encoding="utf-8")
    source = source.split("def orientedFaces",1)[1].split("def faceEdges",1)[0]
    faces = [tuple(int(v)+1 for v in row) for row in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)", source)]
    if vertices.shape != (13,3) or len(faces) != 20:
        raise ValueError("committed cone census")
    return vertices, np.array([(0,)+face for face in faces], dtype=int)


def refine(vertices, tetrahedra):
    """Eight midpoint children; shared faces receive the same four triangles.

    The interior octahedron diagonal is local to a tetrahedron. The proof
    assumes a shape-regular family; finite-run shape numbers are reported
    rather than extrapolating this particular diagonal rule indefinitely.
    """
    output = vertices.tolist()
    midpoint = {}
    for first, second in sorted({tuple(sorted((int(t[i]),int(t[j]))))
                                 for t in tetrahedra for i,j in EDGES}):
        midpoint[first,second] = len(output)
        output.append(((vertices[first]+vertices[second])/2).tolist())
    children = []
    for a,b,c,d in tetrahedra.tolist():
        ab,ac,ad,bc,bd,cd = [midpoint[tuple(sorted(pair))]
                             for pair in ((a,b),(a,c),(a,d),(b,c),(b,d),(c,d))]
        children.extend([(a,ab,ac,ad),(b,ab,bc,bd),(c,ac,bc,cd),(d,ad,bd,cd),
                         (ab,ac,ad,cd),(ab,ac,bc,cd),(ab,ad,bd,cd),(ab,bc,bd,cd)])
    return np.array(output), np.array(children, dtype=int)


def mesh(level):
    vertices, tetrahedra = original_mesh()
    for _ in range(level):
        vertices, tetrahedra = refine(vertices, tetrahedra)
    return vertices, tetrahedra


def quadrature(order=6):
    """Unweighted Gauss--Legendre on a Duffy parameter cube."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = (nodes+1)/2, weights/2
    r,s,t = np.meshgrid(nodes,nodes,nodes,indexing="ij")
    wr,ws,wt = np.meshgrid(weights,weights,weights,indexing="ij")
    lam = np.stack(((1-r)*(1-s)*(1-t), r, (1-r)*s, (1-r)*(1-s)*t),axis=-1).reshape(-1,4)
    weight = (wr*ws*wt*(1-r)**2*(1-s)).ravel()
    return lam, weight


def element_geometry(xyz):
    affine = np.concatenate((np.ones(xyz.shape[:-1]+(1,)),xyz),axis=-1)
    gradient = np.swapaxes(np.linalg.inv(affine)[:,1:,:],1,2)
    volume_jacobian = abs(np.linalg.det(np.swapaxes(xyz[:,1:]-xyz[:,:1],1,2)))
    return gradient, volume_jacobian


def sample_edges(prefix, time, xyz, derivative=None):
    # Degree-two spatial polynomials have exact edge integrals with this rule.
    nodes, weights = np.polynomial.legendre.leggauss(3)
    nodes, weights = (nodes+1)/2, weights/2
    coefficients = []
    for i,j in EDGES:
        delta = xyz[:,j]-xyz[:,i]
        points = xyz[:,i,None,:]+nodes[None,:,None]*delta[:,None,:]
        values = vector(prefix,time,points,derivative)
        coefficients.append(np.einsum("nqc,nc,q->n",values,delta,weights))
    return np.stack(coefficients,axis=-1)


def reconstructed(time, xyz, lam, amplitude=0., omit_dressing_variation=False):
    """Return fields and their full directional derivative, in real amplitude."""
    grad,_ = element_geometry(xyz)
    w = np.stack([lam[None,:,i,None]*grad[:,None,j,:]-lam[None,:,j,None]*grad[:,None,i,:]
                  for i,j in EDGES],axis=2)
    curl_w = np.stack([2*np.cross(grad[:,i],grad[:,j]) for i,j in EDGES],axis=1)
    a,b = sample_edges("A",time,xyz),sample_edges("B",time,xyz)
    at,bt = sample_edges("A",time,xyz,0),sample_edges("B",time,xyz,0)
    a,at = a+amplitude*b,at+amplitude*bt
    A,B = np.einsum("ne,nqec->nqc",a,w),np.einsum("ne,nqec->nqc",b,w)
    At,Bt = np.einsum("ne,nqec->nqc",at,w),np.einsum("ne,nqec->nqc",bt,w)
    magnetic,dmagnetic = np.einsum("ne,nec->nc",a,curl_w),np.einsum("ne,nec->nc",b,curl_w)
    phi = polynomial("phi",time,xyz)+amplitude*polynomial("eta",time,xyz)
    eta = polynomial("eta",time,xyz)
    ph,et = phi@lam.T,eta@lam.T
    E = -At-np.einsum("ni,nic->nc",phi,grad)[:,None,:]
    dE = -Bt-np.einsum("ni,nic->nc",eta,grad)[:,None,:]
    psi = scalar("u","v",time,xyz)+amplitude*scalar("zu","zv",time,xyz)
    psit = scalar("u","v",time,xyz,0)+amplitude*scalar("zu","zv",time,xyz,0)
    zeta,zetat = scalar("zu","zv",time,xyz),scalar("zu","zv",time,xyz,0)
    matrices = []
    for coefficients in (a,b,at,bt):
        matrix = np.zeros((len(xyz),4,4))
        for edge,(i,j) in enumerate(EDGES):
            matrix[:,i,j],matrix[:,j,i] = coefficients[:,edge],-coefficients[:,edge]
        matrices.append(matrix)
    phase,beta,phaset,betat = [np.einsum("nij,qj->nqi",mat,lam) for mat in matrices]
    phasegrad,betagrad = [np.einsum("nij,njc->nic",mat,grad) for mat in matrices[:2]]
    charge = PARAMETERS["charge"]
    unit = np.exp(1j*charge*phase)
    base = lam[None,:,:]*unit
    spatial_base = unit[:,:,:,None]*(grad[:,None,:,:]+1j*charge*lam[None,:,:,None]*phasegrad[:,None,:,:])
    scalar_value = np.einsum("nqi,ni->nq",base,psi)
    scalar_time = np.sum(base*(psit[:,None,:]+1j*charge*phaset*psi[:,None,:]),axis=-1)
    scalar_gradient = np.einsum("nqic,ni->nqc",spatial_base,psi)
    if omit_dressing_variation:
        beta,betat,betagrad = np.zeros_like(beta),np.zeros_like(betat),np.zeros_like(betagrad)
    varied_nodal = zeta[:,None,:]+1j*charge*beta*psi[:,None,:]
    Z = np.sum(base*varied_nodal,axis=-1)
    Zt = np.sum(base*(zetat[:,None,:]+1j*charge*betat*psi[:,None,:]
                     +1j*charge*beta*psit[:,None,:]+1j*charge*phaset*varied_nodal),axis=-1)
    Zgrad = np.sum(spatial_base*varied_nodal[:,:,:,None],axis=2)
    Zgrad += np.sum(base[:,:,:,None]*1j*charge*betagrad[:,None,:,:]*psi[:,None,:,None],axis=2)
    Q = scalar_time+1j*charge*ph*scalar_value
    P = scalar_gradient-1j*charge*A*scalar_value[:,:,None]
    dQ = Zt+1j*charge*et*scalar_value+1j*charge*ph*Z
    dP = Zgrad-1j*charge*B*scalar_value[:,:,None]-1j*charge*A*Z[:,:,None]
    return (E,magnetic[:,None,:],scalar_value,Q,P,dE,dmagnetic[:,None,:],Z,dQ,dP)


def continuum(time, points, amplitude=0.):
    A = vector("A",time,points)+amplitude*vector("B",time,points)
    B = vector("B",time,points)
    At = vector("A",time,points,0)+amplitude*vector("B",time,points,0)
    Bt = vector("B",time,points,0)
    gradient_a = np.stack([vector("A",time,points,j)+amplitude*vector("B",time,points,j) for j in (1,2,3)],axis=-1)
    gradient_b = np.stack([vector("B",time,points,j) for j in (1,2,3)],axis=-1)
    def curl(matrix):
        return np.stack((matrix[...,2,1]-matrix[...,1,2],matrix[...,0,2]-matrix[...,2,0],matrix[...,1,0]-matrix[...,0,1]),axis=-1)
    phi = polynomial("phi",time,points)+amplitude*polynomial("eta",time,points)
    eta = polynomial("eta",time,points)
    gradphi = np.stack([polynomial("phi",time,points,j)+amplitude*polynomial("eta",time,points,j) for j in (1,2,3)],axis=-1)
    gradeta = np.stack([polynomial("eta",time,points,j) for j in (1,2,3)],axis=-1)
    psi = scalar("u","v",time,points)+amplitude*scalar("zu","zv",time,points)
    zeta = scalar("zu","zv",time,points)
    psit = scalar("u","v",time,points,0)+amplitude*scalar("zu","zv",time,points,0)
    zetat = scalar("zu","zv",time,points,0)
    gradpsi = np.stack([scalar("u","v",time,points,j)+amplitude*scalar("zu","zv",time,points,j) for j in (1,2,3)],axis=-1)
    gradzeta = np.stack([scalar("zu","zv",time,points,j) for j in (1,2,3)],axis=-1)
    charge = PARAMETERS["charge"]
    return (-At-gradphi,curl(gradient_a),psi,psit+1j*charge*phi*psi,
            gradpsi-1j*charge*A*psi[...,None],-Bt-gradeta,curl(gradient_b),zeta,
            zetat+1j*charge*eta*psi+1j*charge*phi*zeta,
            gradzeta-1j*charge*B*psi[...,None]-1j*charge*A*zeta[...,None])


def densities(fields):
    E,magnetic,psi,Q,P,dE,dmagnetic,Z,dQ,dP = fields
    norm = lambda v: np.sum(abs(v)**2,axis=-1)
    pair = lambda v,w: np.real(np.sum(np.conj(v)*w,axis=-1))
    mass,g = PARAMETERS["mass_squared"],PARAMETERS["quartic"]
    psi2 = abs(psi)**2
    action = .5*norm(E)-.5*norm(magnetic)+abs(Q)**2-norm(P)-mass*psi2-.5*g*psi2**2
    variation = pair(E,dE)-pair(magnetic,dmagnetic)+2*np.real(np.conj(Q)*dQ)-2*pair(P,dP)
    variation -= 2*(mass+g*psi2)*np.real(np.conj(psi)*Z)
    return action,variation


def evaluate(level, order=6, amplitude=0., omit_dressing_variation=False):
    vertices,tetrahedra = mesh(level)
    lam,weights = quadrature(order)
    times,timeweights = np.polynomial.legendre.leggauss(4)
    times,timeweights = (times+1)/4,timeweights/4
    totals = np.zeros(6)
    for start in range(0,len(tetrahedra),64):
        xyz = vertices[tetrahedra[start:start+64]]
        _,jacobian = element_geometry(xyz)
        points = np.einsum("qi,nic->nqc",lam,xyz)
        weight = jacobian[:,None]*weights[None,:]
        for time,timeweight in zip(times,timeweights,strict=True):
            actual = reconstructed(time,xyz,lam,amplitude,omit_dressing_variation)
            expected = continuum(time,points,amplitude)
            integrands = [*densities(actual),*densities(expected),
                          abs(actual[2]-expected[2])**2,
                          np.sum(abs(actual[4]-expected[4])**2,axis=-1)]
            totals += timeweight*np.array([np.sum(weight*value) for value in integrands])
    _,jac = element_geometry(vertices[tetrahedra])
    diameter = max(np.linalg.norm(vertices[tetrahedra[:,j]]-vertices[tetrahedra[:,i]],axis=-1).max() for i,j in EDGES)
    grad,_ = element_geometry(vertices[tetrahedra])
    # Inradius is 1 / sum_i |grad lambda_i|.
    shapes = np.array([max(np.linalg.norm(vertices[tet[j]]-vertices[tet[i]]) for i,j in EDGES)
                       for tet in tetrahedra])*np.linalg.norm(grad,axis=-1).sum(axis=-1)
    return {"level":level,"vertices":len(vertices),"tetrahedra":len(tetrahedra),
            "diameter":float(diameter),"max_diameter_inradius":float(max(shapes)),
            "volume":float(jac.sum()/6),"action":float(totals[0]),"first_variation":float(totals[1]),
            "continuum_action":float(totals[2]),"continuum_first_variation":float(totals[3]),
            "action_error":float(abs(totals[0]-totals[2])),
            "first_variation_error":float(abs(totals[1]-totals[3])),
            "scalar_l2_error":float(np.sqrt(totals[4])),
            "covariant_gradient_l2_error":float(np.sqrt(totals[5]))}


def build_receipt():
    return {"schema":SCHEMA,"scope":SCOPE,"parameters":PARAMETERS,
            "polynomials":POLYNOMIALS,"quadrature":{"space_duffy_legendre_order":6,"time_legendre_order":4},
            "source_pins":{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in PIN_PATHS},
            "refinements":[evaluate(level) for level in PARAMETERS["levels"]],
            "continuum_trajectory_claimed":False,"observer_history_claimed":False,
            "uniform_bound_certified_by_numerics":False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=OUTPUT)
    args = parser.parse_args()
    receipt = build_receipt()
    args.output.write_bytes(canonical(receipt))
    print(json.dumps(receipt["refinements"],indent=2))


if __name__ == "__main__":
    main()
