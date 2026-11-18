"""Independent geometry, quadrature and action-derivative receipt verifier.

No producer module is imported. Only its declared polynomial coefficients
are read as literal data, after checking every current source pin. Geometry
is reconstructed from the Lean coordinates and face adjacency. The Whitney
field is evaluated through a full antisymmetric matrix, integration uses
weighted Gauss--Jacobi rules, and first variations use five-point differences
of the independently evaluated action. These are numerical implementation
checks; the uniform smooth-window estimate is an analytic paper theorem.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path
import re

import numpy as np
from scipy.special import roots_jacobi
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "runtime/whitney_spatial_consistency_receipt.json"
SCOPE = "MANUFACTURED_SMOOTH_ACTION_VARIATION_REFINEMENT__NO_TRAJECTORY_OR_PHYSICAL_IDENTIFICATION"
PIN_PATHS = ["Lean/Screen/SeamCurrentEdge30Moment.lean",
             "Lean/ObserverPatchHolography/CoreAxioms.lean",
             "Lean/Screen/WhitneySpatialConsistency.lean",
             "code/electromagnetism/whitney_spatial_consistency.py",
             "code/electromagnetism/verify_whitney_spatial_consistency.py",
             "code/electromagnetism/test_whitney_spatial_consistency.py"]
EXPECTED_PARAMETERS = {"charge":0.7,"mass_squared":0.4,"quartic":0.2,
                       "time_interval":[0.0,0.5],"levels":[0,1,2]}
COUNT_KEYS = ("level","vertices","tetrahedra")
NUMERIC_KEYS = ("diameter","max_diameter_inradius","volume","action","first_variation",
                "continuum_action","continuum_first_variation","action_error",
                "first_variation_error","scalar_l2_error","covariant_gradient_l2_error")


def require(condition,message):
    if not condition:
        raise ValueError(message)


def load(path=OUTPUT):
    def pairs(items):
        result={}
        for key,value in items:
            require(key not in result,"duplicate JSON key")
            result[key]=value
        return result
    def constant(_):
        raise ValueError("nonfinite JSON constant")
    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=pairs,parse_constant=constant)


def polynomial_data():
    source = ROOT/"code/electromagnetism/whitney_spatial_consistency.py"
    module = ast.parse(source.read_text(encoding="utf-8"))
    values = [ast.literal_eval(node.value) for node in module.body if isinstance(node,ast.Assign)
              and any(isinstance(target,ast.Name) and target.id == "POLYNOMIALS" for target in node.targets)]
    require(len(values)==1,"one literal manufactured-field table")
    # JSON normalizes tuples exactly as in the written receipt.
    return json.loads(json.dumps(values[0],allow_nan=False))


def independent_meshes():
    source = (ROOT/"Lean/Screen/SeamCurrentEdge30Moment.lean").read_text(encoding="utf-8")
    source = source.split("def portVector",1)[1].split("theorem portVector_positivePort",1)[0]
    golden=(1+np.sqrt(5.))/2
    dictionary={"0":0.,"1":1.,"-1":-1.,"φ":golden,"-φ":-golden}
    boundary=np.array([[dictionary[item.strip()] for item in row.split(",")]
                       for row in re.findall(r"!\[([^\[\]]+)\]",source)])
    require(boundary.shape==(12,3),"twelve source vertices")
    neighbours={tuple(pair) for pair in combinations(range(12),2)
                if abs(np.dot(boundary[pair[0]]-boundary[pair[1]],boundary[pair[0]]-boundary[pair[1]])-4)<1e-12}
    faces={face for face in combinations(range(12),3) if all(pair in neighbours for pair in combinations(face,2))}
    source=(ROOT/"Lean/ObserverPatchHolography/CoreAxioms.lean").read_text(encoding="utf-8")
    source=source.split("def orientedFaces",1)[1].split("def faceEdges",1)[0]
    oriented=[tuple(map(int,row)) for row in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)",source)]
    require(len(neighbours)==30 and len(faces)==20 and {tuple(sorted(f)) for f in oriented}==faces,
            "face table matches geometric adjacency")
    vertices=np.vstack((np.zeros(3),boundary))
    cells=np.array([(0,)+tuple(v+1 for v in f) for f in oriented],dtype=int)
    for level in range(3):
        yield vertices,cells
        if level==2:
            break
        keys=sorted({tuple(sorted(edge)) for cell in cells for edge in combinations(cell.tolist(),2)})
        mid={edge:len(vertices)+i for i,edge in enumerate(keys)}
        new_vertices=np.vstack((vertices,[(vertices[a]+vertices[b])/2 for a,b in keys]))
        new_cells=[]
        for cell in cells.tolist():
            def middle(i,j):
                return mid[tuple(sorted((cell[i],cell[j])))]
            for i in range(4):
                new_cells.append([cell[i]]+[middle(i,j) for j in range(4) if j!=i])
            axis=(middle(0,1),middle(2,3))
            ring=[middle(0,2),middle(0,3),middle(1,3),middle(1,2)]
            for i in range(4):
                new_cells.append([axis[0],axis[1],ring[i],ring[(i+1)%4]])
        # Keep each parent's local ordering identical to its declared scheme
        # on the next refinement; the cyclic construction above establishes
        # the cell sets, then this orientation-independent key selects the
        # declared endpoint/midpoint order.
        ordered=[]
        for a,b,c,d in cells.tolist():
            ab,ac,ad,bc,bd,cd=[mid[tuple(sorted(edge))] for edge in ((a,b),(a,c),(a,d),(b,c),(b,d),(c,d))]
            ordered.extend([(a,ab,ac,ad),(b,ab,bc,bd),(c,ac,bc,cd),(d,ad,bd,cd),
                            (ab,ac,ad,cd),(ab,ac,bc,cd),(ab,ad,bd,cd),(ab,bc,bd,cd)])
        require(Counter(tuple(sorted(c)) for c in new_cells)==Counter(tuple(sorted(c)) for c in ordered),
                "independent octahedral partition")
        vertices,cells=new_vertices,np.array(ordered,dtype=int)


def geometry_summary(vertices,cells,level):
    xyz=vertices[cells]
    rays=xyz[:,1:]-xyz[:,:1]
    volumes=abs(np.linalg.det(rays))/6
    require(np.min(volumes)>0,"positive refined volumes")
    counters=Counter(tuple(sorted(face)) for cell in cells for face in combinations(cell.tolist(),3))
    require(all(value in (1,2) for value in counters.values()),"no duplicate or over-shared faces")
    require(sum(value==1 for value in counters.values())==20*4**level,"conforming boundary face census")
    for start in range(0,len(cells),8):
        if level:
            require(np.ptp(volumes[start:start+8])<1e-12,"equal-volume midpoint children")
    longest=np.max(np.stack([np.linalg.norm(xyz[:,i]-xyz[:,j],axis=-1) for i,j in combinations(range(4),2)]),axis=0)
    surfaces=np.zeros(len(cells))
    for i,j,k in combinations(range(4),3):
        surfaces += np.linalg.norm(np.cross(xyz[:,j]-xyz[:,i],xyz[:,k]-xyz[:,i]),axis=-1)/2
    ratio=longest*surfaces/(3*volumes)
    require(abs(volumes.sum()-20*(3+np.sqrt(5))/6)<1e-11,"same cone volume")
    return {"level":level,"vertices":len(vertices),"tetrahedra":len(cells),
            "diameter":float(longest.max()),"max_diameter_inradius":float(ratio.max()),
            "volume":float(volumes.sum())}


def jacobi_quadrature(order=7):
    nodes=[]
    weights=[]
    for power in (2,1,0):
        root,weight=roots_jacobi(order,power,0)
        nodes.append((root+1)/2)
        weights.append(weight/2**(power+1))
    u,v,w=np.meshgrid(*nodes,indexing="ij")
    first,second,third=np.meshgrid(*weights,indexing="ij")
    lam=np.stack((u,(1-u)*v,(1-u)*(1-v)*w,(1-u)*(1-v)*(1-w)),axis=-1).reshape(-1,4)
    return lam,(first*second*third).ravel()


def symbolic_fields(data):
    time,x,y,z,amplitude=sp.symbols("time x y z amplitude",real=True)
    coordinates=(time,x,y,z)
    symbolic={name:sum(sp.Rational(str(row[0]))*sp.prod(coordinate**power for coordinate,power in zip(coordinates,row[1:],strict=True))
                       for row in rows) for name,rows in data.items()}
    A=[symbolic["A"+axis]+amplitude*symbolic["B"+axis] for axis in "xyz"]
    phi=symbolic["phi"]+amplitude*symbolic["eta"]
    psi=symbolic["u"]+sp.I*symbolic["v"]+amplitude*(symbolic["zu"]+sp.I*symbolic["zv"])
    electric=[-sp.diff(component,time)-sp.diff(phi,coordinate) for component,coordinate in zip(A,(x,y,z),strict=True)]
    magnetic=[sp.diff(A[2],y)-sp.diff(A[1],z),sp.diff(A[0],z)-sp.diff(A[2],x),sp.diff(A[1],x)-sp.diff(A[0],y)]
    e=sp.Rational("0.7")
    Q=sp.diff(psi,time)+sp.I*e*phi*psi
    P=[sp.diff(psi,coordinate)-sp.I*e*component*psi for coordinate,component in zip((x,y,z),A,strict=True)]
    def make(expressions):
        function=sp.lambdify((time,x,y,z,amplitude),expressions,"numpy",cse=True)
        def evaluate(at,points,shift):
            raw=function(at,*np.moveaxis(points,-1,0),shift)
            return np.stack([np.broadcast_to(value,points.shape[:-1]) for value in raw],axis=-1)
        return evaluate
    return {"edges":make(A+[sp.diff(component,time) for component in A]),
            "nodes":make([phi,psi,sp.diff(psi,time)]),
            "continuum":make(electric+magnetic+[psi,Q]+P)}


def density(E,B,psi,Q,P):
    return (.5*np.sum(abs(E)**2,axis=-1)-.5*np.sum(abs(B)**2,axis=-1)+abs(Q)**2
            -np.sum(abs(P)**2,axis=-1)-.4*abs(psi)**2-.1*abs(psi)**4)


def independent_integrals(vertices,cells,evaluators,amplitude=0.):
    lam,weight=jacobi_quadrature()
    time_nodes,time_weights=roots_jacobi(5,0,0)
    time_nodes,time_weights=(time_nodes+1)/4,time_weights/4
    # Simpson's rule integrates the degree-two manufactured edge data exactly.
    line_nodes=np.array([0.,.5,1.])
    line_weights=np.array([1.,4.,1.])/6
    totals=np.zeros(4)
    for start in range(0,len(cells),64):
        xyz=vertices[cells[start:start+64]]
        rays=xyz[:,1:]-xyz[:,:1]
        inverse=np.linalg.inv(rays)
        gradients=np.concatenate((-inverse.sum(axis=-1)[:,:,None],inverse),axis=-1).transpose(0,2,1)
        jacobian=abs(np.linalg.det(rays))
        points=np.einsum("qi,nic->nqc",lam,xyz)
        weights=jacobian[:,None]*weight[None,:]
        for time,timeweight in zip(time_nodes,time_weights,strict=True):
            matrix=np.zeros((len(xyz),4,4))
            velocity=np.zeros_like(matrix)
            for i,j in combinations(range(4),2):
                delta=xyz[:,j]-xyz[:,i]
                edge_points=xyz[:,i,None,:]+line_nodes[None,:,None]*delta[:,None,:]
                raw=evaluators["edges"](time,edge_points,amplitude)
                for target,values in ((matrix,raw[...,:3]),(velocity,raw[...,3:])):
                    value=np.sum(np.sum(values*delta[:,None,:],axis=-1)*line_weights,axis=-1)
                    target[:,i,j],target[:,j,i]=value,-value
            field=np.einsum("qi,nij,njc->nqc",lam,matrix,gradients)
            field_velocity=np.einsum("qi,nij,njc->nqc",lam,velocity,gradients)
            magnetic=np.zeros((len(xyz),3))
            for i,j in combinations(range(4),2):
                magnetic += 2*matrix[:,i,j,None]*np.cross(gradients[:,i],gradients[:,j])
            nodal=evaluators["nodes"](time,xyz,amplitude)
            phi,psi,psit=(nodal[...,i] for i in range(3))
            phi=phi.real
            phi_field=phi@lam.T
            electric=-field_velocity-np.einsum("ni,nic->nc",phi,gradients)[:,None,:]
            theta=np.einsum("nij,qj->nqi",matrix,lam)
            thetat=np.einsum("nij,qj->nqi",velocity,lam)
            gradtheta=np.einsum("nij,njc->nic",matrix,gradients)
            phased=np.exp(.7j*theta)*psi[:,None,:]
            scalar=np.sum(lam[None,:,:]*phased,axis=-1)
            scalar_t=np.sum(lam[None,:,:]*np.exp(.7j*theta)*(psit[:,None,:]+.7j*thetat*psi[:,None,:]),axis=-1)
            scalar_grad=np.einsum("nqi,nic->nqc",phased,gradients)
            scalar_grad += .7j*np.einsum("qi,nqi,nic->nqc",lam,phased,gradtheta)
            Q=scalar_t+.7j*phi_field*scalar
            P=scalar_grad-.7j*field*scalar[:,:,None]
            actual=density(electric,magnetic[:,None,:],scalar,Q,P)
            expected=evaluators["continuum"](time,points,amplitude)
            E,B,ps,Qc,Pc=expected[...,:3],expected[...,3:6],expected[...,6],expected[...,7],expected[...,8:]
            reference=density(E,B,ps,Qc,Pc)
            totals += timeweight*np.array([np.sum(weights*item) for item in
                        (actual,reference,abs(scalar-ps)**2,np.sum(abs(P-Pc)**2,axis=-1))])
    return totals


def verify(receipt):
    require(type(receipt) is dict,"receipt object")
    require(set(receipt)=={"schema","scope","parameters","polynomials","quadrature","source_pins","refinements",
                           "continuum_trajectory_claimed","observer_history_claimed","uniform_bound_certified_by_numerics"},"receipt schema keys")
    require(receipt["schema"]=="oph.whitney_spatial_consistency.v1" and receipt["scope"]==SCOPE,"scope")
    for key in ("continuum_trajectory_claimed","observer_history_claimed","uniform_bound_certified_by_numerics"):
        require(receipt[key] is False,"no numeric/physical promotion: "+key)
    require(json.dumps(receipt["parameters"],sort_keys=True,allow_nan=False)==
            json.dumps(EXPECTED_PARAMETERS,sort_keys=True,allow_nan=False),"fixed campaign parameters")
    require(json.dumps(receipt["quadrature"],sort_keys=True,allow_nan=False)==
            json.dumps({"space_duffy_legendre_order":6,"time_legendre_order":4},sort_keys=True),
            "declared quadrature")
    require(set(receipt["source_pins"])==set(PIN_PATHS),"complete source pin set")
    for path in PIN_PATHS:
        require(receipt["source_pins"][path]==hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),"source pin: "+path)
    # Canonical JSON equality distinguishes integer exponents/orders from
    # floats and booleans; Python numeric equality alone does not.
    require(json.dumps(receipt["polynomials"],sort_keys=True,allow_nan=False)==
            json.dumps(polynomial_data(),sort_keys=True,allow_nan=False),
            "manufactured fields match declared data")
    require(type(receipt["refinements"]) is list and len(receipt["refinements"])==3,"complete scheduled refinements")
    evaluators=symbolic_fields(receipt["polynomials"])
    diagnostics=[]
    for level,((vertices,cells),row) in enumerate(zip(independent_meshes(),receipt["refinements"],strict=True)):
        require(type(row) is dict and set(row)==set(COUNT_KEYS+NUMERIC_KEYS),"refinement schema")
        expected=geometry_summary(vertices,cells,level)
        for key in COUNT_KEYS:
            require(type(row[key]) is int and row[key]==expected[key],"exact refinement count: "+key)
        for key in NUMERIC_KEYS:
            require(type(row[key]) in (int,float) and np.isfinite(row[key]),"finite numeric diagnostic: "+key)
        for key in ("diameter","max_diameter_inradius","volume"):
            require(np.isclose(row[key],expected[key],atol=1e-11,rtol=1e-11),"independent geometry: "+key)
        integrals=independent_integrals(vertices,cells,evaluators)
        step=.002
        minus2,minus1,plus1,plus2=[independent_integrals(vertices,cells,evaluators,shift)[:2]
                                for shift in (-2*step,-step,step,2*step)]
        derivative=(minus2-8*minus1+8*plus1-plus2)/(12*step)
        numerical={"action":integrals[0],"continuum_action":integrals[1],
                   "first_variation":derivative[0],"continuum_first_variation":derivative[1],
                   "action_error":abs(integrals[0]-integrals[1]),
                   "first_variation_error":abs(derivative[0]-derivative[1]),
                   "scalar_l2_error":np.sqrt(integrals[2]),"covariant_gradient_l2_error":np.sqrt(integrals[3])}
        for key,value in numerical.items():
            require(np.isclose(row[key],value,atol=2e-9,rtol=2e-9),"independent quadrature/derivative: "+key)
        diagnostics.append({key:float(value) for key,value in numerical.items()})
    for key in ("action_error","first_variation_error","scalar_l2_error","covariant_gradient_l2_error"):
        require(all(diagnostics[i+1][key]<diagnostics[i][key] for i in range(2)),"recorded refinement decrease: "+key)
    return {"scope":SCOPE,"refinement_levels":[0,1,2],"tetrahedra":[20,160,1280],"vertices":[13,55,309],
            "finite_algebra_lean_declarations":5,"independent_derivative_stencil_points":4,
            "continuum_trajectory_claimed":False,"observer_history_claimed":False,
            "uniform_bound_certified_by_numerics":False,"numeric_diagnostics":diagnostics}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path",type=Path,nargs="?",default=OUTPUT)
    args=parser.parse_args()
    print(json.dumps(verify(load(args.path)),indent=2))


if __name__=="__main__":
    main()
