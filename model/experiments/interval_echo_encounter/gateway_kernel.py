"""Path coefficients and a signed Schur complement for the two-load gateway."""
import hashlib
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import cg

from .gateway_source import HERE, PairSector, candidate, federation, first_paths, geometry, require


def paths():
    results = {}
    for level in range(3):
        board = federation(level)
        s,u,v,b,edges = geometry(board)
        a,c = (first_paths(board,arm,5) for arm in ((s,u),(s,v)))
        incident = [edge for edge in edges if b in edge]
        require(len(incident)==5 and all(a.get(e)==(4,2) and c.get(e)==(5,4) for e in incident),
                "gateway path multiplicities differ")
        require(all(a.get(e,(6,0))[0]>=5 and e not in c for e in edges if b not in e),
                "another receiving seam contributes to the leading contrast")
        results[str(level)] = {"source":[s,u,v],"gateway":b,"receivers":edges,
                               "through_order_five":[[a.get(e),c.get(e)] for e in edges]}
    return results


def kernel():
    protocol = (HERE/"gateway_kernel_protocol.json").read_bytes().replace(b"\r\n",b"\n")
    config = json.loads(protocol)
    require(config["energy_probes"]==[-16,-13,-10] and config["series_terms"]==48,
            "changed kernel probe menu")
    board = federation(0)
    sector = PairSector(board)
    s,u,_,b,edges = geometry(board)
    h = sector.hamiltonian(candidate.KAPPA)
    inside = np.flatnonzero(sector.adjacent)
    outside = np.flatnonzero(~sector.adjacent)
    A = h[inside,:][:,inside]
    coupling = -h[outside,:][:,inside]
    D = h[outside,:][:,outside]
    off = sparse.diags(D.diagonal())-D
    require(min(off.data)>0 and min(coupling.data)>0, "signed paths lost positivity")
    def address(edge):
        return int(np.searchsorted(inside,sector.index(*edge)))
    start, gate = address((s,u)),address((u,b))
    receive = [address(e) for e in edges]
    bare = A[:,[start]].toarray().ravel()
    require(bare[gate]==0 and not np.any(bare[receive]), "bare seam generator crosses the cell")
    rhs = coupling[:,[start]].toarray().ravel()
    rows = {}
    for E in config["energy_probes"]:
        diag = D.diagonal()-E
        inverse = D - E*sparse.eye(len(outside),format="csr")
        y,info = cg(inverse,rhs,atol=0,rtol=1e-13)
        residual = np.linalg.norm(inverse@y-rhs)
        require(info==0 and residual<1e-11,"uncertified Schur solve")
        step = rhs/diag
        summed = step.copy()
        for _ in range(1,config["series_terms"]):
            step = (off@step)/diag
            summed += step
        exact = bare-coupling.T@y
        positive_paths = bare-coupling.T@summed
        rho = np.max(np.asarray(off.sum(axis=1)).ravel()/diag)
        row_sum = float(max(np.asarray(coupling.sum(axis=0)).flat))
        tail = row_sum*max(rhs/diag)*rho**config["series_terms"]/(1-rho)
        # D >= 0 and ||QHP|| <= sqrt(10) give the column's absolute solve error.
        tolerance = np.sqrt(10)*residual/abs(E)+tail+1e-12
        require(np.max(abs(exact-positive_paths))<=tolerance, "independent paths disagree")
        bridge, direct = 1/(4*(6-E)),1/(16*(6-E)**3)
        require(-exact[gate]>=bridge-tolerance and all(-exact[j]>=direct-tolerance
                for j,e in zip(receive,edges,strict=True) if b in e),"gateway bound failed")
        rows[str(E)] = {"adjacent_column":exact.round(11).tolist(), "bridge":round(float(exact[gate]),11),
                        "receiver":exact[receive].round(11).tolist(),
                        "column_norm":round(float(np.linalg.norm(exact)),10),
                        "max_column_disagreement_below_1e_10":bool(tolerance<1e-10)}
        require(tolerance<1e-10,"inaccurate kernel column")
    return {"probe_sha256":hashlib.sha256(protocol).hexdigest(),"paths":paths(),"effective":rows}
