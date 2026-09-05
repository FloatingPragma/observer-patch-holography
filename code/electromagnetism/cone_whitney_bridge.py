"""Finite decoded Maxwell history interpolated onto a supplied tetrahedral cone.

Incidence, extension maps and cochain histories are exact rationals.  Geometric
Gram matrices use analytic barycentric moments evaluated in float64; they are
not interval certificates.  Neither sourced continuum stationarity nor a
source-selected physical geometry is asserted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import sympy as sp

import verify_serial_maxwell_readout as serial_verifier

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "runtime/cone_whitney_bridge_receipt.json"
SERIAL = HERE / "runtime/serial_maxwell_readout_receipt.json"
SCHEMA = "oph.cone_whitney_bridge.v1"
SCOPE = "EXACT_COCHAIN_INTERPOLATION__NUMERIC_GEOMETRIC_ACTION_DEFECT__NO_CONTINUUM_CLOSURE"
PIN_PATHS = [
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "Lean/Screen/SeamCurrentEdge30Moment.lean",
    "Lean/Screen/SerialMaxwellReadout.lean",
    "Lean/Screen/ConeCochainBridge.lean",
    "Lean/Screen/WhitneyTimeBridge.lean",
    "code/electromagnetism/runtime/serial_maxwell_readout_receipt.json",
    "code/electromagnetism/verify_serial_maxwell_readout.py",
    "code/electromagnetism/cone_whitney_bridge.py",
    "code/electromagnetism/verify_cone_whitney_bridge.py",
    "code/electromagnetism/test_cone_whitney_bridge.py",
]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode("ascii")


def exact(matrix):
    """JSON matrix representation; all rational entries are canonical strings."""
    return [[str(v) for v in row] for row in matrix.tolist()]


def vector(matrix):
    return [str(v) for v in matrix]


def numeric(matrix):
    return np.array(matrix.tolist(), dtype=float)


def load_geometry():
    """Parse the actual committed labels, orientations and coordinate table."""
    text = (ROOT / "Lean/Screen/SeamCurrentCarrierQuotient.lean").read_text(encoding="utf-8")
    ends = []
    for name in ("seamLeft", "seamRight"):
        match = re.search(r"def " + name + r"\s*:.*?:=\s*!\[([^\]]+)\]", text, re.S)
        require(match is not None, "missing seam table")
        ends.append([int(x.strip()) for x in match.group(1).split(",")])
    edges = list(zip(*ends, strict=True))
    raw = (ROOT / "Lean/ObserverPatchHolography/CoreAxioms.lean").read_text(encoding="utf-8")
    raw = raw.split("def orientedFaces", 1)[1].split("def faceEdges", 1)[0]
    faces = [tuple(map(int, row)) for row in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)", raw)]
    require(len(edges) == 30 and len(set(edges)) == 30 and len(faces) == 20, "carrier census")
    d, c = sp.zeros(30, 12), sp.zeros(20, 30)
    for e, (u, v) in enumerate(edges):
        require(0 <= u < v < 12, "edge orientation")
        d[e, u], d[e, v] = -1, 1
    for f, (i, j, k) in enumerate(faces):
        for u, v in ((i, j), (j, k), (k, i)):
            c[f, edges.index((min(u, v), max(u, v)))] = 1 if u < v else -1
    raw = (ROOT / "Lean/Screen/SeamCurrentEdge30Moment.lean").read_text(encoding="utf-8")
    raw = raw.split("def portVector", 1)[1].split("theorem portVector_positivePort", 1)[0]
    rows = re.findall(r"!\[([^\[\]]+)\]", raw)
    phi = (1 + sp.sqrt(5)) / 2
    # Accept only the small literal vocabulary used by the committed table.
    values = {"0": sp.Integer(0), "1": sp.Integer(1), "-1": sp.Integer(-1), "φ": phi, "-φ": -phi}
    vertices = [sp.Matrix([values[x.strip()] for x in row.split(",")]) for row in rows]
    require(len(vertices) == 12 and all(len(v) == 3 for v in vertices), "coordinate census")
    require(all(sp.simplify((vertices[v]-vertices[u]).dot(vertices[v]-vertices[u])) == 4
                for u, v in edges), "edge length squared")
    determinants = [sp.simplify(sp.Matrix.hstack(*(vertices[i] for i in face)).det()) for face in faces]
    require(all(det == 3 + sp.sqrt(5) for det in determinants), "oriented tetrahedron geometry")
    require(c*d == sp.zeros(20, 12), "boundary of boundary")
    return edges, faces, d, c, vertices


def cone_maps(d, c):
    """Radial-first edges; boundary-first faces; apex is vertex zero."""
    dc, cc, bc = sp.zeros(42, 13), sp.zeros(50, 42), sp.zeros(20, 50)
    for i in range(12):
        dc[i, 0], dc[i, i+1] = -1, 1
    dc[12:42, 1:13] = d
    cc[:20, 12:42] = c
    cc[20:50, :12], cc[20:50, 12:42] = -d, sp.eye(30)
    bc[:, :20], bc[:, 20:50] = sp.eye(20), -c
    g = (d.T*d + sp.ones(12)/12).inv() - sp.ones(12)/12
    dual = (c*c.T + sp.ones(20)/20).inv() - sp.ones(20)/20
    p0 = sp.ones(1, 12).col_join(sp.zeros(12, 12))/12
    p0[1:13, :] = sp.eye(12)
    p1 = (g*d.T).col_join(sp.eye(30))
    p2 = sp.eye(20).col_join(c.T*dual)
    require(dc*p0 == p1*d, "gauge square")
    require(cc*p1 == p2*c, "curvature square")
    require(bc*p2 == sp.ones(20)/20, "top degree flux obstruction")
    require(cc*dc == sp.zeros(50, 13) and bc*cc == sp.zeros(20, 42), "cone chain")
    return {"D": dc, "C": cc, "B": bc, "P0": p0, "P1": p1, "P2": p2}


def whitney_matrices(vertices, edges, faces, tetrahedra):
    """Analytic moments, not point sampling: int lambda_i lambda_j=V(1+delta)/20.

    One-form coefficients are lambda_i grad(lambda_j)-lambda_j grad(lambda_i).
    Two-form vector proxies are twice the cyclic lambda_i grad_j cross grad_k.
    """
    m1, m2 = np.zeros((42, 42)), np.zeros((50, 50))
    local = []
    for tet in tetrahedra:
        xyz = vertices[list(tet)]
        affine = np.column_stack((np.ones(4), xyz))
        gradients = np.linalg.inv(affine)[1:, :].T
        signed_volume = float(np.linalg.det((xyz[1:] - xyz[0]).T)/6)
        require(signed_volume > 0, "positive geometric volume")
        moments = signed_volume*(np.ones((4, 4))+np.eye(4))/20
        index = {v: i for i, v in enumerate(tet)}
        ei = [e for e, endpoints in enumerate(edges) if set(endpoints) <= set(tet)]
        fi = [f for f, corners in enumerate(faces) if set(corners) <= set(tet)]
        require(len(ei) == 6 and len(fi) == 4, "tetrahedral support")
        t1, t2 = np.zeros((6, 4, 3)), np.zeros((4, 4, 3))
        for e, global_e in enumerate(ei):
            i, j = [index[v] for v in edges[global_e]]
            t1[e, i], t1[e, j] = gradients[j], -gradients[i]
        for f, global_f in enumerate(fi):
            i, j, k = [index[v] for v in faces[global_f]]
            for a, b, q in ((i, j, k), (j, k, i), (k, i, j)):
                t2[f, a] = 2*np.cross(gradients[b], gradients[q])
        m1[np.ix_(ei, ei)] += np.einsum("eic,ij,fjc->ef", t1, moments, t1)
        m2[np.ix_(fi, fi)] += np.einsum("eic,ij,fjc->ef", t2, moments, t2)
        local.append({"vertices": list(tet), "edge_indices": ei, "face_indices": fi,
                      "volume": signed_volume, "gradients": gradients.tolist(),
                      "one_coefficients": t1.tolist(), "two_coefficients": t2.tolist()})
    require(np.linalg.eigvalsh(m1).min() > 0 and np.linalg.eigvalsh(m2).min() > 0, "Gram positivity")
    return m1, m2, local


def decoded_histories(receipt):
    """The upstream verifier authenticates operations before writes are consumed."""
    serial_verifier.verify(receipt)
    histories = []
    for execution in receipt["executions"]:
        state, decoded, origins = {}, {}, []
        for event in execution["events"]:
            for key, value in event["writes"].items():
                state[key] = sp.Rational(value)
                if event["op"] == "decode" and key.startswith("d/"):
                    decoded[key] = sp.Rational(value)
            if event["op"] == "decode":
                origins.append(event["id"])
        require(len(decoded) == 126 and len(origins) == 3, "complete decoded slices")
        a = [sp.Matrix([decoded[f"d/{n}/{12+e}"] for e in range(30)]) for n in range(3)]
        phi = [sp.Matrix([decoded[f"d/{n}/{u}"] for u in range(12)]) for n in range(3)]
        rho = [sp.Matrix([state[f"rho/{n}/{u}"] for u in range(12)]) for n in range(2)]
        current = [sp.Matrix([state[f"J/{n}/{e}"] for e in range(30)]) for n in range(2)]
        histories.append({"gauge": execution["gauge"], "decode_event_ids": origins,
                          "A": a, "phi": phi, "rho": rho, "J": current,
                          "h": state["h"], "public_action": state["field_source_action"]})
    return histories


def polynomial_identity_errors(local, maps):
    """Check every local polynomial coefficient of the Whitney chain maps.

    These are numeric cross-checks of the geometric formulas, separate from
    the exact rational cochain commuting-square certificates.
    """
    dc, cc, bc = [numeric(maps[k]) for k in ("D", "C", "B")]
    derivative0, derivative1, derivative2 = [], [], []
    for t, row in enumerate(local):
        vi, ei, fi = [row[k] for k in ("vertices", "edge_indices", "face_indices")]
        grads = np.array(row["gradients"])
        one, two = [np.array(row[k]) for k in ("one_coefficients", "two_coefficients")]
        dw0 = np.einsum("ev,eic->vic", dc[np.ix_(ei, vi)], one)
        derivative0.append(float(np.max(np.abs(dw0-grads[:, None, :]))))
        curls = np.array([sum((np.cross(grads[i], w[i]) for i in range(4)), np.zeros(3)) for w in one])
        w2c = np.einsum("fe,fic->eic", cc[np.ix_(fi, ei)], two)
        derivative1.append(float(np.max(np.abs(w2c-curls[:, None, :]))))
        divergence = np.einsum("fic,ic->f", two, grads)
        derivative2.append(float(np.max(np.abs(divergence-bc[t, fi]/row["volume"]))))
    result = {"dW0_minus_W1D": max(derivative0), "dW1_minus_W2C": max(derivative1),
              "dW2_minus_W3B": max(derivative2)}
    require(max(result.values()) < 1e-10, "Whitney polynomial commuting identity")
    return result


def execution_metrics(history, d, c, maps, m1, m2, local):
    a, phi, rho, current, h = [history[k] for k in ("A", "phi", "rho", "J", "h")]
    dc, cc, bc, p0, p1, p2 = [maps[k] for k in ("D", "C", "B", "P0", "P1", "P2")]
    electric = [-(a[n+1]-a[n])/h-d*phi[n] for n in range(2)]
    magnetic = [c*v for v in a]
    ac, pc = [p1*v for v in a], [p0*v for v in phi]
    ec, mc = [p1*v for v in electric], [p2*v for v in magnetic]
    rc = [sp.zeros(1, 1).col_join(v) for v in rho]
    jc = [sp.zeros(12, 1).col_join(v) for v in current]
    for n in range(2):
        require(-(ac[n+1]-ac[n])/h-dc*pc[n] == ec[n], "electric square")
        require(mc[n+1]-mc[n]+h*cc*ec[n] == sp.zeros(50, 1), "Faraday")
        require(bc*mc[n] == sp.zeros(20, 1), "magnetic divergence")
    require(bc*mc[2] == sp.zeros(20, 1), "terminal magnetic divergence")
    require(rc[1]-rc[0]+h*dc.T*jc[0] == sp.zeros(13, 1), "impulse continuity")
    terminal_rho = rc[1]-h*dc.T*jc[1]
    source = sum(h*(jc[n].dot(ac[n+1])+rc[n].dot(pc[n])) for n in range(2))
    counting = sum(h*(electric[n].dot(electric[n])/2-magnetic[n+1].dot(magnetic[n+1])/2)
                   for n in range(2)) + source
    require(counting == history["public_action"], "public action binding")
    counting_prism = sum(h*(electric[n].dot(electric[n])/2-
        (magnetic[n].dot(magnetic[n])+magnetic[n].dot(magnetic[n+1])+
         magnetic[n+1].dot(magnetic[n+1]))/6) for n in range(2)) + source
    counting_hat = electric[1]-electric[0]+h*current[0]-h/6*c.T*(magnetic[0]+4*magnetic[1]+magnetic[2])
    # This intentionally wrong current interpolation replaces node impulses
    # by a slab density. Its action error exposes the gauge inconsistency.
    constant_current_defect = sum(h*current[n].dot(a[n]-a[n+1])/2 for n in range(2))
    hf = float(h)
    en, bn = [numeric(v).ravel() for v in ec], [numeric(v).ravel() for v in mc]
    rn, jn = [numeric(v).ravel() for v in rc], [numeric(v).ravel() for v in jc]
    dn, cn = numeric(dc), numeric(cc)
    electric_action = sum(hf*e@m1@e/2 for e in en)
    magnetic_endpoint = sum(-hf*b@m2@b/2 for b in bn[1:])
    magnetic_prism = sum(-hf*(bn[n]@m2@bn[n]+bn[n]@m2@bn[n+1]+bn[n+1]@m2@bn[n+1])/6
                         for n in range(2))
    endpoint_action = electric_action+magnetic_endpoint+float(source)
    prism_action = electric_action+magnetic_prism+float(source)
    gauss = [rn[n]-dn.T@m1@en[n] for n in range(2)]
    ampere_endpoint = m1@(en[1]-en[0])+hf*jn[0]-hf*cn.T@m2@bn[1]
    ampere_hat = m1@(en[1]-en[0])+hf*jn[0]-hf/6*cn.T@m2@(bn[0]+4*bn[1]+bn[2])
    # A degree-one vector polynomial is completely specified by its four
    # barycentric coefficients, retaining all tetrahedra, not a sample norm.
    def coefficients(fields, kind):
        out = []
        for field in fields:
            arr = numeric(field).ravel()
            out.append([np.einsum("e,eic->ic", arr[row[f"{kind}_indices"]],
                                  np.array(row["one_coefficients" if kind == "edge" else "two_coefficients"]))
                        .tolist() for row in local])
        return out
    return {
        "gauge": history["gauge"], "decode_event_ids": history["decode_event_ids"], "h": str(h),
        "source_history": {"A": [vector(v) for v in a], "phi": [vector(v) for v in phi],
                           "rho": [vector(v) for v in rho], "J": [vector(v) for v in current]},
        "cone_cochains": {"A": [vector(v) for v in ac], "phi": [vector(v) for v in pc],
                           "E": [vector(v) for v in ec], "B": [vector(v) for v in mc],
                           "rho_load": [vector(v) for v in rc], "J_load": [vector(v) for v in jc],
                           "terminal_rho": vector(terminal_rho)},
        "whitney_fields": {"E_barycentric_coefficients": coefficients(ec, "edge"),
                            "B_barycentric_coefficients": coefficients(mc, "face")},
        "exact_identities": {"electric_square": True, "Faraday": True,
                             "magnetic_divergence": True, "impulse_continuity": True},
        "action": {"counting": str(counting), "source_pairing": str(source),
                   "volume_endpoint": endpoint_action, "volume_prism": prism_action,
                   "spatial_defect": endpoint_action-float(counting),
                   "temporal_defect": prism_action-endpoint_action,
                   "electric_volume_term": electric_action,
                   "magnetic_endpoint_term": magnetic_endpoint, "magnetic_prism_term": magnetic_prism},
        "counting_temporal_control": {"prism_action": str(counting_prism),
                                      "temporal_defect": str(counting_prism-counting),
                                      "hat_residual": vector(counting_hat),
                                      "hat_residual_squared_norm": str(counting_hat.dot(counting_hat)),
                                      "ordinary_constant_current_action_defect": str(constant_current_defect)},
        "residuals": {"gauss": [v.tolist() for v in gauss],
                      "ampere_endpoint": ampere_endpoint.tolist(), "ampere_hat": ampere_hat.tolist(),
                      "gauss_max_abs": float(np.max(np.abs(gauss))),
                      "ampere_hat_max_abs": float(np.max(np.abs(ampere_hat))),
                      "all_13_scalar_and_42_edge_tests_retained": True},
    }


def build():
    edges, faces, d, c, vertices_exact = load_geometry()
    maps = cone_maps(d, c)
    vertices = np.array([[0.0]*3]+[[float(x) for x in v] for v in vertices_exact])
    cone_edges = [(0, i+1) for i in range(12)]+[(u+1, v+1) for u, v in edges]
    cone_faces = [tuple(i+1 for i in f) for f in faces]+[(0, u+1, v+1) for u, v in edges]
    tetrahedra = [(0, *(i+1 for i in face)) for face in faces]
    m1, m2, local = whitney_matrices(vertices, cone_edges, cone_faces, tetrahedra)
    histories = decoded_histories(serial_verifier.load(SERIAL))
    executions = [execution_metrics(h, d, c, maps, m1, m2, local) for h in histories]
    x, y = executions
    require(x["source_history"]["A"] != y["source_history"]["A"], "nontrivial raw gauge control")
    for key in ("E", "B"):
        require(x["cone_cochains"][key] == y["cone_cochains"][key], "exact gauge-invariant fields")
    require(x["action"] == y["action"] and x["residuals"] == y["residuals"], "gauge-invariant action and residuals")
    chi = [sp.zeros(12, 1), -histories[0]["h"]*(histories[1]["phi"][0]-histories[0]["phi"][0]), sp.zeros(12, 1)]
    for n in range(3):
        require(histories[1]["A"][n]-histories[0]["A"][n] == d*chi[n], "raw gauge potential law")
    for n in range(2):
        require(histories[1]["phi"][n]-histories[0]["phi"][n] == -(chi[n+1]-chi[n])/histories[0]["h"],
                "raw scalar gauge law")
    me = numeric(maps["P1"]).T@m1@numeric(maps["P1"])
    kb = numeric(maps["P2"]*c).T@m2@numeric(maps["P2"]*c)
    metric_payload = {"M1": m1.tolist(), "M2": m2.tolist()}
    missing = [p for p in PIN_PATHS if not (ROOT/p).is_file()]
    require(not missing, "missing implementation pins: " + ", ".join(missing))
    return {
        "schema": SCHEMA, "scope": SCOPE,
        "pins": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PIN_PATHS},
        "numeric_policy": {"exact": "rational incidence, extension maps, decoded histories, cochain identities and source pairing",
                           "metric": "float64 analytic barycentric moments; no interval certificate",
                           "reproduction": "numeric comparisons use tolerances; digests authenticate stored bytes, not cross-platform float identity",
                           "comparison_atol": 1e-10, "comparison_rtol": 1e-10},
        "mesh": {"vertices": vertices.tolist(),
                 "apex_vertex": 0, "port_vertex_map": list(range(1, 13)),
                 "vertices_exact": [["0"]*3]+[vector(v) for v in vertices_exact],
                 "boundary_edges": edges, "boundary_faces": faces, "edges": cone_edges,
                 "faces": cone_faces, "tetrahedra": tetrahedra,
                 "counts": [13, 42, 50, 20], "boundary_edge_squared": "4",
                 "tetra_signed_determinant": "3 + sqrt(5)",
                 "volume": sum(row["volume"] for row in local)},
        "cochain_maps": {"boundary_D": exact(d), "boundary_C": exact(c),
                         **{name: exact(matrix) for name, matrix in maps.items()}},
        "cochain_contract": {"gauge_square": True, "curvature_square": True,
                             "chain_complex": True, "top_degree_defect": "ones(20,20)/20",
                             "magnetic_domain": "zero-total-flux face cochains; B=C A always qualifies",
                             "radial_map_nonzero_entries": sum(v != 0 for v in maps["P1"][:12, :])},
        "whitney": {**metric_payload, "matrix_sha256": hashlib.sha256(canonical(metric_payload)).hexdigest(),
                    "tetrahedra": local, "effective_electric_mass": me.tolist(),
                    "effective_magnetic_stiffness": kb.tolist(),
                    "mass_min_eigenvalues": [float(np.linalg.eigvalsh(m1).min()), float(np.linalg.eigvalsh(m2).min())],
                    "polynomial_identity_errors": polynomial_identity_errors(local, maps),
                    "spatial_formulas": {"one": "w_ij=lambda_i grad(lambda_j)-lambda_j grad(lambda_i)",
                                         "two": "w_ijk=2*(lambda_i grad_j cross grad_k+cyclic)",
                                         "zero": "w_i=lambda_i; physical vector proxies use the supplied Euclidean metric",
                                         "moment": "integral lambda_i lambda_j=V*(1+delta_ij)/20"}},
        "time_source_contract": {"A": "continuous piecewise linear on knots 0,h,2h",
                                 "phi_and_rho": "constant on each open slab",
                                 "E": "constant on each open slab", "B": "continuous piecewise linear",
                                 "current": "boundary-DOF-supported FE covector impulses h*J_n*delta(t-(n+1)*h)",
                                 "endpoint": "terminal rho=rho_last-h*D_cone^T*J_last; retained current includes final endpoint",
                                 "source_action": "sum h*(J_load_n dot A_(n+1)+rho_load_n dot phi_n)",
                                 "gauss_residual": "rho_load-D_cone^T*M1*E; action derivative is h times this",
                                 "ampere_hat_residual": "M1*(E1-E0)+h*J_load_0-h/6*C_cone^T*M2*(B0+4*B1+B2)"},
        "executions": executions,
        "controls": {"raw_gauge_changes": True, "exact_fields_gauge_invariant": True,
                     "gauge_chi": [vector(v) for v in chi], "gauge_endpoints_zero": True,
                     "action_and_residual_gauge_invariant": True,
                     "full_volume_stationarity": False,
                     "uniform_boundary_flux_cone_derivative": ["1"]*20},
        "nonclaims": ["supplied Euclidean embedding, apex, radial Hodge extension and temporal interpolation",
                      "unit constitutive coefficients and numerical length/time units are supplied",
                      "no source-selected spatial map or physical clock; no SI calibration",
                      "nonlocal radial extension is a reconstruction, not a local dynamics selection",
                      "source covectors are not physical charge/current densities without an additional system map",
                      "reported actions are field-plus-source actions; declared Lorentz clock action is not a volume field term",
                      "Whitney conformity is tangential for one-forms and normal for two-forms, not global vector continuity",
                      "full sourced weak Maxwell residuals do not vanish; counting action is not the volume action",
                      "no refinement, convergence, Lorentzian continuum, laboratory measurement or new prediction"],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    packet = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical(packet))
    print(json.dumps({"receipt": str(args.output), "action": packet["executions"][0]["action"],
                      "residual_max": packet["executions"][0]["residuals"]["ampere_hat_max_abs"]}, sort_keys=True))
