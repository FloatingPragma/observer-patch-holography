"""Independent simplex quadrature and full-space variational audit.

No cone producer imports. Incidences are rebuilt from oriented simplices;
Whitney fields are integrated at degree-two tetrahedron quadrature points,
and all free derivatives are evaluated by complex-step differentiation.
The serial instrument is authenticated by its existing independent replay.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import re
from fractions import Fraction

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "runtime/cone_whitney_bridge_receipt.json"
PINS = {
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
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(path=OUTPUT):
    def pairs(items):
        out = {}
        for key, value in items:
            require(key not in out, "duplicate JSON key")
            out[key] = value
        return out
    def constant(_):
        raise ValueError("nonfinite JSON constant")
    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=pairs, parse_constant=constant)


def source_mesh():
    raw = (ROOT / "Lean/Screen/SeamCurrentCarrierQuotient.lean").read_text(encoding="utf-8")
    rows = []
    for name in ("seamLeft", "seamRight"):
        body = raw.split("def " + name, 1)[1].split("![", 1)[1].split("]", 1)[0]
        rows.append([int(v) for v in body.split(",")])
    edges = list(zip(*rows, strict=True))
    raw = (ROOT / "Lean/ObserverPatchHolography/CoreAxioms.lean").read_text(encoding="utf-8")
    body = raw.split("def orientedFaces", 1)[1].split("def faceEdges", 1)[0]
    faces = [tuple(map(int, t)) for t in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)", body)]
    raw = (ROOT / "Lean/Screen/SeamCurrentEdge30Moment.lean").read_text(encoding="utf-8")
    body = raw.split("def portVector", 1)[1].split("theorem portVector_positivePort", 1)[0]
    phi = (1 + sp.sqrt(5)) / 2
    values = {"0": sp.Integer(0), "1": sp.Integer(1), "-1": sp.Integer(-1),
              "φ": phi, "-φ": -phi}
    vertices = [tuple(values[t.strip()] for t in row.split(","))
                for row in re.findall(r"!\[([^\[\]]+)\]", body)]
    require((len(vertices), len(edges), len(faces)) == (12, 30, 20), "source census")
    for u, v in edges:
        require(sp.simplify(sum((vertices[u][k]-vertices[v][k])**2 for k in range(3))) == 4,
                "source edge geometry")
    return vertices, edges, faces


def orientation(simplex, reference):
    require(set(simplex) == set(reference), "simplex support")
    p = [reference.index(v) for v in simplex]
    return (-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))


def coboundary(lower, upper):
    result = np.zeros((len(upper), len(lower)), dtype=np.int64)
    lookup = {frozenset(s): i for i, s in enumerate(lower)}
    require(len(lookup) == len(lower), "duplicate simplex")
    for row, simplex in enumerate(upper):
        for k in range(len(simplex)):
            face = simplex[:k]+simplex[k+1:]
            col = lookup[frozenset(face)]
            result[row, col] = (-1)**k * orientation(face, lower[col])
    return result


def quadrature(vertices, edges, faces, tetrahedra):
    """Four positive degree-two points; no barycentric moment assembly."""
    a, b = (5+3*np.sqrt(5))/20, (5-np.sqrt(5))/20
    barycentric = np.full((4, 4), b)
    np.fill_diagonal(barycentric, a)
    q1, q2, weights, volumes = [], [], [], []
    for tet in tetrahedra:
        xyz = np.asarray([vertices[v] for v in tet], dtype=float)
        affine = np.column_stack((np.ones(4), xyz))
        inverse = np.linalg.inv(affine)
        gradients = inverse[1:, :].T
        volume = abs(np.linalg.det(xyz[1:]-xyz[0])) / 6
        require(volume > 1e-10, "nondegenerate volume")
        volumes.append(volume)
        # Evaluate the globally oriented basis directly at each point.
        for lam in barycentric:
            one, two = np.zeros((len(edges), 3)), np.zeros((len(faces), 3))
            for e, (u, v) in enumerate(edges):
                if u in tet and v in tet:
                    i, j = tet.index(u), tet.index(v)
                    one[e] = lam[i]*gradients[j]-lam[j]*gradients[i]
            for f, (u, v, w) in enumerate(faces):
                if u in tet and v in tet and w in tet:
                    i, j, k = tet.index(u), tet.index(v), tet.index(w)
                    two[f] = 2*(lam[i]*np.cross(gradients[j], gradients[k])
                                -lam[j]*np.cross(gradients[i], gradients[k])
                                +lam[k]*np.cross(gradients[i], gradients[j]))
            q1.append(one)
            q2.append(two)
            weights.append(volume/4)
    q1, q2, weights = np.array(q1), np.array(q2), np.array(weights)
    m1 = np.einsum("q,qic,qjc->ij", weights, q1, q1)
    m2 = np.einsum("q,qic,qjc->ij", weights, q2, q2)
    require(np.linalg.eigvalsh(m1).min() > 1e-8, "positive edge Hodge")
    require(np.linalg.eigvalsh(m2).min() > 1e-8, "positive face Hodge")
    return q1, q2, weights, np.array(volumes), m1, m2


def integrated_action(a, phi, current, rho, h, d, c, q1, q2, weights):
    """Evaluate spacetime fields at independent exact-degree quadrature nodes."""
    value = 0
    for n in range(2):
        e = -(a[n+1]-a[n])/h-d@phi[n]
        efield = np.einsum("i,qic->qc", e, q1)
        value += h/2*np.einsum("q,qc,qc->", weights, efield, efield)
        for s in ((1-1/np.sqrt(3))/2, (1+1/np.sqrt(3))/2):
            b = c@((1-s)*a[n]+s*a[n+1])
            bfield = np.einsum("i,qic->qc", b, q2)
            value -= h/4*np.einsum("q,qc,qc->", weights, bfield, bfield)
        # Impulse current at the right node; scalar load is constant per slab.
        value += h*(current[n]@a[n+1]+rho[n]@phi[n])
    return value


def full_variations(a, phi, current, rho, h, d, c, q1, q2, weights):
    eps = 1e-25
    derivative = []
    for index in range(68):
        av = np.array(a, dtype=complex)
        pv = np.array(phi, dtype=complex)
        if index < 42:
            av[1, index] += 1j*eps
        else:
            pv[(index-42)//13, (index-42)%13] += 1j*eps
        derivative.append(integrated_action(av, pv, current, rho, h, d, c,
                                            q1, q2, weights).imag/eps)
    return np.array(derivative)


def close(actual, expected, name):
    def numeric_types(x):
        return all(numeric_types(y) for y in x) if isinstance(x, list) else type(x) in (int, float)
    require(numeric_types(actual), name + " numeric types")
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    require(a.shape == b.shape and np.isfinite(a).all()
            and np.allclose(a, b, atol=1e-10, rtol=1e-10), name)


def rational_matrix(data, shape):
    require(isinstance(data, list) and len(data) == shape[0], "rational matrix rows")
    require(all(isinstance(row, list) and len(row) == shape[1] for row in data), "rational matrix columns")
    require(all(isinstance(v, str) and re.fullmatch(r"-?\d+(?:/[1-9]\d*)?", v)
                for row in data for v in row), "rational matrix syntax")
    return sp.Matrix([[sp.Rational(v) for v in row] for row in data])


def as_float(m):
    return np.array(m.tolist(), dtype=float)


def authenticated_source(serial_sha):
    # Replay on every call: a stable receipt hash cannot certify that each
    # transitive provider still matches its pin after an earlier verification.
    spec = importlib.util.spec_from_file_location("cone_upstream_serial", HERE / "verify_serial_maxwell_readout.py")
    require(spec is not None and spec.loader is not None, "upstream verifier")
    upstream = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(upstream)
    packet = upstream.load()
    require(hashlib.sha256(upstream.OUTPUT.read_bytes()).hexdigest() == serial_sha, "serial custody")
    upstream.verify(packet)
    out = []
    for execution in packet["executions"]:
        state, decode, ids = {}, {}, []
        for event in execution["events"]:
            state.update(event["writes"])
            if event["op"] == "decode":
                decode.update(event["writes"])
                ids.append(event["id"])
        def row(prefix, n, size):
            return [str(Fraction(state[f"{prefix}/{n}/{i}"])) for i in range(size)]
        out.append({"gauge": execution["gauge"], "h": state["h"], "ids": ids,
                    "action": state["field_source_action"],
                    "A": [[decode[f"d/{n}/{12+e}"] for e in range(30)] for n in range(3)],
                    "phi": [[decode[f"d/{n}/{u}"] for u in range(12)] for n in range(3)],
                    "rho": [row("rho", n, 12) for n in range(2)],
                    "J": [row("J", n, 30) for n in range(2)]})
    return out


def verify(packet):
    require(set(packet) == {"schema", "scope", "pins", "numeric_policy", "mesh",
        "cochain_maps", "cochain_contract", "whitney", "time_source_contract",
        "executions", "controls", "nonclaims"}, "root schema")
    require(packet["schema"] == "oph.cone_whitney_bridge.v1", "schema")
    require(packet["scope"] == "EXACT_COCHAIN_INTERPOLATION__NUMERIC_GEOMETRIC_ACTION_DEFECT__NO_CONTINUUM_CLOSURE", "scope")
    require(set(packet["pins"]) == PINS, "pin census")
    for path, sha in packet["pins"].items():
        require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == sha, "source pin: " + path)
    require(packet["numeric_policy"] == {
        "exact": "rational incidence, extension maps, decoded histories, cochain identities and source pairing",
        "metric": "float64 analytic barycentric moments; no interval certificate",
        "reproduction": "numeric comparisons use tolerances; digests authenticate stored bytes, not cross-platform float identity",
        "comparison_atol": 1e-10, "comparison_rtol": 1e-10}, "numeric policy")
    v, be, bf = source_mesh()
    vertices = np.array([[0.0]*3]+[[float(x) for x in row] for row in v])
    edges = [(0, i+1) for i in range(12)]+[(u+1, w+1) for u, w in be]
    faces = [tuple(i+1 for i in f) for f in bf]+[(0, u+1, w+1) for u, w in be]
    tets = [(0, *(i+1 for i in f)) for f in bf]
    mesh = packet["mesh"]
    require(set(mesh) == {"vertices", "vertices_exact", "boundary_edges", "boundary_faces",
        "edges", "faces", "tetrahedra", "counts", "boundary_edge_squared",
        "tetra_signed_determinant", "volume", "apex_vertex", "port_vertex_map"}, "mesh schema")
    require(type(mesh["apex_vertex"]) is int and mesh["apex_vertex"] == 0
            and mesh["port_vertex_map"] == list(range(1, 13))
            and all(type(i) is int for i in mesh["port_vertex_map"]), "apex and port identification")
    for key, value in {"boundary_edges": be, "boundary_faces": bf, "edges": edges,
                       "faces": faces, "tetrahedra": tets}.items():
        require(mesh[key] == [list(x) for x in value], "mesh " + key)
        require(all(type(i) is int for row in mesh[key] for i in row), "mesh integer identifiers")
    require(mesh["vertices_exact"] == [["0"]*3]+[[str(x) for x in row] for row in v], "exact vertices")
    require(mesh["counts"] == [13, 42, 50, 20] and mesh["boundary_edge_squared"] == "4"
            and mesh["tetra_signed_determinant"] == "3 + sqrt(5)", "geometric classification")
    close(mesh["vertices"], vertices, "vertices")
    d, c, b = [coboundary(lo, hi) for lo, hi in
              [([(i,) for i in range(13)], edges), (edges, faces), (faces, tets)]]
    bd = coboundary([(i,) for i in range(12)], be)
    bc = coboundary(be, bf)
    expected_incidence = {"D": d, "C": c, "B": b, "boundary_D": bd, "boundary_C": bc}
    require(set(packet["cochain_maps"]) == set(expected_incidence) | {"P0", "P1", "P2"}, "map census")
    for key, matrix in expected_incidence.items():
        require(rational_matrix(packet["cochain_maps"][key], matrix.shape) == sp.Matrix(matrix), key)
    p0 = rational_matrix(packet["cochain_maps"]["P0"], (13, 12))
    p1 = rational_matrix(packet["cochain_maps"]["P1"], (42, 30))
    p2 = rational_matrix(packet["cochain_maps"]["P2"], (50, 20))
    ds, cs = sp.Matrix(bd), sp.Matrix(bc)
    # Characterize the maps by boundary trace, mean and Poisson equations;
    # do not reuse the producer's inverse formula.
    require(p0[0, :] == sp.ones(1, 12)/12 and p0[1:, :] == sp.eye(12), "scalar extension")
    radial = p1[:12, :]
    require(p1[12:, :] == sp.eye(30) and sp.ones(1, 12)*radial == sp.zeros(1, 30)
            and ds.T*ds*radial == ds.T, "radial Green characterization")
    require(p2[:20, :] == sp.eye(20) and ds.T*p2[20:, :] == sp.zeros(12, 20)
            and cs*p2[20:, :] == sp.eye(20)-sp.ones(20)/20, "face extension characterization")
    require(sp.Matrix(d)*p0 == p1*ds and sp.Matrix(c)*p1 == p2*cs, "commuting squares")
    require(not np.any(c@d) and not np.any(b@c), "chain identities")
    require(sp.Matrix(b)*p2 == sp.ones(20)/20, "nonzero-flux obstruction")
    require(packet["cochain_contract"] == {"gauge_square": True, "curvature_square": True,
        "chain_complex": True, "top_degree_defect": "ones(20,20)/20",
        "magnetic_domain": "zero-total-flux face cochains; B=C A always qualifies",
        "radial_map_nonzero_entries": 240}, "cochain contract")
    q1, q2, weights, volumes, m1, m2 = quadrature(vertices, edges, faces, tets)
    close(mesh["volume"], sum(volumes), "volume")
    w = packet["whitney"]
    require(set(w) == {"M1", "M2", "matrix_sha256", "tetrahedra", "effective_electric_mass",
                      "effective_magnetic_stiffness", "mass_min_eigenvalues", "spatial_formulas",
                      "polynomial_identity_errors"}, "Whitney schema")
    close(w["M1"], m1, "one-form Gram quadrature")
    close(w["M2"], m2, "two-form Gram quadrature")
    encoded = (json.dumps({"M1": w["M1"], "M2": w["M2"]}, sort_keys=True,
                         separators=(",", ":"), allow_nan=False)+"\n").encode("ascii")
    require(w["matrix_sha256"] == hashlib.sha256(encoded).hexdigest(), "matrix digest")
    pf, cf = as_float(p1), as_float(p2*cs)
    close(w["effective_electric_mass"], pf.T@m1@pf, "effective electric mass")
    close(w["effective_magnetic_stiffness"], cf.T@m2@cf, "effective magnetic stiffness")
    close(w["mass_min_eigenvalues"], [np.linalg.eigvalsh(m1).min(), np.linalg.eigvalsh(m2).min()], "mass positivity")
    require(w["spatial_formulas"] == {"one": "w_ij=lambda_i grad(lambda_j)-lambda_j grad(lambda_i)",
        "two": "w_ijk=2*(lambda_i grad_j cross grad_k+cyclic)",
        "zero": "w_i=lambda_i; physical vector proxies use the supplied Euclidean metric",
        "moment": "integral lambda_i lambda_j=V*(1+delta_ij)/20"}, "spatial formulas")
    require(len(w["tetrahedra"]) == 20, "local polynomial census")
    local_one, local_two = [], []
    for index, (tet, row) in enumerate(zip(tets, w["tetrahedra"], strict=True)):
        require(set(row) == {"vertices", "edge_indices", "face_indices", "volume", "gradients",
                            "one_coefficients", "two_coefficients"}, "local schema")
        ei = [i for i, e in enumerate(edges) if set(e) <= set(tet)]
        fi = [i for i, f in enumerate(faces) if set(f) <= set(tet)]
        require(row["vertices"] == list(tet) and row["edge_indices"] == ei
                and row["face_indices"] == fi, "local support")
        close(row["volume"], volumes[index], "local volume")
        grad = np.linalg.inv(np.column_stack((np.ones(4), vertices[list(tet)])))[1:, :].T
        close(row["gradients"], grad, "barycentric gradient")
        # A linear polynomial is fixed by its values at four quadrature points.
        bary = np.full((4, 4), (5-np.sqrt(5))/20)
        np.fill_diagonal(bary, (5+3*np.sqrt(5))/20)
        one = np.einsum("rq,qec->erc", np.linalg.inv(bary), q1[index*4:index*4+4, ei])
        two = np.einsum("rq,qfc->frc", np.linalg.inv(bary), q2[index*4:index*4+4, fi])
        close(row["one_coefficients"], one, "one-form coefficients")
        close(row["two_coefficients"], two, "two-form coefficients")
        local_one.append((ei, one))
        local_two.append((fi, two))
    return verify_histories(packet, d, c, b, bd, bc, p0, p1, p2, q1, q2,
                            weights, m1, m2, local_one, local_two)


def verify_histories(packet, d, c, b, bd, bc, p0, p1, p2, q1, q2,
                     weights, m1, m2, local_one, local_two):
    require(packet["time_source_contract"] == {
        "A": "continuous piecewise linear on knots 0,h,2h",
        "phi_and_rho": "constant on each open slab", "E": "constant on each open slab",
        "B": "continuous piecewise linear",
        "current": "boundary-DOF-supported FE covector impulses h*J_n*delta(t-(n+1)*h)",
        "endpoint": "terminal rho=rho_last-h*D_cone^T*J_last; retained current includes final endpoint",
        "source_action": "sum h*(J_load_n dot A_(n+1)+rho_load_n dot phi_n)",
        "gauss_residual": "rho_load-D_cone^T*M1*E; action derivative is h times this",
        "ampere_hat_residual": "M1*(E1-E0)+h*J_load_0-h/6*C_cone^T*M2*(B0+4*B1+B2)"}, "time and source contract")
    require(packet["nonclaims"] == [
        "supplied Euclidean embedding, apex, radial Hodge extension and temporal interpolation",
        "unit constitutive coefficients and numerical length/time units are supplied",
        "no source-selected spatial map or physical clock; no SI calibration",
        "nonlocal radial extension is a reconstruction, not a local dynamics selection",
        "source covectors are not physical charge/current densities without an additional system map",
        "reported actions are field-plus-source actions; declared Lorentz clock action is not a volume field term",
        "Whitney conformity is tangential for one-forms and normal for two-forms, not global vector continuity",
        "full sourced weak Maxwell residuals do not vanish; counting action is not the volume action",
        "no refinement, convergence, Lorentzian continuum, laboratory measurement or new prediction"], "interpretation boundary")
    # Differentiate reconstructed polynomial coefficients, including face jumps
    # through the global signed assembly, rather than testing cell interiors alone.
    errors = packet["whitney"]["polynomial_identity_errors"]
    require(set(errors) == {"dW0_minus_W1D", "dW1_minus_W2C", "dW2_minus_W3B"}, "polynomial identity census")
    for value in errors.values():
        require(type(value) in (int, float) and np.isfinite(value) and 0 <= value < 1e-10, "polynomial error bound")
    for t, (row, (ei, one), (fi, two)) in enumerate(zip(
            packet["whitney"]["tetrahedra"], local_one, local_two, strict=True)):
        grad = np.array(row["gradients"])
        local_d = d[np.ix_(ei, row["vertices"])]
        # Values of dW0 and curl(W1) at each barycentric vertex.
        require(np.allclose(np.einsum("eic,ev->ivc", one, local_d),
                            np.broadcast_to(grad, (4, 4, 3)), atol=1e-10, rtol=1e-10), "dW0=W1D")
        curls = sum(np.cross(grad[i], one[:, i, :]) for i in range(4))
        face_curls = np.einsum("fic,fe->eic", two, c[np.ix_(fi, ei)])
        require(np.allclose(face_curls, curls[:, None, :], atol=1e-10, rtol=1e-10), "dW1=W2C")
        divergences = np.einsum("ic,fic->f", grad, two)
        require(np.allclose(divergences, b[t, fi]/row["volume"], atol=1e-10, rtol=1e-10), "dW2=W3B")
    serial_key = "code/electromagnetism/runtime/serial_maxwell_readout_receipt.json"
    sources = authenticated_source(packet["pins"][serial_key])
    require(len(packet["executions"]) == 2, "gauge execution census")
    ds, cs = sp.Matrix(bd), sp.Matrix(bc)
    dc, cc, bb = sp.Matrix(d), sp.Matrix(c), sp.Matrix(b)
    results = []
    for ex, source in zip(packet["executions"], sources, strict=True):
        require(set(ex) == {"gauge", "decode_event_ids", "h", "source_history", "cone_cochains",
            "whitney_fields", "exact_identities", "action", "residuals", "counting_temporal_control"}, "execution schema")
        require(type(ex["gauge"]) is bool and ex["gauge"] == source["gauge"], "gauge selector")
        require(ex["decode_event_ids"] == source["ids"] and all(type(v) is int for v in ex["decode_event_ids"]), "decoded provenance")
        require(ex["h"] == source["h"] == "1/2", "time step")
        require(ex["source_history"] == {k: source[k] for k in ("A", "phi", "rho", "J")}, "actual decoded history")
        a, phi, rho, current = [[sp.Matrix([sp.Rational(v) for v in row])
            for row in source[k]] for k in ("A", "phi", "rho", "J")]
        h = sp.Rational(1, 2)
        e = [-(a[n+1]-a[n])/h-ds*phi[n] for n in range(2)]
        bf = [cs*x for x in a]
        ac, pc, ec, mc = [[p*x for x in seq] for p, seq in
                         ((p1, a), (p0, phi), (p1, e), (p2, bf))]
        rc = [sp.Matrix([0]+list(x)) for x in rho]
        jc = [sp.Matrix([0]*12+list(x)) for x in current]
        values = {"A": ac, "phi": pc, "E": ec, "B": mc, "rho_load": rc, "J_load": jc}
        expected_cochains = {k: [[str(x) for x in row] for row in seq] for k, seq in values.items()}
        expected_cochains["terminal_rho"] = [str(x) for x in rc[1]-h*dc.T*jc[1]]
        require(ex["cone_cochains"] == expected_cochains, "exact cone history")
        for n in range(2):
            require(ec[n] == -(ac[n+1]-ac[n])/h-dc*pc[n], "electric map")
            require(mc[n+1]-mc[n]+h*cc*ec[n] == sp.zeros(50, 1), "Faraday")
        require(all(bb*x == sp.zeros(20, 1) for x in mc), "no magnetic divergence")
        require(rc[1]-rc[0]+h*dc.T*jc[0] == sp.zeros(13, 1), "impulse continuity")
        require(ex["exact_identities"] == {"electric_square": True, "Faraday": True,
            "magnetic_divergence": True, "impulse_continuity": True}, "identity labels")
        anf, pnf, enf, bnf, rnf, jnf = [np.array([as_float(x).ravel() for x in seq])
                                      for seq in (ac, pc, ec, mc, rc, jc)]
        source_action = sum(h*(jc[n].dot(ac[n+1])+rc[n].dot(pc[n])) for n in range(2))
        counting = sp.Rational(source["action"])
        ev = sum(float(h)/2*x@m1@x for x in enf)
        mb = -sum(float(h)/2*x@m2@x for x in bnf[1:])
        endpoint = ev+mb+float(source_action)
        prism = integrated_action(anf, pnf, jnf, rnf, float(h), d, c, q1, q2, weights)
        action = ex["action"]
        require(set(action) == {"counting", "source_pairing", "volume_endpoint", "volume_prism",
            "spatial_defect", "temporal_defect", "electric_volume_term", "magnetic_endpoint_term",
            "magnetic_prism_term"}, "action schema")
        require(action["counting"] == str(counting) and action["source_pairing"] == str(source_action), "exact source action")
        for key, value in {"volume_endpoint": endpoint, "volume_prism": prism,
            "spatial_defect": endpoint-float(counting), "temporal_defect": prism-endpoint,
            "electric_volume_term": ev, "magnetic_endpoint_term": mb,
            "magnetic_prism_term": prism-ev-float(source_action)}.items():
            close(action[key], value, "action " + key)
        derivatives = full_variations(anf, pnf, jnf, rnf, float(h), d, c, q1, q2, weights)
        gauss = derivatives[42:].reshape(2, 13)/float(h)
        hat = derivatives[:42]
        endpoint_res = m1@(enf[1]-enf[0])+float(h)*jnf[0]-float(h)*c.T@m2@bnf[1]
        residuals = ex["residuals"]
        require(set(residuals) == {"gauss", "ampere_endpoint", "ampere_hat", "gauss_max_abs",
            "ampere_hat_max_abs", "all_13_scalar_and_42_edge_tests_retained"}, "residual schema")
        require(residuals["all_13_scalar_and_42_edge_tests_retained"] is True, "variation scope")
        for key, value in {"gauss": gauss, "ampere_hat": hat, "ampere_endpoint": endpoint_res,
            "gauss_max_abs": np.max(np.abs(gauss)), "ampere_hat_max_abs": np.max(np.abs(hat))}.items():
            close(residuals[key], value, "full action derivative " + key)
        require(np.max(np.abs(hat[:12])) > 0.01 and np.max(np.abs(gauss)) > 0.1, "nonzero full volume residual control")
        fields = ex["whitney_fields"]
        require(set(fields) == {"E_barycentric_coefficients", "B_barycentric_coefficients"}, "field schema")
        for key, seq, basis in (("E_barycentric_coefficients", enf, local_one),
                                ("B_barycentric_coefficients", bnf, local_two)):
            expected = [[np.einsum("i,irc->rc", field[indices], poly)
                         for indices, poly in basis] for field in seq]
            close(fields[key], expected, "whole polynomial field " + key)
        # Exact independent counting-metric integration and temporal residual.
        time_prism = source_action + sum(h*e[n].dot(e[n])/2-h*(bf[n].dot(bf[n])+
            bf[n].dot(bf[n+1])+bf[n+1].dot(bf[n+1]))/6 for n in range(2))
        time_hat = e[1]-e[0]+h*current[0]-h/6*cs.T*(bf[0]+4*bf[1]+bf[2])
        constant_extra = -h/2*sum(current[n].dot(a[n+1]-a[n]) for n in range(2))
        require(ex["counting_temporal_control"] == {"prism_action": str(time_prism),
            "temporal_defect": str(time_prism-counting), "hat_residual": [str(x) for x in time_hat],
            "hat_residual_squared_norm": str(time_hat.dot(time_hat)),
            "ordinary_constant_current_action_defect": str(constant_extra)}, "counting temporal control")
        results.append({"a": a, "phi": phi, "ec": ec, "mc": mc, "action": prism,
                        "gauss": gauss, "hat": hat})
    chi = [sp.zeros(12, 1), -sp.Rational(1, 2)*(results[1]["phi"][0]-results[0]["phi"][0]), sp.zeros(12, 1)]
    for n in range(3):
        require(results[1]["a"][n]-results[0]["a"][n] == ds*chi[n], "actual vector gauge law")
    for n in range(2):
        require(results[1]["phi"][n]-results[0]["phi"][n] == -2*(chi[n+1]-chi[n]), "actual scalar gauge law")
    require(packet["controls"] == {"raw_gauge_changes": True, "exact_fields_gauge_invariant": True,
        "gauge_chi": [[str(x) for x in row] for row in chi], "gauge_endpoints_zero": True,
        "action_and_residual_gauge_invariant": True, "full_volume_stationarity": False,
        "uniform_boundary_flux_cone_derivative": ["1"]*20}, "control contract")
    require(chi[1] != sp.zeros(12, 1), "nontrivial gauge")
    require(results[0]["ec"] == results[1]["ec"] and results[0]["mc"] == results[1]["mc"], "exact field gauge invariance")
    for key in ("action", "gauss", "hat"):
        require(np.allclose(results[0][key], results[1][key], atol=1e-10, rtol=1e-10), "gauge " + key)
    return {"tetrahedra": 20, "gauge_histories": 2, "field_variations_per_history": 68,
            "volume_action": float(results[0]["action"]),
            "gauss_max_abs": float(np.max(np.abs(results[0]["gauss"]))),
            "ampere_hat_max_abs": float(np.max(np.abs(results[0]["hat"])))}


if __name__ == "__main__":
    print(json.dumps(verify(load()), sort_keys=True))
