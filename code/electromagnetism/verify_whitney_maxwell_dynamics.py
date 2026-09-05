"""Independent geometric action, exact stability and serial provenance audit.

The dynamics producer is never imported. Source simplices and spacetime
quadrature come from the previous independent cone verifier. Linear solves
are checked through equations, including every scalar and radial variation.
The exact local bound is reconstructed from source coordinates and its LDL
identity is checked symbolically; rational intervals certify pivot signs.
"""
from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import importlib.util
import itertools
import json
from math import isqrt
from pathlib import Path
import re

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
_geometry_spec = importlib.util.spec_from_file_location("whitney_dynamics_geometry", HERE / "verify_cone_whitney_bridge.py")
if _geometry_spec is None or _geometry_spec.loader is None:
    raise RuntimeError("independent geometric verifier unavailable")
geometry = importlib.util.module_from_spec(_geometry_spec)
_geometry_spec.loader.exec_module(geometry)
OUTPUT = HERE / "runtime/whitney_maxwell_dynamics_receipt.json"
PINS = {
    "Lean/Screen/WhitneyMaxwellDynamics.lean",
    "code/electromagnetism/runtime/cone_whitney_bridge_receipt.json",
    "code/electromagnetism/verify_cone_whitney_bridge.py",
    "code/electromagnetism/whitney_maxwell_dynamics.py",
    "code/electromagnetism/verify_whitney_maxwell_dynamics.py",
    "code/electromagnetism/test_whitney_maxwell_dynamics.py",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(path=OUTPUT):
    return geometry.load(path)


def rational(value):
    require(isinstance(value, str) and re.fullmatch(r"-?\d+(?:/[1-9]\d*)?", value), "rational syntax")
    result = Q(value)
    require(str(result) == value, "canonical rational")
    require(np.isfinite(float(result)), "finite rational register")
    return result


def floats(rows):
    if isinstance(rows, list):
        return [floats(x) for x in rows]
    return float(rational(rows))


def float64_encoding(value):
    require(Q.from_float(float(value)) == value, "exact float64 output encoding")


def close(actual, expected, name, tol=1e-9):
    def typed(x):
        return all(typed(y) for y in x) if isinstance(x, list) else type(x) in (int, float)
    require(typed(actual), name + " numeric types")
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
            and np.allclose(a, b, atol=tol, rtol=tol), name)


def verify_parent_summary(recorded, fresh):
    """Compare exact census and numerical diagnostics at their own precision.

    The caller must obtain ``fresh`` by fully replaying the parent. Its
    float64 quadrature diagnostics are not bitwise portable between BLAS
    platforms; the parent's own absolute/relative tolerance is 1e-10.
    """
    counts = {"tetrahedra", "gauge_histories", "field_variations_per_history"}
    diagnostics = {"volume_action", "gauss_max_abs", "ampere_hat_max_abs"}
    require(isinstance(recorded, dict) and set(recorded) == set(fresh)
            == counts | diagnostics, "parent replay summary schema")
    for key in counts:
        require(type(recorded[key]) is int and recorded[key] == fresh[key],
                "parent replay exact census " + key)
    for key in diagnostics:
        close(recorded[key], fresh[key], "parent replay numerical " + key,
              tol=1e-10)


def exact_mv(matrix, vector):
    return [sum((int(a)*v for a, v in zip(row, vector, strict=True)), Q(0)) for row in matrix]


def certify_stability(packet, vertices, edges, tets, mass, stiffness):
    require(set(packet) == {"field", "normalization", "local_edges", "gradient_gram",
        "mass_over_volume", "stiffness_over_volume", "bindings", "ldl",
        "global_consequence", "field_bounds"}, "stability schema")
    require(packet["field"] == "a+b*sqrt(5), encoded as [a,b] rational strings", "exact field")
    require(packet["normalization"] == "each local matrix divided by its positive tetrahedron volume", "local normalization")
    pairs = list(itertools.combinations(range(4), 2))
    require(packet["local_edges"] == [list(x) for x in pairs], "local edges")
    require(all(type(x) is int for pair in packet["local_edges"] for x in pair), "local edge integer types")
    xyz = sp.Matrix([vertices[v] for v in tets[0]])
    affine = sp.ones(4, 1).row_join(xyz)
    grad = affine.inv()[1:, :].T.applyfunc(sp.simplify)
    gram = (grad*grad.T).applyfunc(sp.simplify)
    local_m, local_k = sp.zeros(6), sp.zeros(6)
    for e, (i, j) in enumerate(pairs):
        for f, (k, l) in enumerate(pairs):
            local_m[e, f] = sp.simplify(((1+(i == k))*gram[j, l]-(1+(i == l))*gram[j, k]
                -(1+(j == k))*gram[i, l]+(1+(j == l))*gram[i, k])/20)
            local_k[e, f] = sp.simplify(4*(gram[i, k]*gram[j, l]-gram[i, l]*gram[j, k]))

    def scalar(value):
        require(isinstance(value, list) and len(value) == 2, "quadratic field encoding")
        return sp.Rational(rational(value[0]))+sp.sqrt(5)*sp.Rational(rational(value[1]))

    def matrix(value, n):
        require(isinstance(value, list) and len(value) == n and all(isinstance(r, list) and len(r) == n for r in value), "exact matrix shape")
        return sp.Matrix([[scalar(x) for x in row] for row in value])

    def equal(a, b, name):
        require(a.shape == b.shape and all(sp.expand(x) == 0 for x in a-b), name)

    equal(matrix(packet["gradient_gram"], 4), gram, "source-coordinate gradient Gram")
    equal(matrix(packet["mass_over_volume"], 6), local_m, "exact local mass")
    equal(matrix(packet["stiffness_over_volume"], 6), local_k, "exact local curl energy")
    require(set(packet["ldl"]) == {"24", "48"}, "bound census")
    denominator = 10**50
    floor_root = isqrt(5*denominator**2)
    low, high = Q(floor_root, denominator), Q(floor_root+1, denominator)
    require(low*low < 5 < high*high, "rational square-root enclosure")
    for bound, cert in packet["ldl"].items():
        require(set(cert) == {"lower", "diagonal", "positive_pivots"}, "LDL schema")
        lower = matrix(cert["lower"], 6)
        require(isinstance(cert["diagonal"], list) and len(cert["diagonal"]) == 6, "pivot census")
        diagonal = sp.diag(*[scalar(x) for x in cert["diagonal"]])
        require(all(lower[i, i] == 1 and all(lower[i, j] == 0 for j in range(i+1, 6)) for i in range(6)), "unit lower triangular")
        equal(lower*diagonal*lower.T, int(bound)*local_m-local_k, "exact LDL identity")
        for aa, bb in cert["diagonal"]:
            a, b = rational(aa), rational(bb)
            require(min(a+b*low, a+b*high) > 0, "strict exact pivot sign")
        require(cert["positive_pivots"] is True, "pivot assertion")
    require(isinstance(packet["bindings"], list) and len(packet["bindings"]) == 20, "local binding census")
    assembled_m, assembled_k = np.zeros((42, 42)), np.zeros((42, 42))
    seen = set()
    for tet, binding in zip(tets, packet["bindings"], strict=True):
        require(set(binding) == {"vertices", "edge_indices", "edge_signs"}, "binding schema")
        require(binding["vertices"] == list(tet), "tetrahedron identity")
        require(all(type(x) is int for x in binding["vertices"]), "tetrahedron integer types")
        txyz = sp.Matrix([vertices[v] for v in tet])
        rays = txyz[1:, :]-sp.ones(3, 1)*txyz[0, :]
        first_rays = xyz[1:, :]-sp.ones(3, 1)*xyz[0, :]
        equal((rays*rays.T).applyfunc(sp.simplify), (first_rays*first_rays.T).applyfunc(sp.simplify), "exact simplex congruence")
        volume = sp.simplify(abs(rays.det())/6)
        require(volume > 0, "positive simplex volume")
        ids, signs = [], []
        for i, j in pairs:
            pair = (tet[i], tet[j])
            sign = 1 if pair in edges else -1
            ids.append(edges.index(pair if sign == 1 else pair[::-1])); signs.append(sign)
        require(all(type(x) is int for x in binding["edge_indices"]+binding["edge_signs"]), "binding integer types")
        require(binding["edge_indices"] == ids and binding["edge_signs"] == signs, "signed edge assembly")
        seen.update(ids)
        orientation = np.outer(signs, signs)
        assembled_m[np.ix_(ids, ids)] += float(volume)*np.array(local_m, dtype=float)*orientation
        assembled_k[np.ix_(ids, ids)] += float(volume)*np.array(local_k, dtype=float)*orientation
    require(seen == set(range(42)), "global positive definiteness coverage")
    close(assembled_m.tolist(), mass, "local-to-quadrature mass binding")
    close(assembled_k.tolist(), stiffness, "local-to-quadrature stiffness binding")
    require(packet["global_consequence"] == "K<24M; h=1/2 implies M-h^2*K/12>M/2", "global strict bound")
    require(packet["field_bounds"] == "zero-current intervals: E^T M E<=4H and B^T M2 B<=16H; no potential bound", "field coercivity scope")


class Replay:
    def __init__(self, events):
        require(isinstance(events, list), "event list")
        self.events, self.cursor, self.state, self.writers = events, 0, {}, {}

    def step(self, op, args, keys, exact=None, numeric=None):
        require(self.cursor < len(self.events), "missing event")
        e = self.events[self.cursor]
        require(set(e) == {"id", "op", "args", "reads", "parents", "writes"}, "event schema")
        require(type(e["id"]) is int and e["id"] == self.cursor, "event identity")
        require(e["op"] == op and e["args"] == args and all(type(x) is int for x in e["args"]), "operation program")
        require(set(e["reads"]) == set(keys), "read support")
        parents = set()
        for key in keys:
            require(key in self.state, "unwritten register")
            row = e["reads"][key]
            require(set(row) == {"writer", "value"} and type(row["writer"]) is int, "read witness schema")
            require(row["writer"] == self.writers[key] and rational(row["value"]) == self.state[key], "authenticated read-after-write")
            parents.add(row["writer"])
        require(e["parents"] == sorted(parents) and all(type(x) is int for x in e["parents"]), "exact parent support")
        writes = {key: rational(value) for key, value in e["writes"].items()}
        values = {key: self.state[key] for key in keys}
        if exact is not None:
            require(writes == exact(values), "exact " + op + " operation")
        else:
            require(numeric is not None, "missing semantic checker")
            numeric(values, writes)
        self.state.update(writes)
        self.writers.update({key: self.cursor for key in writes})
        self.cursor += 1
        return self.cursor-1


def window(a, phi, rho, j, d, c, q1, q2, weights, m, m2):
    e = -2*np.diff(a, axis=0)-phi[:2]@d.T
    b = a@c.T
    derivative = geometry.full_variations(a, phi, j, rho, 0.5, d, c, q1, q2, weights)
    gauss, ampere = derivative[42:].reshape(2, 13)/0.5, derivative[:42]
    return {"E": e.tolist(), "B": b.tolist(),
        "action": float(geometry.integrated_action(a, phi, j, rho, 0.5, d, c, q1, q2, weights)),
        "gauss": gauss.tolist(), "ampere_hat": ampere.tolist(),
        "gauss_max_abs": float(np.max(np.abs(gauss))), "ampere_max_abs": float(np.max(np.abs(ampere))),
        "faraday_max_abs": float(np.max(np.abs(np.diff(b, axis=0)+0.5*e@c.T)))}


def replay_execution(ex, parent, edges, d, c, q1, q2, weights, m, m2, k):
    require(set(ex) == {"gauge", "events", "decode_event_ids", "decoded", "projection", "metrics",
        "instrumented_slices", "probe_cycles", "writable_slots"}, "execution schema")
    require(type(ex["gauge"]) is bool, "gauge Boolean")
    require(all(type(ex[key]) is int and ex[key] == value for key, value in
        (("instrumented_slices", 3), ("probe_cycles", 252), ("writable_slots", 55))), "instrument scope")
    incoming = parent["executions"][0]["cone_cochains"]
    chi = [Q((u*u+3*u) % 11-5, 7) if ex["gauge"] else Q(0) for u in range(13)]
    phis = [[-2*x for x in chi], [2*x for x in chi], [Q(0)]*13]
    initial = {"h": Q(1, 2), "clock": Q(0)}
    initial.update({f"inherited_A0/{i}": rational(x) for i, x in enumerate(incoming["A"][0])})
    initial.update({f"inherited_E0/{i}": rational(x) for i, x in enumerate(incoming["E"][0])})
    initial.update({f"rho/0/{i}": rational(x) for i, x in enumerate(incoming["rho_load"][0])})
    initial.update({f"J/{n}/{i}": rational(x) for n in range(2) for i, x in enumerate(incoming["J_load"][n])})
    initial.update({f"phi_input/{n}/{u}": x for n in range(3) for u, x in enumerate(phis[n])})
    r = Replay(ex["events"])
    r.step("inputs", [], [], exact=lambda _: initial)
    def projection(v, out):
        require(set(out) == {f"projection_z/{u}" for u in range(13)} | {f"initial_E/{e}" for e in range(42)}, "projection writes")
        for value in out.values():
            float64_encoding(value)
        old = np.array([float(v[f"inherited_E0/{e}"]) for e in range(42)])
        rho = np.array([float(v[f"rho/0/{u}"]) for u in range(13)])
        z = np.array([float(out[f"projection_z/{u}"]) for u in range(13)])
        new = np.array([float(out[f"initial_E/{e}"]) for e in range(42)])
        close((d.T@m@d@z).tolist(), rho-d.T@m@old, "Gauss projection equation")
        require(abs(z.sum()) < 1e-9, "mean-zero projection")
        close(new.tolist(), old+d@z, "gradient correction")
        require(np.linalg.norm(new-old) > 0.1, "new initial condition is explicit")
    r.step("gauss_project", [], [f"inherited_E0/{e}" for e in range(42)]+[f"rho/0/{u}" for u in range(13)], numeric=projection)
    r.step("source_next", [], ["h"]+[f"rho/0/{u}" for u in range(13)]+[f"J/0/{e}" for e in range(42)],
        exact=lambda v: {f"rho/1/{u}": v[f"rho/0/{u}"]-v["h"]*x for u, x in enumerate(exact_mv(d.T, [v[f"J/0/{e}"] for e in range(42)]))})
    ids = []
    for n in range(3):
        if n == 0:
            r.step("seed_zero", [], [f"inherited_A0/{e}" for e in range(42)]+[f"phi_input/0/{u}" for u in range(13)],
                exact=lambda v: {**{f"x/{u}": v[f"phi_input/0/{u}"] for u in range(13)}, **{f"x/{13+e}": v[f"inherited_A0/{e}"] for e in range(42)}})
        elif n == 1:
            def first(v):
                grad = exact_mv(d, [v[f"d/0/{u}"] for u in range(13)])
                return {**{f"x/{u}": v[f"phi_input/1/{u}"] for u in range(13)},
                    **{f"x/{13+e}": v[f"d/0/{13+e}"]-v["h"]*(v[f"initial_E/{e}"]+grad[e]) for e in range(42)}}
            r.step("seed_one", [], ["h"]+[f"d/0/{s}" for s in range(55)]+[f"initial_E/{e}" for e in range(42)]+[f"phi_input/1/{u}" for u in range(13)], exact=first)
        else:
            def advance(v, out):
                require(set(out) == {f"x/{s}" for s in range(55)}, "advance output domain")
                require(all(out[f"x/{u}"] == v[f"phi_input/2/{u}"] for u in range(13)), "advance scalar output")
                for e in range(42):
                    float64_encoding(out[f"x/{13+e}"])
                a = np.array([[float(v[f"d/{t}/{13+e}"]) for e in range(42)] for t in (0, 1)]+
                    [[float(out[f"x/{13+e}"]) for e in range(42)]])
                phi = np.array([[float(v[f"d/{t}/{u}"]) for u in range(13)] for t in (0, 1)]+[[float(out[f"x/{u}"]) for u in range(13)]])
                j = np.array([[float(v[f"J/0/{e}"]) for e in range(42)], [0.0]*42])
                residual = geometry.full_variations(a, phi, j, np.zeros((2, 13)), float(v["h"]), d, c, q1, q2, weights)[:42]
                require(np.max(np.abs(residual)) < 1e-9, "all 42 independent variational advance equations")
            r.step("advance", [], ["h"]+[f"d/{t}/{s}" for t in (0, 1) for s in range(55)]+[f"J/0/{e}" for e in range(42)]+[f"phi_input/2/{u}" for u in range(13)], numeric=advance)
        for u in range(13):
            r.step("baseline", [n, u], [f"x/{u}"], exact=lambda v, u=u: {f"baseline/{n}/{u}": v[f"x/{u}"]})
        before = [r.state[f"x/{s}"] for s in range(55)]
        for e, endpoints in enumerate(edges):
            for side, port in enumerate(endpoints):
                x, y, base, response = f"x/{port}", f"x/{13+e}", f"baseline/{n}/{port}", f"response/{n}/{e}/{side}"
                r.step("probe", [n, e, side], [x, y], exact=lambda v: {x: (v[x]+v[y])/2, y: (v[x]+v[y])/2})
                r.step("response", [n, e, side], [x], exact=lambda v: {response: v[x]})
                r.step("feedback", [n, e, side], [base, response, "clock"], exact=lambda v: {x: v[base], y: 2*v[response]-v[base], "clock": v["clock"]+1})
                require([r.state[f"x/{s}"] for s in range(55)] == before, "every probe restores every coordinate")
        def decode(v):
            out = {f"d/{n}/{u}": v[f"baseline/{n}/{u}"] for u in range(13)}
            for e, (left, right) in enumerate(edges):
                x = 2*v[f"response/{n}/{e}/0"]-v[f"baseline/{n}/{left}"]
                require(x == 2*v[f"response/{n}/{e}/1"]-v[f"baseline/{n}/{right}"], "two-endpoint decode")
                out[f"d/{n}/{13+e}"] = x
            return out
        ids.append(r.step("decode", [n], [f"baseline/{n}/{u}" for u in range(13)]+[f"response/{n}/{e}/{s}" for e in range(42) for s in range(2)], exact=decode))
    data = [[str(r.state[f"d/{n}/{s}"]) for s in range(55)] for n in range(3)]
    require(ex["decoded"] == data and ex["decode_event_ids"] == ids, "decoded-record binding")
    arr = np.array(floats(data))
    rho = np.array([[float(r.state[f"rho/{n}/{u}"]) for u in range(13)] for n in range(2)])
    current = np.array([[float(r.state[f"J/{n}/{e}"]) for e in range(42)] for n in range(2)])
    metrics = window(arr[:, 13:], arr[:, :13], rho, current, d, c, q1, q2, weights, m, m2)
    require(set(ex["metrics"]) == set(metrics), "metric census")
    for key, value in metrics.items():
        close(ex["metrics"][key], value, "full action " + key)
    require(max(metrics["gauss_max_abs"], metrics["ampere_max_abs"], metrics["faraday_max_abs"]) < 1e-9, "68-component stationarity and Faraday")
    def public(v, out):
        expected = {"public/action": metrics["action"], "public/gauss_max": metrics["gauss_max_abs"], "public/ampere_max": metrics["ampere_max_abs"], "public/cycles": 252,
            **{f"public/E/{n}/{e}": x for n, row in enumerate(metrics["E"]) for e, x in enumerate(row)},
            **{f"public/B/{n}/{f}": x for n, row in enumerate(metrics["B"]) for f, x in enumerate(row)}}
        require(set(out) == set(expected) and out["public/cycles"] == v["clock"] == 252, "public record census")
        for key in expected:
            float64_encoding(out[key])
            close(float(out[key]), expected[key], "public " + key)
    r.step("public", [], ["h", "clock"]+[f"d/{n}/{s}" for n in range(3) for s in range(55)]+[f"rho/{n}/{u}" for n in range(2) for u in range(13)]+[f"J/{n}/{e}" for n in range(2) for e in range(42)], numeric=public)
    require(r.cursor == len(r.events), "unaccounted event")
    require(ex["projection"] == {"inherited_E0": [str(r.state[f"inherited_E0/{e}"]) for e in range(42)], "projected_E0": [str(r.state[f"initial_E/{e}"]) for e in range(42)], "z": [str(r.state[f"projection_z/{u}"]) for u in range(13)]}, "projection provenance")
    return arr, metrics


def verify(packet):
    """Validate all custody before physics; no verification-result cache."""
    require(set(packet) == {"schema", "scope", "pins", "parent_replay", "numeric_policy", "dynamics", "stability_certificate", "executions", "continuation", "scalar_controls", "nonclaims"}, "root schema")
    require(packet["schema"] == "oph.whitney_maxwell_dynamics.v1", "schema")
    require(packet["scope"] == "NUMERIC_STATIONARY_FINITE_VOLUME_MAXWELL__EXACT_CLASSICAL_READBACK__EXACT_LOCAL_STABILITY_BOUND", "scope")
    require(set(packet["pins"]) == PINS, "pin census")
    for path, sha in packet["pins"].items():
        require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == sha, "provider pin " + path)
    parent = geometry.load()
    verify_parent_summary(packet["parent_replay"], geometry.verify(parent))
    require(packet["numeric_policy"] == {"evolution": "float64 solves and mass matrices; no interval or exact-trajectory certificate", "registers": "exact rational encodings of finite numeric results; exact pair-average restoration", "stability": "exact Q(sqrt(5)) LDL certificate, independent of float eigensolvers", "atol": 1e-9, "rtol": 1e-9}, "numeric policy")
    v, be, bf = geometry.source_mesh()
    vertices = [(0, 0, 0)]+v
    edges = [(0, u+1) for u in range(12)]+[(u+1, w+1) for u, w in be]
    faces = [tuple(u+1 for u in f) for f in bf]+[(0, u+1, w+1) for u, w in be]
    tets = [(0, *(u+1 for u in f)) for f in bf]
    d = geometry.coboundary([(u,) for u in range(13)], edges)
    c = geometry.coboundary(edges, faces)
    q1, q2, weights, _, m, m2 = geometry.quadrature(vertices, edges, faces, tets)
    k = c.T@m2@c
    certify_stability(packet["stability_certificate"], vertices, edges, tets, m, k)
    require(isinstance(packet["executions"], list) and len(packet["executions"]) == 2 and [x["gauge"] for x in packet["executions"]] == [False, True], "gauge history census")
    results = [replay_execution(ex, parent, edges, d, c, q1, q2, weights, m, m2, k) for ex in packet["executions"]]
    for key in ("E", "B", "action"):
        close(results[0][1][key], results[1][1][key], "gauge invariance " + key)
    # Continuation and interpretation checks are kept below in the same verifier.
    verify_continuation(packet["continuation"], results[0][0], parent, d, c, m, m2, k)
    verify_contracts(packet)
    return {"gauge_histories": 2, "events_per_history": len(packet["executions"][0]["events"]),
        "instrumented_slices_per_history": 3, "probe_cycles_per_history": 252,
        "field_variations_per_history": 68, "continuation_slabs": 64,
        "exact_global_stiffness_bound": 24, "gauss_max_abs": results[0][1]["gauss_max_abs"],
        "ampere_max_abs": results[0][1]["ampere_max_abs"]}


def verify_continuation(tail, initial, parent, d, c, m, m2, k):
    require(set(tail) == {"scope", "slabs", "h", "A", "rho_load", "J_load", "energy", "source_work",
        "electric_mass_norm_squared", "magnetic_mass_norm_squared", "gauss_max_abs", "ampere_max_abs",
        "faraday_max_abs", "work_identity_max_abs", "source_free_energy_drift", "radial_potential_max_abs",
        "static_extension_defect_max_abs", "source_free_electric_bound_margin", "source_free_magnetic_bound_margin"}, "continuation schema")
    require(tail["scope"] == "uninstrumented numeric continuation from the first two decoded temporal-gauge slices", "continuation custody scope")
    require(type(tail["slabs"]) is int and tail["slabs"] == 64 and type(tail["h"]) is float and tail["h"] == 0.5, "continuation grid")
    for name, shape in (("A", (65, 42)), ("rho_load", (64, 13)), ("J_load", (64, 42))):
        close(tail[name], np.zeros(shape)+np.asarray(tail[name]), "finite shaped " + name)
    a, rho, j = [np.array(tail[key], dtype=float) for key in ("A", "rho_load", "J_load")]
    close(a[:2].tolist(), initial[:2, 13:], "decoded continuation seed")
    # Compare its first evolved slice to the instrumented advance as well.
    close(a[2].tolist(), initial[2, 13:], "same first advance")
    parent_sources = parent["executions"][0]["cone_cochains"]
    expected_j = np.zeros((64, 42))
    expected_j[0] = floats(parent_sources["J_load"][0])
    close(j.tolist(), expected_j, "prescribed current continuation")
    close(rho[0].tolist(), floats(parent_sources["rho_load"][0]), "initial charge continuation")
    require(np.max(np.abs(np.diff(rho, axis=0)+0.5*j[:-1]@d)) < 1e-9, "continuity at every interior knot")
    e, b = -2*np.diff(a, axis=0), a@c.T
    # Equivalent cross-endpoint energy supplies a different evaluation path
    # from the producer's corrected-kinetic/midpoint expression.
    f = m+k/24
    energy = np.array([x@f@x/2+y@m2@z/2 for x, y, z in zip(e, b[:-1], b[1:], strict=True)])
    work = -0.25*np.sum((e[:-1]+e[1:])*j[:-1], axis=1)
    electric_norm = np.einsum("ni,ij,nj->n", e, m, e)
    magnetic_norm = np.einsum("ni,ij,nj->n", b, m2, b)
    gauss = rho-e@m@d
    # Write Ampere as a second difference in potentials rather than reuse
    # the producer's magnetic averaging expression.
    acceleration = 4*(a[2:]-2*a[1:-1]+a[:-2])
    ampere = -(acceleration@f+a[1:-1]@k-j[:-1])/2
    faraday = np.diff(b, axis=0)+0.5*e@c.T
    lift = np.array(floats(parent["cochain_maps"]["P1"]))
    radial_defect = a[:, :12]-a[:, 12:]@lift[:12, :].T
    expected = {"energy": energy, "source_work": work,
        "electric_mass_norm_squared": electric_norm, "magnetic_mass_norm_squared": magnetic_norm,
        "gauss_max_abs": np.max(np.abs(gauss)), "ampere_max_abs": np.max(np.abs(ampere)),
        "faraday_max_abs": np.max(np.abs(faraday)),
        "work_identity_max_abs": np.max(np.abs(np.diff(energy)-work)),
        "source_free_energy_drift": np.max(np.abs(energy[1:]-energy[1])),
        "radial_potential_max_abs": np.max(np.abs(a[:, :12])),
        "static_extension_defect_max_abs": np.max(np.abs(radial_defect)),
        "source_free_electric_bound_margin": np.min(4*energy[1]-electric_norm[1:]),
        "source_free_magnetic_bound_margin": np.min(16*energy[1]-magnetic_norm[1:])}
    for key, value in expected.items():
        close(tail[key], value, "continuation " + key)
    for key in ("gauss_max_abs", "ampere_max_abs", "faraday_max_abs", "work_identity_max_abs", "source_free_energy_drift"):
        require(expected[key] < 1e-9, "small continuation residual " + key)
    require(min(expected["source_free_electric_bound_margin"], expected["source_free_magnetic_bound_margin"]) > 0, "positive field control margins")
    require(expected["static_extension_defect_max_abs"] > 0.1, "independent radial dynamics")


def verify_contracts(packet):
    require(packet["dynamics"] == {"h": "1/2", "F": "M+h^2*K/6", "K": "C^T*M2*C",
        "recurrence": "F*A[n+1]=(2M-2h^2*K/3)*A[n]-F*A[n-1]+h^2*J[n-1]-h*M*D*(phi[n]-phi[n-1])",
        "projection": "E0_new=E0_old+D*(D^T*M*D+ones(13,13)/13)^(-1)*(rho0-D^T*M*E0_old)",
        "energy": "H=1/2*E^T*(M-h^2*K/12)*E+1/2*midpoint^T*K*midpoint, E=-(A[n+1]-A[n])/h-D*phi[n]",
        "work": "H[n+1]-H[n]=-h/2*(E[n]+E[n+1])^T*J[n]",
        "source": "prescribed inherited covectors; current is h*J[n] delta at t[n+1]; conventional physical charge has opposite load sign"}, "action and source contract")
    require(packet["nonclaims"] == [
        "this projected and evolved history replaces, rather than repairs or relabels, the historical serial/cone history",
        "finite-element stationarity for all 42 edge and 13 scalar coordinates, not arbitrary smooth-test continuum stationarity",
        "three slices are instrumented; the 64-slab continuation has no serial execution trace",
        "sources are prescribed; no coupled charged-path or matter stationarity",
        "supplied geometry, constitutive coefficients, clock step, projection and writable classical coordinate ports",
        "no source selection, spatial refinement limit, laboratory measurement, SI calibration or new empirical prediction",
        "field stability does not bound gauge potentials; gradient histories may shear"], "interpretation boundaries")
    controls = packet["scalar_controls"]
    require(set(controls) == {"scope", "histories"} and controls["scope"] == "scalar recurrence controls; lambda=48,64 are not claimed mesh eigenvalues", "scalar control scope")
    require(isinstance(controls["histories"], list) and len(controls["histories"]) == 3, "scalar control census")
    for lam, row in zip((1, 48, 64), controls["histories"], strict=True):
        require(set(row) == {"lambda", "h", "h2lambda", "characteristic_half_trace", "A"}, "scalar control schema")
        require(rational(row["lambda"]) == lam and rational(row["h"]) == Q(1, 2), "scalar parameters")
        z = Q(lam, 4)
        ratio = (1-z/3)/(1+z/6)
        require(rational(row["h2lambda"]) == z and rational(row["characteristic_half_trace"]) == ratio, "scalar characteristic trace")
        require(isinstance(row["A"], list) and len(row["A"]) == 17, "scalar trajectory census")
        a = [rational(x) for x in row["A"]]
        require(a[:2] == [0, 1] and all(a[n+1]-2*ratio*a[n]+a[n-1] == 0 for n in range(1, 16)), "exact scalar recurrence")
        if lam == 48:
            require(a == [Q((-1)**(n+1)*n) for n in range(17)], "critical Jordan growth")
        if lam == 64:
            require(abs(a[-1]) > 1000 and ratio < -1, "supercritical growth control")


if __name__ == "__main__":
    print(json.dumps(verify(load()), sort_keys=True))
