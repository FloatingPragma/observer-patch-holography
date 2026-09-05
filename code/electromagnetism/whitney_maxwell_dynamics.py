"""Stationary finite-volume Whitney history with exact classical readback.

The evolution and metric projection use float64 linear solves. Register values
are exact rational encodings of their finite numeric results; pair-average
feedback is exact in those registers. A separate Q(sqrt(5)) local certificate
proves a strict global stiffness bound for the supplied geometric mesh.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import sympy as sp

import verify_cone_whitney_bridge as parent_verifier

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PARENT = HERE / "runtime/cone_whitney_bridge_receipt.json"
OUTPUT = HERE / "runtime/whitney_maxwell_dynamics_receipt.json"
SCHEMA = "oph.whitney_maxwell_dynamics.v1"
PIN_PATHS = [
    "Lean/Screen/WhitneyMaxwellDynamics.lean",
    "code/electromagnetism/runtime/cone_whitney_bridge_receipt.json",
    "code/electromagnetism/verify_cone_whitney_bridge.py",
    "code/electromagnetism/whitney_maxwell_dynamics.py",
    "code/electromagnetism/verify_whitney_maxwell_dynamics.py",
    "code/electromagnetism/test_whitney_maxwell_dynamics.py",
]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n").encode("ascii")


def qfloat(value):
    value = float(value)
    require(np.isfinite(value), "nonfinite solver result")
    return Q.from_float(value)


@dataclass(frozen=True)
class Quadratic:
    """Exact a+b*sqrt(5); all arithmetic and positivity tests are rational."""
    a: Q = Q(0)
    b: Q = Q(0)

    @staticmethod
    def cast(x):
        return x if isinstance(x, Quadratic) else Quadratic(Q(x), Q(0))

    def __add__(self, other):
        other = self.cast(other)
        return Quadratic(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Quadratic(-self.a, -self.b)

    def __sub__(self, other):
        return self + -self.cast(other)

    def __rsub__(self, other):
        return self.cast(other) + -self

    def __mul__(self, other):
        other = self.cast(other)
        return Quadratic(self.a*other.a+5*self.b*other.b,
                         self.a*other.b+self.b*other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        norm = other.a**2-5*other.b**2
        require(norm != 0, "zero quadratic divisor")
        return self*Quadratic(other.a/norm, -other.b/norm)

    def __rtruediv__(self, other):
        return self.cast(other)/self

    def positive(self):
        if self.b == 0:
            return self.a > 0
        if self.a >= 0 and self.b > 0:
            return True
        if self.a <= 0 and self.b < 0:
            return False
        if self.a > 0:
            return self.a**2 > 5*self.b**2
        return 5*self.b**2 > self.a**2

    def encoded(self):
        return [str(self.a), str(self.b)]


def exact_ldl(matrix):
    n = len(matrix)
    lower = [[Quadratic.cast(int(i == j)) for j in range(n)] for i in range(n)]
    diagonal = []
    for i in range(n):
        pivot = matrix[i][i]-sum(lower[i][k]*lower[i][k]*diagonal[k] for k in range(i))
        require(pivot.positive(), "nonpositive exact LDL pivot")
        diagonal.append(pivot)
        for j in range(i+1, n):
            lower[j][i] = (matrix[j][i]-sum(lower[j][k]*lower[i][k]*diagonal[k]
                                          for k in range(i)))/pivot
    for i in range(n):
        for j in range(n):
            require(sum(lower[i][k]*diagonal[k]*lower[j][k] for k in range(n)) == matrix[i][j], "LDL identity")
    return {"lower": [[x.encoded() for x in row] for row in lower],
            "diagonal": [x.encoded() for x in diagonal], "positive_pivots": True}


def local_stability_certificate(parent):
    """The positive volume factor is divided out of both local forms."""
    phi = Quadratic(Q(1, 2), Q(1, 2))
    # The three rays from the central apex have Gram matrix 2I+phi*11^T.
    inv = [[Quadratic.cast(Q(int(i == j), 2))-phi/(2*(2+3*phi))
            for j in range(3)] for i in range(3)]
    g = [[Quadratic() for _ in range(4)] for _ in range(4)]
    for i in range(3):
        for j in range(3):
            g[i+1][j+1] = inv[i][j]
    for i in range(1, 4):
        g[0][i] = g[i][0] = -sum(g[j][i] for j in range(1, 4))
    g[0][0] = sum(g[i][j] for i in range(1, 4) for j in range(1, 4))
    pairs = list(itertools.combinations(range(4), 2))
    mass, stiffness = [], []
    for i, j in pairs:
        mr, kr = [], []
        for k, l in pairs:
            mr.append(((1+int(i == k))*g[j][l]-(1+int(i == l))*g[j][k]
                       -(1+int(j == k))*g[i][l]+(1+int(j == l))*g[i][k])/20)
            kr.append(4*(g[i][k]*g[j][l]-g[i][l]*g[j][k]))
        mass.append(mr)
        stiffness.append(kr)
    bindings = []
    edges = [tuple(e) for e in parent["mesh"]["edges"]]
    symbolic = [sp.Matrix([sp.sympify(x) for x in row]) for row in parent["mesh"]["vertices_exact"]]
    for tet in parent["mesh"]["tetrahedra"]:
        rays = [symbolic[tet[i]]-symbolic[tet[0]] for i in range(1, 4)]
        for i in range(3):
            for j in range(3):
                require(sp.simplify(rays[i].dot(rays[j])-(2*int(i == j)+(1+sp.sqrt(5))/2)) == 0,
                        "exact tetrahedral congruence")
        ids, signs = [], []
        for i, j in pairs:
            edge = (tet[i], tet[j])
            if edge in edges:
                ids.append(edges.index(edge)); signs.append(1)
            else:
                ids.append(edges.index(edge[::-1])); signs.append(-1)
        bindings.append({"vertices": tet, "edge_indices": ids, "edge_signs": signs})
    certificates = {}
    for bound in (24, 48):
        form = [[bound*mass[i][j]-stiffness[i][j] for j in range(6)] for i in range(6)]
        certificates[str(bound)] = exact_ldl(form)
    return {"field": "a+b*sqrt(5), encoded as [a,b] rational strings",
            "normalization": "each local matrix divided by its positive tetrahedron volume",
            "local_edges": [list(p) for p in pairs],
            "gradient_gram": [[x.encoded() for x in row] for row in g],
            "mass_over_volume": [[x.encoded() for x in row] for row in mass],
            "stiffness_over_volume": [[x.encoded() for x in row] for row in stiffness],
            "bindings": bindings, "ldl": certificates,
            "global_consequence": "K<24M; h=1/2 implies M-h^2*K/12>M/2",
            "field_bounds": "zero-current intervals: E^T M E<=4H and B^T M2 B<=16H; no potential bound"}


def load_parent():
    packet = parent_verifier.load(PARENT)
    # Always replay; the parent verifier freshly checks all transitive source
    # bytes and authenticates the original serial execution again.
    result = parent_verifier.verify(packet)
    d = np.array(packet["cochain_maps"]["D"], dtype=float)
    c = np.array(packet["cochain_maps"]["C"], dtype=float)
    m = np.array(packet["whitney"]["M1"], dtype=float)
    m2 = np.array(packet["whitney"]["M2"], dtype=float)
    return packet, d, c, m, m2, c.T@m2@c, result


def exact_incidence_mv(matrix, values):
    return [sum((int(x)*v for x, v in zip(row, values, strict=True)), Q(0)) for row in matrix]


class Recorder:
    def __init__(self):
        self.state, self.writer, self.events = {}, {}, []

    def emit(self, op, args, keys, calculate):
        require(len(keys) == len(set(keys)), "duplicate read key")
        reads = {key: {"writer": self.writer[key], "value": str(self.state[key])} for key in keys}
        writes = {key: Q(value) for key, value in calculate({k: self.state[k] for k in keys}).items()}
        event = {"id": len(self.events), "op": op, "args": args, "reads": reads,
                 "parents": sorted({row["writer"] for row in reads.values()}),
                 "writes": {key: str(value) for key, value in writes.items()}}
        self.events.append(event)
        self.state.update(writes)
        self.writer.update({key: event["id"] for key in writes})


def window_metrics(a, phi, rho, current, h, d, c, m, m2):
    electric = np.array([-(a[n+1]-a[n])/h-d@phi[n] for n in range(2)])
    magnetic = np.array([c@v for v in a])
    source = sum(h*(current[n]@a[n+1]+rho[n]@phi[n]) for n in range(2))
    action = source+sum(h*electric[n]@m@electric[n]/2-h/6*(magnetic[n]@m2@magnetic[n]+
        magnetic[n]@m2@magnetic[n+1]+magnetic[n+1]@m2@magnetic[n+1]) for n in range(2))
    gauss = np.array([rho[n]-d.T@m@electric[n] for n in range(2)])
    ampere = m@(electric[1]-electric[0])+h*current[0]-h/6*c.T@m2@(magnetic[0]+4*magnetic[1]+magnetic[2])
    return {"E": electric.tolist(), "B": magnetic.tolist(), "action": float(action),
            "gauss": gauss.tolist(), "ampere_hat": ampere.tolist(),
            "gauss_max_abs": float(np.max(np.abs(gauss))), "ampere_max_abs": float(np.max(np.abs(ampere))),
            "faraday_max_abs": float(np.max(np.abs(magnetic[1:]-magnetic[:-1]+h*electric@c.T)))}


def instrument(parent, d, c, m, m2, k, gauge=False):
    h = Q(1, 2)
    incoming = parent["executions"][0]["cone_cochains"]
    chi = [Q((u*u+3*u) % 11-5, 7) if gauge else Q(0) for u in range(13)]
    phi = [[-v/h for v in chi], [v/h for v in chi], [Q(0)]*13]
    rec = Recorder()
    inputs = {"h": h, "clock": Q(0)}
    inputs.update({f"inherited_A0/{e}": Q(v) for e, v in enumerate(incoming["A"][0])})
    inputs.update({f"inherited_E0/{e}": Q(v) for e, v in enumerate(incoming["E"][0])})
    inputs.update({f"rho/0/{u}": Q(v) for u, v in enumerate(incoming["rho_load"][0])})
    inputs.update({f"J/{n}/{e}": Q(v) for n in range(2) for e, v in enumerate(incoming["J_load"][n])})
    inputs.update({f"phi_input/{n}/{u}": v for n in range(3) for u, v in enumerate(phi[n])})
    rec.emit("inputs", [], [], lambda _: inputs)

    def project(values):
        old = np.array([float(values[f"inherited_E0/{e}"]) for e in range(42)])
        rho = np.array([float(values[f"rho/0/{u}"]) for u in range(13)])
        defect = rho-d.T@m@old
        require(abs(defect.sum()) < 1e-10, "compatible Gauss load")
        z = np.linalg.solve(d.T@m@d+np.ones((13, 13))/13, defect)
        out = {f"projection_z/{u}": qfloat(v) for u, v in enumerate(z)}
        out.update({f"initial_E/{e}": qfloat(v) for e, v in enumerate(old+d@z)})
        return out
    rec.emit("gauss_project", [], [f"inherited_E0/{e}" for e in range(42)]+
             [f"rho/0/{u}" for u in range(13)], project)
    def source_next(values):
        boundary = exact_incidence_mv(d.T, [values[f"J/0/{e}"] for e in range(42)])
        return {f"rho/1/{u}": values[f"rho/0/{u}"]-values["h"]*boundary[u] for u in range(13)}
    rec.emit("source_next", [], ["h"]+[f"rho/0/{u}" for u in range(13)]+[f"J/0/{e}" for e in range(42)], source_next)
    decoded_ids = []
    edges = parent["mesh"]["edges"]
    for n in range(3):
        if n == 0:
            keys = [f"inherited_A0/{e}" for e in range(42)]+[f"phi_input/0/{u}" for u in range(13)]
            rec.emit("seed_zero", [], keys, lambda v: {
                **{f"x/{u}": v[f"phi_input/0/{u}"] for u in range(13)},
                **{f"x/{13+e}": v[f"inherited_A0/{e}"] for e in range(42)}})
        elif n == 1:
            keys = ["h"]+[f"d/0/{s}" for s in range(55)]+[f"initial_E/{e}" for e in range(42)]+[f"phi_input/1/{u}" for u in range(13)]
            def seed_one(v):
                grad = exact_incidence_mv(d, [v[f"d/0/{u}"] for u in range(13)])
                return {**{f"x/{u}": v[f"phi_input/1/{u}"] for u in range(13)},
                        **{f"x/{13+e}": v[f"d/0/{13+e}"]-v["h"]*(v[f"initial_E/{e}"]+grad[e]) for e in range(42)}}
            rec.emit("seed_one", [], keys, seed_one)
        else:
            keys = ["h"]+[f"d/{t}/{s}" for t in (0, 1) for s in range(55)]+[f"J/0/{e}" for e in range(42)]+[f"phi_input/2/{u}" for u in range(13)]
            def advance(v):
                aa = [np.array([float(v[f"d/{t}/{13+e}"]) for e in range(42)]) for t in (0, 1)]
                dp = np.array([float(v[f"d/1/{u}"]-v[f"d/0/{u}"]) for u in range(13)])
                j = np.array([float(v[f"J/0/{e}"]) for e in range(42)])
                f = m+float(v["h"]**2)/6*k
                b = 2*m-2*float(v["h"]**2)/3*k
                out = np.linalg.solve(f, b@aa[1]-f@aa[0]+float(v["h"]**2)*j-float(v["h"])*m@d@dp)
                return {**{f"x/{u}": v[f"phi_input/2/{u}"] for u in range(13)},
                        **{f"x/{13+e}": qfloat(vv) for e, vv in enumerate(out)}}
            rec.emit("advance", [], keys, advance)
        for u in range(13):
            rec.emit("baseline", [n, u], [f"x/{u}"], lambda v, u=u, n=n: {f"baseline/{n}/{u}": v[f"x/{u}"]})
        for e, endpoints in enumerate(edges):
            for side, port in enumerate(endpoints):
                xkey, akey = f"x/{port}", f"x/{13+e}"
                response = f"response/{n}/{e}/{side}"
                baseline = f"baseline/{n}/{port}"
                rec.emit("probe", [n, e, side], [xkey, akey], lambda v, x=xkey, a=akey: {x: (v[x]+v[a])/2, a: (v[x]+v[a])/2})
                rec.emit("response", [n, e, side], [xkey], lambda v, x=xkey, r=response: {r: v[x]})
                rec.emit("feedback", [n, e, side], [baseline, response, "clock"], lambda v, base=baseline, r=response, x=xkey, a=akey: {x: v[base], a: 2*v[r]-v[base], "clock": v["clock"]+1})
        keys = [f"baseline/{n}/{u}" for u in range(13)]+[f"response/{n}/{e}/{side}" for e in range(42) for side in range(2)]
        def decode(v, n=n):
            out = {f"d/{n}/{u}": v[f"baseline/{n}/{u}"] for u in range(13)}
            for e, (left, right) in enumerate(edges):
                value = 2*v[f"response/{n}/{e}/0"]-v[f"baseline/{n}/{left}"]
                require(value == 2*v[f"response/{n}/{e}/1"]-v[f"baseline/{n}/{right}"], "endpoint consistency")
                out[f"d/{n}/{13+e}"] = value
            return out
        rec.emit("decode", [n], keys, decode)
        decoded_ids.append(len(rec.events)-1)
        require(all(rec.state[f"x/{s}"] == rec.state[f"d/{n}/{s}"] for s in range(55)), "exact serial restoration")
    keys = ["h", "clock"]+[f"d/{n}/{s}" for n in range(3) for s in range(55)]+[f"rho/{n}/{u}" for n in range(2) for u in range(13)]+[f"J/{n}/{e}" for n in range(2) for e in range(42)]
    def public(v):
        a = np.array([[float(v[f"d/{n}/{13+e}"]) for e in range(42)] for n in range(3)])
        ph = np.array([[float(v[f"d/{n}/{u}"]) for u in range(13)] for n in range(3)])
        rho = np.array([[float(v[f"rho/{n}/{u}"]) for u in range(13)] for n in range(2)])
        j = np.array([[float(v[f"J/{n}/{e}"]) for e in range(42)] for n in range(2)])
        metrics = window_metrics(a, ph, rho, j, float(v["h"]), d, c, m, m2)
        return {"public/action": qfloat(metrics["action"]), "public/gauss_max": qfloat(metrics["gauss_max_abs"]),
                "public/ampere_max": qfloat(metrics["ampere_max_abs"]), "public/cycles": v["clock"],
                **{f"public/E/{n}/{e}": qfloat(x) for n, row in enumerate(metrics["E"]) for e, x in enumerate(row)},
                **{f"public/B/{n}/{ff}": qfloat(x) for n, row in enumerate(metrics["B"]) for ff, x in enumerate(row)}}
    rec.emit("public", [], keys, public)
    s = rec.state
    a = np.array([[float(s[f"d/{n}/{13+e}"]) for e in range(42)] for n in range(3)])
    ph = np.array([[float(s[f"d/{n}/{u}"]) for u in range(13)] for n in range(3)])
    rho = np.array([[float(s[f"rho/{n}/{u}"]) for u in range(13)] for n in range(2)])
    current = np.array([[float(s[f"J/{n}/{e}"]) for e in range(42)] for n in range(2)])
    projection = {"inherited_E0": [str(s[f"inherited_E0/{e}"]) for e in range(42)],
                  "projected_E0": [str(s[f"initial_E/{e}"]) for e in range(42)],
                  "z": [str(s[f"projection_z/{u}"]) for u in range(13)]}
    return {"gauge": gauge, "events": rec.events, "decode_event_ids": decoded_ids,
            "decoded": [[str(s[f"d/{n}/{slot}"]) for slot in range(55)] for n in range(3)],
            "projection": projection, "metrics": window_metrics(a, ph, rho, current, float(h), d, c, m, m2),
            "instrumented_slices": 3, "probe_cycles": 252, "writable_slots": 55}


def continuation(execution, parent, d, c, m, m2, k, slabs=64):
    h = 0.5
    data = np.array([[float(Q(x)) for x in row] for row in execution["decoded"]])
    a = [data[0, 13:], data[1, 13:]]
    j = np.zeros((slabs, 42))
    j[0] = np.array([float(Q(x)) for x in parent["executions"][0]["cone_cochains"]["J_load"][0]])
    rho = [np.array([float(Q(x)) for x in parent["executions"][0]["cone_cochains"]["rho_load"][0]])]
    for n in range(slabs-1):
        rho.append(rho[-1]-h*d.T@j[n])
    f, b = m+h*h*k/6, 2*m-2*h*h*k/3
    for n in range(1, slabs):
        a.append(np.linalg.solve(f, b@a[n]-f@a[n-1]+h*h*j[n-1]))
    a, rho = np.array(a), np.array(rho)
    e, magnetic = -(a[1:]-a[:-1])/h, a@c.T
    eff = m-h*h*k/12
    midpoint = (a[1:]+a[:-1])/2
    energy = np.einsum("ni,ij,nj->n", e, eff, e)/2+np.einsum("ni,ij,nj->n", midpoint, k, midpoint)/2
    work = -h/2*np.einsum("ni,ni->n", e[:-1]+e[1:], j[:-1])
    gauss = rho-e@m@d
    ampere = (e[1:]-e[:-1])@m+h*j[:-1]-h/6*(magnetic[:-2]+4*magnetic[1:-1]+magnetic[2:])@m2@c
    emass = np.einsum("ni,ij,nj->n", e, m, e)
    bmass = np.einsum("ni,ij,nj->n", magnetic, m2, magnetic)
    lift = np.array([[float(Q(x)) for x in row] for row in parent["cochain_maps"]["P1"]])
    radial_defect = a[:,:12]-a[:,12:]@lift[:12,:].T
    return {"scope": "uninstrumented numeric continuation from the first two decoded temporal-gauge slices",
            "slabs": slabs, "h": h, "A": a.tolist(), "rho_load": rho.tolist(), "J_load": j.tolist(),
            "energy": energy.tolist(), "source_work": work.tolist(),
            "electric_mass_norm_squared": emass.tolist(), "magnetic_mass_norm_squared": bmass.tolist(),
            "gauss_max_abs": float(np.max(np.abs(gauss))), "ampere_max_abs": float(np.max(np.abs(ampere))),
            "faraday_max_abs": float(np.max(np.abs(magnetic[1:]-magnetic[:-1]+h*e@c.T))),
            "work_identity_max_abs": float(np.max(np.abs(energy[1:]-energy[:-1]-work))),
            "source_free_energy_drift": float(np.max(np.abs(energy[1:]-energy[1]))),
            "radial_potential_max_abs": float(np.max(np.abs(a[:,:12]))),
            "static_extension_defect_max_abs": float(np.max(np.abs(radial_defect))),
            "source_free_electric_bound_margin": float(np.min(4*energy[1]-emass[1:])),
            "source_free_magnetic_bound_margin": float(np.min(16*energy[1]-bmass[1:]))}


def scalar_controls():
    result = []
    for stiffness in (1, 48, 64):
        h, lam = Q(1, 2), Q(stiffness)
        f, b = 1+h*h*lam/6, 2-2*h*h*lam/3
        a = [Q(0), Q(1)]
        for _ in range(15):
            a.append((b*a[-1]-f*a[-2])/f)
        result.append({"lambda": str(lam), "h": str(h), "h2lambda": str(h*h*lam),
                       "characteristic_half_trace": str(b/(2*f)), "A": [str(v) for v in a]})
    require(result[1]["A"] == [str(Q((-1)**(n+1)*n)) for n in range(17)], "boundary generalized mode")
    return {"scope": "scalar recurrence controls; lambda=48,64 are not claimed mesh eigenvalues",
            "histories": result}


def build():
    parent, d, c, m, m2, k, replay = load_parent()
    stability = local_stability_certificate(parent)
    executions = [instrument(parent, d, c, m, m2, k, gauge) for gauge in (False, True)]
    tail = continuation(executions[0], parent, d, c, m, m2, k)
    for ex in executions:
        require(ex["metrics"]["gauss_max_abs"] < 1e-10 and ex["metrics"]["ampere_max_abs"] < 1e-10, "stationary two-slab history")
    for key in ("E", "B", "action", "gauss", "ampere_hat"):
        require(np.allclose(executions[0]["metrics"][key], executions[1]["metrics"][key], atol=1e-10, rtol=1e-10), "numerical gauge covariance")
    require(tail["gauss_max_abs"] < 1e-9 and tail["ampere_max_abs"] < 1e-9 and tail["work_identity_max_abs"] < 1e-9,
            "continuation equations")
    missing = [p for p in PIN_PATHS if not (ROOT/p).is_file()]
    require(not missing, "missing implementation pins: "+", ".join(missing))
    packet = {"schema": SCHEMA,
        "scope": "NUMERIC_STATIONARY_FINITE_VOLUME_MAXWELL__EXACT_CLASSICAL_READBACK__EXACT_LOCAL_STABILITY_BOUND",
        "pins": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PIN_PATHS},
        "parent_replay": replay,
        "numeric_policy": {"evolution": "float64 solves and mass matrices; no interval or exact-trajectory certificate",
                           "registers": "exact rational encodings of finite numeric results; exact pair-average restoration",
                           "stability": "exact Q(sqrt(5)) LDL certificate, independent of float eigensolvers",
                           "atol": 1e-9, "rtol": 1e-9},
        "dynamics": {"h": "1/2", "F": "M+h^2*K/6", "K": "C^T*M2*C",
                     "recurrence": "F*A[n+1]=(2M-2h^2*K/3)*A[n]-F*A[n-1]+h^2*J[n-1]-h*M*D*(phi[n]-phi[n-1])",
                     "projection": "E0_new=E0_old+D*(D^T*M*D+ones(13,13)/13)^(-1)*(rho0-D^T*M*E0_old)",
                     "energy": "H=1/2*E^T*(M-h^2*K/12)*E+1/2*midpoint^T*K*midpoint, E=-(A[n+1]-A[n])/h-D*phi[n]",
                     "work": "H[n+1]-H[n]=-h/2*(E[n]+E[n+1])^T*J[n]",
                     "source": "prescribed inherited covectors; current is h*J[n] delta at t[n+1]; conventional physical charge has opposite load sign"},
        "stability_certificate": stability, "executions": executions, "continuation": tail,
        "scalar_controls": scalar_controls(),
        "nonclaims": ["this projected and evolved history replaces, rather than repairs or relabels, the historical serial/cone history",
                      "finite-element stationarity for all 42 edge and 13 scalar coordinates, not arbitrary smooth-test continuum stationarity",
                      "three slices are instrumented; the 64-slab continuation has no serial execution trace",
                      "sources are prescribed; no coupled charged-path or matter stationarity",
                      "supplied geometry, constitutive coefficients, clock step, projection and writable classical coordinate ports",
                      "no source selection, spatial refinement limit, laboratory measurement, SI calibration or new empirical prediction",
                      "field stability does not bound gauge potentials; gradient histories may shear"]}
    # Round-trip ensures that the returned API and written receipt use the same
    # JSON types. Finite-value validation is enforced by allow_nan=False.
    return json.loads(canonical(packet))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    packet = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical(packet))
    print(json.dumps({"receipt": str(args.output), "events_per_history": len(packet["executions"][0]["events"]),
                      "action": packet["executions"][0]["metrics"]["action"],
                      "gauss_max": packet["continuation"]["gauss_max_abs"],
                      "ampere_max": packet["continuation"]["ampere_max_abs"]}, sort_keys=True))
