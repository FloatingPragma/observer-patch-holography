"""Independent SymPy matrix and versioned-register replay; no producer imports.

Rebuild the incidence matrices from oriented triangles, differentiate the
finite action, enumerate every admissible two-step path, and enforce the
instrument program's actual read sets. Parent lists never create reads.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "runtime/serial_maxwell_readout_receipt.json"
PATHS = {
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "Lean/Screen/NeutralPairJointStationaryWitness.lean",
    "Lean/Screen/SerialMaxwellReadout.lean",
    "code/electromagnetism/local_face_maxwell_action.py",
    "code/electromagnetism/serial_maxwell_readout.py",
    "code/electromagnetism/verify_serial_maxwell_readout.py",
    "code/electromagnetism/test_serial_maxwell_readout.py",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(path=OUTPUT):
    def pairs(items):
        result = {}
        for k, v in items:
            require(k not in result, "duplicate key")
            result[k] = v
        return result
    def constant(_):
        raise ValueError("nonfinite JSON")
    return json.loads(path.read_text(), object_pairs_hook=pairs, parse_constant=constant)


def carrier():
    seams = (ROOT / "Lean/Screen/SeamCurrentCarrierQuotient.lean").read_text()
    def row(name):
        segment = seams.split("def " + name, 1)[1].split("![", 1)[1].split("]", 1)[0]
        result = [int(t) for t in re.findall(r"\d+", segment)]
        require(len(result) == 30, "seam count")
        return result
    edges = list(zip(row("seamLeft"), row("seamRight")))
    raw = (ROOT / "Lean/ObserverPatchHolography/CoreAxioms.lean").read_text()
    raw = raw.split("def orientedFaces", 1)[1].split("def faceEdges", 1)[0]
    triangles = [tuple(int(n) for n in t) for t in re.findall(r"\((\d+),\s*(\d+),\s*(\d+)\)", raw)]
    require(len(triangles) == 20 and len(set(edges)) == 30, "carrier census")
    d = sp.zeros(30, 12)
    for i, (u, v) in enumerate(edges):
        d[i, u], d[i, v] = -1, 1
    c = sp.zeros(20, 30)
    for f, tri in enumerate(triangles):
        for u, v in zip(tri, tri[1:]+tri[:1]):
            if (u, v) in edges:
                c[f, edges.index((u, v))] += 1
            else:
                c[f, edges.index((v, u))] -= 1
    require(c*d == sp.zeros(20, 12), "boundary of boundary")
    return edges, d, c


def reference():
    raw = (ROOT / "Lean/Screen/NeutralPairJointStationaryWitness.lean").read_text()
    denominator = raw.split("def potentialDenominator : ℤ :=", 1)[1].split()[0]
    require(denominator == "9192", "reference denominator")
    raw = raw.split("def jointPotentialZ", 1)[1].split("/-- Real seam", 1)[0]
    # Scan the inner vectors separately from the producer's regular expression.
    rows = []
    for part in raw.split("![")[2:]:
        values = [sp.Rational(int(t), 9192) for t in re.findall(r"-?\d+", part.split("]", 1)[0])]
        require(len(values) == 30, "reference arity")
        rows.append(sp.Matrix(values))
    require(len(rows) == 3, "reference slices")
    return rows


def fields_and_action(d, c, h, a, phi, rho, current):
    electric = [-(a[n+1]-a[n])/h-d*phi[n] for n in range(2)]
    magnetic = [c*a[n] for n in range(3)]
    action = sum(h*(electric[n].dot(electric[n])/2-magnetic[n+1].dot(magnetic[n+1])/2
                   +current[n].dot(a[n+1])+rho[n].dot(phi[n])) for n in range(2))
    return electric, magnetic, action


def physics(d, c, edges, h, a, phi, rho, current):
    electric, magnetic, action = fields_and_action(d, c, h, a, phi, rho, current)
    for n in range(2):
        require(d.T*electric[n] == rho[n], "Gauss")
        require(magnetic[n+1]-magnetic[n]+h*c*electric[n] == sp.zeros(20, 1), "Faraday")
    require(electric[1]-electric[0]-h*(c.T*magnetic[1]-current[0]) == sp.zeros(30, 1), "Ampere")
    require(rho[1]-rho[0]+h*d.T*current[0] == sp.zeros(12, 1), "continuity")
    # Independently differentiate the action in every free field coordinate;
    # do not substitute an asserted Euler--Lagrange residual for this check.
    variables = sp.symbols("v:54")
    av = [a[0], sp.Matrix(variables[:30]), a[2]]
    pv = [sp.Matrix(variables[30:42]), sp.Matrix(variables[42:])]
    varied = fields_and_action(d, c, h, av, pv, rho, current)[2]
    substitution = dict(zip(variables, list(a[1])+list(phi[0])+list(phi[1])))
    require(all(sp.diff(varied, v).subs(substitution) == 0 for v in variables), "field action derivative")
    # All admissible intermediate ports, including the two endpoint choices.
    differences = []
    for q, edge in ((1, 0), (-1, 29)):
        start, end = edges[edge]
        def admissible(u, v):
            return u == v or (u, v) in edges or (v, u) in edges
        def hop(n, u, v):
            coupling = h*phi[n][u]
            if (u, v) in edges:
                coupling -= a[n+1][edges.index((u, v))]
            elif (v, u) in edges:
                coupling += a[n+1][edges.index((v, u))]
            return q*coupling - (4 if u != v else 0)
        base = hop(0, start, end)+hop(1, end, end)
        for middle in range(12):
            if admissible(start, middle) and admissible(middle, end):
                differences.append(hop(0, start, middle)+hop(1, middle, end)-base)
    require(len(differences) == 8 and all(v == 0 for v in differences), "complete path stationarity")
    return electric, magnetic, action


def replay(execution):
    require(set(execution) == {"gauge", "events"} and type(execution["gauge"]) is bool, "execution schema")
    edges, d, c = carrier()
    h = sp.Rational(1, 2)
    a_ref = reference()
    chi = sp.Matrix([sp.Rational((u*u+3*u) % 11-5, 7) if execution["gauge"] else 0 for u in range(12)])
    initial = [a_ref[0], a_ref[1]+d*chi]
    p_ref = [-chi/h, chi/h, sp.zeros(12, 1)]
    state, writer = {}, {}
    index = 0
    events = execution["events"]

    def step(op, args, keys, evaluate):
        nonlocal index
        require(index < len(events), "truncated execution")
        event = events[index]
        require(type(event.get("id")) is int and all(type(a) is int for a in event.get("args", [])), "integer identity")
        require(all(type(p) is int for p in event.get("parents", [])), "integer parent identity")
        require(all(type(r.get("writer")) is int for r in event.get("reads", {}).values()), "integer writer identity")
        reads = {k: {"value": str(state[k]), "writer": writer[k]} for k in keys}
        # Values are fetched from the replay state using the opcode's support,
        # never from the receipt's parent list or its alleged read values.
        out = {k: sp.Rational(v) for k, v in evaluate({k: state[k] for k in keys}).items()}
        expected = {"id": index, "op": op, "args": args, "reads": reads,
                    "parents": sorted({writer[k] for k in keys}),
                    "writes": {k: str(v) for k, v in out.items()}}
        require(event == expected, f"semantic/version replay mismatch at {index} ({op})")
        for k, v in out.items():
            state[k], writer[k] = v, index
        index += 1

    inputs = {"clock": 0, "h": h, "q/0": 1, "q/1": -1, "tau/0": 3, "tau/1": 3}
    for j, e in enumerate((0, 29)):
        for n in range(3):
            inputs[f"port/{j}/{n}"] = edges[e][int(n != 0)]
    step("inputs", [], [], lambda _: inputs)
    sourcekeys = ["h", "q/0", "q/1"]+[f"port/{j}/{n}" for j in range(2) for n in range(3)]

    def source(v):
        out = {}
        for n in range(2):
            charge = sp.zeros(12, 1)
            flow = sp.zeros(30, 1)
            for j in range(2):
                u, w = int(v[f"port/{j}/{n}"]), int(v[f"port/{j}/{n+1}"])
                charge[u] += v[f"q/{j}"]
                if u != w:
                    edge = edges.index((min(u, w), max(u, w)))
                    flow[edge] -= v[f"q/{j}"]*sp.sign(w-u)/v["h"]
            out.update({f"rho/{n}/{u}": charge[u] for u in range(12)})
            out.update({f"J/{n}/{e}": flow[e] for e in range(30)})
        return out
    step("sources", [], sourcekeys, source)
    for n in range(3):
        if n < 2:
            packet = list(p_ref[n])+list(initial[n])
            step("seed", [n], [], lambda _, packet=packet: {f"x/{s}": packet[s] for s in range(42)})
        else:
            keys = [f"d/{t}/{s}" for t in range(2) for s in range(42)]+[f"J/0/{e}" for e in range(30)]+["h"]
            def advance(v):
                aa = [sp.Matrix([v[f"d/{t}/{12+e}"] for e in range(30)]) for t in range(2)]
                pp = [sp.Matrix([v[f"d/{t}/{u}"] for u in range(12)]) for t in range(2)]
                ee = -(aa[1]-aa[0])/v["h"]-d*pp[0]
                # Advance E by Ampere, then recover A from the electric
                # definition: distinct from the producer's second difference.
                ee += v["h"]*(c.T*c*aa[1]-sp.Matrix([v[f"J/0/{e}"] for e in range(30)]))
                result = aa[1]-v["h"]*(ee+d*pp[1])
                return {**{f"x/{u}": 0 for u in range(12)}, **{f"x/{12+e}": result[e] for e in range(30)}}
            step("advance", [], keys, advance)
        before = [state[f"x/{s}"] for s in range(42)]
        for u in range(12):
            step("baseline", [n, u], [f"x/{u}"], lambda v, u=u: {f"b/{n}/{u}": v[f"x/{u}"]})
        for e, pair in enumerate(edges):
            for side, u in enumerate(pair):
                x, y = f"x/{u}", f"x/{12+e}"
                base, record = f"b/{n}/{u}", f"r/{n}/{e}/{side}"
                step("probe", [n, e, side], [x, y], lambda v: {k: (v[x]+v[y])/2 for k in (x, y)})
                step("response", [n, e, side], [x], lambda v: {record: v[x]})
                step("feedback", [n, e, side], [base, record, "clock"], lambda v: {
                    x: v[base], y: 2*v[record]-v[base], "clock": v["clock"]+1})
                require([state[f"x/{s}"] for s in range(42)] == before, "same-state restoration")
        keys = [f"b/{n}/{u}" for u in range(12)]+[f"r/{n}/{e}/{s}" for e in range(30) for s in range(2)]
        def decode(v):
            out = {f"d/{n}/{u}": v[f"b/{n}/{u}"] for u in range(12)}
            for e, (u, w) in enumerate(edges):
                # Invert the 2x2 baseline/response observation system.
                solution = sp.Matrix([[1, 0], [sp.Rational(1, 2), sp.Rational(1, 2)]]).inv()*sp.Matrix(
                    [v[f"b/{n}/{u}"], v[f"r/{n}/{e}/0"]])
                require((v[f"b/{n}/{w}"]+solution[1])/2 == v[f"r/{n}/{e}/1"], "second endpoint check")
                out[f"d/{n}/{12+e}"] = solution[1]
            return out
        step("decode", [n], keys, decode)
        require([state[f"d/{n}/{s}"] for s in range(42)] == before, "42-slot decoder")

    publickeys = [f"d/{n}/{s}" for n in range(3) for s in range(12, 42)]
    publickeys += [f"d/{n}/{u}" for n in range(2) for u in range(12)]
    publickeys += [f"rho/{n}/{u}" for n in range(2) for u in range(12)]
    publickeys += [f"J/{n}/{e}" for n in range(2) for e in range(30)]+["h", "clock"]
    publickeys += [f"tau/{j}" for j in range(2)]+[f"port/{j}/{n}" for j in range(2) for n in range(3)]
    aa = [sp.Matrix([state[f"d/{n}/{12+e}"] for e in range(30)]) for n in range(3)]
    pp = [sp.Matrix([state[f"d/{n}/{u}"] for u in range(12)]) for n in range(2)]
    rho = [sp.Matrix([state[f"rho/{n}/{u}"] for u in range(12)]) for n in range(2)]
    current = [sp.Matrix([state[f"J/{n}/{e}"] for e in range(30)]) for n in range(2)]
    require(aa == [a_ref[0], a_ref[1]+d*chi, a_ref[2]], "reference solution agreement")
    electric, magnetic, action = physics(d, c, edges, h, aa, pp, rho, current)

    def public(v):
        # Recompute from the declared public operation's actual fetched values.
        a = [sp.Matrix([v[f"d/{n}/{12+e}"] for e in range(30)]) for n in range(3)]
        p = [sp.Matrix([v[f"d/{n}/{u}"] for u in range(12)]) for n in range(2)]
        r = [sp.Matrix([v[f"rho/{n}/{u}"] for u in range(12)]) for n in range(2)]
        j = [sp.Matrix([v[f"J/{n}/{e}"] for e in range(30)]) for n in range(2)]
        ee, bb, value = fields_and_action(d, c, v["h"], a, p, r, j)
        clock_value = sum(2*v[f"tau/{j}"]**2 - sum(
            4 for n in range(2) if v[f"port/{j}/{n}"] != v[f"port/{j}/{n+1}"]) for j in range(2))
        return {**{f"E/{n}/{e}": ee[n][e] for n in range(2) for e in range(30)},
                **{f"B/{n+1}/{f}": bb[n+1][f] for n in range(2) for f in range(20)},
                **{f"public_rho/{n}/{u}": r[n][u] for n in range(2) for u in range(12)},
                **{f"public_J/{n}/{e}": j[n][e] for n in range(2) for e in range(30)},
                "field_source_action": value, "declared_lorentz_clock_action": clock_value,
                "coupled_action": value+clock_value, "completed_probe_cycles": v["clock"]}
    step("public", [], publickeys, public)
    require(index == len(events) == 585, "complete event census")
    require(state["clock"] == 180, "completed feedback cycles")
    return {"electric": electric, "magnetic": magnetic, "action": action,
            "events": index, "cycles": 180, "path_variations": 8,
            "field_variations": 54, "raw": aa+pp}


def verify(receipt):
    require(set(receipt) == {"schema", "scope", "assumptions", "pins", "executions"}, "root schema")
    require(receipt["schema"] == "oph.serial_maxwell_readout.v1", "schema")
    require(receipt["scope"] == "EXACT_CLASSICAL_SOFTWARE_INSTRUMENT__PHYSICAL_ATTACHMENT_OPEN", "scope")
    require(receipt["assumptions"] == ["writable classical scalar ports and exact retained records",
        "declared potential typing, carrier incidence, action and h=1/2",
        "declared two initial slices and neutral charged port paths",
        "probe-cycle counter is not spacetime or SI time"], "assumptions")
    require(set(receipt["pins"]) == PATHS, "custody census")
    for path, digest in receipt["pins"].items():
        require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest, "source digest")
    require(len(receipt["executions"]) == 2 and [r["gauge"] for r in receipt["executions"]] == [False, True], "gauge control census")
    results = [replay(execution) for execution in receipt["executions"]]
    require(results[0]["raw"] != results[1]["raw"], "nontrivial gauge control")
    for field in ("electric", "magnetic", "action"):
        require(results[0][field] == results[1][field], "public gauge invariance")
    return {k: v for k, v in results[0].items() if k not in ("electric", "magnetic", "raw")}


if __name__ == "__main__":
    print(verify(load()))
