"""Exact record-assisted Maxwell instrument, with an authenticated execution trace.

Three potential slices occupy the same 42 writable scalar slots in succession.
Each probe restores its two slots using retained records. The third slice is
advanced from decoded records, not loaded from the reference solution table.
This bounded software observer has local state, ports, readback, feedback and
a public evidence bundle; it supplies no laboratory or spacetime calibration.
"""
from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import re

from local_face_maxwell_action import parse_vector, parse_faces, build_incidence

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "runtime/serial_maxwell_readout_receipt.json"
PROVIDERS = [
    "Lean/Screen/SeamCurrentCarrierQuotient.lean",
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "Lean/Screen/NeutralPairJointStationaryWitness.lean",
    "Lean/Screen/SerialMaxwellReadout.lean",
    "code/electromagnetism/local_face_maxwell_action.py",
    "code/electromagnetism/serial_maxwell_readout.py",
    "code/electromagnetism/verify_serial_maxwell_readout.py",
    "code/electromagnetism/test_serial_maxwell_readout.py",
]


def geometry():
    text = (ROOT / PROVIDERS[0]).read_text()
    left, right = [parse_vector(text, name, 30) for name in ("seamLeft", "seamRight")]
    faces = parse_faces((ROOT / PROVIDERS[1]).read_text())
    b, c = build_incidence(left, right, faces)
    return left, right, b, c


def reference():
    text = (ROOT / PROVIDERS[2]).read_text().split("def jointPotentialZ", 1)[1]
    rows = re.findall(r"!\[([^\[\]]+)\]", text.split("/-- Real seam", 1)[0])
    return [[Q(int(v), 9192) for v in row.split(",")] for row in rows]


def mv(matrix, vector):
    return [sum((Q(x) * y for x, y in zip(row, vector, strict=True)), Q(0))
            for row in matrix]


def dot(x, y):
    return sum((a*b for a, b in zip(x, y, strict=True)), Q(0))


class Recorder:
    def __init__(self):
        self.state, self.writer, self.events = {}, {}, []

    def emit(self, op, args, keys, calculate):
        # Read versions and values are captured before calculating the write.
        reads = {key: {"writer": self.writer[key], "value": str(self.state[key])}
                 for key in keys}
        values = {key: self.state[key] for key in keys}
        writes = {key: Q(value) for key, value in calculate(values).items()}
        index = len(self.events)
        self.events.append({"id": index, "op": op, "args": args, "reads": reads,
                            "parents": sorted({v["writer"] for v in reads.values()}),
                            "writes": {key: str(value) for key, value in writes.items()}})
        self.state.update(writes)
        self.writer.update({key: index for key in writes})


def execute(gauge=False):
    left, right, boundary, curl = geometry()
    transpose_c = list(zip(*curl))
    initial = reference()[:2]  # The third reference row is never an execution input.
    chi = [Q((u*u+3*u) % 11 - 5, 7) if gauge else Q(0) for u in range(12)]
    phi = [[-2*v for v in chi], [2*v for v in chi], [Q(0)]*12]
    initial[1] = [v + chi[r]-chi[l] for v, l, r in zip(initial[1], left, right)]
    rec = Recorder()
    inputs = {"clock": Q(0), "h": Q(1, 2), "q/0": Q(1), "q/1": Q(-1),
              "tau/0": Q(3), "tau/1": Q(3)}
    for j, edge in enumerate((0, 29)):
        for n in range(3):
            inputs[f"port/{j}/{n}"] = Q(left[edge] if n == 0 else right[edge])
    rec.emit("inputs", [], [], lambda _: inputs)

    source_keys = ["h", "q/0", "q/1"] + [f"port/{j}/{n}" for j in range(2) for n in range(3)]

    def sources(v):
        out = {}
        for n in range(2):
            for u in range(12):
                out[f"rho/{n}/{u}"] = sum(v[f"q/{j}"] for j in range(2)
                                             if v[f"port/{j}/{n}"] == u)
            for e, (l, r) in enumerate(zip(left, right)):
                flux = Q(0)
                for j in range(2):
                    start, end = v[f"port/{j}/{n}"], v[f"port/{j}/{n+1}"]
                    # Positive charge moving left-to-right has J=-q/h in
                    # the existing boundary/continuity convention.
                    flux += v[f"q/{j}"] / v["h"] * (
                        int(start == r and end == l) - int(start == l and end == r))
                out[f"J/{n}/{e}"] = flux
        return out
    rec.emit("sources", [], source_keys, sources)

    for n in range(3):
        if n < 2:
            rec.emit("seed", [n], [], lambda _, n=n: {
                **{f"x/{u}": phi[n][u] for u in range(12)},
                **{f"x/{12+e}": initial[n][e] for e in range(30)}})
        else:
            keys = [f"d/{t}/{s}" for t in (0, 1) for s in range(42)] + [
                f"J/0/{e}" for e in range(30)] + ["h"]

            def advance(v):
                h = v["h"]
                a0, a1 = [[v[f"d/{t}/{12+e}"] for e in range(30)] for t in (0, 1)]
                dp = [v[f"d/1/{u}"]-v[f"d/0/{u}"] for u in range(12)]
                ka = mv(transpose_c, mv(curl, a1))
                return {**{f"x/{u}": Q(0) for u in range(12)}, **{
                    f"x/{12+e}": 2*a1[e]-a0[e]-h*h*ka[e]+h*h*v[f"J/0/{e}"]
                    - h*(dp[right[e]]-dp[left[e]]) for e in range(30)}}
            rec.emit("advance", [], keys, advance)

        for u in range(12):
            rec.emit("baseline", [n, u], [f"x/{u}"],
                     lambda v, u=u, n=n: {f"b/{n}/{u}": v[f"x/{u}"]})
        for e in range(30):
            for side, parent in enumerate((left[e], right[e])):
                x, y = f"x/{parent}", f"x/{12+e}"
                b, r = f"b/{n}/{parent}", f"r/{n}/{e}/{side}"
                rec.emit("probe", [n, e, side], [x, y],
                         lambda v, x=x, y=y: {x: (v[x]+v[y])/2, y: (v[x]+v[y])/2})
                rec.emit("response", [n, e, side], [x], lambda v, x=x, r=r: {r: v[x]})
                rec.emit("feedback", [n, e, side], [b, r, "clock"],
                         lambda v, b=b, r=r, x=x, y=y: {
                             x: v[b], y: 2*v[r]-v[b], "clock": v["clock"]+1})
        keys = [f"b/{n}/{u}" for u in range(12)] + [f"r/{n}/{e}/{s}" for e in range(30) for s in range(2)]

        def decode(v, n=n):
            out = {f"d/{n}/{u}": v[f"b/{n}/{u}"] for u in range(12)}
            for e in range(30):
                a = 2*v[f"r/{n}/{e}/0"]-v[f"b/{n}/{left[e]}"]
                other = 2*v[f"r/{n}/{e}/1"]-v[f"b/{n}/{right[e]}"]
                if a != other:
                    raise ValueError("inconsistent endpoint readouts")
                out[f"d/{n}/{12+e}"] = a
            return out
        rec.emit("decode", [n], keys, decode)

    keys = [f"d/{n}/{s}" for n in range(3) for s in range(12, 42)]
    keys += [f"d/{n}/{u}" for n in range(2) for u in range(12)]
    keys += [f"rho/{n}/{u}" for n in range(2) for u in range(12)]
    keys += [f"J/{n}/{e}" for n in range(2) for e in range(30)] + ["h", "clock"]
    keys += [f"tau/{j}" for j in range(2)] + [f"port/{j}/{n}" for j in range(2) for n in range(3)]

    def public(v):
        a = [[v[f"d/{n}/{12+e}"] for e in range(30)] for n in range(3)]
        out, action = {}, Q(0)
        for n in range(2):
            efield = [-(a[n+1][e]-a[n][e])/v["h"]
                      -(v[f"d/{n}/{right[e]}"]-v[f"d/{n}/{left[e]}"]) for e in range(30)]
            bfield = mv(curl, a[n+1])
            out.update({f"E/{n}/{e}": x for e, x in enumerate(efield)})
            out.update({f"B/{n+1}/{f}": x for f, x in enumerate(bfield)})
            action += v["h"]*(dot(efield, efield)/2-dot(bfield, bfield)/2
                      +dot([v[f"J/{n}/{e}"] for e in range(30)], a[n+1])
                      +dot([v[f"rho/{n}/{u}"] for u in range(12)],
                           [v[f"d/{n}/{u}"] for u in range(12)]))
        out["field_source_action"] = action
        clock_action = sum(v[f"tau/{j}"]**2 - 4*int(v[f"port/{j}/{n+1}"] != v[f"port/{j}/{n}"])
                           for j in range(2) for n in range(2))
        out["declared_lorentz_clock_action"] = clock_action
        out["coupled_action"] = action + clock_action
        out.update({f"public_rho/{n}/{u}": v[f"rho/{n}/{u}"] for n in range(2) for u in range(12)})
        out.update({f"public_J/{n}/{e}": v[f"J/{n}/{e}"] for n in range(2) for e in range(30)})
        out["completed_probe_cycles"] = v["clock"]
        return out
    rec.emit("public", [], keys, public)
    return {"gauge": gauge, "events": rec.events}


def build():
    return {"schema": "oph.serial_maxwell_readout.v1",
            "scope": "EXACT_CLASSICAL_SOFTWARE_INSTRUMENT__PHYSICAL_ATTACHMENT_OPEN",
            "assumptions": ["writable classical scalar ports and exact retained records",
                            "declared potential typing, carrier incidence, action and h=1/2",
                            "declared two initial slices and neutral charged port paths",
                            "probe-cycle counter is not spacetime or SI time"],
            "pins": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in PROVIDERS},
            "executions": [execute(False), execute(True)]}


if __name__ == "__main__":
    OUTPUT.write_text(json.dumps(build(), sort_keys=True, separators=(",", ":")) + "\n")
    print(OUTPUT)
