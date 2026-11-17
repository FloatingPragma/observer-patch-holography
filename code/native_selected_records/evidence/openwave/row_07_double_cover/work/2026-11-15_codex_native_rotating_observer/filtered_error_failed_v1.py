"""Complete probability audit under correlated preparation/filter phase error."""
import cmath
import json
from model import (CALLS, EDGES, I, add, density, eraser, fine, mm, scaled,
                   star, trace, transpose)


def adj(a):
    return tuple(tuple(complex(x).conjugate() for x in row) for row in transpose(a))


def cp(a, rho):
    return mm(mm(a, rho), adj(a))


def prob(rho):
    x = trace(rho)
    assert abs(complex(x).imag) < 1e-12 * max(abs(x.real), 1e-300)
    return float(x.real)


def phase_ops(edges, delta):
    target, blank = eraser(fine(edges, "native")), eraser(fine(edges, "native", True))
    return tuple(scaled(add(b, a, sign * cmath.exp(1j * delta)), (4 * d) ** -.5)
                 for (a, d), (b, db) in zip(target, blank) for sign in (1, -1))


def diagnostic(delta, turns):
    initial, warm = density((1, 1, 1)), fine(EDGES, "native")
    ops = phase_ops((EDGES[0],), delta)
    prepared, records = scaled(I, 0), []
    rejected = 0.0
    for w, m in enumerate(warm):
        rho = scaled(cp(m, initial), 1 / 6)
        outputs = [cp(a, rho) for a in ops]
        ps = [prob(x) for x in outputs]
        records.append({"kind": "warmup_and_preparation", "warmup_flag": w,
                        "path_eraser_probabilities": ps, "accept_flag": 1})
        prepared = add(prepared, outputs[1])
        rejected += sum(ps) - ps[1]
    rho = prepared
    for k in range(3):
        outputs = [cp(a, rho) for a in phase_ops(EDGES, delta)]
        ps = [prob(x) for x in outputs]
        records.append({"kind": "optical_filter", "index": k, "reach": prob(rho),
                        "path_eraser_probabilities": ps, "accept_flag": 1})
        assert abs(sum(ps) / prob(rho) - 1) < 1e-12
        rejected += sum(ps) - ps[1]
        rho = outputs[1]
    q = density((1, -1, 0))
    prep_fidelity = prob(mm(q, rho)) / prob(rho)
    left = right = cross = scaled(rho, .5)
    vertex = 2
    for hop in range(3 * turns):
        timeout = (prob(left) + prob(right)) / 256
        records.append({"kind": "marker_timeout", "probability": timeout})
        rejected += timeout
        left, right, cross = (scaled(x, 255 / 256) for x in (left, right, cross))
        vertex = (vertex + 1) % 3
        r_ops, b_ops = star(vertex, "native")[1], star(2, "native")[1]
        for event in range(4):
            flags, states = [], []
            for (a, d), (b, db) in zip(r_ops, b_ops):
                ll, rr = scaled(cp(a, left), 1 / d), scaled(cp(b, right), 1 / d)
                cc = scaled(mm(mm(a, cross), adj(b)), 1 / d)
                flags.append(prob(ll) + prob(rr))
                states.append((ll, rr, cc))
            reach = prob(left) + prob(right)
            assert abs(sum(flags) / reach - 1) < 1e-11
            records.append({"kind": "target", "vertex": vertex, "reach": reach,
                            "flag_probabilities": flags})
            rejected += sum(flags[1:])
            left, right, cross = states[0]
    for k in range(6 * turns):
        fail = prob(right) / 2
        records.append({"kind": "reference_coin", "index": k, "failure_probability": fail})
        rejected += fail
        right, cross = scaled(right, .5), scaled(cross, 2 ** -.5)
    success = prob(left) + prob(right)
    visibility = 2 * trace(cross).real / success
    fidelity = prob(mm(q, left)) / prob(left)
    normalization_error = abs(1 - rejected - success)
    assert normalization_error < 1e-12
    return {"delta_radians": delta, "turns": turns, "visibility": visibility,
            "filtered_preparation_contrast_fidelity": prep_fidelity,
            "rotated_contrast_fidelity": fidelity, "success": success,
            "normalization_error": normalization_error, "records": records}


cases = [diagnostic(d, t) for d in (0, 1e-5, -1e-5, 1e-3, -1e-3, 1e-2, -1e-2, 1e-1, -1e-1)
         for t in (1, 2)]
failures = [c for c in cases if abs(c["delta_radians"]) <= .01 and not
            (c["rotated_contrast_fidelity"] > .999 and
             (c["visibility"] < -.99 if c["turns"] == 1 else c["visibility"] > .97))]
print(json.dumps({"candidate": "A-native-incidence-rotation-filtered-v2", "cases": cases,
                  "gate_failures": [{k: c[k] for k in ("delta_radians", "turns", "visibility",
                                                       "rotated_contrast_fidelity")} for c in failures],
                  "literal_endpoint_calls": CALLS}, indent=2))
if failures:
    raise SystemExit(1)
