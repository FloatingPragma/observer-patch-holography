"""Stable amplitude-ensemble full audit of correlated analyzer phase errors.

No density-matrix cancellation is used for tiny rejected probabilities.
"""
import cmath
import json
from model import CALLS, EDGES, add, eraser, fine, mv, scaled, star


def norm(v):
    return sum(abs(x) ** 2 for x in v)


def total(ensemble):
    return sum(norm(v) for v in ensemble)


def apply(a, ensemble):
    return [mv(a, v) for v in ensemble]


def phase_ops(edges, delta):
    target, blank = eraser(fine(edges, "native")), eraser(fine(edges, "native", True))
    return tuple(scaled(add(b, a, sign * cmath.exp(1j * delta)), (4 * d) ** -.5)
                 for (a, d), (b, db) in zip(target, blank) for sign in (1, -1))


def contrast_fidelity(ensemble):
    return sum(abs(v[0] - v[1]) ** 2 / 2 for v in ensemble) / total(ensemble)


def diagnostic(delta, turns):
    initial = (3 ** -.5,) * 3
    warm = fine(EDGES, "native")
    ops = phase_ops((EDGES[0],), delta)
    prepared, records, rejected = [], [], 0.0
    for w, m in enumerate(warm):
        ensemble = [tuple(x / 6 ** .5 for x in mv(m, initial))]
        outputs = [apply(a, ensemble) for a in ops]
        ps = [total(x) for x in outputs]
        records.append({"kind": "warmup_and_preparation", "warmup_flag": w,
                        "path_eraser_probabilities": ps, "accept_flag": 1})
        prepared.extend(outputs[1])
        rejected += sum(ps) - ps[1]
    ensemble = prepared
    for k in range(3):
        outputs = [apply(a, ensemble) for a in phase_ops(EDGES, delta)]
        ps = [total(x) for x in outputs]
        records.append({"kind": "optical_filter", "index": k, "reach": total(ensemble),
                        "path_eraser_probabilities": ps, "accept_flag": 1})
        assert abs(sum(ps) / total(ensemble) - 1) < 1e-12
        rejected += sum(ps) - ps[1]
        ensemble = outputs[1]
    prep_fidelity = contrast_fidelity(ensemble)
    left = right = [tuple(x / 2 ** .5 for x in v) for v in ensemble]
    vertex = 2
    for hop in range(3 * turns):
        timeout = (total(left) + total(right)) / 256
        records.append({"kind": "marker_timeout", "probability": timeout})
        rejected += timeout
        left, right = [[tuple(x * (255 / 256) ** .5 for x in v) for v in arm]
                       for arm in (left, right)]
        vertex = (vertex + 1) % 3
        r_ops, b_ops = star(vertex, "native")[1], star(2, "native")[1]
        for event in range(4):
            flags, states = [], []
            for (a, d), (b, db) in zip(r_ops, b_ops):
                ll, rr = apply(scaled(a, d ** -.5), left), apply(scaled(b, d ** -.5), right)
                flags.append(total(ll) + total(rr))
                states.append((ll, rr))
            reach = total(left) + total(right)
            assert abs(sum(flags) / reach - 1) < 1e-11
            records.append({"kind": "target", "vertex": vertex, "reach": reach,
                            "flag_probabilities": flags})
            rejected += sum(flags[1:])
            left, right = states[0]
    for k in range(6 * turns):
        fail = total(right) / 2
        records.append({"kind": "reference_coin", "index": k, "failure_probability": fail})
        rejected += fail
        right = [tuple(x / 2 ** .5 for x in v) for v in right]
    success = total(left) + total(right)
    cross = sum(sum(complex(x).conjugate() * y for x, y in zip(a, b))
                for a, b in zip(left, right))
    visibility = 2 * cross.real / success
    fidelity = contrast_fidelity(left)
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
