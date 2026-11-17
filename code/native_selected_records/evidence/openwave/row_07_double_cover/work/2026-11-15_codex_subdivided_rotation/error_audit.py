"""Complete amplitude-ensemble phase-error audit with periodic native filters."""
import cmath
import json
from math import sqrt
from native import CALLS, EDGES, add, dot, eraser, fine, mv, scale, source


def norm(v):
    return sum(abs(x) ** 2 for x in v)


def total(ensemble):
    return sum(norm(v) for v in ensemble)


def act(a, ensemble):
    return [mv(a, v) for v in ensemble]


def phase_ops(edges, delta):
    target, blank = eraser(fine(edges, "native")), eraser(fine(edges, "native", True))
    return [scale(add(b, a, sign * cmath.exp(1j * delta)), (4 * d) ** -.5)
            for (a, d), (b, db) in zip(target, blank) for sign in (1, -1)]


def density(ensemble):
    return [[sum(v[i] * complex(v[j]).conjugate() for v in ensemble) for j in range(3)] for i in range(3)]


def encode(ensemble):
    return [[[complex(x).real, complex(x).imag] for x in v] for v in ensemble]


def audit(m, turns, delta, systems):
    prep_ops = phase_ops((EDGES[0],), delta)
    records, accepted = [], []
    for index, branch in enumerate(fine(EDGES, "native")):
        v = tuple(x / sqrt(18) for x in mv(branch, (1, 1, 1)))
        outputs = [mv(a, v) for a in prep_ops]
        records.append({"kind": "preparation", "warmup_flag": index,
                        "joint_flag_probabilities": [norm(x) for x in outputs]})
        accepted.append(outputs[1])
    # The committed uniform-input warmup is diagonal in f,g,c. This checked
    # three-mode decomposition compresses18 branches without discarding them.
    t = abs(1 - cmath.exp(1j * delta)) ** 2
    weights = [1 / 54, t / 54, 5 * t / 108]
    basis = [(1 / sqrt(2), -1 / sqrt(2), 0),
             (1 / sqrt(6), 1 / sqrt(6), -2 / sqrt(6)), (1 / sqrt(3),) * 3]
    ensemble = [tuple(sqrt(w) * x for x in v) for w, v in zip(weights, basis)]
    old, new = density(accepted), density(ensemble)
    assert max(abs(old[i][j] - new[i][j]) for i in range(3) for j in range(3)) < 1e-15
    rejected = 1 - total(ensemble)
    left = right = [tuple(x / sqrt(2) for x in v) for v in ensemble]
    filters = phase_ops(EDGES, delta)
    for step in range(1, 3 * m * turns + 1):
        k = step % m
        rot = (2 * m + step) % (3 * m)
        ref = (2 * m + min(k, m - k)) % (3 * m)
        for marker in ("rotating", "reference"):
            reach = total(left) + total(right)
            fail = reach / 256
            rejected += fail
            records.append({"kind": "marker", "arm": marker, "step": step, "reach": reach, "timeout": fail})
            left, right = [[tuple(sqrt(255 / 256) * x for x in v) for v in arm] for arm in (left, right)]
        operations = [("target", systems[rot], systems[ref])] * 4 + [("plane_filter", filters, filters)]
        for event, (kind, ops, refs) in enumerate(operations):
            states = [(act(a, left), act(b, right)) for a, b in zip(ops, refs)]
            ps = [total(a) + total(b) for a, b in states]
            reach = total(left) + total(right)
            assert abs(sum(ps) / reach - 1) < 2e-11
            accept = 0 if kind == "target" else 1
            rejected += sum(ps) - ps[accept]
            records.append({"kind": kind, "step": step, "event": event, "reach": reach,
                            "flag_probabilities": ps, "accept_flag": accept,
                            "input_amplitudes": [encode(left), encode(right)]})
            left, right = states[accept]
    success = total(left) + total(right)
    cross = sum(sum(complex(x).conjugate() * y for x, y in zip(a, b)) for a, b in zip(left, right))
    visibility = 2 * cross.real / success
    fidelity = sum(abs(v[0] - v[1]) ** 2 / 2 for v in left) / total(left)
    assert abs(1 - rejected - success) < 1e-12
    return {"M": m, "turns": turns, "delta_radians": delta, "visibility": visibility,
            "rotated_contrast_fidelity": fidelity, "success": success,
            "normalization_error": abs(1 - rejected - success), "records": records}


cases = []
for m in (4, 8):
    systems = [[scale(a, d ** -.5) for a, d in source(slot, m, "native", True)[1]] for slot in range(3 * m)]
    for delta in (0, .01, -.01):
        for turns in (1, 2):
            cases.append(audit(m, turns, delta, systems))
    if m == 8:
        cases += [audit(m, 2, delta, systems) for delta in (.1, -.1)]
failures = [c for c in cases if abs(c["delta_radians"]) <= .01 and not
            (c["rotated_contrast_fidelity"] > .9999 and
             (c["visibility"] < -.99 if c["turns"] == 1 else c["visibility"] > .97))]
print(json.dumps({"candidate": "A-native-subdivided-observer-rotation-v3", "cases": cases,
                  "gate_failures": [{k: c[k] for k in ("M", "turns", "delta_radians", "visibility",
                                                       "rotated_contrast_fidelity")} for c in failures],
                  "literal_endpoint_calls": CALLS,
                  "compression": "uniform-input warmup18 branches checked against exact3-mode mixture; all subsequent outcomes retained"}, indent=2))
if failures:
    raise SystemExit(1)
