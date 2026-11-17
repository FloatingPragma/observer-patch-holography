"""Complete exact outcome ledger for the actual subdivided observer process."""
import json
from fractions import Fraction as F
from math import atan2, degrees, sqrt
from native import CALLS, EDGES, I, add, dot, endpoint, eraser, fine, logprob, mm, mv, scale, source


def encode_ops(ops):
    return [{"numerator": [[str(x) for x in row] for row in a], "squared_denominator": d}
            for a, d in ops]


def contrast_ops(edges, law):
    target, blank = eraser(fine(edges, law)), eraser(fine(edges, law, True))
    return tuple((add(b, a, sign), 4 * d)
                 for (a, d), (b, db) in zip(target, blank) for sign in (1, -1))


def preparation(law):
    ops, warm = contrast_ops((EDGES[0],), law), fine(EDGES, law)
    input_vector = (F(1), F(0), F(0))
    records, accepted = [], [[F(0)] * 3 for _ in range(3)]
    for w, m in enumerate(warm):
        v = mv(m, input_vector)
        probs = []
        for k, (a, d) in enumerate(ops):
            out = mv(a, v)
            probs.append(dot(out, out) / (6 * d))
            if k == 1:
                for i in range(3):
                    for j in range(3):
                        accepted[i][j] += out[i] * out[j] / (6 * d)
        records.append({"warmup_flag": w, "joint_path_eraser_probabilities": list(map(str, probs))})
    probability = sum(accepted[i][i] for i in range(3))
    f = (F(1), F(-1), F(0))
    assert all(accepted[i][j] == probability * f[i] * f[j] / 2 for i in range(3) for j in range(3))
    assert sum(F(p) for row in records for p in row["joint_path_eraser_probabilities"]) == 1
    return probability, {"records": records, "accepted_probability": str(probability),
                         "accepted_density": [[str(x / probability) for x in row] for row in accepted],
                         "kraus": encode_ops(ops)}


def hop(old, desired, law, cap=8):
    a, b = sorted((old, desired))
    active, histories = [((), (), F(1))], []
    for attempt in range(cap):
        nxt = []
        for coins, currents, p in active:
            for coin in (0, 1):
                before = int(old == a)
                left, right = endpoint(before, 1 - before, coin, law)
                current = before - left
                new_coins, new_currents, q = (*coins, coin), (*currents, current), p / 2
                if current:
                    decoded = b if current == 1 else a
                    assert decoded == desired
                    histories.append({"coins": new_coins, "currents": new_currents,
                                      "probability": str(q), "decoded_position": decoded, "outcome": "hop"})
                elif attempt + 1 == cap:
                    histories.append({"coins": new_coins, "currents": new_currents,
                                      "probability": str(q), "decoded_position": old, "outcome": "timeout"})
                else:
                    nxt.append((new_coins, new_currents, q))
        active = nxt
    success = sum(F(r["probability"]) for r in histories if r["outcome"] == "hop")
    assert sum(F(r["probability"]) for r in histories) == 1
    decoded = next(r["decoded_position"] for r in histories if r["outcome"] == "hop")
    return success, decoded, {"address": [a, b], "histories": histories}


def xy(slot, m):
    vertices = ((1.0, 0.0), (-.5, sqrt(3) / 2), (-.5, -sqrt(3) / 2))
    v, k = divmod(slot % (3 * m), m)
    return tuple((1 - k / m) * vertices[v][i] + k / m * vertices[(v + 1) % 3][i] for i in (0, 1))


def angle(old, new, m):
    a, b = xy(old, m), xy(new, m)
    return degrees(atan2(a[0] * b[1] - a[1] * b[0], dot(a, b)))


def execute(m, turns, direction, law, systems, plane):
    reach, prep = preparation(law)
    rejected = 1 - reach
    a = b = (F(1), F(-1), F(0))
    start = rot = ref = 2 * m
    records, angles, reference_angles = [], [], []
    for step in range(1, 3 * m * turns + 1):
        k = step % m
        requested = ((start + direction * step) % (3 * m),
                     (start + direction * min(k, m - k)) % (3 * m))
        decoded = []
        for arm, old, new in zip(("rotating", "reference"), (rot, ref), requested):
            success, position, record = hop(old, new, law)
            timeout = reach * (1 - success)
            rejected += timeout
            rotation = angle(old, position, m)
            (angles if arm == "rotating" else reference_angles).append(rotation)
            records.append({"kind": "marker", "arm": arm, "step": step,
                            "reach": str(reach), "timeout": str(timeout),
                            "angle_degrees": rotation, "decoded_coordinate": xy(position, m), **record})
            reach *= success
            decoded.append(position)
        rot, ref = decoded
        p, rops = systems[rot]
        q, bops = systems[ref]
        for event in range(4):
            norm = dot(a, a) + dot(b, b)
            probs = [F(dot(mv(g, a), mv(g, a)) + dot(mv(f, b), mv(f, b)), d * norm)
                     for (g, d), (f, dd) in zip(rops, bops)]
            assert sum(probs) == 1
            failure = reach * (1 - probs[0])
            rejected += failure
            records.append({"kind": "target", "step": step, "event": event,
                            "rotating_slot": rot, "reference_slot": ref,
                            "reach": str(reach), "first_failure": str(failure),
                            "flag_probabilities": list(map(str, probs)),
                            "audit_input_vectors": [[str(x) for x in v] for v in (a, b)]})
            reach *= probs[0]
            a, b = mv(p, a), mv(q, b)
        norm = dot(a, a) + dot(b, b)
        probs = [F(dot(mv(g, a), mv(g, a)) + dot(mv(g, b), mv(g, b)), d * norm) for g, d in plane]
        assert sum(probs) == 1
        assert probs[1] == F(1, 48 if law == "native" else 12)
        failure = reach * (1 - probs[1])
        rejected += failure
        records.append({"kind": "plane_filter", "step": step, "reach": str(reach),
                        "first_failure": str(failure), "flag_probabilities": list(map(str, probs)),
                        "audit_input_vectors": [[str(x) for x in v] for v in (a, b)]})
        reach *= probs[1]
    norm, ab = dot(a, a) + dot(b, b), dot(a, b)
    visibility = 2 * ab / norm
    overlap = (1 if ab >= 0 else -1) * sqrt(float(ab * ab / (dot(a, a) * dot(b, b))))
    assert rejected + reach == 1
    assert abs(sum(angles) - direction * 360 * turns) < 1e-10
    assert abs(sum(reference_angles)) < 1e-10
    assert max(map(abs, angles + reference_angles)) < 90
    return {"M": m, "h": 4, "turns": turns, "direction": direction, "law": law,
            "visibility": float(visibility), "exact_visibility": str(visibility), "ray_overlap": overlap,
            "success": str(reach), "log10_success": logprob(reach), "all_rejections": str(rejected),
            "conditional_path_plus": str((1 + visibility) / 2), "conditional_path_minus": str((1 - visibility) / 2),
            "frame_turn_degrees": sum(angles), "reference_turn_degrees": sum(reference_angles),
            "maximum_frame_step_degrees": max(map(abs, angles)), "preparation": prep,
            "final_arm_vectors": [[str(x) for x in v] for v in (a, b)], "records": records}


if __name__ == "__main__":
    cases, catalogs = [], {}
    for m in (4, 8):
        for law in ("native", "swap"):
            systems = [source(i, m, law, True) for i in range(3 * m)]
            plane = contrast_ops(EDGES, law)
            catalogs[f"{m}:{law}"] = {"target": [encode_ops(ops) for p, ops in systems],
                                      "plane_filter": encode_ops(plane),
                                      "invalid_chart": "add |invalid><invalid| as a separate absorbing flag"}
            for turns in (1, 2):
                for direction in (1, -1):
                    cases.append(execute(m, turns, direction, law, systems, plane))
    failures = []
    for c in cases:
        okay = (c["visibility"] < -.99 if c["turns"] == 1 else c["visibility"] > .97)
        okay = okay and abs(c["ray_overlap"]) > .9999
        if c["law"] == "swap":
            okay = F(c["exact_visibility"]) == 1
        if not okay:
            failures.append({k: c[k] for k in ("M", "turns", "direction", "law", "visibility")})
    print(json.dumps({"candidate": "A-native-subdivided-observer-rotation-v3", "cases": cases,
                      "instrument_catalog": catalogs, "gate_failures": failures,
                      "literal_endpoint_calls": CALLS,
                      "compression": "fine Kraus matrices + each input ray and all flag probabilities reconstruct every output"}, indent=2))
    if failures:
        raise SystemExit(1)
