"""Complete stationary-reference reproduction of the known M8 pilot cases."""
import json
import math
from fractions import Fraction as F
from pathlib import Path
from native import CALLS, EDGES, dot, logprob, mv, source
from confirm import angle, contrast_ops, encode_ops, hop, preparation, xy


def integer_ops(ops):
    return [([list(map(int, row)) for row in a], d) for a, d in ops]


def probabilities(left, right, a, b):
    common = math.lcm(*(x.denominator for x in (*a, *b)))
    aa, bb = [int(x * common) for x in a], [int(x * common) for x in b]
    norm = dot(aa, aa) + dot(bb, bb)
    out = []
    for (g, d), (h, dd) in zip(left, right):
        av, bv = mv(g, aa), mv(h, bb)
        out.append(F(dot(av, av) + dot(bv, bv), d * norm))
    assert sum(out) == 1
    return out


base = Path(__file__).resolve().parent
prior = json.loads((base / "CONFIRM_RECEIPT.json").read_text())
pilot = json.loads((base / "PILOT_RECEIPT.json").read_text())["cases"]
cases, catalogs = [], {}
for law in ("native", "swap"):
    systems = [source(i, 8, law, True) for i in range(24)]
    plane = contrast_ops(EDGES, law)
    catalog = {"target": [encode_ops(ops) for p, ops in systems], "plane_filter": encode_ops(plane),
               "invalid_chart": "orthogonal absorbing projector"}
    assert catalog["target"] == prior["instrument_catalog"][f"8:{law}"]["target"]
    assert catalog["plane_filter"] == prior["instrument_catalog"][f"8:{law}"]["plane_filter"]
    catalogs[f"8:{law}"] = catalog
    integer = [integer_ops(ops) for p, ops in systems]
    plane_integer = integer_ops(plane)
    for turns in (1, 2):
        reach, prep = preparation(law)
        rejected, records = 1 - reach, []
        a = b = (F(1), F(-1), F(0))
        rot, ref = 16, 16
        angles = []
        for step in range(1, 24 * turns + 1):
            success, decoded, tape = hop(rot, (16 + step) % 24, law)
            timeout = reach * (1 - success)
            rejected += timeout
            delta = angle(rot, decoded, 8)
            angles.append(delta)
            records.append({"kind": "marker", "arm": "rotating", "step": step, "reach": str(reach),
                            "timeout": str(timeout), "angle_degrees": delta,
                            "decoded_coordinate": xy(decoded, 8), **tape})
            reach *= success
            rot = decoded
            for event in range(4):
                ps = probabilities(integer[rot], integer[ref], a, b)
                failure = reach * (1 - ps[0])
                records.append({"kind": "target", "step": step, "event": event,
                                "rotating_slot": rot, "reference_slot": ref, "reach": str(reach),
                                "first_failure": str(failure), "flag_probabilities": list(map(str, ps)),
                                "audit_input_vectors": [[str(x) for x in v] for v in (a, b)]})
                rejected += failure
                reach *= ps[0]
                a, b = mv(systems[rot][0], a), mv(systems[ref][0], b)
            ps = probabilities(plane_integer, plane_integer, a, b)
            failure = reach * (1 - ps[1])
            records.append({"kind": "plane_filter", "step": step, "reach": str(reach),
                            "first_failure": str(failure), "flag_probabilities": list(map(str, ps)),
                            "audit_input_vectors": [[str(x) for x in v] for v in (a, b)]})
            rejected += failure
            reach *= ps[1]
        norm, ab = dot(a, a) + dot(b, b), dot(a, b)
        visibility = 2 * ab / norm
        overlap = (1 if ab >= 0 else -1) * math.sqrt(float(ab * ab / (dot(a, a) * dot(b, b))))
        expected = next(c for c in pilot if c["M"] == 8 and c["h"] == 4 and c["turns"] == turns
                        and c["law"] == law and c["direction"] == 1 and c["reference"] == "fixed")
        assert visibility == F(expected["exact_visibility"])
        assert rejected + reach == 1
        expected_cost = F(2)
        for row in records:
            unit = (F(255, 128) if law == "native" else F(1)) if row["kind"] == "marker" else (5 if row["kind"] == "target" else 1)
            expected_cost += F(row["reach"]) * unit
        cases.append({"M": 8, "h": 4, "turns": turns, "direction": 1, "law": law,
                      "reference": "stationary marker and source coupling", "visibility": float(visibility),
                      "exact_visibility": str(visibility), "ray_overlap": overlap, "success": str(reach),
                      "log10_success": logprob(reach), "all_rejections": str(rejected), "preparation": prep,
                      "conditional_path_plus": str((1 + visibility) / 2), "conditional_path_minus": str((1 - visibility) / 2),
                      "frame_turn_degrees": sum(angles), "reference_turn_degrees": 0,
                      "maximum_frame_step_degrees": max(map(abs, angles)),
                      "expected_endpoint_opportunities_per_attempt": str(expected_cost),
                      "log10_expected_endpoint_opportunities_per_herald": logprob(expected_cost / reach),
                      "final_arm_vectors": [[str(x) for x in v] for v in (a, b)], "records": records})
failures = []
for c in cases:
    okay = c["visibility"] > 0
    if c["law"] == "native":
        okay = c["visibility"] < 0 and c["ray_overlap"] < -.99999 if c["turns"] == 1 else c["visibility"] > 0 and c["ray_overlap"] > .99999
    if not okay:
        failures.append({k: c[k] for k in ("law", "turns", "visibility", "ray_overlap")})
print(json.dumps({"candidate": "A-native-subdivided-observer-rotation-v3", "cases": cases,
                  "instrument_catalog": catalogs, "gate_failures": failures, "literal_endpoint_calls": CALLS,
                  "kind": "full-outcome stationary-reference reproduction of known pilot values"}, indent=2))
if failures:
    raise SystemExit(1)
