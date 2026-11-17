"""V2: complete native optical-plane filtering before the rotating apparatus."""
import json
from fractions import Fraction as F
from model import (CALLS, EDGES, I, add, apply, complete, density, eraser, fine,
                   logprob, number, scaled, trace)
from process import execute


def filter_ops(law):
    target = eraser(fine(EDGES, law))
    blank = eraser(fine(EDGES, law, vacuum=True))
    ops = tuple((add(b, a, sign), 4 * d)
                for (a, d), (b, db) in zip(target, blank) for sign in (1, -1))
    assert complete(ops)
    return ops


def extend(original):
    law = original["law"]
    ops = filter_ops(law)
    prep = original["preparation"]
    rho = tuple(tuple(F(x) for x in row) for row in prep["accepted_density"])
    reach = F(prep["success_probability"])
    rows, multiplier = [], F(1)
    for index in range(3):
        outputs = [scaled(apply(a, rho), F(1, d)) for a, d in ops]
        probs = [trace(part) for part in outputs]
        assert sum(probs) == 1
        assert probs[1] == F(1, 48 if law == "native" else 12)
        rows.append({"filter_index": index, "reach_probability": str(reach),
                     "accept_flag": 1, "conditional_probabilities": list(map(str, probs)),
                     "unnormalized_output_densities": [[[str(x) for x in row] for row in part]
                                                         for part in outputs],
                     "first_failure_probability": str(reach * (1 - probs[1]))})
        rho = scaled(outputs[1], 1 / probs[1])
        reach *= probs[1]
        multiplier *= probs[1]
    # Exact composition: each accepted filter is scalar on the admitted
    # optical plane, so the subsequent original full instrument is unchanged
    # except for its unconditional reach. Its full failure records remain.
    for record in original["records"]:
        for key in ("reach_probability", "timeout_probability", "failure_probability",
                    "first_failure_probability"):
            if key in record:
                record[key] = str(F(record[key]) * multiplier)
    for key in ("success", "raw_success"):
        original[key] = number(F(original[key]["exact"]) * multiplier)
    success = F(original["success"]["exact"])
    original["all_rejections"] = number(1 - success)
    original["expected_preparations_per_success"] = number(1 / success)
    original["log10_success"] = logprob(success)
    original["optical_filters"] = rows
    original["filter_success_multiplier"] = str(multiplier)
    original["optical_filter_kraus"] = [{"matrix": [[str(x) for x in r] for r in a],
                                          "squared_denominator": d} for a, d in ops]
    return original


if __name__ == "__main__":
    cases = [extend(execute(h=4, turns=t, law=law, direction=d, start=s))
             for s in range(3) for d in (1, -1) for t in (1, 2) for law in ("native", "swap")]
    cases += [extend(execute(h=6, turns=t, direction=d)) for d in (1, -1) for t in (1, 2)]
    print(json.dumps({"candidate": "A-native-incidence-rotation-filtered-v2", "cases": cases,
                      "literal_endpoint_calls": CALLS,
                      "compression": "complete branch instruments composed exactly; accepted plane filters are scalar"}, indent=2))
