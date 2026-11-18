"""Independent finite count-clock checks; the continuum proof is analytic.

The clock producer is never imported. The parent causet verifier independently
reconstructs the selected source graph, record execution and interval counts.
Rational fourth-power comparisons check every exported clock bound directly.
"""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import types

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "source_count_clock_receipt.json"
PIN_PATHS = (
    "code/causal_refinement/source_count_clock.py",
    "code/causal_refinement/verify_source_count_clock.py",
    "code/causal_refinement/test_source_count_clock.py",
    "Lean/Time/SourceCountClock.lean",
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact(actual, expected):
    require(type(actual) is type(expected), "wrong type")
    if type(expected) is dict:
        require(actual.keys() == expected.keys(), "wrong fields")
        for key in expected:
            exact(actual[key], expected[key])
    elif type(expected) is list:
        require(len(actual) == len(expected), "wrong list length")
        for a, b in zip(actual, expected):
            exact(a, b)
    else:
        require(actual == expected, "wrong value")


def pairs(items):
    result = {}
    for key, value in items:
        require(key not in result, "duplicate JSON key")
        result[key] = value
    return result


def forbidden(token):
    raise ValueError("noninteger JSON token")


def load(path=OUTPUT):
    data = Path(path).read_bytes()
    require(len(data) <= 30_000, "receipt size")
    return json.loads(data.decode("utf-8"), object_pairs_hook=pairs,
                      parse_float=forbidden, parse_constant=forbidden)


def fresh(path):
    module = types.ModuleType("independent_clock_parent")
    module.__file__ = str(path)
    exec(compile(path.read_bytes(), str(path), "exec"), module.__dict__)
    return module


def fraction(value):
    require(type(value) is str and len(value) < 1000, "rational string required")
    result = Fraction(value)
    require(str(result) == value, "noncanonical rational")
    return result


def bound_pair(pair):
    require(type(pair) is list and len(pair) == 2, "root bound pair")
    lo, hi = map(fraction, pair)
    require(0 <= lo <= hi, "positive ordered root bounds")
    for endpoint in (lo, hi):
        require((endpoint*(1 << 80)).denominator == 1, "dyadic grid")
    return lo, hi


def check_floor(endpoint, value):
    require(endpoint**4 <= value < (endpoint+Fraction(1, 1 << 80))**4,
            "incorrect fourth-root floor")


def check_ceil(endpoint, value):
    require(value <= endpoint**4 and
            (endpoint == 0 or (endpoint-Fraction(1, 1 << 80))**4 < value),
            "incorrect fourth-root ceiling")


def verify(receipt):
    require(type(receipt) is dict, "receipt object")
    expected_keys = {"schema", "source_pins", "source_receipt_sha256", "source_q",
                     "authenticated_events", "source_trace_sha256", "interval_counts",
                     "continuum_diamond_inside_source_window", "dyadic_bits",
                     "conditional_error_controls", "ratio_to_first_interval_bounds", "scope"}
    require(receipt.keys() == expected_keys, "receipt fields")
    exact(receipt["schema"], "source-count-clock-v1")
    exact(receipt["source_pins"], {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in PIN_PATHS})
    parent_path = HERE / "source_net_causet_receipt.json"
    exact(receipt["source_receipt_sha256"], sha256(parent_path.read_bytes()).hexdigest())
    parent = fresh(HERE / "verify_source_net_causet.py")
    packet = parent.load()
    exact(packet["source_pins"], {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in parent.PIN_PATHS})
    rows = [r for r in packet["levels"] if type(r.get("q")) is int and r["q"] == 5]
    require(len(rows) == 1, "selected source level")
    row = rows[0]
    parent.verify_level(row)
    counts = [r["inclusive_event_count"] for r in row["center_intervals"]]
    exact(receipt["source_q"], 5)
    exact(receipt["authenticated_events"], row["forward_execution"]["event_count"])
    exact(receipt["source_trace_sha256"], row["forward_execution"]["audit_trace_sha256"])
    exact(receipt["interval_counts"], counts)
    exact(receipt["continuum_diamond_inside_source_window"],
          [r["continuum_diamond_inside_cube"] for r in row["center_intervals"]])
    exact(receipt["dyadic_bits"], 80)
    root_rows = receipt["ratio_to_first_interval_bounds"]
    require(type(root_rows) is list and len(root_rows) == len(counts), "clock rows")
    for pair, count in zip(root_rows, counts):
        lo, hi = bound_pair(pair)
        value = Fraction(count, counts[0])
        check_floor(lo, value)
        check_ceil(hi, value)
    controls = receipt["conditional_error_controls"]
    require(type(controls) is list and len(controls) == 3, "conditional controls")
    for control, (n,m,en,em) in zip(controls, [(16,81,1,2),(1,16,2,1),(16,81,0,0)]):
        require(type(control) is dict and control.keys() ==
                {"counts", "count_error_bounds", "proper_time_ratio_bounds"}, "control fields")
        exact(control["counts"], [n, m])
        exact(control["count_error_bounds"], [str(en), str(em)])
        lo, hi = bound_pair(control["proper_time_ratio_bounds"])
        check_floor(lo, Fraction(max(0,n-en), m+em))
        check_ceil(hi, Fraction(n+en, m-em))
    exact(receipt["scope"], {
        "timestamp_density_or_spatial_coordinate_used_by_clock_readout": False,
        "selected_source_level_freshly_replayed": True,
        "source_commitment_checked": True,
        "external_signature_verified": False,
        "error_control_values_are_supplied_arithmetic_tests": True,
        "proper_time_ratio_limit_is_analytic": True,
        "finite_decoder_control_only": True,
        "finite_clock_accuracy_certified": False,
        "native_read_law_selected": False,
        "physical_clock_identified": False,
    })
    return {"schema": "source-count-clock-v1", "source_q": 5,
            "authenticated_events": row["forward_execution"]["event_count"],
            "interval_counts": counts, "exact_clock_bounds": len(counts),
            "conditional_error_controls": len(controls),
            "proper_time_ratio_limit_is_analytic": True,
            "finite_clock_accuracy_certified": False,
            "physical_clock_identified": False}


if __name__ == "__main__":
    print(json.dumps(verify(load()), sort_keys=True))
