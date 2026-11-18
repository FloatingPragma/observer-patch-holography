"""Exact retrospective count-clock readout; no timestamps or density input."""
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import types

HERE = Path(__file__).resolve().parent
RER = HERE.parents[1]
OUTPUT = HERE / "source_count_clock_receipt.json"
PIN_PATHS = (
    "code/causal_refinement/source_count_clock.py",
    "code/causal_refinement/verify_source_count_clock.py",
    "code/causal_refinement/test_source_count_clock.py",
    "Lean/Time/SourceCountClock.lean",
)


def positive_int(value):
    if type(value) is not int or value <= 0:
        raise ValueError("positive integer required")
    return value


def rational(value):
    if type(value) not in (int, str, Fraction):
        raise ValueError("exact rational required")
    return Fraction(value)


def fourth_root_bounds(value, bits=80):
    x = rational(value)
    if x < 0 or type(bits) is not int or not 1 <= bits <= 256:
        raise ValueError("nonnegative rational and bounded integer precision required")
    grid = 1 << bits
    k = isqrt(isqrt((x.numerator * grid**4) // x.denominator))
    lo = Fraction(k, grid)
    hi = lo if lo**4 == x else Fraction(k+1, grid)
    assert lo**4 <= x <= hi**4
    return lo, hi


def count_clock(numerator, denominator, bits=80):
    return fourth_root_bounds(Fraction(positive_int(numerator), positive_int(denominator)), bits)


def certified_ratio_bounds(numerator, denominator, numerator_error, denominator_error, bits=80):
    """Conditional continuum ratio bound: errors are independently supplied count units."""
    n, d = positive_int(numerator), positive_int(denominator)
    en, ed = rational(numerator_error), rational(denominator_error)
    if en < 0 or ed < 0 or d <= ed:
        raise ValueError("nonnegative errors and a strictly positive denominator bound required")
    lo = fourth_root_bounds(max(Fraction(0), n-en)/(d+ed), bits)[0]
    hi = fourth_root_bounds((n+en)/(d-ed), bits)[1]
    return lo, hi


def canonical(material):
    return (json.dumps(material, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False)+"\n").encode("ascii")


def authenticated_order(trace, expected_digest):
    """Check records against a source-pinned digest, not an external signature.

    Reads name immutable versions, so a valid old snapshot remains readable
    after a later write. Numerical event labels have no temporal meaning here.
    """
    if (type(expected_digest) is not str or len(expected_digest) != 64 or
            any(c not in "0123456789abcdef" for c in expected_digest)):
        raise ValueError("authenticated trace digest required")
    state, latest, ancestors = {}, {}, {}
    chain = "0"*64
    for material in trace:
        if type(material) is not list or len(material) != 3:
            raise ValueError("event material shape")
        label, reads, write = material
        if (type(label) is not list or len(label) != 2 or
                any(type(x) is not int or x < 0 for x in label)):
            raise ValueError("event identity")
        event = tuple(label)
        if event in ancestors or type(reads) is not list:
            raise ValueError("duplicate event or malformed reads")
        past, resources = set(), set()
        for read in reads:
            if (type(read) is not list or len(read) != 4 or
                    type(read[0]) is not int or type(read[1]) is not int or
                    type(read[3]) is not int):
                raise ValueError("read shape")
            resource, version, writer, value = read
            if resource in resources or type(writer) is not list or len(writer) != 2:
                raise ValueError("duplicate read or writer shape")
            if any(type(x) is not int for x in writer):
                raise ValueError("writer identity")
            parent = tuple(writer)
            if state.get((resource, version)) != (parent, value) or parent not in ancestors:
                raise ValueError("read is not the actual register writer/version/value")
            resources.add(resource)
            past.add(parent)
            past.update(ancestors[parent])
        if (type(write) is not list or len(write) != 4 or
                any(type(write[i]) is not int for i in (0, 1, 3)) or
                type(write[2]) is not list or len(write[2]) != 2 or
                any(type(x) is not int for x in write[2]) or write[2] != label):
            raise ValueError("write shape or writer identity")
        resource, version, _, value = write
        if resource < 0 or version <= 0:
            raise ValueError("write resource/version")
        if resource in latest and version != latest[resource]+1:
            raise ValueError("write did not advance the actual register version")
        if resource not in latest and version != 1:
            raise ValueError("first register version")
        state[(resource, version)] = (event, value)
        latest[resource] = version
        ancestors[event] = frozenset(past)
        chain = sha256(bytes.fromhex(chain)+canonical(material)).hexdigest()
    if chain != expected_digest:
        raise ValueError("trace differs from its authenticated source commitment")
    return ancestors


def interval_count(ancestors, lower, upper):
    anchors = []
    for value in (lower, upper):
        if (type(value) not in (tuple, list) or len(value) != 2 or
                any(type(x) is not int or x < 0 for x in value)):
            raise ValueError("exact nonnegative integer anchor pair required")
        anchors.append(tuple(value))
    lower, upper = anchors
    if lower not in ancestors or upper not in ancestors:
        raise ValueError("missing interval anchor")
    if lower != upper and lower not in ancestors[upper]:
        raise ValueError("interval anchors not ordered")
    return sum(e == lower or lower in ancestors[e]
               for e in ancestors[upper] | {upper})


def fresh_module(path):
    """Compile authenticated current source bytes, bypassing stale bytecode caches."""
    module = types.ModuleType("count_clock_source_"+path.stem)
    module.__file__ = str(path)
    exec(compile(path.read_bytes(), str(path), "exec"), module.__dict__)
    return module


def source_fixture():
    """Replay the retained q=5 finite law; this is a decoder control, not clock convergence."""
    directory = RER / "code/causal_refinement"
    verifier = fresh_module(directory / "verify_source_net_causet.py")
    packet = verifier.load()
    for name, digest in packet["source_pins"].items():
        if sha256((RER/name).read_bytes()).hexdigest() != digest:
            raise ValueError("stale parent source pin: "+name)
    row = next(r for r in packet["levels"] if r["q"] == 5)
    verifier.verify_level(row)
    producer = fresh_module(directory / "source_net_causet.py")
    _, _, _, _, _, neighbors = producer.geometry(row["fibonacci_index"])
    trace, previous = [], []
    for layer in range(row["layer_steps"]+1):
        current = []
        for site in range(len(neighbors)):
            reads = [] if layer == 0 else [[r, layer, [layer-1, r], previous[r]] for r in neighbors[site]]
            value = 1+site if layer == 0 else 1+sum(r[3] for r in reads)
            current.append(value)
            trace.append([[layer, site], reads, [site, layer+1, [layer, site], value]])
        previous = current
    order = authenticated_order(trace, row["forward_execution"]["audit_trace_sha256"])
    center = row["intervention_source_id"]
    counts = [interval_count(order, (0, center), (k, center)) for k in range(1, row["layer_steps"]+1)]
    expected = [x["inclusive_event_count"] for x in row["center_intervals"]]
    if counts != expected:
        raise ValueError("independent ancestry count disagrees with parent distance-based count")
    return trace, order, row, counts


def report():
    _, order, row, counts = source_fixture()
    return {
        "source_receipt_sha256": sha256((RER/"code/causal_refinement/source_net_causet_receipt.json").read_bytes()).hexdigest(),
        "schema": "source-count-clock-v1",
        "source_pins": {name: sha256((RER/name).read_bytes()).hexdigest() for name in PIN_PATHS},
        "source_q": 5,
        "authenticated_events": len(order),
        "source_trace_sha256": row["forward_execution"]["audit_trace_sha256"],
        "interval_counts": counts,
        "continuum_diamond_inside_source_window": [r["continuum_diamond_inside_cube"] for r in row["center_intervals"]],
        "dyadic_bits": 80,
        "conditional_error_controls": [
            {"counts": [n, m], "count_error_bounds": [str(en), str(em)],
             "proper_time_ratio_bounds": [str(v) for v in certified_ratio_bounds(n,m,en,em)]}
            for n,m,en,em in [(16,81,1,2),(1,16,2,1),(16,81,0,0)]
        ],
        "ratio_to_first_interval_bounds": [[str(v) for v in count_clock(n, counts[0])] for n in counts],
        "scope": {"timestamp_density_or_spatial_coordinate_used_by_clock_readout": False,
                  "selected_source_level_freshly_replayed": True,
                  "source_commitment_checked": True,
                  "external_signature_verified": False,
                  "error_control_values_are_supplied_arithmetic_tests": True,
                  "proper_time_ratio_limit_is_analytic": True,
                  "finite_decoder_control_only": True,
                  "finite_clock_accuracy_certified": False,
                  "native_read_law_selected": False, "physical_clock_identified": False},
    }


if __name__ == "__main__":
    OUTPUT.write_bytes((json.dumps(report(), indent=2, sort_keys=True)+"\n").encode("utf-8"))
    print(OUTPUT)
