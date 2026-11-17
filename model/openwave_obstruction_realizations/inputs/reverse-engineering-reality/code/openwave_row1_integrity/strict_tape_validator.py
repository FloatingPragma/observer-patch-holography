#!/usr/bin/env python3
"""Strict native-tape validator for claimed readouts of the two-load route.

Pinned object (not imported, not modified):
    openwave-tuning/next/ROW1_NATIVE/routes/multiload/multiload_readout.py
        SHA-256 442d697747ee4adbbcafea5c0a5676e848cbb94173c652042893964145d94eef
    openwave-tuning/next/ROW1_NATIVE/routes/multiload/MULTILOAD_READOUT.json
        SHA-256 b70871dbefcb731cb7b956a7e8a8b89c722fbfcc49781f571c92a167b73f8e1f
    pointer policy of routes/auditor_bridge/native_cycle_meter.py
        SHA-256 4db8d707879c416ca5002b782a70e07aa71d4d9a7fe7111b76a64828c1a05329

A claimed readout is the dictionary returned by ``multiload_readout.execute`` together with the inputs of
that call (``with_inputs``).  ``validate(result)`` returns ``(ok, reasons)``; ``ok`` is True only if every
rule below holds.  Each reason starts with the name of the rule that failed.

    R0  the claim is complete and well typed, and the preparation is admitted (twelve nonnegative
        integers of total two)
    R1  the source tape replays under the native rule from the stated initial loads, event by event
        (seam, coin, current), to the stated final loads
    R2  the regional class equals the one rebuilt from the emitted currents only (twelve counters modulo 6,
        zero start); class 3 is refused
    R3  the pointer tape replays under the native rule on the pointer carrier from its stated initial port:
        discovery scan, routing to port 0, the requested cycle traversal, in this order, each attempt a
        native event with its coin and current, the cursor advancing only after a moving current, and the
        tape ends exactly where the policy ends
    R4  the claimed winding equals the signed count of cycle steps in that replayed tape divided by six, it
        equals the traversal requested by the rebuilt class, and the pointer's final loads are the replayed
        ones
    R5  the fine-history probability equals 2^(-attempts), the attempt counts add up, and the number of
        pointer attempts is within the stated cap (30 + 3 + 12) * retries
    R6  after a timeout no winding and no charge are reported, and the stated stage is the replayed one
    R7  the stated conditional success probability equals p_reset(x) * (1 - 2^-R)^(6|w|)

Standard library only; plain integers and fractions.Fraction.  Nothing is read from disk.
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction

__all__ = ["validate", "validate_pointer_tape", "with_inputs", "native", "SEAMS", "CYCLE", "PINS",
           "route_to_zero", "reset_success", "attempt_cap"]

PINS = {
    "routes/multiload/multiload_readout.py": "442d697747ee4adbbcafea5c0a5676e848cbb94173c652042893964145d94eef",
    "routes/multiload/MULTILOAD_READOUT.json": "b70871dbefcb731cb7b956a7e8a8b89c722fbfcc49781f571c92a167b73f8e1f",
    "routes/auditor_bridge/native_cycle_meter.py": "4db8d707879c416ca5002b782a70e07aa71d4d9a7fe7111b76a64828c1a05329",
}

# the thirty seams of the twelve-port carrier in the order of the scan (run.py compares this tuple with the
# push kernel and with the route's own EDGES)
SEAMS = ((0, 1), (0, 5), (0, 7), (0, 10), (0, 11), (1, 5), (1, 7), (1, 8), (1, 9), (2, 3), (2, 4), (2, 6),
         (2, 10), (2, 11), (3, 4), (3, 6), (3, 8), (3, 9), (4, 5), (4, 9), (4, 11), (5, 9), (5, 11), (6, 7),
         (6, 8), (6, 10), (7, 8), (7, 10), (8, 9), (10, 11))
SEAM_SET = frozenset(SEAMS)
CYCLE = (0, 1, 8, 6, 2, 10)
PORTS = 12
SUCCESS_REASON = "complete_closed_pointer_word"
TIMEOUTS = ("discovery_timeout", "routing_timeout", "measurement_timeout")


def native(loads, seam, coin):
    """d = n_a - n_b; j = (d - (2c - 1)(|d| mod 2)) // 2; n_a -= j; n_b += j.  Returns (loads, j)."""
    a, b = seam
    d = loads[a] - loads[b]
    j = (d - (2 * coin - 1) * (abs(d) % 2)) // 2
    out = list(loads)
    out[a] -= j
    out[b] += j
    return tuple(out), j


def route_to_zero(start):
    """Shortest path from start to port 0 with the tie-breaking of the pinned meter: breadth-first search
    from start, neighbours in seam order, first discovery fixes the predecessor."""
    previous = {start: None}
    todo = deque([start])
    while todo:
        v = todo.popleft()
        if v == 0:
            break
        for a, b in SEAMS:
            if v != a and v != b:
                continue
            w = a if b == v else b
            if w not in previous:
                previous[w] = v
                todo.append(w)
    path = []
    current = 0
    while current != start:
        path.append(current)
        current = previous[current]
    return tuple(reversed(path))


def reset_success(port, retries):
    """p_x = sum_k 2^(-Rk) (1 - 2^-R)^(1 + d(y_k, 0)) over the five incident seams in scan order."""
    hold = Fraction(1, 2 ** retries)
    total = Fraction(0)
    k = 0
    for a, b in SEAMS:
        if port != a and port != b:
            continue
        y = b if port == a else a
        total += hold ** k * (1 - hold) ** (1 + len(route_to_zero(y)))
        k += 1
    return total


def attempt_cap(retries):
    """The cap stated by the route: (30 + 3 + 12) attempts blocks of `retries` (720 at retries = 16)."""
    return (30 + 3 + 12) * retries


def with_inputs(result, G, schedule, source_coins, region, pointer_initial, retries=16):
    """The claimed readout: the dictionary returned by multiload_readout.execute plus the inputs of the call."""
    claim = dict(result)
    claim["inputs"] = {"G": list(G), "schedule": list(schedule), "source_coins": list(source_coins),
                       "region": list(region), "pointer_initial": pointer_initial, "retries": retries}
    return claim


def _is_int(x):
    return type(x) is int


def _decode(residue):
    if residue == 3:
        return None
    return residue if residue <= 2 else residue - 6


def validate_pointer_tape(meter, pointer_initial, commanded, retries):
    """Rules R3 to R6 for one pointer result.  Returns (reasons, replay) where replay holds the outcome of the
    replayed policy: {'outcome', 'signed_steps', 'final_loads', 'attempts'}."""
    reasons = []
    if not isinstance(meter, dict):
        return ["R0: pointer result is not a record"], None
    for key in ("success", "reason", "records", "final_loads", "retained_signed_cycle_steps", "winding",
                "fine_history_probability"):
        if key not in meter:
            reasons.append(f"R0: pointer result lacks '{key}'")
    recs = meter.get("records")
    if not isinstance(recs, list):
        reasons.append("R0: pointer tape is not a list of records")
        recs = []
    if not (_is_int(pointer_initial) and 0 <= pointer_initial < PORTS):
        return reasons + ["R0: stated pointer initial port is not a port"], None
    if not (_is_int(retries) and 1 <= retries <= 128):
        return reasons + ["R0: stated retries outside 1..128"], None
    if not (_is_int(commanded) and -2 <= commanded <= 2):
        return reasons + ["R0: requested winding outside the two-load interval"], None

    loads = tuple(int(v == pointer_initial) for v in range(PORTS))
    state = {"i": 0, "loads": loads, "signed": 0, "broken": False}

    def take(stage, seam, direction):
        """Consume one record that must be an attempt of `stage` on `seam`.  Returns the native current, or
        None if the tape has ended or the record does not belong to the policy at this point."""
        i = state["i"]
        if i >= len(recs):
            return None
        rec = recs[i]
        if not isinstance(rec, dict) or any(k not in rec for k in ("stage", "seam", "coin", "current")):
            reasons.append(f"R3: pointer record {i} is not a complete (stage, seam, coin, current) record")
            state["broken"] = True
            return None
        rseam = tuple(rec["seam"]) if isinstance(rec["seam"], (list, tuple)) else rec["seam"]
        if rseam not in SEAM_SET:
            reasons.append(f"R3: pointer record {i}: {rec['seam']} is not a seam of the pointer carrier")
            state["broken"] = True
            return None
        if rec["stage"] != stage or rseam != seam:
            reasons.append(f"R3: pointer record {i}: policy requires stage '{stage}' on seam {list(seam)}, "
                           f"tape has '{rec['stage']}' on {list(rseam)}")
            state["broken"] = True
            return None
        if not _is_int(rec["coin"]) or rec["coin"] not in (0, 1):
            reasons.append(f"R3: pointer record {i}: coin {rec['coin']!r} is not a native coin")
            state["broken"] = True
            return None
        if rec.get("requested_cycle_direction") != direction:
            reasons.append(f"R3: pointer record {i}: requested cycle direction {rec.get('requested_cycle_direction')!r}, "
                           f"policy has {direction!r}")
        before = state["loads"]
        after, j = native(before, seam, rec["coin"])
        if not _is_int(rec["current"]) or rec["current"] != j:
            reasons.append(f"R3: pointer record {i}: current {rec['current']!r} is not the native current {j} "
                           f"of seam {list(seam)} with coin {rec['coin']}")
            state["broken"] = True
            return None
        if j and stage == "measure":
            p = before.index(1)
            q = after.index(1)
            kp, kq = CYCLE.index(p), CYCLE.index(q)
            state["signed"] += 1 if (kp + 1) % 6 == kq else -1
        state["loads"] = after
        state["i"] = i + 1
        return j

    outcome = None
    known = None
    for seam in SEAMS:                                    # discovery scan
        for _ in range(retries):
            j = take("discover", seam, None)
            if j is None:
                outcome = "incomplete"
                break
            if j:
                known = seam[1] if j > 0 else seam[0]
                break
        if outcome or known is not None:
            break
    if outcome is None and known is None:
        outcome = "discovery_timeout"
    if outcome is None:                                   # routing to port 0
        for dest in route_to_zero(known):
            seam = (min(known, dest), max(known, dest))
            moved = False
            for _ in range(retries):
                j = take("reset_route", seam, None)
                if j is None:
                    outcome = "incomplete"
                    break
                if j:
                    moved = True
                    known = dest
                    break
            if outcome:
                break
            if not moved:
                outcome = "routing_timeout"
                break
    if outcome is None:                                   # the requested traversal: 6|w| steps
        cursor = 0
        direction = 1 if commanded > 0 else -1
        for _ in range(6 * abs(commanded)):
            nxt = (cursor + direction) % 6
            seam = (min(CYCLE[cursor], CYCLE[nxt]), max(CYCLE[cursor], CYCLE[nxt]))
            moved = False
            for _ in range(retries):
                j = take("measure", seam, direction)
                if j is None:
                    outcome = "incomplete"
                    break
                if j:
                    moved = True
                    cursor = nxt
                    break
            if outcome:
                break
            if not moved:
                outcome = "measurement_timeout"
                break
        if outcome is None:
            outcome = "success"
    if outcome == "incomplete" and not state["broken"]:
        reasons.append(f"R3: pointer tape ends after {len(recs)} attempts before the policy terminates")
    if outcome != "incomplete" and state["i"] != len(recs):
        reasons.append(f"R3: pointer tape has {len(recs) - state['i']} records after the policy terminated")
    replay = {"outcome": outcome, "signed_steps": state["signed"], "final_loads": list(state["loads"]),
              "attempts": state["i"]}

    claimed_success = meter.get("success")
    winding = meter.get("winding")
    if claimed_success is True:
        if outcome != "success":
            reasons.append(f"R4: success is claimed, the replayed tape gives '{outcome}'")
        if not _is_int(winding):
            reasons.append("R4: success is claimed without an integer winding")
        elif state["signed"] != 6 * winding:
            reasons.append(f"R4: claimed winding {winding}, the replayed tape holds {state['signed']} signed cycle "
                           f"steps ({Fraction(state['signed'], 6)} turns)")
        if _is_int(winding) and winding != commanded:
            reasons.append(f"R4: claimed winding {winding} is not the requested traversal {commanded}")
        if meter.get("reason") != SUCCESS_REASON:
            reasons.append(f"R4: success is claimed with reason {meter.get('reason')!r}")
    else:
        if claimed_success is not False:
            reasons.append("R0: 'success' is not a boolean")
        if winding is not None:
            reasons.append(f"R6: winding {winding!r} is reported without a completed traversal")
        if outcome == "success":
            reasons.append("R6: a timeout is claimed, the replayed tape completes the traversal")
        elif outcome in TIMEOUTS and meter.get("reason") != outcome:
            reasons.append(f"R6: stated stage {meter.get('reason')!r}, replayed stage '{outcome}'")
    if "retained_signed_cycle_steps" in meter and meter["retained_signed_cycle_steps"] != state["signed"]:
        reasons.append(f"R4: retained signed cycle steps {meter['retained_signed_cycle_steps']!r}, replayed "
                       f"{state['signed']}")
    if "final_loads" in meter and list(meter["final_loads"]) != list(state["loads"]) and not state["broken"] \
            and outcome != "incomplete":
        reasons.append("R4: stated final pointer loads are not the replayed ones")
    if "final_loads" in meter and (state["broken"] or outcome == "incomplete"):
        fl = meter["final_loads"]
        if not (isinstance(fl, list) and len(fl) == PORTS and all(_is_int(v) for v in fl) and sum(fl) == 1
                and min(fl) == 0):
            reasons.append("R4: stated final pointer loads are not one unit on one port")
    if "fine_history_probability" in meter and meter["fine_history_probability"] != str(Fraction(1, 2 ** len(recs))):
        reasons.append(f"R5: pointer fine-history probability {meter['fine_history_probability']} is not "
                       f"2^-{len(recs)}")
    if len(recs) > attempt_cap(retries):
        reasons.append(f"R5: {len(recs)} pointer attempts exceed the stated cap {attempt_cap(retries)}")
    return reasons, replay


def validate(result, inputs=None):
    """Accept a claimed readout of multiload_readout.execute only if rules R0 to R7 hold.

    result : the dictionary returned by execute, with the inputs of the call under result['inputs']
             (see with_inputs) or passed as `inputs`:
             {'G', 'schedule', 'source_coins', 'region', 'pointer_initial', 'retries'}
    returns: (ok, reasons)   reasons is a list of strings, empty when ok
    """
    reasons = []
    if not isinstance(result, dict):
        return False, ["R0: the claim is not a record"]
    inp = inputs if inputs is not None else result.get("inputs")
    if not isinstance(inp, dict) or any(k not in inp for k in ("G", "region", "pointer_initial", "retries")):
        return False, ["R0: the inputs of the call (G, region, pointer_initial, retries) are not stated"]
    for key in ("source_records", "counter_class", "meter", "decoded_spatial_degree", "decoded_relative_charge",
                "combined_attempts", "fine_joint_probability", "audit_only"):
        if key not in result:
            reasons.append(f"R0: the claim lacks '{key}'")
    G = inp["G"]
    if not (isinstance(G, (list, tuple)) and len(G) == PORTS and all(_is_int(v) and v >= 0 for v in G)
            and sum(G) == 2):
        return False, reasons + ["R0: preparation is not twelve nonnegative integers of total two"]
    G = tuple(G)
    region = inp["region"]
    if not (isinstance(region, (list, tuple)) and len(set(region)) == len(region)
            and all(_is_int(v) and 0 <= v < PORTS for v in region)):
        return False, reasons + ["R0: region is not a set of ports"]
    retries = inp["retries"]
    audit = result.get("audit_only") if isinstance(result.get("audit_only"), dict) else {}
    if "initial_G" in audit and tuple(audit["initial_G"]) != G:
        reasons.append("R1: stated initial loads differ from the preparation of the call")

    # R1: source tape under the native rule; R2: counters from the emitted currents only
    records = result.get("source_records")
    if not isinstance(records, list):
        reasons.append("R0: source tape is not a list of records")
        records = []
    schedule = inp.get("schedule")
    coins = inp.get("source_coins")
    if schedule is not None and len(schedule) != len(records):
        reasons.append(f"R1: source tape has {len(records)} records for a schedule of {len(schedule)} events")
    loads = G
    counters = [0] * PORTS
    for k, rec in enumerate(records):
        if not isinstance(rec, dict) or any(key not in rec for key in ("seam", "coin", "current")):
            reasons.append(f"R1: source record {k} is not a complete (seam, coin, current) record")
            continue
        e, c, j = rec["seam"], rec["coin"], rec["current"]
        if not (_is_int(e) and 0 <= e < len(SEAMS)) or not (_is_int(c) and c in (0, 1)) or not _is_int(j):
            reasons.append(f"R1: source record {k} is outside the native alphabet")
            continue
        if schedule is not None and k < len(schedule) and schedule[k] != e:
            reasons.append(f"R1: source record {k} is on seam {e}, the schedule has {schedule[k]}")
        if coins is not None and k < len(coins) and coins[k] != c:
            reasons.append(f"R1: source record {k} has coin {c}, the call has {coins[k]}")
        loads, jn = native(loads, SEAMS[e], c)
        if jn != j:
            reasons.append(f"R1: source record {k}: current {j} is not the native current {jn}")
        a, b = SEAMS[e]
        counters[a] = (counters[a] + j) % 6               # the emitted current, nothing else
        counters[b] = (counters[b] - j) % 6
    if "final_n" in audit:
        if list(audit["final_n"]) != list(loads):
            reasons.append("R1: the source tape does not replay to the stated final loads")
        elif any(counters[v] != (G[v] - loads[v]) % 6 for v in range(PORTS)):
            reasons.append("R2: counters rebuilt from the currents differ from G minus the final loads modulo 6")
    else:
        reasons.append("R0: final source loads are not stated")
    rebuilt = sum(counters[v] for v in region) % 6
    if result.get("counter_class") != rebuilt or not _is_int(result.get("counter_class")):
        reasons.append(f"R2: stated regional class {result.get('counter_class')!r}, rebuilt from the currents {rebuilt}")
    if "counters" in result and list(result["counters"]) != counters:
        reasons.append("R2: stated counters differ from the counters rebuilt from the currents")
    if "expected_w" in audit and "final_n" in audit and len(audit["final_n"]) == PORTS:
        if audit["expected_w"] != sum(G[v] - audit["final_n"][v] for v in region):
            reasons.append("R1: stated expected winding is not G(A) minus the stated final load of the region")
    commanded = _decode(rebuilt)
    if commanded is None:
        return False, reasons + ["R2: rebuilt class 3 lies outside the two-load interval; no readout exists"]

    # R3, R4, R6 and the pointer part of R5
    meter = result.get("meter")
    pointer_reasons, replay = validate_pointer_tape(meter, inp["pointer_initial"], commanded, retries)
    reasons += pointer_reasons
    if replay is None:
        return False, reasons
    meter = meter if isinstance(meter, dict) else {}
    success = meter.get("success") is True
    winding = meter.get("winding")
    want_degree = winding if success else None
    want_charge = -winding if success and _is_int(winding) else None
    if result.get("decoded_spatial_degree") != want_degree:
        rule = "R4" if success else "R6"
        reasons.append(f"{rule}: decoded spatial degree {result.get('decoded_spatial_degree')!r}, pointer result "
                       f"gives {want_degree!r}")
    if result.get("decoded_relative_charge") != want_charge:
        rule = "R4" if success else "R6"
        reasons.append(f"{rule}: decoded relative charge {result.get('decoded_relative_charge')!r}, pointer result "
                       f"gives {want_charge!r}")

    # R5: attempts and fine weights
    n_pointer = len(meter.get("records", [])) if isinstance(meter.get("records"), list) else 0
    total = len(records) + n_pointer
    if result.get("combined_attempts") != total:
        reasons.append(f"R5: combined attempts {result.get('combined_attempts')!r}, the tapes hold {total}")
    if result.get("fine_joint_probability") != str(Fraction(1, 2 ** total)):
        reasons.append(f"R5: joint fine-history probability {result.get('fine_joint_probability')} is not 2^-{total}")

    # R7: the stated conditional law
    if "conditional_meter_success" in result:
        p = reset_success(inp["pointer_initial"], retries) * (1 - Fraction(1, 2 ** retries)) ** (6 * abs(commanded))
        if result["conditional_meter_success"] != str(p) or result.get("conditional_meter_timeout") != str(1 - p):
            reasons.append("R7: stated conditional success or timeout probability is not p_reset(x)(1-2^-R)^(6|w|)")
    return (not reasons), reasons
