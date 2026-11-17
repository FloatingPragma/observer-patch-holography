"""Recorded marker controller and complete, balanced rotating-observer process.

Imports frozen pilot primitives only; every added operation is exposed here.
"""
from fractions import Fraction as F
from math import sqrt
from model import (CALLS, EDGES, I, add, apply, complete, density, dot, endpoint,
                   eraser, fine, logprob, mm, mv, number, scaled, star, trace,
                   transpose)


def marker_hop(vertex, nxt, cap, law):
    """Enumerate stopped histories, using only addressed current to update reader."""
    edge = tuple(sorted((vertex, nxt)))
    active = [((), (), F(1))]
    histories = []
    for attempt in range(cap):
        continuing = []
        for coins, currents, probability in active:
            for coin in (0, 1):
                state = [int(x == vertex) for x in range(3)]
                a, b = edge
                old = state[a]
                state[a], state[b] = endpoint(state[a], state[b], coin, law)
                current = old - state[a]
                record = (*currents, current)
                tape = (*coins, coin)
                p = probability / 2
                if current:
                    # A reader knows its initial address, not simulator state.
                    read_position = b if current == 1 else a
                    assert read_position == nxt == state.index(1)
                    histories.append({"coins": tape, "currents": record,
                                      "probability": str(p), "outcome": "hop",
                                      "decoded_position": read_position})
                elif attempt + 1 == cap:
                    assert state.index(1) == vertex
                    histories.append({"coins": tape, "currents": record,
                                      "probability": str(p), "outcome": "timeout",
                                      "decoded_position": vertex})
                else:
                    assert state.index(1) == vertex
                    continuing.append((tape, record, p))
        active = continuing
    assert sum(F(h["probability"]) for h in histories) == 1
    success = sum(F(h["probability"]) for h in histories if h["outcome"] == "hop")
    return success, {"address": edge, "histories": histories}


def preparation(law, rho, start=2, full=False):
    leaves = tuple(v for v in range(3) if v != start)
    warm = fine(EDGES, law)
    edge = eraser(fine((leaves,), law))
    ref = eraser(fine((leaves,), law, vacuum=True))
    ops = tuple((add(ar, ae, sign), 4 * de)
                for (ae, de), (ar, dr) in zip(edge, ref) for sign in (1, -1))
    assert complete(ops)
    f = tuple(int(i == leaves[0]) - int(i == leaves[1]) for i in range(3))
    output, records = scaled(I, 0), []
    for w, m in enumerate(warm):
        mixed = scaled(apply(m, rho), F(1, 6))
        for k, (a, d) in enumerate(ops):
            part = scaled(apply(a, mixed), F(1, d))
            row = {"warmup_flag": w, "eraser_flag": k // 2,
                   "path": "+" if k % 2 == 0 else "-", "probability": str(trace(part))}
            if full:
                row["unnormalized_output_density"] = [[str(q) for q in r] for r in part]
            records.append(row)
            if k == 1:
                output = add(output, part)
    probability = trace(output)
    assert sum(F(r["probability"]) for r in records) == 1
    assert scaled(output, 1 / probability) == density(f)
    return probability, {"input_density": [[str(q) for q in r] for r in rho],
                         "success_probability": str(probability),
                         "accepted_density": [[str(q) for q in r] for r in density(f)],
                         "records": records,
                         "preparation_kraus": [{"matrix": [[str(q) for q in r] for r in a],
                                                 "squared_denominator": d} for a, d in ops]}


def execute(h=4, turns=1, law="native", direction=1, cap=8, start=2,
            rho=None, balanced=True, full=True):
    rho = density((1, 0, 0)) if rho is None else rho
    preparation_p, preparation_record = preparation(law, rho, start, full)
    systems = tuple(star(v, law) for v in range(3))
    assert all(complete(ops) for p, ops in systems)
    leaves = tuple(v for v in range(3) if v != start)
    a = b = tuple(F(int(i == leaves[0]) - int(i == leaves[1])) for i in range(3))
    vertex, reach = start, preparation_p
    records, rejected, stage = [], 1 - reach, 0
    for hop in range(3 * turns):
        nxt = (vertex + direction) % 3
        success, tape = marker_hop(vertex, nxt, cap, law)
        timeout = reach * (1 - success)
        rejected += timeout
        records.append({"kind": "marker", "hop": hop,
                        "angle_degrees": 120 * direction,
                        "reach_probability": str(reach), "timeout_probability": str(timeout),
                        **tape})
        reach *= success
        # Every successful native tape has the same decoded address.
        vertex = next(r["decoded_position"] for r in tape["histories"] if r["outcome"] == "hop")
        p, ops = systems[vertex]
        ref, ref_ops = systems[start]
        for _ in range(h):
            norm = dot(a, a) + dot(b, b)
            probs, outputs = [], []
            for (g, d), (gref, dr) in zip(ops, ref_ops):
                aa, bb = mv(g, a), mv(gref, b)
                q = F(dot(aa, aa) + dot(bb, bb), d * norm)
                probs.append(q)
                outputs.append([str(x) for x in (*aa, *bb)])
            assert sum(probs) == 1
            failure = reach * (1 - probs[0])
            rejected += failure
            row = {"kind": "source", "stage": stage, "decoded_vertex": vertex,
                   "reach_probability": str(reach),
                   "conditional_flags": [str(q) for q in probs],
                   "failure_probability": str(failure)}
            if full:
                row["flag_output_numerators"] = outputs
                row["flag_squared_denominators"] = [d for g, d in ops]
            records.append(row)
            reach *= probs[0]
            a, b = mv(p, a), mv(ref, b)
            stage += 1
    norm = dot(a, a) + dot(b, b)
    raw_visibility = 2 * dot(a, b) / norm
    raw_success = reach
    if balanced:
        # Six native vacuum coins per turn: known projection-limit attenuation.
        # Coherent erase on rotating arm; computational coin read on reference.
        # Success B0=diag(1,1/sqrt(2)); reject B1=diag(0,1/sqrt(2)).
        original_reach, nr, nb = reach, dot(a, a), dot(b, b)
        for k in range(1, 6 * turns + 1):
            assert all(endpoint(0, 0, coin, law) == (0, 0) for coin in (0, 1))
            first_failure = original_reach * nb / norm / 2 ** k
            records.append({"kind": "reference_coin", "index": k,
                            "reach_probability": str(reach),
                            "first_failure_probability": str(first_failure),
                            "failure_output": "reference arm only, unchanged source ray",
                            "native_coin_probabilities": ["1/2", "1/2"]})
            reach -= first_failure
            rejected += first_failure
        b = tuple(x / 2 ** (3 * turns) for x in b)
    norm = dot(a, a) + dot(b, b)
    visibility = 2 * dot(a, b) / norm
    overlap = float(dot(a, b)) / sqrt(float(dot(a, a) * dot(b, b)))
    plus, minus = (1 + visibility) / 2, (1 - visibility) / 2
    assert rejected + reach == 1
    assert vertex == start
    return {"h": h, "turns": turns, "law": law, "direction": direction,
            "cap": cap, "start": start, "balanced": balanced,
            "visibility": number(visibility), "raw_visibility": number(raw_visibility),
            "ray_overlap": overlap, "conditional_path_plus": number(plus),
            "conditional_path_minus": number(minus), "success": number(reach),
            "raw_success": number(raw_success), "log10_success": logprob(reach),
            "expected_preparations_per_success": number(1 / reach),
            "final_arm_vectors": [[str(q) for q in a], [str(q) for q in b]],
            "all_rejections": number(rejected), "normalization_exact": True,
            "preparation": preparation_record if full else None,
            "records": records if full else None}


def committed_visibility(h, turns, start=2, direction=1):
    """All target fine flags committed before eraser. No amplitude summation."""
    leaves = [v for v in range(3) if v != start]
    f = tuple(int(i == leaves[0]) - int(i == leaves[1]) for i in range(3))
    cross = scaled(density(f), F(1, 2))
    vertex = start
    for hop in range(3 * turns):
        vertex = (vertex + direction) % 3
        rot = fine(tuple(e for e in EDGES if vertex in e), "native")
        ref = fine(tuple(e for e in EDGES if start in e), "native")
        for _ in range(h):
            total = scaled(I, 0)
            for a, b in zip(rot, ref):
                total = add(total, mm(mm(a, cross), transpose(b)))
            cross = scaled(total, F(1, 4))
    return {"visibility": number(2 * trace(cross)),
            "eraser_success_given_prep_and_marker": str(F(1, 12 ** (3 * h * turns))),
            "all_committed_flags_retained": True}
