"""M10.4 common-state bridge verification, not the full common-controller archive.

Native fine instruments -> actual cooling/counter density -> complete choice
 erasure -> two-arm subdivided rotation -> complete re-encoding -> Ramsey read.
The full packet additionally checks global integer embeddings, current-only
entry/routing, all prefix records and native work. No file/network/process access.
"""

import json
import math
from fractions import Fraction
import numpy as np

TRI = ((0, 1), (1, 2), (0, 2))
FLAG_CHECKS = 0
CALLS = {"native": 0, "swap_hold": 0, "swap": 0}
MAX_PROBABILITY_DEFECT = 0.0


def native(q, edge, coin, law="native"):
    CALLS[law] += 1
    a, b = edge
    d = q[a] - q[b]
    j = (
        (d - (2 * coin - 1) * (abs(d) % 2)) // 2
        if law == "native"
        else d if law == "swap" else coin * d
    )
    out = list(q)
    out[a] -= j
    out[b] += j
    return tuple(out), [a, b, coin, j]


def fine(edges, law, vacuum=False):
    branches = []
    columns = []
    for e, edge in enumerate(edges):
        for c in (0, 1):
            matrices = np.zeros((3, 3, 3))
            for p in range(3):
                q = tuple(int(x == p) for x in range(3))
                out, rec = native((0, 0, 0) if vacuum else q, edge, c, law)
                j = rec[3]
                y = p if vacuum else out.index(1)
                matrices[j + 1, y, p] = 1
                columns.append([e, c, p, j, y])
            branches.extend(matrices)
    return np.asarray(branches) / math.sqrt(2 * len(edges)), columns


def helmert(branches):
    count = len(branches)
    out = [sum(branches) / math.sqrt(count)]
    prefix = branches[0].copy()
    for k in range(1, count):
        out.append((prefix - k * branches[k]) / math.sqrt(k * (k + 1)))
        prefix += branches[k]
    return np.asarray(out)


def instrument(edges, law):
    return helmert(fine(edges, law)[0])


def contrast(edges, law):
    target = instrument(edges, law)
    blank = helmert(fine(edges, law, True)[0])
    return np.asarray([(b + sign * a) / 2 for a, b in zip(target, blank) for sign in (1, -1)])


def routing(slot, m):
    v, k = divmod(slot % (3 * m), m)
    w = (v + 1) % 3
    return [e for e in TRI if v in e] * (m - k) + [e for e in TRI if w in e] * k


def channel(rho, ops, accept):
    outs = np.einsum("kxa,ab,kyb->kxy", ops, rho, ops.conj())
    probabilities = np.trace(outs, axis1=1, axis2=2).real
    selected = sum(outs[k] for k in accept)
    p = float(np.trace(selected).real)
    return selected / p, p, probabilities


def clock(k):
    edges = ((0, 1), (0, 1), (0, 2), (0, 2), (1, 2))
    R = np.eye(33)
    columns = []
    for p in range(3):
        for e, edge in enumerate(edges):
            for c in (0, 1):
                q = tuple(int(x == p) for x in range(3))
                out, rec = native(q, edge, c)
                j = rec[3]
                y = out.index(1)
                cc = c ^ int(j != 0)
                i = 11 * p + 2 * e + c
                o = 11 * y + 2 * e + cc
                R[:, i] = 0
                R[o, i] = 1
                columns.append([p, e, c, j, y, cc])
    active = 0
    for word in range(4**k):
        for shift in range(2*k):
            out, record = native((0, 0), (0, 1), (word >> shift) & 1)
            assert out == (0, 0) and record[3] == 0
        active += word == 0
    assert active == 1
    eps = math.sqrt(active / 4**k)
    a = np.r_[np.full(10, eps / math.sqrt(10)), math.sqrt(1 - eps * eps)]
    A = np.kron(np.eye(3), a[:, None])
    P = A.T @ R @ A
    C = np.kron(np.eye(3), 2 * np.outer(a, a) - np.eye(11))
    W = R @ C
    v = a.copy()
    v[0] -= 1
    H = np.eye(11) - 2 * np.outer(v, v) / (v @ v)
    decoders = np.asarray([np.kron(np.eye(3), H[j : j + 1]) for j in range(11)])
    Q = 256 * 2**k
    waves = []
    live = A.astype(complex)
    for _ in range(Q):
        waves.append(live.copy())
        live = W @ live
    ops = np.fft.fft(np.asarray(waves), axis=0) / Q
    return {
        "R": R,
        "A": A,
        "P": P,
        "W": W,
        "a": a,
        "decode": decoders,
        "qpe": ops,
        "Q": Q,
        "epsilon": eps,
        "columns": columns,
        "edges": edges,
    }


def hop(old, desired, ring, graph_size, cap=8):
    edge = tuple(sorted((old, desired)))
    physical = tuple(ring[x] for x in edge)
    live = [([], [], 1.0)]
    history = []
    for attempt in range(cap):
        nxt = []
        for coins, currents, p in live:
            for c in (0, 1):
                q = [0] * graph_size
                q[ring[old]] = 1
                out, rec = native(q, physical, c)
                j = rec[3]
                co = coins + [c]
                cu = currents + [j]
                prob = p / 2
                if j or attempt == cap - 1:
                    decoded = edge[1] if j == 1 else edge[0] if j == -1 else old
                    history.append(
                        {
                            "coins": co,
                            "currents": cu,
                            "probability": prob,
                            "decoded": decoded,
                            "outcome": "hop" if j else "timeout",
                        }
                    )
                else:
                    nxt.append((co, cu, prob))
        live = nxt
    success = sum(r["probability"] for r in history if r["outcome"] == "hop")
    return success, {"edge": list(edge), "physical_edge": list(physical), "histories": history}


def apply_pair(rho, rot, ref, index):
    # All maps have one COMMON eraser flag; neither arm is measured separately.
    probs = np.zeros(len(rot))
    outs = []
    for a, b in zip(rot, ref):
        op = np.zeros((6, 6), complex)
        op[:3, :3] = a
        op[3:, 3:] = b
        out = op @ rho @ op.conj().T
        outs.append(out)
        probs[len(outs) - 1] = np.trace(out).real
    p = float(probs[index])
    return outs[index] / p, p, probs


def ramsey(rho, W, stride):
    unit = np.linalg.matrix_power(W, stride)
    power = np.eye(len(W), dtype=complex)
    records = []
    for n in range(17):
        z = np.trace(rho @ power)
        records.append(
            {"time": n, "x_plus": float((1 + z.real) / 2), "y_plus": float((1 + z.imag) / 2)}
        )
        power = unit @ power
    return records


def encode(value):
    array = np.asarray(value)
    return {"real": array.real.tolist(), "imag": array.imag.tolist()}


def reencode(rho, cl):
    blocks = rho.reshape(2, 3, 2, 3)
    ops = cl["qpe"]
    pos = ops[1 : cl["Q"] // 2]
    output = np.einsum("kxa,iajb,kyb->ixjy", pos, blocks, pos.conj(), optimize=True).reshape(
        66, 66
    )
    gram = np.einsum("kxa,kxb->kab", ops.conj(), ops)
    marginal = blocks[0, :, 0, :] + blocks[1, :, 1, :]
    probs = np.einsum("kab,ba->k", gram, marginal).real
    p = float(np.trace(output).real)
    return output / p, p, probs


def selector_success(count):
    """All native blank words, with retained coherent word and capped retry."""
    width = (count - 1).bit_length()
    base = 2**width
    for word in range(base):
        for shift in range(width):
            out, record = native((0, 0), (0, 1), (word >> shift) & 1)
            assert out == (0, 0) and record[3] == 0
    bad = Fraction(base - count, base)
    success = sum(bad**t * Fraction(count, base) for t in range(8))
    assert success + bad**8 == 1
    return success


def check_probabilities(probs):
    global MAX_PROBABILITY_DEFECT, FLAG_CHECKS
    FLAG_CHECKS += len(probs)
    defect = abs(float(sum(probs))-1)
    MAX_PROBABILITY_DEFECT = max(MAX_PROBABILITY_DEFECT, defect)
    assert defect < 1e-10


def interval(records, eta):
    """Record-only unsigned frequency bounds; eta is a separate certificate."""
    r = 2 * records[1]["x_plus"] - 1
    lower_cos = max(-1., (r - eta) / (1 - eta))
    upper_cos = min(1., (r + eta) / (1 - eta))
    low, high = math.acos(upper_cos), math.acos(lower_cos)
    return {"fixed_eta": eta, "lower": low, "upper": high, "width": high-low}


def run_member(refinement):
    cl = clock(refinement)
    weighted = cl["edges"]
    rho = np.diag([0.0, 1.0, 0.0])
    logprep = 0.0
    for ops, index in [(contrast(weighted, "native"), 1)] + [(instrument(weighted, "native"), 0)] * 8:
        rho, p, probs = channel(rho, ops, [index])
        check_probabilities(probs)
        logprep += math.log10(p)
    rho, p, probs = channel(rho, cl["qpe"], range(1, cl["Q"] // 2))
    logprep += math.log10(p)
    check_probabilities(probs)
    rho, p, probs = channel(rho, cl["decode"], [0])
    logprep += math.log10(p)
    check_probabilities(probs)
    permutation = np.eye(3)[[1, 2, 0]]
    position = permutation @ rho @ permutation.T
    results = []
    s3 = selector_success(3)
    s5 = selector_success(5)
    for m in (4, 8):
        assert selector_success(2 * m) == 1
        for law in ("native", "swap"):
            sources = [instrument(routing(slot, m), law) for slot in range(3 * m)]
            plane = contrast(TRI, law)
            for turns in (1, 2):
                for direction in (1, -1):
                    rho = np.kron(np.ones((2, 2)) / 2, position)
                    reach = logprep
                    rot = ref = 2 * m
                    for step in range(1, 3 * m * turns + 1):
                        k = step % m
                        wanted = (
                            (2 * m + direction * step) % (3 * m),
                            2 * m,
                        )
                        for old, new in zip((rot, ref), wanted):
                            if old == new:
                                continue
                            p, hist = hop(old, new, list(range(3 * m)), 3 * m)
                            assert abs(sum(h["probability"] for h in hist["histories"]) - 1) < 1e-12
                            reach += math.log10(p)
                        rot, ref = wanted
                        for _ in range(4):
                            rho, p, probs = apply_pair(rho, sources[rot], sources[ref], 0)
                            check_probabilities(probs)
                            reach += math.log10(p)
                        rho, p, probs = apply_pair(rho, plane, plane, 1)
                        check_probabilities(probs)
                        reach += math.log10(p)
                    change = np.kron(np.eye(2), permutation.T)
                    rho = change @ rho @ change.T
                    final, p, probs = reencode(rho, cl)
                    reach += math.log10(p)
                    check_probabilities(probs)
                    blocks = final.reshape(2, 33, 2, 33)
                    visibility = 2 * np.trace(blocks[0, :, 1, :]).real
                    sign = -1 if law == "native" and turns == 1 else 1
                    selected_path_probability = (1+sign*visibility)/2
                    assert 0 < selected_path_probability <= 1+1e-12
                    signed_normalized_coherence = float(np.trace(blocks[0,:,1,:]).real /
                        math.sqrt(np.trace(blocks[0,:,0,:]).real*np.trace(blocks[1,:,1,:]).real))
                    arc = (
                        blocks[0, :, 0, :]
                        + blocks[1, :, 1, :]
                        + sign * (blocks[0, :, 1, :] + blocks[1, :, 0, :])
                    ) / (1 + sign * visibility)
                    records = ramsey(arc, cl["W"], 2**refinement)
                    values, vectors = np.linalg.eigh(cl["P"])
                    value, vector = values[1], vectors[:, 1]
                    u = cl["A"] @ vector
                    w = (cl["R"] @ u - value*u) / math.sqrt(1-value*value)
                    plus, minus = (u-1j*w)/math.sqrt(2), (u+1j*w)/math.sqrt(2)
                    populations = [float(np.vdot(v, arc @ v).real) for v in (plus, minus)]
                    impurity = 1-sum(populations)
                    pair_impurity = 1-sum(p*p for p in populations)
                    single_interval = interval(records, .002)
                    pair_records = []
                    for record in records:
                        z = 2*record["x_plus"]-1 + 1j*(2*record["y_plus"]-1)
                        pair_records.append({"time":record["time"], "x_plus":float((1+(z*z).real)/2),
                                             "y_plus":float((1+(z*z).imag)/2)})
                    pair_interval = interval(pair_records, .01)
                    assert impurity < .002 and pair_impurity < .01
                    expected = 2**(refinement+1) * math.asin(2**(-refinement) / math.sqrt(5))
                    assert single_interval["lower"] <= expected <= single_interval["upper"]
                    assert pair_interval["lower"] <= 2*expected <= pair_interval["upper"]
                    heldout = max(abs(2*r["x_plus"]-1-math.cos(r["time"]*expected)) for r in records[2:])
                    assert heldout <= .004+1e-10
                    assert (
                        (visibility < 0 if turns == 1 else visibility > 0)
                        if law == "native"
                        else visibility > 0
                    )
                    factor = s5**11 * s3 ** (3 * m * turns)
                    reach += math.log10(factor.numerator) - math.log10(factor.denominator)
                    results.append(
                        {
                            "k": refinement,
                            "M": m,
                            "law": law,
                            "turns": turns,
                            "direction": direction,
                            "visibility": float(visibility),
                            "source_prediction": expected,
                            "stable_single_interval": single_interval,
                            "stable_pair_interval": pair_interval,
                            "single_off_band": impurity,
                            "pair_off_band": pair_impurity,
                            "single_ramsey": records,
                            "pair_ramsey": pair_records,
                            "heldout_real_moment_error": heldout,
                            "log10_interference_success_given_classical_prefix": reach,
                            "selected_clock_path_probability": selected_path_probability,
                            "log10_clock_success_given_classical_prefix": reach+math.log10(selected_path_probability),
                            "log10_independent_pair_clock_success_given_prefixes": 2*(reach+math.log10(selected_path_probability)),
                            "signed_normalized_coherence": signed_normalized_coherence,
                        }
                    )
    return results

results = [case for k in range(7) for case in run_member(k)]
print(json.dumps({
    "candidate": "independent stationary-reference verification of common source",
    "reference": "stationary",
    "flag_probabilities_checked": FLAG_CHECKS,
    "verification_only": True,
    "cases": results,
    "literal_endpoint_calls": CALLS,
    "maximum_probability_sum_defect": MAX_PROBABILITY_DEFECT,
    "complete_outcome_probabilities_checked": True,
    "control_scope": "Always-swap rotation-stage controls; charge and clock preparation stay native.",
    "pair_scope": "Actual returned density tensor square; independently admitted global pair source is checked in the sealed full packet.",
    "full_current_only_entry_and_global_graph_replayed": False,
    "physical_award": False,
}, indent=2))
