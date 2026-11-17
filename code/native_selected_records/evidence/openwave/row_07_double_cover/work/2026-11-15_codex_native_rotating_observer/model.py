"""Complete finite native-current instruments. Stdlib, exact rational transfers."""
from fractions import Fraction as F
from math import log10, sqrt

CALLS = {"native": 0, "swap": 0}
EDGES = ((0, 1), (1, 2), (0, 2))
I = tuple(tuple(F(i == j) for j in range(3)) for i in range(3))


def endpoint(a, b, coin, law):
    CALLS[law] += 1
    if law == "native":
        lo = (a + b) // 2
        hi = a + b - lo
        return (lo, hi) if coin == 0 else (hi, lo)
    if law == "swap":
        return b, a
    raise ValueError(law)


def mv(a, v):
    return tuple(sum(x * y for x, y in zip(row, v)) for row in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def add(a, b, scale=1):
    return tuple(tuple(x + scale * y for x, y in zip(ar, br))
                 for ar, br in zip(a, b))


def scaled(a, s):
    return tuple(tuple(s * x for x in row) for row in a)


def transpose(a):
    return tuple(zip(*a))


def mm(a, b):
    return tuple(tuple(dot(row, col) for col in transpose(b)) for row in a)


def fine(edges, law, vacuum=False):
    """Integer branch matrices; each has amplitude 1/sqrt(2*len(edges))."""
    out = []
    for a, b in edges:
        for coin in (0, 1):
            bins = [[[0] * 3 for _ in range(3)] for _ in range(3)]
            for x in range(3):
                state = [int(j == x) for j in range(3)]
                old = state[a]
                if vacuum:
                    endpoint(0, 0, coin, law)
                    current, y = 0, x
                else:
                    state[a], state[b] = endpoint(state[a], state[b], coin, law)
                    current, y = old - state[a], state.index(1)
                bins[current + 1][y][x] = 1
            out.extend(tuple(tuple(row) for row in m) for m in bins)
    return tuple(out)


def eraser(branches):
    """Helmert outputs (integer numerator, denominator in squared amplitude)."""
    z = scaled(I, 0)
    prefix = z
    rows = []
    branch_count = len(branches) // 3
    for k, m in enumerate(branches):
        if k:
            rows.append((add(prefix, m, -k), branch_count * k * (k + 1)))
        prefix = add(prefix, m)
    return ((prefix, branch_count * len(branches)), *rows)


def star(vertex, law):
    edges = tuple(e for e in EDGES if vertex in e)
    branches = fine(edges, law)
    p = scaled(eraser(branches)[0][0], F(1, 4))
    return p, eraser(branches)


def complete(ops):
    total = scaled(I, 0)
    for a, d in ops:
        total = add(total, scaled(mm(transpose(a), a), F(1, d)))
    return total == I


def density(v):
    return tuple(tuple(F(x * y, dot(v, v)) for y in v) for x in v)


def trace(a):
    return sum(a[i][i] for i in range(3))


def apply(a, rho):
    return mm(mm(a, rho), transpose(a))


def prep(law, rho):
    warm = fine(EDGES, law)
    edge = eraser(fine((EDGES[0],), law))
    ref = eraser(fine((EDGES[0],), law, vacuum=True))
    ops = tuple((add(ar, ae, sign), 4 * de)
                for (ae, de), (ar, dr) in zip(edge, ref) for sign in (1, -1))
    assert complete(ops)
    output, records = scaled(I, 0), []
    for w, m in enumerate(warm):
        mixed = scaled(apply(m, rho), F(1, 6))
        for k, (a, d) in enumerate(ops):
            part = scaled(apply(a, mixed), F(1, d))
            p = trace(part)
            records.append({"warmup_flag": w, "eraser_flag": k // 2,
                            "path": "+" if k % 2 == 0 else "-", "probability": str(p)})
            if k == 1:
                output = add(output, part)
    probability = trace(output)
    assert sum(F(r["probability"]) for r in records) == 1
    assert scaled(output, 1 / probability) == density((1, -1, 0))
    return probability, records


def number(q):
    return {"exact": str(q), "float": float(q)}


def logprob(p):
    return log10(p.numerator) - log10(p.denominator) if p else None


def run(h, turns, law="native", direction=1, cap=8, full=False):
    """Readout consumes current-confirmed marker program and final path flags."""
    systems = tuple(star(v, law) for v in range(3))
    assert all(complete(ops) for p, ops in systems)
    a = b = (F(1), F(-1), F(0))
    vertex, reach, source_success = 2, F(1), F(1)
    marker_success = 1 - F(1, 2 ** cap) if law == "native" else F(1)
    records, rejected, stage = [], F(0), 0
    for hop in range(3 * turns):
        nxt = (vertex + direction) % 3
        timeout = reach * (1 - marker_success)
        rejected += timeout
        records.append({"kind": "marker", "from": vertex, "to": nxt,
                        "angle_degrees": 120 * direction,
                        "hop_attempt_probabilities": [str(F(1, 2 ** k)) if law == "native"
                                                      else str(F(k == 1)) for k in range(1, cap + 1)],
                        "timeout_probability": str(timeout)})
        reach *= marker_success
        vertex = nxt
        p, ops = systems[vertex]
        ref, ref_ops = systems[2]
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
            row = {"kind": "source", "stage": stage, "vertex": vertex,
                   "reach_probability": str(reach),
                   "conditional_flags": [str(q) for q in probs],
                   "failure_probability": str(failure)}
            if full:
                row["flag_output_numerators"] = outputs
                row["flag_squared_denominators"] = [d for g, d in ops]
            records.append(row)
            reach *= probs[0]
            source_success *= probs[0]
            a, b = mv(p, a), mv(ref, b)
            stage += 1
    norm = dot(a, a) + dot(b, b)
    visibility = 2 * dot(a, b) / norm
    overlap = float(dot(a, b)) / sqrt(float(dot(a, a) * dot(b, b)))
    plus, minus = (1 + visibility) / 2, (1 - visibility) / 2
    assert rejected + reach == 1
    assert vertex == 2
    return {"h": h, "turns": turns, "law": law, "direction": direction,
            "cap": cap, "visibility": number(visibility), "ray_overlap": overlap,
            "conditional_path_plus": number(plus), "conditional_path_minus": number(minus),
            "success_after_preparation": number(reach), "source_success": number(source_success),
            "log10_success_after_preparation": logprob(reach),
            "all_first_rejections": number(rejected), "normalization_exact": True,
            "final_arm_vectors": [[str(x) for x in a], [str(x) for x in b]],
            "records": records if full else None}
