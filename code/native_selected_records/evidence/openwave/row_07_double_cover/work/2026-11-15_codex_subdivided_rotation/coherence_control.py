"""Exact committed-record channel control on the same fine observer program."""
import json
from fractions import Fraction as F
from native import CALLS, I, dot, fine, logprob, mm, routing


def superoperator(left, right):
    return [[sum(a[i][k] * b[j][l] for a, b in zip(left, right))
             for k in range(3) for l in range(3)] for i in range(3) for j in range(3)]


def apply(operator, vector, denominator):
    return [dot(row, vector) / denominator for row in operator]


def trace(v):
    return v[0] + v[4] + v[8]


def plane(v):
    p = [[F(i == j) - F(1, 3) for j in range(3)] for i in range(3)]
    a = [v[3 * i:3 * i + 3] for i in range(3)]
    out = mm(mm(p, a), p)
    return [x / 48 for row in out for x in row]


results = []
for m in (4, 8):
    branches = [fine(routing(slot, m), "native") for slot in range(3 * m)]
    cache = {}
    for turns in (1, 2):
        f = (1, -1, 0)
        left = right = cross = [F(f[i] * f[j], 4) for i in range(3) for j in range(3)]
        prep = F(5, 144)
        left, right, cross = ([prep * x for x in v] for v in (left, right, cross))
        rejected, records = 1 - prep, []
        for step in range(1, 3 * m * turns + 1):
            k = step % m
            rot, ref = (2 * m + step) % (3 * m), (2 * m + min(k, m - k)) % (3 * m)
            for arm in ("rotating", "reference"):
                reach = trace(left) + trace(right)
                rejected += reach / 256
                records.append({"kind": "marker", "arm": arm, "reach": str(reach), "timeout": str(reach / 256)})
                left, right, cross = ([F(255, 256) * x for x in v] for v in (left, right, cross))
            for pair in ((rot, rot), (ref, ref), (rot, ref)):
                if pair not in cache:
                    cache[pair] = superoperator(branches[pair[0]], branches[pair[1]])
            for event in range(4):
                reach = trace(left) + trace(right)
                incoming = [list(map(str, v)) for v in (left, right, cross)]
                left = apply(cache[rot, rot], left, 48 * m * m)
                right = apply(cache[ref, ref], right, 48 * m * m)
                cross = apply(cache[rot, ref], cross, 48 * m * m)
                success = trace(left) + trace(right)
                assert success == reach / (12 * m)
                rejected += reach - success
                records.append({"kind": "committed_target", "rotating_slot": rot, "reference_slot": ref,
                                "input_joint_density_blocks": incoming, "reach": str(reach),
                                "first_rejection": str(reach - success), "success": str(success)})
            reach = trace(left) + trace(right)
            incoming = [list(map(str, v)) for v in (left, right, cross)]
            left, right, cross = (plane(v) for v in (left, right, cross))
            success = trace(left) + trace(right)
            rejected += reach - success
            records.append({"kind": "coherent_plane_filter", "reach": str(reach),
                            "input_joint_density_blocks": incoming, "first_rejection": str(reach - success)})
        success = trace(left) + trace(right)
        visibility = 2 * trace(cross) / success
        assert rejected + success == 1
        results.append({"M": m, "turns": turns, "visibility": float(visibility),
                        "exact_visibility": str(visibility), "success": str(success),
                        "log10_success": logprob(success), "all_rejections": str(rejected), "records": records})
print(json.dumps({"candidate": "A-native-subdivided-observer-rotation-v3", "diagnostic": "commit target fine flags",
                  "cases": results, "literal_endpoint_calls": CALLS,
                  "full_outcome_rule": "K(r,mu)=Helmert[mu,r]*diag(Mrot_r,Mref_r)/sqrt(4M); everyr,mu retained",
                  "compression": "exact9x9 channel sums, all input density blocks stored; no history sampling"}, indent=2))
