"""Source dependence of the reduced carrier Hamiltonian, memory and dark-port variance."""

import argparse
from collections import deque
import json
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv

from ...action_response import measure as action
from ...action_response.four_ticks import raw_effect, volume_ratio, parse_interval

HERE = Path(__file__).resolve().parent
PAIRS = [(0, 39), (1, 38)]
IDENTITY = sp.eye(4)


def read(matrices, source):
    def quad(e):
        values = dict(reference=raw_effect(e, source[0], sp.Rational("6e-7"), 0))
        values.update(action.source_jets(e, source))
        for field, key, err in (
            ("raw_o", "ordinary", "12e-7"),
            ("raw_p", "unpair", "18e-7"),
        ):
            values[field] = raw_effect(e, source[1][key], sp.Rational(err), 0)
        return values

    re = quad([(m + m.H) / 2 for m in matrices])
    im = quad([(m - m.H) / (2 * sp.I) for m in matrices])
    return {k: iv.mpc(re[k], im[k]) for k in re}


def product(a, b):
    return {
        k: (
            x * b[k] if k == "reference" else x * b["reference"] + a["reference"] * b[k]
        )
        for k, x in a.items()
    }


def serialize(vals, rho):
    out = {k: str(vals[k]) for k in ("reference", "ordinary", "pairing")}
    if 0 not in vals["ordinary"]:
        direct = vals["pairing"] / vals["ordinary"]
        shared = 1 / rho + (vals["raw_p"] - vals["raw_o"] / rho) / vals["ordinary"]
        ratio = iv.mpf([max(direct.a, shared.a), min(direct.b, shared.b)])
        out["collar_over_ordinary"] = str(ratio)
    return out


def blank_path(board, route, start, end):
    excluded = set(tuple(sorted(x)) for x in (route[:2], route[2:]))
    stream = deque([(start, [start])])
    seen = {start}
    while stream:
        p, path = stream.popleft()
        for v in sorted(board.neighbours[p]):
            if v not in seen and tuple(sorted((p, v))) not in excluded:
                if v == end:
                    return path + [int(v)]
                seen.add(v)
                stream.append((v, path + [int(v)]))
    raise ValueError("No head-preserving blank preparation route")


def effects(board, route, q):
    native = action.native
    values = {
        n: [] for n in ("left", "right", "third", "dark", "even_mean", "odd_mean")
    }
    z = sp.diag(1, -1)
    for parity in (0, 1):
        twice = native.generator(board, route, q, bool(parity))

        def block(p, r):
            return sp.Matrix(twice[4 * p : 4 * p + 4, 4 * r : 4 * r + 4].toarray()) / 2

        left, right = block(1, 0), block(39, 38)
        values["left"].append(left)
        values["right"].append(right)
        assert block(38, 1) == -IDENTITY / 2 and block(39, 0) == sp.zeros(4)
        # The only shortest three-hop route is 0,1,38,39. Verify it against
        # powers of the complete generator, without dropping diagonal energies.
        x = np.zeros((twice.shape[0], 4), dtype=np.int64)
        x[:4] = np.eye(4, dtype=int)
        assert int(abs(twice).sum(axis=1).max()) ** 3 < 2**63
        for _ in range(3):
            assert not np.any(x[156:160])
            x = twice @ x
        third = sp.Matrix(x[156:160]) / 8
        assert third == -right * left / 2
        values["third"].append(third)
        cs = -sp.Rational(parity, 4) * IDENTITY
        tag = q * (
            (2 + sp.Rational(parity, 4)) * sp.kronecker_product(z, sp.eye(2))
            + 2 * sp.kronecker_product(sp.eye(2), z)
        )
        ds = []
        for idx, (p, r) in enumerate(PAIRS):
            delta = block(p, p) - block(r, r) + block(p, r) - block(r, p)
            assert delta == cs + (-1) ** idx * tag
            ds.append(delta)
        values["dark"].append(sum((u.H * u / 8 for u in ds), sp.zeros(4)))
        values["even_mean"].append(cs)
        values["odd_mean"].append(tag)
    return values


def run():
    iv.dps = 34
    native = action.native
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    source = action.front.inputs()
    rho = volume_ratio(source)
    panels = {}
    equivalence = {}
    cert = json.loads(
        (HERE.parent.parent / "full_algebra/certificate.json").read_text()
    )
    entropy = parse_interval(
        cert["signed_transfers"]["pairing_entropy"]
    ) / parse_interval(cert["signed_transfers"]["ordinary_entropy"])
    for q in (-1, 1):
        ops = effects(board, route, q)
        vals = {k: read(ms, source) for k, ms in ops.items()}
        # Every collar change preserves each tagged left marginal and the
        # right marginal summed over left parity. The matching rotation also
        # fixes that right marginal. Ordinary preserves the right marginal.
        for name in ("left", "right", "even_mean", "odd_mean"):
            assert 0 in vals[name]["pairing"].real and 0 in vals[name]["pairing"].imag
            vals[name]["pairing"] = vals[name]["raw_p"] = iv.mpc(0)
        for name in ("left", "right", "even_mean", "odd_mean"):
            vals[name]["ordinary"] = vals[name]["raw_o"]
        vals["right"]["ordinary"] = vals["right"]["raw_o"] = iv.mpc(0)
        left, right = vals["left"], vals["right"]
        lr = product(left, right)
        # Sigma_{2,39,0} = M3 - M2*M1 - M1*M2 + M1^3 on this route.
        # Both subtracted products and the cubic are (-1/2) l*r.
        memory = {k: vals["third"][k] + lr[k] / 2 for k in lr}
        mean1 = product(vals["even_mean"], vals["even_mean"])
        mean2 = product(vals["odd_mean"], vals["odd_mean"])
        variance = {k: vals["dark"][k] - (mean1[k] + mean2[k]) / 4 for k in lr}
        joint = {
            "instantaneous_left_hop_real": (left, "real"),
            "instantaneous_left_hop_imaginary": (left, "imag"),
            "third_moment_real": (vals["third"], "real"),
            "third_moment_imaginary": (vals["third"], "imag"),
            "memory_real": (memory, "real"),
            "memory_imaginary": (memory, "imag"),
            "dark_probability_quadratic": (vals["dark"], "real"),
            "dark_incoherent_quadratic": (variance, "real"),
        }
        panel = {
            k: serialize({i: getattr(x, comp) for i, x in row.items()}, rho)
            for k, (row, comp) in joint.items()
        }
        ref = memory["reference"]
        # This transfer is invariant under an arbitrary rephasing of the
        # remote material ports, including a source-dependent phase choice.
        modulus = {
            k: (
                abs(ref)
                if k == "reference"
                else (ref.real * d.real + ref.imag * d.imag) / abs(ref)
            )
            for k, d in memory.items()
        }
        panel["memory_modulus"] = serialize(modulus, rho)
        weight = read([x.H * x for x in ops["third"]], source)["reference"].real
        elastic = abs(vals["third"]["reference"]) ** 2
        panel["local_sharp_early_arrival"] = dict(
            inclusive_weight=str(weight),
            coherent_weight=str(elastic),
            coherent_fraction=str(elastic / weight),
        )
        equivalence[q] = ops
        panels[str(q)] = panel
    for name in ("third", "left", "right"):
        assert equivalence[-1][name] == [x.H for x in equivalence[1][name]]
    assert (
        equivalence[1]["dark"]
        == equivalence[-1]["dark"]
        == [sp.diag(4, 0, 0, 4), sp.diag(145, 1, 1, 145) / 32]
    )
    ratios = [
        parse_interval(panels["1"][k]["collar_over_ordinary"])
        for k in ("memory_modulus", "dark_incoherent_quadratic")
    ]
    assert ratios[0].b < ratios[1].a and ratios[1].b < entropy.a
    return dict(
        routes=[blank_path(board, route, p, r) for p, r in PAIRS],
        clock="native exchange time",
        entropy_collar_over_ordinary=str(entropy),
        charges=panels,
        certified_modulus_and_dark_ratios_disjoint=bool(ratios[0].b < ratios[1].a),
        certified_modulus_and_entropy_ratios_disjoint=bool(ratios[0].b < entropy.a),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    text = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(text)
    else:
        print(text, end="")
