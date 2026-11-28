"""Compare physical two-port modular references with independently prepared energy pilots."""

import argparse
import json

import numpy as np
from scipy.linalg import eigh

from .measure import HERE, check, layout, sector


EQUAL = np.diag([0.5, 0, 0, 0.5])


def parity(sign):
    v = np.array([0.0, 1, sign, 0]) / np.sqrt(2)
    return np.outer(v, v)


EVEN, ODD = parity(1), parity(-1)
H = 8 * EVEN + 9 * ODD


def local_rho(words, rho, edge):
    """Partial trace with sharp record ordering (00,01,10,11)."""
    a, b = edge
    groups = {}
    for i, word in enumerate(words):
        rest = tuple(x for x in word if x not in edge)
        groups.setdefault(rest, []).append((2 * (a in word) + (b in word), i))
    reduced = np.zeros((4, 4))
    for group in groups.values():
        ix, full = map(list, zip(*group, strict=True))
        reduced[np.ix_(ix, ix)] += rho[np.ix_(full, full)]
    return reduced


def calibrate(ref, delta):
    vals, v = eigh(ref)
    if vals.min() <= 1e-14:
        return dict(singular=True, minimum_weight=float(vals.min()))
    K = (v * (-np.log(vals))) @ v.T

    def trace(a, b):
        return float(np.sum(a.T * b))

    pilots = [trace(K, q) for q in (EQUAL, EVEN, ODD)]
    beta8 = (pilots[1] - pilots[0]) / 8
    beta1 = pilots[2] - pilots[1]
    delta_e = trace(H, delta)
    delta_k = trace(K, delta)
    return dict(
        singular=False,
        minimum_weight=float(vals.min()),
        beta_equal_to_even=beta8,
        beta_even_to_odd=beta1,
        held_energy=delta_e,
        held_modular=delta_k,
        residual_nats=delta_k - beta8 * delta_e,
    )


def partial_traces(sec, indices, edges):
    w, v = eigh(sec.matrix / 2)
    blocks = abs(w[:, None] - w[None, :]) < 1e-10
    weights = np.bincount(indices, minlength=len(sec.words)) / len(indices)
    rho = v @ (((v.T * weights) @ v) * blocks) @ v.T
    return np.array([local_rho(sec.words, rho, e) for e in edges])


def run():
    ports, edges, _, pairs = layout(1)
    # Use the original head order on each sorted pair. The full regional
    # reference retains the actual imbalance of the two single-head marginals.
    singles = sector(len(ports), edges, 1, 8)
    doubles = sector(len(ports), edges, 2, 8)
    components = [np.tile(np.diag([1.0, 0, 0, 0]), (len(edges), 1, 1))]
    for slot in (1, 0):
        components.append(
            partial_traces(
                singles, [singles.words.index((pair[slot],)) for pair in pairs], edges
            )
        )
    components.append(
        partial_traces(doubles, [doubles.words.index(pair) for pair in pairs], edges)
    )
    j = json.loads(HERE.parent.joinpath("formed_heads.json").read_text())
    out = []
    for s in j["source"]:
        if s["orientation"] != "matching":
            continue
        original = np.array(s["formed_sharp_pairs"]).reshape(2, 2)
        independent = np.outer(original.sum(axis=1), original.sum(axis=0))
        ref, joint = [
            np.tensordot(m.ravel(), components, axes=1) for m in (independent, original)
        ]
        result = []
        for edge, a, b in zip(edges, ref, joint, strict=True):
            np.testing.assert_allclose([a.trace(), b.trace()], 1, atol=3e-14)
            result.append(
                dict(
                    edge=edge,
                    formed=calibrate(b, a - b),
                    independent=calibrate(a, a - b),
                )
            )

        def masked(q):
            return np.diag((q.ravel() + q.ravel()[::-1]) / 2)

        held = calibrate(masked(original), masked(independent) - masked(original))
        np.testing.assert_allclose(held["beta_even_to_odd"], 0, atol=5e-15)
        np.testing.assert_allclose(
            held["held_modular"], s["modular_change_fixed_formed"]
        )
        # An independently specified Gibbs state is the positive control for
        # a common thermometer. A number chemical potential cancels for all pilots.
        energies, vec = eigh(H)
        gibbs = vec @ np.diag(np.exp(-energies / 8)) @ vec.T
        gibbs /= gibbs.trace()
        control = calibrate(gibbs, masked(independent) - masked(original))
        np.testing.assert_allclose(
            [control["beta_equal_to_even"], control["beta_even_to_odd"]], [1 / 8] * 2
        )
        np.testing.assert_allclose(control["residual_nats"], 0, atol=1e-14)
        out.append(
            dict(
                source_sign=s["source_sign"],
                held_masked=held,
                thermal_control=control,
                free_contacts=result,
            )
        )
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = run()
    path = HERE / "modular.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2) + "\n")
    else:
        check(json.loads(path.read_text()), json.loads(json.dumps(data)))
