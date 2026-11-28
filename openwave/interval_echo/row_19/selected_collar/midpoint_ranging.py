"""Test the necessary midpoint range ratios at every native port of the six-station survey."""

import argparse
import json
from hashlib import sha256
from itertools import combinations
from math import log, sqrt

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import finite_dose_curvature as curve

lab = curve.lab
N = 2**28


def panels(h, start, targets):
    """Launch from one port, retaining all four head amplitudes and all requested targets."""
    src = np.zeros((h.shape[0], 4))
    src[4 * start : 4 * start + 4] = np.eye(4)
    u = -0.25j * h
    t = expm_multiply(u, src, start=0, stop=11, num=12, traceA=u.diagonal().sum())
    return t[1:].reshape(11, -1, 4, 4)[:, targets]


def first(forward, diagonal):
    hits = []
    for i, U in enumerate(forward[:10]):
        hits.append(U - sum(diagonal[i - j - 1] @ hits[j] for j in range(i)))
    return hits


def echo(out, home, away):
    # The integer generator is real symmetric, so reversing a block transposes it.
    go = first(out, away)
    back = first(np.swapaxes(out, -1, -2), home)
    effects = []
    for tick in range(2, 12):
        total = np.zeros((4, 4), complex)
        for a in range(1, tick):
            amp = back[tick - a - 1] @ go[a - 1]
            total += amp.conj().T @ amp
        effects.append(total)
    return np.array(effects)


def lower_range(freq, err):
    # P(completed by tick k) <= P(first-leg hit by k-1)
    # <= 4 sum_{s<=k-1} P(free position at s equals the addressed relay).
    cumulative = np.minimum(1, 4 * np.cumsum(freq + err, axis=0))
    return 3 - cumulative.sum(axis=0) / 4


def run():
    source, states, *_ = curve.trials()
    arms = {k: states[k] for k in ("reference", "held")}
    native = lab.six.survey.native
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    ports = (*route, 5, 6)
    candidates = np.arange(board.n)
    rng = np.random.default_rng(28112084)
    free = {k: np.zeros((6, 10, board.n)) for k in arms}
    archives = {}
    # Save only four columns per sender per flight sector.
    for odd in (0, 1):
        for charge in (-1, 1):
            h = native.generator(board, route, charge, bool(odd))
            assert (h - h.T).nnz == 0
            samples = [panels(h, p, candidates) for p in ports]
            for i, arr in enumerate(samples):
                for k, x in arms.items():
                    free[k][i] += (
                        np.einsum(
                            "tnai,ij,tnaj->tn",
                            arr[:10],
                            x[odd],
                            arr[:10].conj(),
                            optimize=True,
                        ).real
                        / 2
                    )
            archives[odd, charge] = h, samples
    # One 99% family for both free-position banks and every possible echo bank.
    count_tests = 2 * 6 * board.n * 11
    ferr = sqrt(log(2 * count_tests / 0.01) / (2 * N)) + 1e-11
    lower, raw = {}, {}
    needed = set()
    for k, bank in free.items():
        assert np.max(abs(bank.sum(2) - 1)) < 1e-11
        hist = np.empty(bank.shape, dtype=int)
        for i, row in enumerate(bank):
            for t, p in enumerate(row):
                order = np.argsort(p, kind="stable")
                draws = rng.multinomial(N, np.maximum(p[order], 0))
                hist[i, t, order] = draws
        raw[k] = hist.tolist()
        lower[k] = np.stack([lower_range(row / N, ferr) for row in hist])
        for i, p in enumerate(ports):
            needed.update(
                (i, int(m)) for m in candidates if m != p and lower[k][i, m] < 2.5
            )
    # Include all original station pairs explicitly for the denominators.
    near = sorted(
        needed | {(i, p) for i, home in enumerate(ports) for p in ports if p != home}
    )
    probs = {k: np.zeros((len(near), 10)) for k in arms}
    for (odd, charge), (h, samples) in archives.items():
        diagonals = {p: samples[i][:, p] for i, p in enumerate(ports)}
        for m in sorted(set(m for _, m in near) - set(ports)):
            diagonals[m] = panels(h, m, [m])[:, 0]
        for j, (i, m) in enumerate(near):
            eff = echo(samples[i][:, m], diagonals[ports[i]], diagonals[m])
            for k, x in arms.items():
                probs[k][j] += np.einsum("tij,ji->t", eff, x[odd]).real / 2
    reports = {}
    for k, arr in probs.items():
        records = {}
        upper = np.full((6, board.n), 3.0)
        for (i, m), p in zip(near, arr, strict=True):
            law = np.r_[p, 1 - p.sum()]
            assert law.min() > -1e-13
            hits = rng.multinomial(N, np.maximum(law, 0))
            value = hits @ curve.reg.SCORE / N
            lo, hi = value - 2.5 * ferr, value + 2.5 * ferr
            lower[k][i, m] = max(lo, lower[k][i, m])
            upper[i, m] = min(hi, 3)
            records[f"{ports[i]}-{m}"] = hits.tolist()
        latency_upper = [upper[i, candidates != p].min() for i, p in enumerate(ports)]
        latency_safe = [lower[k][i, candidates != p].min() for i, p in enumerate(ports)]
        edges = []
        for i, j in combinations(range(6), 2):
            a, b = ports[i], ports[j]
            available = (candidates != a) & (candidates != b)
            ratios = np.stack(
                [
                    2 * lower[k][i, available] / upper[i, b],
                    2 * lower[k][j, available] / upper[j, a],
                ]
            )
            single = ratios.max(0)
            best = int(np.argmin(single))
            # Every convex combination lower bounds the larger of the two legs.
            # A grid provides rigorous dual witnesses without assuming an optimal weight.
            weights = np.linspace(0, 1, 129)
            bounds = (
                weights[:, None] * ratios[0] + (1 - weights[:, None]) * ratios[1]
            ).min(1)
            best_weight = int(bounds.argmax())
            apparent = 2 * lower[k][[i, j]][:, available]
            allowed = np.array(
                [upper[i, b] + latency_upper[i], upper[j, a] + latency_upper[j]]
            )
            residual = apparent - allowed[:, None]
            affine = (
                weights[:, None] * residual[0] + (1 - weights[:, None]) * residual[1]
            ).min(1)
            ai = int(affine.argmax())
            witnesses = []
            for m in candidates[available]:
                need_lo = [
                    2 * lower[k][i, m] - upper[i, b],
                    2 * lower[k][j, m] - upper[j, a],
                ]
                need_hi = [
                    2 * upper[i, m] - lower[k][i, b],
                    2 * upper[j, m] - lower[k][j, a],
                ]
                if min(need_lo) > 0 and all(
                    x < latency_safe[z] for x, z in zip(need_hi, (i, j), strict=True)
                ):
                    witnesses.append(int(m))
            edges.append(
                {
                    "ports": [a, b],
                    "affine_excess_lower": float(affine[ai] - 1e-12),
                    "affine_dual_weight": float(weights[ai]),
                    "individually_admissible_affine_sites": witnesses,
                    "best_single_port_by_bound": int(candidates[available][best]),
                    "single_ratio_lower": float(single[best] - 1e-12),
                    "mixture_ratio_lower": float(bounds[best_weight] - 1e-12),
                    "mixture_dual_weight": float(weights[best_weight]),
                    "minimum_required_absolute_log_error": float(
                        log(bounds[best_weight]) / 2 - 1e-12
                    ),
                }
            )
        reports[k] = {
            "edges": edges,
            "echo_counts": records,
            "sender_latency_upper": list(latency_upper),
            "sender_safe_latency": list(latency_safe),
        }
    return {
        "ports": list(ports),
        "native_ports": board.n,
        "trials_per_setting": N,
        "seed": 28112084,
        "position_counts": raw,
        "targeted_extra_echoes": len(needed),
        "directed_echoes": len(near),
        "results": reports,
        "inputs": {
            p: sha256((lab.HERE / p).read_bytes()).hexdigest()
            for p in (
                "MIDPOINT_RANGING.md",
                "MIDPOINT_LATENCY.md",
                "finite_pairing_dose.json",
                "ordinary_stationarity.json",
            )
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = run()
    path = lab.HERE / "midpoint_ranging.json"
    if cli.parse_args().write:
        path.write_text(
            json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)
    for k, r in result["results"].items():
        print(
            k,
            result["directed_echoes"],
            min(p["mixture_ratio_lower"] for p in r["edges"]),
        )
