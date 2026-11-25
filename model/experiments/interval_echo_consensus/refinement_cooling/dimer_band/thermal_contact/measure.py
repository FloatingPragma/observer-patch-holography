"""Independent body thermometers, routed thermal contacts, and complete joint records."""

import argparse
from collections import Counter
from functools import cache
from itertools import product
import json
from math import ceil, log2, pi
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import federation
from ..measure import gates, layouts, outcomes, pointer
from ..surface import branches

HERE = Path(__file__).resolve().parent
BODIES = {"A": (0,), "B": (1, 2), "C": (3, 4, 5)}
CALIBRATION_SHOTS, SHOTS = 2**18, 2**19


def proposals(owner, exposure):
    """Unconditional lower, upper, refusal from the enrolled native instruments."""
    if owner in BODIES["B"]:
        lower, upper, retry = branches("native")[0]
        scale = sum(retry**i for i in range(4))
        return np.array([lower * scale, upper * scale, retry**4])
    _, decoder, leak = gates("native")
    raw = sum(outcomes(decoder, leak, i, exposure) for i in (0, 1)) / 2
    measured = np.r_[raw[-2:], raw[:-2].sum()]
    # Sharp half contacts have dyadic populations. Certify their exact grid
    # before sampling so zero-probability events consume no platform-specific draws.
    denominator = 2 ** (exposure + 1)
    dyadic = np.rint(denominator * measured) / denominator
    np.testing.assert_allclose(measured, dyadic, atol=2e-13, rtol=0)
    return dyadic


def source(owners, exposure, copy):
    local = [proposals(i, exposure if i in BODIES["C"] else 1) for i in owners]
    words = list(product(range(2), repeat=len(owners)))
    law = Counter()
    for w in words:
        p = float(np.prod([q[x] for q, x in zip(local, w, strict=True)]))
        law[w] += p * (0.5 if copy else 1)
        if copy:
            altered = list(w)
            altered[owners.index(2)] = w[owners.index(1)]
            law[tuple(altered)] += p / 2
    # A refusal ends this enrollment, and identifies the first refusing owner.
    weight, refused = 1.0, []
    for low, high, fail in local:
        refused.append(weight * fail)
        weight *= low + high
    probs = np.array([law[x] for x in words] + refused)
    np.testing.assert_allclose(probs.sum(), 1, atol=5e-13, rtol=0)
    return words, probs / probs.sum()


def route_certificate():
    board = federation(0)
    graph = nx.Graph()
    graph.add_edges_from(board.seams)
    ports = []
    for i in range(6):
        private = tuple(np.flatnonzero(board.degree[12 * i : 12 * i + 12] == 6))
        a, b, _, _ = layouts_for_routes()[private]
        ports.append((12 * i + a, 12 * i + b))
    exchange, closure = native_law.local_term_parts()
    # This is a column check of the entire blank-collar action, including its
    # phases, so the material route also transports energy superpositions.
    perm = [native_law.SRC.index((b, a)) for a, b in native_law.SRC]
    expected = np.eye(native_law.DIM)[:, perm]
    for generator, target in [
        (exchange + KAPPA * closure, expected),
        (KAPPA * closure, np.eye(native_law.DIM)[:, :9]),
    ]:
        np.testing.assert_allclose(
            expm(-1j * pi * generator)[:, :9], target, atol=2e-13, rtol=0
        )
    records = []
    for i in range(6):
        for j in range(i + 1, 6):
            path_pair = []
            labels = list(range(board.n))
            for x, y in zip(ports[i], ports[j], strict=True):
                path = nx.shortest_path(graph, x, y)
                edges = list(zip(path, path[1:]))
                path_pair.append(path)
                for a, b in edges + edges[:-1][::-1]:
                    labels[a], labels[b] = labels[b], labels[a]
            wanted = list(range(board.n))
            for a, b in zip(ports[i], ports[j], strict=True):
                wanted[a], wanted[b] = wanted[b], wanted[a]
            if wanted != labels:
                raise ValueError("A relay register retained a change")
            records.append(dict(pair=[i, j], paths=path_pair))
    return records


layouts_for_routes = cache(layouts)


def certify_sign_controls():
    h, decoder, _ = gates("native")
    free = expm(-1j * pi * h)[:2, :2]
    sign_basis = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    np.testing.assert_allclose(
        sign_basis.T @ free @ sign_basis, np.diag([1, -1]), atol=3e-13
    )
    # A private load-squared exposure of pi/(4*kappa) on rail B flips sign.
    private_pulse = np.diag(
        np.exp(-1j * pi / (4 * KAPPA) * KAPPA * np.array([0, 2]) ** 2)
    )
    np.testing.assert_allclose(
        sign_basis.T @ private_pulse @ sign_basis, [[0, 1], [1, 0]], atol=2e-13
    )
    np.testing.assert_allclose(
        abs(decoder[:2, :2] @ sign_basis) ** 2, np.eye(2), atol=2e-13
    )
    np.testing.assert_allclose(
        abs(expm(-1j * pi / 2 * h)[:2, 0]) ** 2, [0.5, 0.5], atol=2e-13
    )


def menu(left, right):
    names = [("idle",)] + [("phase", i) for i in range(left + right)]
    names += [("swap", i, j) for i in range(left) for j in range(left, left + right)]
    return names + [("idle",)] * (2 ** ceil(log2(len(names))) - len(names))


def apply_offer(word, tag, control):
    w = list(word)
    if tag[0] == "swap" and control == "native":
        _, i, j = tag
        w[i], w[j] = w[j], w[i]
    return tuple(w)


def acquire(rng, first, second, exposure=1, copy=False, control="native"):
    owners = BODIES[first] + BODIES[second]
    words, probs = source(owners, exposure, copy)
    tags = menu(len(BODIES[first]), len(BODIES[second]))
    # Select the contact setting only on successfully prepared pairs.
    census = rng.multinomial(SHOTS, probs)
    records = []
    for w, n in zip(words, census[: len(words)], strict=True):
        allocated = rng.multinomial(n, np.full(len(tags), 1 / len(tags)))
        for k, count in enumerate(allocated):
            new = apply_offer(w, tags[k], control)
            records.append([k, word_id(w), word_id(new), int(count)])
    return dict(
        first=first,
        second=second,
        exposure=exposure,
        copy=copy,
        control=control,
        owners=owners,
        tags=tags,
        first_refusals=census[len(words) :].tolist(),
        records=records,
    )


def word_id(bits):
    return sum(x << i for i, x in enumerate(bits))


def coherence(rng, mode):
    # A sharp physical rail is coherent in energy. All input factors are
    # identical; a complete register SWAP preserves this whole state.
    h = gates("native")[0]
    tags = menu(1, 3)
    counts = rng.multinomial(SHOTS, np.full(len(tags), 1 / len(tags)))
    banks = []
    for count, tag in zip(counts, tags, strict=True):
        if tag[0] == "phase" and mode == "native":
            p = float(abs(expm(-1j * pi * h)[1, 0]) ** 2)
            np.testing.assert_allclose(p, 1, atol=2e-13, rtol=0)
            moved = int(rng.binomial(count, 1.0))
            banks.append([int(count) - moved, 1 << tag[1], moved])
        else:
            banks.append([int(count), 0, 0])
    return dict(mode=mode, tags=tags, banks=banks)


def run():
    pointer()
    certify_sign_controls()
    routes = route_certificate()
    rng = np.random.default_rng(893014)
    calibrations = []
    for body, exposure in (("A", 1), ("B", 1), ("C", 1), ("C", 0), ("C", 3)):
        owners = BODIES[body]
        words, p = source(owners, exposure, False)
        n = rng.multinomial(CALIBRATION_SHOTS, p)
        calibrations.append(
            dict(
                body=body,
                exposure=exposure,
                owners=owners,
                words=[
                    [word_id(w), int(k)]
                    for w, k in zip(words, n[: len(words)], strict=True)
                ],
                first_refusals=n[len(words) :].tolist(),
            )
        )
    banks = [acquire(rng, a, b) for a, b in (("A", "B"), ("B", "C"), ("A", "C"))]
    banks.extend(
        [
            acquire(rng, "A", "C", 0),
            acquire(rng, "A", "C", 3),
            acquire(rng, "A", "C", 0, control="exchange_off"),
            acquire(rng, "A", "B", copy=True),
            acquire(rng, "A", "B", copy=True, control="exchange_off"),
        ]
    )
    return dict(
        shots=SHOTS,
        calibration_shots=CALIBRATION_SHOTS,
        routes=routes,
        calibrations=calibrations,
        contacts=banks,
        coherence=[coherence(rng, mode) for mode in ("native", "parked_phase")],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.dumps(run(), separators=(",", ":")) + "\n"
    receipt = HERE / "records.json"
    if args.write:
        receipt.write_text(out)
    elif receipt.read_text() != out:
        raise SystemExit("Thermal contact receipt has changed")
