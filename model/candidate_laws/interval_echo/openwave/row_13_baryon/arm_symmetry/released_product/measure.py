"""Released three-local-read instrument with full-host numerical uncertainty."""

import argparse
import json
from itertools import combinations, product
from math import factorial
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import beta

from .. import measure as source
from ...open_release.flight_check import neighbourhood
from model.candidate_laws.interval_echo.candidate import KAPPA
from model.executable_model import carrier
from model.experiments.interval_echo_gauge.trimer_transport import (
    local_bases,
    measure as colour,
)

HERE = Path(__file__).resolve().parent
WORDS = sorted(product(range(3), repeat=3), key=lambda w: (len(set(w)) == 3, w))
STARTS = 2**32
PANELS = 19
COMPARISONS = PANELS * 21
READ_ERROR = 1e-10


def masks(rails):
    edges = set(carrier.cell().seams)
    counts = [
        sum(
            tuple(sorted((int(rails[i, w[i]]), int(rails[j, w[j]])))) in edges
            for i, j in combinations(range(3), 2)
        )
        for w in WORDS
    ]
    support = {counts[i] for i, w in enumerate(WORDS) if len(set(w)) == 3}
    if len(support) != 1:
        raise ValueError(
            "This aperture requires a noncommuting connectivity measurement"
        )
    return counts, support.pop() >= 2


def compress(level, kappa, sec, initial, snapshots):
    space, lost = neighbourhood(level)
    h = space.hamiltonian(kappa) + sparse.diags(lost / 2)
    home = np.flatnonzero(np.all(space.tri < 12, axis=1))
    if [tuple(w) for w in space.tri[home]] != sec.states:
        raise ValueError("Home basis ordering differs")
    state = np.zeros(space.dim, complex)
    state[home] = initial
    # Center the spectrum in [-6*kappa,18]. A quarter-step uses eight
    # centered Taylor substeps; total truncation + sparse rounding < 1e-9.
    midpoint = (18 - 6 * kappa) / 2
    step = (-1j / 32) * (h - sparse.eye(space.dim, format="csr") * midpoint)
    norm = max(np.abs(step).sum(axis=1).A1)
    if norm >= 2.1:
        raise ValueError("Taylor and rounding certificate no longer applies")
    tail = 32 * np.exp(norm) * norm**25 / factorial(25)
    if tail >= 1e-11:
        raise ValueError("Insufficient polynomial degree")
    out = {}
    for tick in range(1, 33):
        term = state
        new = state.copy()
        for j in range(1, 25):
            term = (step @ term) / j
            new += term
        state = new
        if tick / 32 in snapshots:
            out[tick / 32] = state[home].copy()
    return out


def cumulative_bank(signal, err, source_range, rng):
    masses = np.asarray([signal[w] for w in WORDS])
    prefix = np.sqrt(np.cumsum(masses))
    lo, hi = source_range
    # Each prefix is a projective effect after a unitary, or a positive
    # contraction with a forgotten dephasing read. Its square root contracts.
    lower = list(np.maximum(0, lo * np.maximum(0, prefix - err) ** 2 - 1e-14))
    upper = list(np.minimum(1, hi * (prefix + err) ** 2 + 1e-14))
    bounds = np.array([*zip(lower, upper), (lo, hi)])
    cuts = np.unique(np.r_[0.0, bounds.ravel(), 1.0])
    bins = rng.multinomial(STARTS, np.diff(cuts))
    table = {}
    for left, right, n in zip(cuts[:-1], cuts[1:], bins, strict=True):
        if n:
            u = (left + right) / 2
            a = int(np.count_nonzero(bounds[:, 1] < u))
            b = int(np.count_nonzero(bounds[:, 0] < u))
            table[a, b] = table.get((a, b), 0) + int(n)
    return [[a, b, n] for (a, b), n in table.items()]


def event(bank, which):
    good = set(which)
    low = high = 0
    for a, b, n in bank:
        covered = good.intersection(range(a, b + 1))
        high += n if covered else 0
        low += n if len(covered) == b - a + 1 else 0
    return low, high


def interval(klo, khi, nlo, nhi):
    tail = 0.01 / (2 * COMPARISONS)
    return (
        0.0 if klo == 0 else float(beta.ppf(tail, klo, nhi - klo + 1)),
        1.0 if khi >= nlo else float(beta.isf(tail, khi + 1, nlo - khi)),
    )


def decode(banks, connected):
    eq = []
    grouped = []
    admissions = []
    for b in banks:
        if sum(x[2] for x in b) != STARTS:
            raise ValueError("Missing original launches")
        n = event(b, range(27))
        if n[0] == 0:
            raise ValueError("No certain aperture admissions")
        admissions.append(n)
        for i, j in combinations(range(3), 2):
            k = event(b, [ix for ix, w in enumerate(WORDS) if w[i] == w[j]])
            eq.append(interval(*k, *n))
        repeat = event(b, [ix for ix, w in enumerate(WORDS) if len(set(w)) < 3])
        triple = event(b, [ix for ix, w in enumerate(WORDS) if len(set(w)) == 1])
        grouped.append(
            np.array(interval(*repeat, *n)) + 2 * np.array(interval(*triple, *n))
        )
    q = interval(*admissions[0], STARTS, STARTS)
    separate = np.sum(eq, axis=0)
    together = np.sum(grouped, axis=0)
    sums = max(separate[0], together[0]), min(separate[1], together[1])
    if sums[0] > sums[1]:
        raise ValueError("Simultaneous bounds for the same statistic are inconsistent")
    f = [max(0.0, 1 - sums[1] / 3 - READ_ERROR), min(1.0, 1 - sums[0] / 6 + READ_ERROR)]
    joint = f if connected else [0.0, 0.0]
    return dict(
        admissions=admissions,
        separate_pair_sum=separate.tolist(),
        grouped_pair_sum=together.tolist(),
        fidelity=f,
        connected_singlet=joint,
        joint_per_million=[1e6 * a * b for a, b in zip(q, joint, strict=True)],
    )


def observe(sec, rails, vector, err, meta, case, frames, rng):
    raw = colour.tensor(sec, rails, vector)
    counts, connected = masks(rails)
    repeat = colour.native_word("AB")
    changed = colour.native_word("BA")
    bounds = meta["original_start_source_lower"], meta["original_start_source_upper"]
    banks = []
    for frame in frames:
        a, b = (
            (frame @ repeat, frame @ changed) if case == "mismatch" else (frame, frame)
        )
        probs = colour.detector(raw, a, b, b, case == "dephased")
        banks.append(cumulative_bank(probs, err, bounds, rng))
    return dict(
        connected_words=sum(x >= 2 for x in counts),
        determinant_connected=connected,
        bins=banks,
        **decode(banks, connected),
    )


def run():
    record = json.loads((source.HERE / "results.json").read_text())
    sec, _, arms, _, gates = source.apparatus()
    if KAPPA != 8:
        raise ValueError(
            "The archived finite-source and Duhamel certificates require contrast eight"
        )
    h = sec.hamiltonian(KAPPA, False) + np.diag(gates / 2)
    eig, vec = np.linalg.eigh(h)
    ix = np.argmin(abs(eig - source.TARGET))
    start = vec[:, ix]
    meta = record["source"]
    residual = np.linalg.norm(h @ start - eig[ix] * start)
    if (residual + 1e-11) / meta["gap_lower"] > 1e-8:
        raise ValueError("Fresh start vector is outside source tolerance")
    frames = [local_bases.compile_basis(f) for f in local_bases.bases()]
    rng = np.random.default_rng(13110271)
    panels = []
    for level in (0, 1, 2):
        cells = 20 * 4**level
        for channel in ("native", "exchange") if not level else ("native",):
            times = (0.25, 1.0) if not level and channel == "native" else (1.0,)
            held = compress(
                level, KAPPA if channel == "native" else 0, sec, start, times
            )
            for t in times:
                ref = record["open_boards"][str(cells)][channel][str(t)]
                for omit in range(4) if channel == "native" else (3,):
                    rails = np.delete(arms, omit, axis=1)
                    q = np.linalg.norm(colour.tensor(sec, rails, held[t])) ** 2
                    if abs(q - ref["colours"][omit]["triplet"]) > 1e-10:
                        raise ValueError("Independent endpoint replay disagrees")
                    cases = (
                        ("native", "dephased", "mismatch")
                        if (level == 0 and channel == "native" and t == 1 and omit == 3)
                        else ("native",)
                    )
                    for case in cases:
                        obs = observe(
                            sec,
                            rails,
                            held[t],
                            ref["home_vector_error_upper"] + 5e-8,
                            meta,
                            case,
                            frames,
                            rng,
                        )
                        panels.append(
                            dict(
                                cells=cells,
                                flight=channel,
                                wait=t,
                                omit=omit,
                                read=case,
                                **obs,
                            )
                        )
    if len(panels) != PANELS:
        raise ValueError("Confidence family needs updating")
    return dict(starts_per_frame=STARTS, words=WORDS, panels=panels)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    data = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "records.json"
    if p.parse_args().write:
        path.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    else:
        saved = json.loads(path.read_bytes())
        if saved["words"] != data["words"] or saved["starts_per_frame"] != STARTS:
            raise ValueError("Changed enrollment")
        for result, expected in zip(data["panels"], saved["panels"], strict=True):
            for key in result:
                if key in {
                    "fidelity",
                    "connected_singlet",
                    "joint_per_million",
                    "separate_pair_sum",
                    "grouped_pair_sum",
                }:
                    np.testing.assert_allclose(
                        result[key], expected[key], atol=1e-9, rtol=0
                    )
                elif result[key] != expected[key]:
                    raise ValueError(f"Changed release record {key}")
    for row in data["panels"]:
        print(
            row["cells"],
            row["flight"],
            row["wait"],
            row["omit"],
            row["read"],
            row["fidelity"],
            row["joint_per_million"],
        )


if __name__ == "__main__":
    main()
