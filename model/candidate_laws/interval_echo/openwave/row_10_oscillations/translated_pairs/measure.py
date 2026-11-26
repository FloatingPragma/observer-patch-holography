"""Sharp neutral standards on three source bonds of the commissioned release ring."""

import argparse
import json
import math
from itertools import combinations, permutations
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
    LAW,
    blocks,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release.geometry import (
    layout,
)

HERE = Path(__file__).resolve().parent
TIMES = (64, 128, 256, 512)
EXPOSURE = 2**24
SEEDS = {8: 26112610, 7: 26112611}
COMPARISONS = 2048


def generator(ring, kappa, blocked=False):
    """Assemble integer 8H from the native host-subtracted blank-seam terms."""
    n = len(ring)
    words = tuple(permutations(range(n), 2))
    index = {w: i for i, w in enumerate(words)}
    seams = [(i, (i + 1) % n) for i in range(n)]
    if blocked:
        src_owner = ring[0] // 12
        seams = [
            (a, b)
            for a, b in seams
            if (ring[a] // 12 == src_owner) == (ring[b] // 12 == src_owner)
        ]
    block = blocks(kappa)["blank"]
    host = LAW.SRC.index((1, 1))
    h = np.zeros((len(words), len(words)), dtype=int)
    for col, (vac, exc) in enumerate(words):
        loads = np.ones(n, dtype=int)
        loads[vac], loads[exc] = 0, 2
        for a, b in seams:
            word = (loads[a], loads[b])
            j = LAW.SRC.index(word)
            h[col, col] += block[j, j] - block[host, host]
            if word[0] != word[1]:
                swap = {a: b, b: a}
                other = (swap.get(vac, vac), swap.get(exc, exc))
                h[index[other], col] += block[LAW.SRC.index(word[::-1]), j]
    np.testing.assert_array_equal(h, h.T)
    if not blocked:
        shift = [index[((a + 1) % n, (b + 1) % n)] for a, b in words]
        np.testing.assert_array_equal(h, h[np.ix_(shift, shift)])
        # Reflect the ring about source bond zero and interchange both signs.
        mirror = [index[((1 - b) % n, (1 - a) % n)] for a, b in words]
        np.testing.assert_array_equal(h, h[np.ix_(mirror, mirror)])
        for i in range(n):
            j = (-i) % n
            for src_word, target_word in (
                ((i, (i + 1) % n), (j, (j + 1) % n)),
                (((i + 1) % n, i), ((j + 1) % n, j)),
            ):
                assert mirror[index[src_word]] == index[target_word]
    return words, h


def calibrate_echo(ring):
    """Both signs are written in the canonical endpoint convention on each seam."""
    echo = LAW.echo()
    for bond, a in enumerate(ring):
        b = ring[(bond + 1) % len(ring)]
        for x, y in LAW.SRC:
            canonical = (x, y) if a < b else (y, x)
            j = LAW.SRC.index(canonical)
            marked = echo[:, j].reshape(11, 9)
            sign = 1 if a < b else -1
            pointer = [False, *(sign * d == x - y for _, d in LAW.LABELS)]
            np.testing.assert_allclose(
                np.linalg.norm(marked[pointer]) ** 2, 1, atol=1e-13
            )


def read_probabilities(ring, strength, blocked):
    words, h = generator(ring, strength, blocked)
    ev, vectors = eigh(h / 8, driver="evd")
    idx = {word: i for i, word in enumerate(words)}
    source = [idx[(s, s + 1)] for s in range(3)]
    out = {}
    for time in TIMES:
        amplitudes = (vectors * np.exp(-1j * ev * time)) @ vectors[source].T
        probabilities = abs(amplitudes) ** 2
        np.testing.assert_allclose(probabilities.sum(axis=0), 1, rtol=0, atol=5e-13)
        for bond in range(5, 11):
            edge = (bond, (bond + 1) % len(ring))
            for src in range(3):
                first, second = probabilities[[idx[edge], idx[edge[::-1]]], src]
                out[src, bond, time] = [first, second, 1 - first - second]
    return out


def enrollment(fixed_bond):
    for strength in (candidate.KAPPA, 4):
        for time in TIMES:
            for bond in range(fixed_bond - 2, fixed_bond + 1):
                yield strength, time, "calibration", 0, bond
            for source in (0, 1, 2):
                for setting, bond in (
                    ("fixed", fixed_bond),
                    ("translated", fixed_bond + source),
                    ("blocked", fixed_bond),
                    ("host", fixed_bond),
                ):
                    yield strength, time, setting, source, bond


def collect(fixed_bond):
    ring = layout()[0]
    if len(ring) != 14 or len({p // 12 for p in ring[:4]}) != 1:
        raise ValueError("The registered source needs four consecutive local ports.")
    calibrate_echo(ring)
    rng = np.random.default_rng(SEEDS[fixed_bond])
    probabilities = {
        (s, blocked): read_probabilities(ring, s, blocked)
        for s in (candidate.KAPPA, 4)
        for blocked in (False, True)
    }
    rows = []
    for strength, time, setting, src, bond in enrollment(fixed_bond):
        weights = (
            [0, 0, 1]
            if setting == "host"
            else probabilities[strength, setting == "blocked"][src, bond, time]
        )
        rows.append(
            [
                strength,
                time,
                setting,
                src,
                bond,
                rng.multinomial(EXPOSURE, weights).tolist(),
            ]
        )
    return dict(
        ring=list(ring),
        exposure=EXPOSURE,
        seed=SEEDS[fixed_bond],
        fixed_bond=fixed_bond,
        outcomes=["vacancy_excess", "excess_vacancy", "other"],
        banks=rows,
    )


def decode(record):
    fixed_bond = record["fixed_bond"]
    if (
        fixed_bond not in SEEDS
        or record["seed"] != SEEDS[fixed_bond]
        or record["ring"] != list(layout()[0])
    ):
        raise ValueError("unexpected enrollment or surveyed ring")
    if record["exposure"] != EXPOSURE or record["outcomes"] != [
        "vacancy_excess",
        "excess_vacancy",
        "other",
    ]:
        raise ValueError("unexpected exposure or contrast labels")
    data = {}
    for *key, counts in record["banks"]:
        key = tuple(key)
        if (
            key in data
            or len(counts) != 3
            or any(type(v) is not int or v < 0 for v in counts)
            or sum(counts) != EXPOSURE
        ):
            raise ValueError("invalid or duplicate terminal bank")
        data[key] = np.array(counts) / EXPOSURE
    if set(data) != set(enrollment(fixed_bond)):
        raise ValueError("missing or unregistered detector bank")
    width = math.sqrt(math.log(2 * COMPARISONS / 0.01) / (2 * EXPOSURE))
    reports = {}
    for strength in (candidate.KAPPA, 4):
        covered = {"fixed": 0, "translated": 0}
        maximum_residual = 0.0
        for time in TIMES:
            for src in range(3):
                for mode, bond, distance in (
                    ("fixed", fixed_bond, fixed_bond - src),
                    ("translated", fixed_bond + src, fixed_bond),
                ):
                    pred = data[strength, time, "calibration", 0, distance]
                    held = data[strength, time, mode, src, bond]
                    residuals = abs(pred - held)
                    maximum_residual = max(maximum_residual, float(residuals.max()))
                    covered[mode] += int(sum(residuals <= 2 * width))
        resolved = []
        for a, b in combinations(range(3), 2):
            differences = [
                abs(
                    data[strength, time, "fixed", a, fixed_bond][:2]
                    - data[strength, time, "fixed", b, fixed_bond][:2]
                ).max()
                for time in TIMES
            ]
            resolved.append(
                dict(
                    sources=[a, b],
                    largest_signed_bin_difference=float(max(differences)),
                    distinguishable=bool(max(differences) > 2 * width),
                )
            )
        sham_events = sum(
            int(round((1 - data[strength, time, mode, src, fixed_bond][2]) * EXPOSURE))
            for time in TIMES
            for mode in ("blocked", "host")
            for src in range(3)
        )
        reports[str(int(strength))] = dict(
            covered=covered,
            per_setting=36,
            source_comparisons=resolved,
            largest_calibration_residual=maximum_residual,
            sham_events=sham_events,
        )
    return dict(simultaneous_half_width=width, strengths=reports)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path, help="write the simulated count banks")
    parser.add_argument("--receiver", type=int, choices=tuple(SEEDS), default=8)
    args = parser.parse_args()
    if args.acquire:
        record = collect(args.receiver)
        args.acquire.write_text(json.dumps(record, separators=(",", ":")) + "\n")
    else:
        name = "counts.json" if args.receiver == 8 else "second_counts.json"
        record = json.loads((HERE / name).read_text())
    print(json.dumps(decode(record), indent=2))


if __name__ == "__main__":
    main()
