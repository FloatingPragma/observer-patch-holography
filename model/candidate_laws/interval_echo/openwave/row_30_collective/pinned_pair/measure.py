"""Read the released two-load source through its original seam detectors."""

import argparse
import json
from math import log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from ...row_04 import open_pair_check as independent
from ...row_04 import open_pair_release as pair
from .... import candidate


HERE = Path(__file__).resolve().parent
PHASES = (0, pi / 3, -pi / 3, pi / 2, pi)
WAITS = (2, 3, 5, 8, 13, 21)
SHOTS = 2**27
SEED = 30041923


def tip_kick(angle):
    """Exchange onto the native degree-five phase reference and back."""
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    swap = expm(-1j * pi * (exchange + 8 * contrast))
    # In the (156,166) order, X carries 10 and Z carries 11.
    columns = np.eye(law.DIM)[:, [law.SRC.index((1, 0)), law.SRC.index((1, 1))]]
    dwell = (angle % (2 * pi)) / 8
    after_swap = swap @ columns
    after_swap[:law.N_SRC] *= np.array(
        [np.exp(-8j * dwell * n**2) for n, _ in law.SRC]
    )[:, None]
    result = swap @ after_swap
    phases = columns.T @ result
    if not np.allclose(result, columns @ np.diag([1, np.exp(-1j * angle)]),
                       atol=2e-13, rtol=0):
        raise ArithmeticError("Native impulse does not return the two source words")
    return phases


def acquire():
    exchange, contrast, launch, home, detector, _ = independent.board_generator()
    density, preparation = pair.source()
    source = density[np.ix_([0, 2], [0, 2])]
    exhausted = preparation["exhaustion"]
    if abs(source.trace() + exhausted - 1) > 1e-13:
        raise ArithmeticError("Two-word preparation lost outcomes")
    kicks = [tip_kick(angle) for angle in PHASES]
    rng = np.random.default_rng(SEED)
    panels = []
    for kappa in (8, 0):
        flights = independent.flight_amplitudes(exchange, contrast, kappa, WAITS, launch, detector)
        for wait, amplitudes in zip(WAITS, flights, strict=True):
            for arm, prepared in (("coherent", source), ("dephased", np.diag(source.diagonal()))):
                banks = []
                for kick in kicks:
                    evolution = amplitudes @ kick
                    home_probs = np.einsum("ia,ab,ib->i", evolution, prepared, evolution.conj()).real
                    probabilities = np.r_[home_probs, 1 - exhausted - home_probs.sum(), exhausted]
                    if probabilities.min() < -2e-13:
                        raise ArithmeticError("Negative terminal probability")
                    probabilities = np.maximum(0, probabilities).round(13)
                    probabilities /= probabilities.sum()
                    banks.append(rng.multinomial(SHOTS, probabilities).tolist())
                panels.append(dict(kappa=kappa, wait=wait, arm=arm, banks=banks))
    return dict(seed=SEED, shots=SHOTS, phases_in_pi=[0, "1/3", "-1/3", "1/2", 1],
                outcomes=[*map(list, home), "accepted_quiet", "exhausted"], panels=panels)


def bright_weight(angle):
    return (365 - 364 * np.cos(angle)) / 730


def read(records):
    home = [list(edge) for edge in independent.board_generator()[3]]
    outcomes = [*home, "accepted_quiet", "exhausted"]
    if (outcomes != records["outcomes"] or records["shots"] != SHOTS
            or records["phases_in_pi"] != [0, "1/3", "-1/3", "1/2", 1]):
        raise ValueError("Unexpected receiver addresses, phase order or bank size")
    eligible = [i for i, label in enumerate(outcomes)
                if label not in ([156, 157], [156, 166], "exhausted")]
    allowance = sqrt(log(2 * 8000 / .01) / (2 * SHOTS))
    seen, response, worst = set(), {}, 0.
    for row in records["panels"]:
        key = (row["kappa"], row["wait"], row["arm"])
        counts = np.asarray(row["banks"])
        if (key in seen or counts.shape != (5, 32) or counts.dtype.kind not in "iu"
                or np.any(counts < 0) or np.any(counts.sum(axis=1) != SHOTS)):
            raise ValueError("Duplicate panel or incomplete integer bank")
        seen.add(key)
        selected = counts[:, eligible] / SHOTS
        all_scores = np.column_stack([selected, selected.sum(axis=1)])
        for index, angle in enumerate(PHASES[:-1]):
            ratio = bright_weight(angle) / bright_weight(pi) if key[2] == "coherent" else 1.
            residual = abs(all_scores[index] - ratio * all_scores[-1])
            worst = max(worst, float(residual.max() / (allowance * (1 + ratio))))
        response[key] = selected.sum(axis=1)
    expected = {(k, w, a) for k in (8, 0) for w in WAITS for a in ("coherent", "dephased")}
    if seen != expected:
        raise ValueError("Unexpected acquisition schedule")
    if worst > 1:
        raise ArithmeticError("Phase-weight prediction rejected")
    return dict(
        banks=len(seen) * len(PHASES), simultaneous_half_width=allowance,
        largest_error_in_allowances=worst,
        native_pi=[float(response[8, w, "coherent"][-1]) for w in WAITS],
        exchange_only_pi=[float(response[0, w, "coherent"][-1]) for w in WAITS],
        native_zero=[float(response[8, w, "coherent"][0]) for w in WAITS],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="create new acquisition at this path")
    args = parser.parse_args()
    records = acquire()
    result = read(records)
    if args.output:
        with args.output.open("x", encoding="utf-8") as f:
            json.dump(records, f, indent=1, allow_nan=False)
            f.write("\n")
    elif records != json.loads((HERE / "counts.json").read_text()):
        raise ArithmeticError("Acquisition differs from recorded counts")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
