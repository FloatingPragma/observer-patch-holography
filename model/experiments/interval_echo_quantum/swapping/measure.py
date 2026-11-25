"""Acquire central records, independent late settings and both terminal signs."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_quantum.experiment import sample
from . import circuit, compiler

HERE = Path(__file__).resolve().parent
SHOTS = 2**22
SEED = 738502
SIGNS = np.array([1, -1, 1])
CHSH = np.array([[1, 1], [1, -1]])
CHECKER = np.array([[1, -1], [-1, 1]])


def acquire(prob, coins, rng):
    table = np.zeros(prob.shape, dtype=int)
    central = prob.sum(axis=(4, 5))[:, :, 0, 0]
    syndromes = sample(SHOTS, central, rng)
    for b, d in product(range(2), repeat=2):
        for x, n in enumerate(sample(syndromes[b, d], coins[0], rng)):
            for z, m in enumerate(sample(n, coins[1], rng)):
                table[b, d, x, z] = sample(m, prob[b, d, x, z] / central[b, d], rng)
    return table


def radius(n):
    return np.sqrt(2 * np.log(2 * 2048 / 0.01) / n)


def score(table, mode, menu):
    table = np.asarray(table)
    if mode not in circuit.MODES or menu not in circuit.MENU:
        raise ValueError("Unknown protocol")
    if (
        not np.issubdtype(table.dtype, np.integer)
        or table.shape != (2, 2, 2, 2, 3, 3)
        or table.min() < 0
        or table.sum() != SHOTS
    ):
        raise ValueError("Expected one complete nonnegative count bank")
    n = table.sum(axis=(4, 5))
    if np.any(n == 0):
        raise ValueError("A required syndrome and setting was never acquired")
    signed = np.einsum("bdxzac,a,c->bdxz", table, SIGNS, SIGNS)
    pooled_n = n.sum(axis=(0, 1))
    zsign = np.array([[1, -1], [1, -1]])
    xsign = np.array([[1, 1], [-1, -1]])
    if mode == "separable_bob":
        zsign = CHECKER

    def expectation(frame):
        return (signed * frame[:, :, None, None]).sum(axis=(0, 1)) / pooled_n

    ez, ex = expectation(zsign), expectation(xsign)
    outcome = dict(
        launches=int(n.sum()),
        syndromes=n.sum(axis=(2, 3)).tolist(),
        outer_quiets=[
            int(table[:, :, :, :, 2, :].sum()),
            int(table[:, :, :, :, :, 2].sum()),
        ],
        settings=pooled_n.tolist(),
    )
    if menu == "chsh":
        frame = np.stack((zsign, xsign), axis=2)[:, :, :, None]
        witness = float((CHSH * np.stack((ez[0], ex[1]))).sum())
        outcome.update(
            S=witness,
            lower_99=witness - radius(pooled_n).sum(),
            erased_S=float((CHSH * signed.sum(axis=(0, 1)) / pooled_n).sum()),
            conditional_S=(CHSH * frame * signed / n).sum(axis=(2, 3)).tolist(),
            conditional_lower_99=(CHSH * frame * signed / n - radius(n))
            .sum(axis=(2, 3))
            .tolist(),
        )
    else:
        i, j = float(ez.mean()), float((CHECKER * ex).mean())
        vals = np.abs([i, j])
        error = radius(pooled_n).mean()
        outcome.update(
            I=i,
            J=j,
            network_score=float(np.sqrt(vals).sum()),
            lower_99=float(np.sqrt(np.maximum(vals - error, 0)).sum()),
            upper_99=float(np.sqrt(vals + error).sum()),
        )
    # Each cell here is a +/-1 observable. The radius on a binary outcome
    # probability is half the signed radius, also for quiet.
    local = []
    for axis in (4, 5):
        marg = table.sum(axis=(0, 1, axis)) / pooled_n[:, :, None]
        remote_axis = 0 if axis == 4 else 1
        gap = abs(np.diff(marg, axis=remote_axis)).max()
        allowance = (radius(pooled_n) / 2).sum(axis=remote_axis).max()
        local.append(
            {
                "max_remote_rate_difference": float(gap),
                "joint_allowance": float(allowance),
            }
        )
    outcome["no_signalling_C_A"] = local
    return outcome


def run(compiled=False):
    circuit.check_readers()
    line, routes = circuit.apparatus()
    coins = [circuit.settings_probability(line[k]) for k in (0, -1)]
    rng = np.random.default_rng(SEED)
    phase = compiler.native_phase if compiled else circuit.phase
    banks = {}
    for mode, menu in product(circuit.MODES, circuit.MENU):
        prob = circuit.distributions(line, routes, mode, menu, phase)
        table = acquire(prob, coins, rng)
        banks[f"{mode}/{menu}"] = {
            "counts": table.tolist(),
            "score": score(table, mode, menu),
        }
    return {
        "owners": line,
        "swap_routes": [[list(e) for e in r] for r in routes],
        "shots_per_bank": SHOTS,
        "seed": SEED,
        "axes": [
            "central_left_bit",
            "central_right_bit",
            "setting_A",
            "setting_C",
            "pointer_A(+,-,quiet)",
            "pointer_C(+,-,quiet)",
        ],
        "banks": banks,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument(
        "--compiled", action="store_true", help="Use the ten-contact phase compiler"
    )
    args = parser.parse_args()
    receipt = run(args.compiled)
    if args.compiled:
        print("Gate residual bounds", compiler.phase_bound())
    path = HERE / "records.json"
    if args.write:
        path.write_text(json.dumps(receipt, indent=1) + "\n")
    elif json.loads(path.read_text()) != receipt:
        raise ValueError("Record replay changed")
    for name, bank in receipt["banks"].items():
        s = bank["score"]
        print(
            name, round(s.get("S", s.get("network_score")), 6), round(s["lower_99"], 6)
        )


if __name__ == "__main__":
    main()
