"""Interrogate the registered antipodal port controls on a single and a pair."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import currents
from model.executable_model import cell
from ..__main__ import ordered_word
from .certificate import pair_bright

HERE = Path(__file__).parent
SHOTS = 131072
STEPS = 2048


def acquisition(pair, bright_rate, cycles):
    car = cell()
    p = currents.SOURCE_PORT
    anti = currents._geometry(car)[0][p]
    sector, sources, _, _, bands, _, _, _, operators, _ = currents._service(
        car, pair, 1.0
    )
    h = np.array(operators)
    basis = np.column_stack(bands)
    axis = bands[0].T @ sources[:, p]
    axis /= np.linalg.norm(axis)
    bright = np.outer(axis, axis)
    d = h[p] + h[anti]
    loads = 1 + int(pair)
    rate = bright_rate if pair else 0.5
    plane = (loads / 2 - rate) / 2
    target = np.zeros((6, 6))
    target[:3, :3] = plane * np.eye(3) + (rate - plane) * bright
    np.testing.assert_allclose(d, target, atol=1e-12, rtol=0)
    np.testing.assert_allclose(h[p] @ h[anti], h[anti] @ h[p], atol=1e-12)
    proposals = {
        "sham": (0, 0),
        "phase_anchor": (0, np.pi / loads),
        "native": (np.pi, 0),
        "determinant_lock": (np.pi, np.pi * (bright_rate - 0.5)),
        "fixed_line": (np.pi, -np.pi * bright_rate / 2),
        "plane_flip": (4 * np.pi * cycles, 0),
    }
    rng = np.random.default_rng(66193 + loads)
    out = {}
    for arm in (*proposals, "zero_occupation"):
        time, calibration = proposals["plane_flip" if arm == "zero_occupation" else arm]
        gen = 1j * h.imag if arm == "zero_occupation" else h
        word = expm(1j * time * gen[anti]) @ expm(1j * time * gen[p])
        if calibration:
            word = ordered_word(gen, np.ones(12), calibration, STEPS) @ word
        ideal = expm(1j * (time * (gen[p] + gen[anti]) + calibration * gen.sum(axis=0)))
        err = float(np.linalg.norm(word - ideal, 2))
        if err > 1e-5:
            raise ValueError("finite pulse failed the ideal group word")
        whole = np.eye(len(sources)) + basis @ (word - np.eye(6)) @ basis.T
        arms = {}
        for origin in (p, min(i for i in range(12) if i not in (p, anti))):
            prob = abs(whole @ sources[:, origin]) ** 2
            if pair:
                prob = np.r_[prob[sector.seam_states], prob[~sector.adjacent].sum()]
            np.testing.assert_allclose(prob.sum(), 1, atol=1e-11, rtol=0)
            counts = rng.multinomial(SHOTS, prob / prob.sum())
            event = (
                np.any(np.isin(car.seams, [p, anti]), axis=1).nonzero()[0]
                if pair
                else np.array([p, anti])
            )
            arms[str(origin)] = dict(
                counts=counts.tolist(),
                probabilities=prob.round(11).tolist(),
                return_count=int(counts[event].sum()),
            )
        out[arm] = dict(operator_error=round(err, 10), launches=arms)
    return out


def cosine_interval(arms):
    # Six return means (two sectors), simultaneous two-sided 99% Hoeffding.
    width = np.sqrt(np.log(12 / 0.01) / (2 * SHOTS))
    # Account for the symmetric product error in the ideal anchor probability:
    # L<=12*(pi/N)*1.4, two-sided state distance <= unitary operator error.
    n = len(arms["sham"]["launches"]["1"]["counts"])
    loads = 1 if n == 12 else 2
    length = 12 * (np.pi / loads) * (1.4 if loads == 1 else 0.81)
    anchor_bias = np.exp(length / STEPS) * length**3 / (3 * STEPS**2)
    bounds = []
    for arm in ("native", "sham", "phase_anchor"):
        f = arms[arm]["launches"]["1"]["return_count"] / SHOTS
        radius = width + (anchor_bias if arm == "phase_anchor" else 0)
        bounds.append((max(0.0, f - radius), min(1.0, f + radius)))
    if bounds[1][0] <= bounds[2][1]:
        raise ValueError("phase reference has no resolved contrast")
    estimates = [
        1 + 2 * (reading - sham) / (sham - anchor)
        for reading, sham, anchor in product(*bounds)
    ]
    return [float(min(estimates)), float(max(estimates))]


def run():
    lower, upper = pair_bright()
    # The 12 eigenvalues sum to 2I on this triplet, so trace(D)=1.
    plane = (1 - (lower + upper) / 2) / 2
    cycles = next(
        n for n in range(1, 513) if abs(float(2 * n * plane % 1) - 0.5) < 0.003
    )
    for n in range(1, cycles + 1):
        turns = [2 * n * (1 - x) / 2 for x in (lower, upper)]
        distances = [abs(float(x % 1) - 0.5) for x in turns]
        if (n < cycles and min(distances) <= 0.003) or (
            n == cycles and max(distances) >= 0.003
        ):
            raise ValueError("spectral interval does not certify cycle selection")
    arms = [
        acquisition(pair, float((lower + upper) / 2), cycles) for pair in (False, True)
    ]
    intervals = [cosine_interval(a) for a in arms]
    if intervals[1][1] >= intervals[0][0]:
        raise ValueError("the two phase readings do not distinguish the proposed lock")
    return dict(
        shots=SHOTS,
        symmetric_slices=STEPS,
        flip_cycles=cycles,
        pair_bright_enclosure=[0.719927, 0.719929],
        single=arms[0],
        pair=arms[1],
        registered_cosine_intervals=np.round(intervals, 7).tolist(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    data = run()
    if args.record:
        with args.record.open("x") as f:
            f.write(json.dumps(data, indent=1) + "\n")
    elif data != json.loads((HERE / "records.json").read_text()):
        raise ValueError("observer acquisition differs from the deposited banks")
    print(
        "Intervals for the colour determinant and the pair bright phase:",
        data["registered_cosine_intervals"],
    )
