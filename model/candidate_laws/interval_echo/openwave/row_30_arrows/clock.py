"""Joint fuel-clock and transported matter readings inside the same finite patch."""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import collar
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from . import phase
from .experiment import HERE, canonical, sector


SHOTS = 1_048_576
WORDS = list(itertools.product(range(2), repeat=3))


def hands():
    car = federation(0)
    order = sorted(range(car.n), key=lambda p: (p // 12, (-p) % 12))
    _, sites, info = sector(car, order)
    assert info["omitted_states"] == 0
    reachable = set(sites.tolist())
    rank = np.argsort(order)
    edges = [
        (a, b)
        for a, b in car.seams
        if a // 12 == b // 12 == 6
        and a not in reachable
        and b not in reachable
        and rank[a] > rank[b]
    ]
    return (
        next(
            (e, f) for e, f in itertools.combinations(edges, 2) if set(e).isdisjoint(f)
        ),
        info,
    )


def instrument():
    law = candidate.native_law
    k = candidate.KAPPA
    dim = law.DIM
    half = np.eye(dim)
    zero = law.SRC.index((0, 0))
    blank = half[:, zero]
    signed = [
        half[:, (1 + law.LABELS.index((0, d))) * law.N_SRC + zero] for d in (1, -1)
    ]
    even = (signed[0] + signed[1]) / np.sqrt(2)
    odd = (signed[0] - signed[1]) / np.sqrt(2)
    h0, h1 = law.local_term_parts()
    h = h0 + k * h1
    free = np.kron(np.eye(2), h)
    # Basis: blank/fuel present; sign-even/fuel absent; sign-odd/fuel absent.
    v = np.zeros((2 * dim, 3))
    v[dim:, 0] = blank
    v[:dim, 1:] = np.column_stack((even, odd))
    drift = np.diag([0.0, k, k + 1])
    leakage = np.linalg.norm(free @ v - v @ drift)
    blind = np.linalg.norm(
        law.echo() @ np.column_stack(signed) - np.column_stack(signed)
    )
    pulses = []
    for i, (parity, angle) in enumerate(((1, np.pi / 2), (-1, np.pi / 4)), 1):
        # Same fuel-to-record couplings as the registered collar reader, defined on all 9 load inputs.
        transition = np.zeros((dim, dim))
        for load in range(law.N_SRC):
            for d, sign in ((1, 1), (-1, parity)):
                dst = (1 + law.LABELS.index((0, d))) * law.N_SRC + load
                transition[dst, load] = sign / np.sqrt(2)
        z = np.zeros_like(h)
        full = free + collar.G * np.block([[z, transition], [transition.T, z]])
        sub = drift.copy()
        sub[i, 0] = sub[0, i] = collar.G
        leakage = max(leakage, np.linalg.norm(full @ v - v @ sub))
        pulses.append(expm(-1j * sub * angle / collar.G))
    if max(leakage, blind) > 1e-10:
        raise ArithmeticError("the empty-endpoint sign or fuel-sector check failed")
    return (pulses[1] @ pulses[0])[0], np.diag(drift), leakage, blind


def angle_error(a, b):
    return float(np.angle(np.exp(1j * (a - b))))


def run():
    edges, inventory = hands()
    effect, levels, leak, blind = instrument()
    matter = phase.run()
    assert matter == json.loads((HERE / "phase.json").read_text())

    def probability(t, mode="positive"):
        if mode == "stationary":
            initial = np.array([0, 1, 0])
        else:
            initial = np.array([0, 1, 1 if mode == "positive" else -1]) / np.sqrt(2)
        return float(abs(effect @ (np.exp(-1j * levels * t) * initial)) ** 2)

    calibration = []
    for i, t in enumerate(np.arange(4) * np.pi / 2):
        p = probability(t)
        clicks = int(np.random.default_rng([221119, 0, i]).binomial(SHOTS, p))
        calibration.append(
            dict(advance=t, probability=p, clicks=clicks, quiet=SHOTS - clicks)
        )
    observations = np.array([r["clicks"] / SHOTS for r in calibration])
    offset = observations.mean()
    cos_coeff, sin_coeff = (observations[:2] - observations[2:]) / 2
    design = np.array([[cos_coeff, sin_coeff], [sin_coeff, -cos_coeff]])
    visibility = np.hypot(cos_coeff, sin_coeff)
    if visibility < 0.20:
        raise ArithmeticError("the calibrated fuel reader cannot resolve the hand")
    cases = []
    for i, t in enumerate(phase.TIMES):
        (witness,) = [
            r
            for r in matter["trials"]
            if r["orientation"] == "reversed_cell_labels"
            and r["receiver"] == 78
            and r["t"] == t
        ]
        for j, mode in enumerate(("positive", "negative", "stationary")):
            terminal = {}
            for m, mixer in enumerate(("plus", "minus")):
                p = [
                    probability(t, mode),
                    probability(t + np.pi / 2, mode),
                    witness["arms"]["coherent"]["settings"][mixer]["probability"],
                ]
                law = np.prod(np.where(np.asarray(WORDS), p, 1 - np.asarray(p)), axis=1)
                counts = np.random.default_rng([221119, 1, i, j, m]).multinomial(
                    SHOTS, law
                )
                obs = counts @ np.asarray(WORDS) / SHOTS
                quadrature = np.linalg.solve(design, obs[:2] - offset)
                angle = float(np.arctan2(quadrature[1], quadrature[0]) % (2 * np.pi))
                terminal[mixer] = dict(
                    probabilities=law.tolist(),
                    joint_counts=counts.tolist(),
                    marginal_click_rates=obs.tolist(),
                    hand_vector_norm=float(np.linalg.norm(quadrature)),
                    phase=angle if mode != "stationary" else None,
                    circular_error=(
                        angle_error(angle, t + int(mode == "negative") * np.pi)
                        if mode != "stationary"
                        else None
                    ),
                )
            contrast = (
                terminal["plus"]["marginal_click_rates"][2]
                - terminal["minus"]["marginal_click_rates"][2]
            )
            expected_matter = witness["arms"]["coherent"]["contrast"]
            cases.append(
                dict(
                    exposure=t,
                    preparation=mode,
                    settings=terminal,
                    recorded_matter_contrast=contrast,
                    prior_matter_contrast=expected_matter,
                    difference=contrast - expected_matter,
                )
            )
    angle_errors = [
        abs(s["circular_error"])
        for c in cases
        for s in c["settings"].values()
        if s["phase"] is not None
    ]
    stopped = np.array(
        [
            s["marginal_click_rates"][:2]
            for c in cases
            if c["preparation"] == "stationary"
            for s in c["settings"].values()
        ]
    )
    constraints = dict(
        max_held_circular_error_radians=max(angle_errors),
        max_stopped_hand_range=float(np.ptp(stopped, axis=0).max()),
        max_absolute_matter_contrast_difference=max(
            abs(c["difference"]) for c in cases
        ),
    )
    if (
        max(angle_errors) >= 0.025
        or constraints["max_stopped_hand_range"] >= 0.01
        or constraints["max_absolute_matter_contrast_difference"] >= 0.01
    ):
        raise ArithmeticError(f"held clock or matter comparison failed: {constraints}")
    return canonical(
        dict(
            scope="one terminal two-bit clock per trial in cell 6, time modulo 2*pi",
            active_seams_in_flight=630,
            reachable_states=inventory["reachable_states"],
            clock_seams=[list(e) for e in edges],
            matter_mixer=[82, 78],
            parity_energies=[float(levels[1]), float(levels[2])],
            hand_angular_frequency=float(levels[2] - levels[1]),
            fuel_gain=collar.G,
            complete_column_leakage=leak,
            blind_echo_identity_error=blind,
            blind_fuel_probability_with_both_converters_off=0.0,
            shots_per_setting=SHOTS,
            word_order=[list(w) for w in WORDS],
            fields=["cosine hand fuel", "quarter-advanced hand fuel", "matter click"],
            calibration=calibration,
            offset=offset,
            cos_coefficient=cos_coeff,
            sin_coefficient=sin_coeff,
            visibility=visibility,
            cases=cases,
            checks=constraints,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="create a fresh archive receipt")
    output = parser.parse_args().output
    result = run()
    if output:
        with output.open("x") as f:
            json.dump(result, f, allow_nan=False, indent=2)
            f.write("\n")
    elif result != json.loads((HERE / "clock.json").read_text()):
        raise SystemExit("The finite-patch clock receipt changed")
    print(result["clock_seams"], result["checks"])
