"""Count all material endpoints and compare to the independent guide Green flux."""

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from .. import measure as parent
from . import receiver

HERE = Path(__file__).resolve().parent
PLAN = json.loads((HERE / "protocol.json").read_text())


@dataclass(frozen=True)
class Setting:
    name: str
    face: tuple
    probe: int
    source: tuple
    blind: bool = False
    stiffness: int = 8
    short: bool = False
    parked: bool = False
    equal_degree: bool = False
    same_sign: bool = False


def settings():
    faces = [tuple(face) for face in PLAN["faces"]]
    for face in faces:
        for probe in (1, 2):
            yield Setting("matched", face, probe, (2, 2, 2))
    face = faces[2]
    for source in ((1, 1, 1), (1, 2, 1), (2, 1, 2)):
        for probe in (1, 2):
            yield Setting("changed source", face, probe, source)
    for source in ((1, 2, 1), (2, 1, 2)):
        for probe in (1, 2):
            yield Setting("blind", face, probe, source, blind=True)
    for probe in (1, 2):
        for control in ("short", "parked", "equal_degree", "same_sign"):
            yield Setting(control, face, probe, (1, 2, 1), **{control: True})
        yield Setting("changed stiffness", face, probe, (1, 2, 1), stiffness=4)
        yield Setting("empty guide", face, probe, (1, 0, 1), blind=True)


def table(setting):
    time = receiver.RELEASE / (2 if setting.short else 1)
    return np.array(
        [
            [
                receiver.right_probability(
                    setting.probe,
                    setting.stiffness,
                    word,
                    setting.source,
                    setting.blind,
                    sign,
                    time,
                    setting.parked,
                    setting.equal_degree,
                )
                for word in receiver.WORDS
            ]
            for sign in ((1, 1) if setting.same_sign else (-1, 1))
        ]
    )


def signed_readout(setting):
    """Analytic prediction from the ternary phases before solving the guide field."""
    time = receiver.RELEASE / (2 if setting.short else 1)
    detuning = 0 if setting.equal_degree else 0.5
    amplitude = (
        detuning * (1 - np.cos(time * np.hypot(1, detuning))) / (2 * (1 + detuning**2))
    )
    if setting.parked:
        amplitude = 0
    result = np.zeros((10, 10, 10))
    for sign, words in zip(
        (1 if setting.same_sign else -1, 1), parent.addresses(setting.face), strict=True
    ):
        if setting.stiffness == 4:
            weight = np.ones_like(words)
        else:
            weight = (words == 7) * (all(setting.source) / (8 if setting.blind else 1))
        result += sign * amplitude / 2 * weight
    return result


def visits_before_drain():
    f = np.zeros((10, 10, 10))
    for corner in parent.green.CORNERS:
        f[tuple(v - 1 for v in corner)] = 1 / 8
    return parent.green.poisson(f, (1 / 8,) * 3)[1:-1, 1:-1, 1:-1]


def count_panel(rng, appointments, setting):
    # Marginalize the independently chosen trits and standards before the binary
    # draw. The outcome law is identical to drawing each intermediate switch.
    probs = table(setting)
    near = rng.binomial(appointments, 0.5)
    counts = []
    for weights, indices, values in zip(
        (near, appointments - near), parent.addresses(setting.face), probs, strict=True
    ):
        cohort = np.array(
            [weights[indices == i].sum() for i in range(8)], dtype=np.int64
        )
        rights = rng.binomial(cohort, np.rint(values * 2**32) / 2**32)
        counts.append(dict(starts=cohort.tolist(), right=rights.tolist()))
    return counts


def bound(reward, u, panels, releases=None, alpha=0.01):
    # Bound the second moment of an entire independent, deadline-stopped source
    # walk. A left or right read centered at 1/2 has square exactly 1/4.
    half_grid = 2**-33
    a = abs(reward) + half_grid
    remaining = parent.green.poisson(a, (1 / 8,) * 3)[1:-1, 1:-1, 1:-1]
    moment = np.sum(u * (0.25 + 2 * a * remaining))
    x = np.log(2 * panels / alpha)
    n = PLAN["releases"] if releases is None else releases
    h = PLAN["maximum_appointments"]
    spectral = 1 - 0.75 * (1 - np.cos(np.pi / 11))
    tail = 0.5 * np.sqrt(u.size) * spectral**h / (1 - spectral)
    return float(
        np.sqrt(2 * moment * x / n)
        + 2 * h * x / (3 * n)
        + 120 * half_grid
        + tail
        + 1e-11
    )


def calibrations():
    """Separate ordinary-body and dimer banks use the registered shutters."""
    rng = np.random.default_rng(61493)
    shots = 2**28
    records = {}
    for name, duration in [
        ("trit subpulse", np.pi / 3),
        ("full release", receiver.RELEASE),
        ("half release", receiver.RELEASE / 2),
    ]:
        u = receiver.original.rails.candidate.native_law
        exchange, contrast = u.local_term_parts()
        rail = [u.SRC.index((1, 0)), u.SRC.index((0, 1))]
        rotation = expm(-1j * duration * (exchange + 8 * contrast)[np.ix_(rail, rail)])
        probability = abs(rotation[1, 0]) ** 2
        records[name] = int(rng.binomial(shots, round(probability, 12)))
    for name, flags in [
        ("ordinary full match", (1, 1, 1)),
        ("ordinary single miss", (1, 0, 1)),
    ]:
        prob = receiver.right_probability(
            1, 8, flags, (2, 2, 2), False, 1, receiver.RELEASE, False, False
        )
        records[name] = int(rng.binomial(shots, round(prob, 12)))
    return dict(
        starts_per_bank=shots,
        terminal_counts=records,
        trit_subpulses=[0, 4, 8],
        clock_branch=[0, "pi"],
    )


def run():
    if parent.N != PLAN["releases"] or parent.H != PLAN["maximum_appointments"]:
        raise ValueError(
            "guide acquisition must match the frozen release count and deadline"
        )
    rng = np.random.default_rng(PLAN["seed"])
    actual, ticks = parent.visits(rng)
    u = visits_before_drain()
    total = int(actual.sum())
    configs = list(settings())
    panels = []
    for setting in configs:
        arms = count_panel(rng, actual, setting)
        count = sum(sum(row["right"]) for row in arms)
        shift = (count - total / 2) / parent.N
        reward = signed_readout(setting)
        expected = float(np.sum(u * reward))
        error = bound(reward, u, len(configs))
        if abs(shift - expected) > error:
            raise ArithmeticError(
                f"{setting} failed the held comparison: {shift} vs {expected} +/- {error}"
            )
        panels.append(
            dict(
                setting=asdict(setting),
                arms=arms,
                mean_per_release=shift,
                predicted=expected,
                simultaneous99_half_width=error,
            )
        )
    return dict(
        protocol="protocol.json",
        calibration=calibrations(),
        source_releases=parent.N,
        appointments_per_panel=total,
        last_drain=ticks,
        native_edge_and_pads=receiver.apparatus()[:4],
        equal_degree_control=receiver.apparatus(True)[:4],
        panels=panels,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    path = HERE / "results.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        previous = json.loads(path.read_text())
        if any(
            new["arms"] != old["arms"]
            for new, old in zip(result["panels"], previous["panels"], strict=True)
        ):
            raise ValueError("seeded right-port acquisitions differ")
