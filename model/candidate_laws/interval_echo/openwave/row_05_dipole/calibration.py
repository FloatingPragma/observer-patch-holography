"""Measure the local clock, material charge and external twist on independent probes."""

import argparse
import json
import numpy as np
from scipy.linalg import expm
from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell
from model.experiments.interval_echo_encounter.conjugation import free_seam
from .experiment import HERE, CHARGE, FIELDS, portable, source

SHOTS = 2**26


def read(value, rng):
    p = np.round(np.clip((1 + value) / 2, 0, 1), 13)
    k = int(rng.binomial(SHOTS, p))
    return [k, SHOTS - k]


def decode(counts):
    return (counts[0] - counts[1]) / SHOTS


def run():
    rng = np.random.default_rng(510244)
    h = free_seam(candidate.KAPPA)[
        np.ix_(
            [law.SRC.index(e) for e in ((0, 1), (1, 0))],
            [law.SRC.index(e) for e in ((0, 1), (1, 0))],
        )
    ]
    clock = []
    for t in range(257):
        a, b = expm(-1j * h * t)[:, 0]
        clock.append(
            [read(abs(a) ** 2 - abs(b) ** 2, rng), read(2 * (a.conj() * b).imag, rng)]
        )
    angles = np.unwrap([np.arctan2(decode(y), decode(x)) for x, y in clock])
    angles -= angles[0]

    # Sign calibration from the explicit signed outcome on a unit reference.
    def signed(a):
        v = law.echo()[:, law.SRC.index((a, 1))]
        blocks = (abs(v) ** 2).reshape(-1, law.N_SRC).sum(axis=1)
        signed_prob = np.zeros(5)
        for weight, (_, sign) in zip(blocks[1:], law.LABELS, strict=True):
            signed_prob[int(sign) + 2] += weight
        return signed_prob

    positive_or_negative = int(np.argmax(signed(0)) - 2)
    assert positive_or_negative == CHARGE and signed(1)[2] > 1 - 1e-12
    primary = json.loads((HERE / "results.json").read_text())
    plus, areas, _, _ = source(tuple(primary["owner_zero_gates"]))
    charge_counts = rng.multinomial(SHOTS, (abs(plus) ** 2).round(12)).tolist()
    by_edge = dict(zip(cell().seams, areas, strict=True))

    def area(a, b):
        return by_edge[tuple(sorted((a, b)))] * (1 if a < b else -1)

    face = max(
        cell().faces(),
        key=lambda f: abs(sum(area(a, b) for a, b in zip(f, (f[1], f[2], f[0])))),
    )
    walk = list(zip(face, (face[1], face[2], face[0])))
    triangle = sum(area(a, b) for a, b in walk)
    fields = {}
    for B in FIELDS:
        tours = []
        for cycle in (walk, [(b, a) for a, b in walk[::-1]]):
            amp = 1.0 + 0j
            for a, b in cycle:
                step = h.astype(complex)
                step[1, 0] *= np.exp(1j * CHARGE * B * area(a, b))
                step[0, 1] = step[1, 0].conj()
                amp *= expm(-1j * np.pi * step)[1, 0]
            tours.append(amp)
        z = tours[0] * tours[1].conjugate()
        counts = [read(z.real, rng), read(z.imag, rng)]
        measured = np.arctan2(decode(counts[1]), decode(counts[0])) / (
            2 * positive_or_negative * triangle
        )
        fields[str(B)] = dict(counts=counts, measured=measured)
    slope = primary["current"]["plus"]["moment"] - primary["current"]["minus"]["moment"]
    target = next(
        r for r in primary["precession"] if r["owners"] == 320 and r["kappa"] == 8
    )
    holds = []
    for r in target["rows"]:
        for obs in r["read"]:
            x, y = np.array(obs["counts_xy"])[:, :2] @ np.array([1, -1])
            measured = np.arctan2(y, x)
            t = obs["time"]
            # JSON encodes floating fields and the literal zero distinctly.
            B = r["field"]
            bhat = fields[str(B)]["measured"]
            holds.append(
                dict(
                    time=t, field=B, phase=measured, forecast=-slope * bhat * angles[t]
                )
            )
    return dict(
        shots=SHOTS,
        seed=510244,
        clock_counts_z_y=clock,
        clock_read=angles.tolist(),
        signed_echo_bins={str(a): signed(a).tolist() for a in (0, 1, 2)},
        vacancy_address_counts=charge_counts,
        charge=positive_or_negative,
        face=face,
        oriented_area=triangle,
        fields=fields,
        held=holds,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    data = portable(run())
    path = HERE / "calibration.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(data, indent=1) + "\n")
    elif json.loads(path.read_text()) != data:
        raise ValueError("Local instrument calibrations differ")
    print(
        "Verified an independent material clock, signed charge and opposite-tour flux calibration."
    )


if __name__ == "__main__":
    main()
