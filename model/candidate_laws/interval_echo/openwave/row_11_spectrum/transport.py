"""Interfere an ordinary charge impulse and read its freely moving material density."""

import argparse
import json

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from .measure import HERE, sector, symmetry_basis, word_table, portable
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from ..quark_winding.belt_walls import belt


SHOTS = 2**21
TIMES = (0, 16, 64, 256, 1024, 4096, 16384)


def observe(words, ground, excited, gap, received_gap, rng):
    n = words.shape[1]
    ports = np.arange(n)
    phase = -2j * np.pi * ports / n
    ramp = np.exp((words.astype(float) - 1) @ phase)
    amp = np.vdot(excited, ramp * ground)
    accepted = float(abs(amp) ** 2)
    kick_counts = int(rng.binomial(SHOTS, np.round(accepted, 10)))
    success = (1 + accepted) / 4
    panels = []
    for t in TIMES:
        state = (
            ground + amp * np.exp(-1j * round(float(gap), 9) * t) * excited
        ) / np.sqrt(1 + accepted)
        p = abs(state) ** 2
        read = []
        means = []
        for j in ports:
            tagged = np.bincount(words[:, j], weights=p, minlength=3) * success
            prob = np.r_[tagged, 1 - success]
            counts = rng.multinomial(SHOTS, np.round(prob, 10))
            read.append(counts.tolist())
            means.append((counts[2] - counts[0]) / sum(counts[:3]))
        z = np.exp(phase) @ means
        if t == 0:
            initial = z
        prediction = initial * np.exp(-1j * received_gap * t)
        panels.append(
            dict(
                time=t,
                port_012_and_failure=read,
                charge_harmonic=[z.real, z.imag],
                forecast=[prediction.real, prediction.imag],
                error=abs(z - prediction),
            )
        )
    return dict(
        kick_acceptance=accepted,
        kick_passes=kick_counts,
        paired_filter_rate=success,
        panels=panels,
    )


def run():
    cal = json.loads((HERE / "results.json").read_text())["belt"]
    rng = np.random.default_rng(110425)
    n = 10
    board = belt.ring(n)
    book = word_table(n)
    group = [np.roll(np.arange(n), -i) for i in range(n)]
    res = []
    for p in range(5):
        words, ids, kinetic, cost, _ = sector(board, p, book)
        for k in (8, 4):
            modes, energy, line_gaps = [], [], []
            for m in (0, 1):
                u = symmetry_basis(words, ids, book[0], group, m)
                h = (
                    u.conj().T @ (kinetic + sparse.diags(k * cost, dtype=float)) @ u
                ).toarray()
                if len(h) == 1:
                    e = np.array([h[0, 0].real])
                    v = np.eye(1)
                else:
                    e, v = eigh(h, subset_by_index=(0, 1))
                modes.append(np.asarray(u @ v[:, 0]).ravel())
                energy.append(e[0])
                line_gaps.append(float(e[1] - e[0]) if len(e) > 1 else None)
            if any(g is not None and g < 1e-8 for g in line_gaps):
                raise ValueError("Ideal scalar passband contains multiple modes")
            observed = cal[str(k)]["compositions"][p]["recoils"][0]["read"]["energy"]
            read = observe(words, *modes, energy[1] - energy[0], observed, rng)
            res.append(dict(p=p, k=k, separation_to_next=line_gaps, **read))
    min_draws = min(
        sum(port[:3])
        for c in res
        for row in c["panels"]
        for port in row["port_012_and_failure"]
    )
    # Split 1% error between the two experiments; union-bound both quadratures.
    angle_error = np.arcsin(np.sqrt(2) * np.sqrt(2 * np.log(2 * 580 / 0.005) / 262144))
    count_radius = np.sqrt(2) * np.sqrt(2 * np.log(2 * 140 / 0.005) * 5 / min_draws)
    for c in res:
        initial = np.linalg.norm(c["panels"][0]["charge_harmonic"])
        for r in c["panels"]:
            bound = 2 * count_radius + (initial + count_radius) * r["time"] * (
                angle_error / 65536 + 1e-9
            )
            r["simultaneous_error_bound"] = bound
            if r["error"] > bound:
                raise ValueError(
                    "Held signed pattern disagrees with its recoil frequency"
                )
    return dict(
        shots=SHOTS,
        seed=110425,
        min_draws=min_draws,
        charge_read_radius=count_radius,
        phase_read_radius=angle_error,
        cases=res,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    p = HERE / "transport.json"
    d = json.loads(json.dumps(run(), allow_nan=False))
    if parser.parse_args().write:
        p.write_text(json.dumps(d, indent=1) + "\n")
    else:
        compare(portable(d), portable(json.loads(p.read_text())), "transport")
    print(
        "Verified the kick, preparation failures and signed port records for all ten cases."
    )


if __name__ == "__main__":
    main()
