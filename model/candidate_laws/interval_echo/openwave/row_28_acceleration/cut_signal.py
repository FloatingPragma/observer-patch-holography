"""Transmit a bit across a material cut using one receiving patch's port record."""

import argparse
import json

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components

from model.candidate_laws.interval_echo import candidate

from .experiment import HERE, layout
from .signals import adjacency

LAW = candidate.native_law
SHOTS = 4096
SEED = 280222
DURATIONS = (0, 1, 2)  # in units of pi/64


def deposited_collar(sign, coin):
    initial = (2, 0) if sign == 2 else (0, 2)
    collar = (
        LAW.echo()[:, LAW.SRC.index(initial)]
        .reshape(-1, 9)[:, LAW.SRC.index((1, 1))]
        .copy()
    )
    if coin is not None:
        collar.fill(0)
        collar[1 + LAW.LABELS.index((coin, sign))] = 1

    return collar


def receive(collar, bit, direction, contact, half_swap):
    """Apply one physical readout, retaining joint conditional memory amplitudes."""
    # The receiver owns the addressed port and the rail. Keep both number sectors
    # in the propagation and trace out the sender and the collar at the end.
    wave = np.zeros((11, 3, 3, 3), dtype=complex)
    wave[:, bit, :, :] = collar[:, None, None] * half_swap[
        :, LAW.SRC.index((0, 1))
    ].reshape(3, 3)
    if direction == "right":
        wave = wave.swapaxes(1, 2)
    wave = (contact @ wave.reshape(99, 3)).reshape(11, 3, 3, 3)
    if direction == "right":
        wave = wave.swapaxes(1, 2)
    recovered = (wave.reshape(33, 9) @ half_swap.conj()).reshape(11, 3, 3, 3)
    return recovered


def port_distribution(sign, coin, bit, direction, contact, half_swap):
    recovered = receive(
        deposited_collar(sign, coin), bit, direction, contact, half_swap
    )
    return (abs(recovered) ** 2).sum(axis=(0, 1, 3))


def run():
    n, edges, *_ = layout(1)
    with np.load(HERE / "records.npz", allow_pickle=False) as bank:
        closed = bank[f"epoch_{SEED}"] > 0
        signs = bank[f"sign_{SEED}"]
    _, groups = connected_components(adjacency(n, edges[~closed]), directed=False)
    component = groups == groups[265]
    cut = np.flatnonzero(component[edges[:, 0]] != component[edges[:, 1]])
    assert component.sum() == 468 and len(cut) == 25 and closed[cut].all()
    # Check that each endpoint can park its rail locally on the same side of the cut.
    ports_with_local_rail = {int(v) for ep in edges[~closed] for v in ep}
    assert all(int(v) in ports_with_local_rail for e in cut for v in edges[e])
    exchange, contrast = LAW.local_term_parts()
    rng = np.random.default_rng(281117)
    readings = []
    for kappa in (candidate.KAPPA, 0):
        half = expm(-0.5j * np.pi * LAW.seam_term(kappa))
        for duration in DURATIONS:
            contact = expm(-1j * np.pi * duration / 64 * (exchange + kappa * contrast))
            for sign in (-2, 2):
                for direction in ("left", "right"):
                    for coin in (None, 0, 1):
                        by_bit = [
                            port_distribution(sign, coin, bit, direction, contact, half)
                            for bit in (0, 1)
                        ]
                        # Exact written-sector identity at zero duration or sender bit zero.
                        for bit, p in enumerate(by_bit):
                            assert abs(p.sum() - 1) < 1e-12 and p[2] < 1e-25
                            expected = (
                                1
                                - (
                                    np.cos(4 * kappa * np.pi * duration / 64)
                                    + np.cos((4 * kappa + 1) * np.pi * duration / 64)
                                )
                                / 2
                            ) / 2
                            assert (
                                abs(p[1] - bit * expected / (1 if coin is None else 2))
                                < 2e-13
                            )
                        for edge in cut:
                            if signs[edge] != sign:
                                continue
                            counts = [int(rng.binomial(SHOTS, p[1])) for p in by_bit]
                            readings.append(
                                [
                                    int(edge),
                                    direction,
                                    coin,
                                    int(kappa),
                                    duration,
                                    *counts,
                                ]
                            )
    return dict(
        board_cells=80,
        seed=SEED,
        station=265,
        ports_in_component=int(component.sum()),
        endpoints={str(i): edges[i].tolist() for i in cut},
        shots_per_bit=SHOTS,
        columns=[
            "seam",
            "sending_endpoint",
            "read_coin_or_null",
            "kappa",
            "exposure_pi_over_64",
            "clicks_bit_zero",
            "clicks_bit_one",
        ],
        records=readings,
    )


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", action="store_true")
    args = cli.parse_args()
    path = HERE / "cut_signal.json"
    observed = run()
    if args.record:
        with path.open("x") as output:
            json.dump(observed, output, separators=(",", ":"))
            output.write("\n")
    else:
        assert observed == json.loads(path.read_text())
        print(
            "Verified every sender setting and receiver count across all 25 material cuts."
        )


if __name__ == "__main__":
    main()
