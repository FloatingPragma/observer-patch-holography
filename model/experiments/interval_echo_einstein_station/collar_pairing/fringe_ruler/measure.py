"""A registered four-tick two-emitter interferometer on the formed 80-cell board."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.clock import _dimer_pulse
from ..ranging import measure as radar

HERE = Path(__file__).resolve().parent
PHASES = (0, 2, 4, 6, 1, 7)  # in pi/4; the last two quadratures are held
STOCKS = ("product", "half", "formed", "ordinary", "scrambled")
SHOTS = 2**48
SEED = 823517


def board_setup():
    source = radar.formation
    board = source.federation(1)
    ports = list(source.region(board))
    cell = ports[0] // 12
    emit = min(
        (a, b)
        for a, b in board.seams
        if a // 12 == cell == b // 12
        and board.degree[a] == board.degree[b] == 5
        and a not in ports
        and b not in ports
    )
    return board, ports, emit


def prepare(charge, phase):
    """Isolate the bare seam, then the second rail's private term on unit background."""
    law = radar.formation.emission.LAW
    exchange, contact = law.local_term_parts()
    h = exchange + 8 * contact
    eye = np.eye(law.DIM)[
        :, [law.SRC.index(word) for word in ((1 + charge, 1), (1, 1 + charge))]
    ]
    dimer = eye.T @ h @ eye
    np.testing.assert_allclose(h @ eye, eye @ dimer, atol=1e-14)
    np.testing.assert_allclose(dimer, [[8.5, -0.5], [-0.5, 8.5]], atol=1e-14)
    tick = np.pi / 96
    raw = expm(-1j * 48 * tick * dimer)[:, 0]
    contrast = 8 * ((1 + charge) ** 2 - 1)
    dwell = ((np.pi / 2 - phase * np.pi / 4) * np.sign(contrast) % (2 * np.pi)) / abs(
        contrast
    )
    pulses = int(round(dwell / tick))
    assert 0 <= pulses <= 21 and abs(pulses * tick - dwell) < 1e-14
    raw[1] *= np.exp(-1j * contrast * pulses * tick)
    return raw


def screens(charge, odd, board, ports, emit):
    gen = radar.generator(board, ports, charge, 8, odd)
    indices = [4 * v + i for v in emit for i in range(4)]
    src = np.zeros((gen.shape[0], 8))
    src[indices] = np.eye(8)
    flight = -4j * gen
    waves = expm_multiply(flight, src, traceA=flight.diagonal().sum())
    return waves.reshape(board.n, 4, 2, 4)[ports]


def frequencies(waves, coin, rails, scramble=False):
    def read(packet):
        return np.einsum("yoi,ij,yoj->y", packet, coin, packet.conj()).real

    first, second = (waves[:, :, j] * rails[j] for j in (0, 1))
    return read(first) + read(second) if scramble else read(first + second)


def populations():
    board, ports, emit = board_setup()
    rho, tau, _ = radar.prepared_states(board, ports, 1)
    panels = {}
    for charge in (-1, 1):
        even, odd = (
            screens(charge, sector, board, ports, emit) for sector in (False, True)
        )
        panels[charge] = {}
        for phase in PHASES:
            rails = prepare(charge, phase)
            for name, coin in zip(
                STOCKS, (tau, (rho + tau) / 2, rho, tau, rho), strict=True
            ):
                p = frequencies(even, coin, rails, name == "scrambled")
                if name == "ordinary":
                    p = (p + frequencies(odd, coin, rails)) / 2
                panels[charge][name, phase] = p
    return ports, emit, panels


def acquire():
    ports, emit, pops = populations()
    rng = np.random.default_rng(SEED)
    raw = {}
    for q, panels in pops.items():
        raw[str(q)] = {}
        for stock in STOCKS:
            rows = []
            for p in PHASES:
                probs = panels[stock, p]
                cells = rng.multinomial(SHOTS, [*probs, 1 - probs.sum()])
                rows.append(cells.tolist())
            raw[str(q)][stock] = rows
    hand = abs(_dimer_pulse(1)[1, 0]) ** 2
    return dict(
        ports=ports,
        emit=list(emit),
        shots=SHOTS,
        seed=SEED,
        phases_in_pi_over_four=list(PHASES),
        four_equal_ticks={
            "one_tick_starts": SHOTS,
            "one_tick_hits": int(rng.binomial(SHOTS, hand)),
        },
        counts=raw,
    )


def main():
    from .readout import examine

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    file = HERE / "results.json"
    if args.write:
        receipt = acquire()
        receipt["readout"] = examine(receipt)
        with file.open("x") as fp:
            fp.write(json.dumps(receipt, indent=2) + "\n")
    else:
        receipt = json.loads(file.read_text())
        verify(receipt)
        print(json.dumps(examine(receipt), indent=2))


def verify(receipt):
    from .readout import radius

    ports, emit, live = populations()
    assert ports == receipt["ports"] and list(emit) == receipt["emit"]
    assert receipt["shots"] == SHOTS and receipt["phases_in_pi_over_four"] == list(
        PHASES
    )
    for q in (-1, 1):
        for stock, rows in receipt["counts"][str(q)].items():
            for phase, counts in zip(PHASES, rows, strict=True):
                assert len(counts) == 5 and min(counts) >= 0 and sum(counts) == SHOTS
                p = live[q][stock, phase]
                seen = np.asarray(counts) / SHOTS
                assert np.all(abs(np.r_[p, 1 - p.sum()] - seen) <= radius(seen, SHOTS))
    clocks = receipt["four_equal_ticks"]
    seen = clocks["one_tick_hits"] / clocks["one_tick_starts"]
    assert abs(seen - abs(_dimer_pulse(1)[1, 0]) ** 2) < radius(seen, SHOTS)


if __name__ == "__main__":
    main()
