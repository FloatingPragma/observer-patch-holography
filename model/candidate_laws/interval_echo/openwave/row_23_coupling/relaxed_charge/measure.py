"""Static phase and address records in a bath shared by two ternary sign probes."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import diags
from scipy.sparse.csgraph import connected_components
from model.executable_model import EMITTER_SEAM, federation

from ..free_bath import native
from . import readback as read

OUTPUT = Path(__file__).with_name("records.json.gz")


def occupied(h):
    """The registered selector accepts a uniformly populated unresolved lowest line."""
    n = h.shape[0]
    block = h.toarray()
    take = min(n, 4)
    while True:
        vals, vectors = eigh(block, subset_by_index=(0, take - 1), check_finite=False)
        keep = vals <= vals[0] + read.LINE_WIDTH
        if take == n or not keep.all():
            break
        take = min(2 * take, n)
    if abs(vals[0]) * read.TIMES[0] + read.ANGLE_ERROR >= np.pi:
        raise ArithmeticError("First clock angle is ambiguous")
    vec = vectors[:, keep]
    if np.max(abs(h @ vec - vec * vals[keep])) > 2e-10:
        raise ArithmeticError("Invalid stationary line")
    return vals[keep], np.mean(abs(vec) ** 2, axis=1)


def quantize(p):
    p = np.asarray(p, dtype=float)
    if p.min() < -1e-12 or abs(p.sum() - 1) > 2e-10:
        raise ArithmeticError("Invalid terminal probability")
    ticks = np.rint(read.GRID * p.clip(0, 1)).astype(np.int64)
    ticks[p.argmax()] += read.GRID - ticks.sum()
    return ticks / read.GRID


def clock_bank(lines, rng):
    energies = np.rint(lines * 2**32) / 2**32
    amplitude = np.exp(1j * np.outer(read.TIMES, energies)).mean(axis=1)
    probs = (1 + np.column_stack((amplitude.real, amplitude.imag))) / 2
    plus = np.rint(probs.clip(0, 1) * read.GRID) / read.GRID
    return rng.binomial(read.SHOTS, plus).tolist()


def acquire(h, rng, *, address):
    vals, p = occupied(h)
    out = dict(
        clock=clock_bank(vals, rng),
        multiplicity=len(vals),
        dimension=len(p),
        energy_lines=[round(float(v), 10) for v in vals],
    )
    if address:
        out["address"] = rng.multinomial(read.SHOTS, quantize(p)).tolist()
    return out


def calibration(k, rng):
    e = native.contacts(k).diagonal().reshape(3, 3)
    gain, counts = {}, {}
    for t in (-1, 1):
        # Four matched holds compare adding a vacancy beside a charged or neutral guard.
        tau = np.arange(256) * (4 * np.pi / 256)
        amp = np.ones(len(tau), complex)
        for mobile, pin, sign in [(0, t + 1, 1), (1, t + 1, -1), (0, 1, -1), (1, 1, 1)]:
            amp *= np.exp(-1j * tau * sign * e[mobile, pin])
        q = quantize(abs(np.fft.ifft(amp)) ** 2)
        bins = rng.multinomial(read.SHOTS, q)
        peak = int(bins.argmax())
        gain[str(t)] = (peak if peak < 128 else peak - 256) / 2
        counts[str(t)] = bins.tolist()
    direct_counts, bins = [], []
    for adjacent in (False, True):
        tau = 2 * np.pi * np.arange(256) / 256
        amp = np.ones(len(tau), complex)
        for a, b, sign in ((2, 2, 1), (2, 0, -1), (0, 2, -1), (0, 0, 1)):
            energy = e[a, 0] + e[a, 2] + e[b, 1] + (e[a, b] if adjacent else 0)
            amp *= np.exp(-1j * tau * sign * energy)
        record = rng.multinomial(read.SHOTS, quantize(abs(np.fft.ifft(amp)) ** 2))
        peak = int(record.argmax())
        bins.append(peak if peak < 128 else peak - 256)
        direct_counts.append(record.tolist())
    return dict(
        gain=gain,
        local_clock_counts=counts,
        direct_four_sign_bins=bins,
        direct_counts=direct_counts,
    )


def scan(board, home, k, cal, rng, shutter_rng):
    d = native.distances(board, home)
    h, free = native.vacancy_operator(board, k, {home: 0})
    source = acquire(h, rng, address=True)
    source["radius"] = read.radius(source["address"], d[free])
    shells = sorted({1, 2, 3, 4, int(d.max() // 2), int(d.max())})
    rows = []
    for shell in shells:
        y = int(np.flatnonzero(d == shell)[0])
        half = native.distances(board, y)
        full, cut, shutter, neutral = {}, {}, {}, {}
        for s in (-1, 1):
            for t in (0, -1, 1):
                g, free = native.vacancy_operator(board, k, {home: s + 1, y: t + 1})
                if connected_components(g, directed=False, return_labels=False) != 1:
                    raise ArithmeticError(
                        "The intact stationary stock has disconnected supports"
                    )
                key = f"{s},{t}"
                full[key] = acquire(g, rng, address=True)
                if t == 0:
                    neutral[s] = g
                else:
                    sides = (d[free] <= half[free], d[free] > half[free])
                    cut[key] = [
                        acquire(g[side][:, side], rng, address=False) for side in sides
                    ]
                    if shell >= 3:
                        where = d <= half
                        lost = np.array(
                            [
                                sum(where[u] != where[v] for u in board.neighbours[v])
                                for v in free
                            ]
                        )
                        e = native.contacts(k).diagonal().reshape(3, 3)
                        changed = g - diags(lost * (e[0, 1] - e[1, 1]))
                        shutter[key] = [
                            acquire(changed[side][:, side], shutter_rng, address=False)
                            for side in sides
                        ]
        mask = np.isin(free, board.neighbours[y])
        direct = cal["direct_four_sign_bins"][int(shell == 1)] / 4
        summary = read.report(
            full, cut, {"neighbors": mask}, cal["gain"], direct, neutral
        )
        if shutter:
            sides = {
                key: [read.clock_energy(bank["clock"]) for bank in banks]
                for key, banks in shutter.items()
            }
            summary["shutter_one_vacancy"] = read.signed(
                {key: min(pair) for key, pair in sides.items()}
            )
            summary["shutter_two_vacancies"] = read.signed(
                {key: sum(pair) for key, pair in sides.items()}
            )
        rows.append(
            dict(
                distance=shell,
                port=y,
                intact=full,
                cut=cut,
                shutter=shutter,
                readings=summary,
            )
        )
    return dict(owners=board.cells, source=home, source_alone=source, shell_banks=rows)


def run():
    rng = np.random.Generator(np.random.PCG64(23042917))
    cal = {str(k): calibration(k, rng) for k in (8, 0)}
    shutter_rng = np.random.Generator(np.random.PCG64(23042917).jumped())
    boards = {}
    for level in (0, 1):
        board = federation(level)
        native.location_instrument(board)
        for k in (8, 0):
            boards[f"{board.cells},{k}"] = scan(
                board, EMITTER_SEAM[level][0], k, cal[str(k)], rng, shutter_rng
            )
    return dict(
        shots=read.SHOTS,
        clock_times=read.TIMES.tolist(),
        calibration=cal,
        boards=boards,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    # Floats only annotate the receipt; all raw terminal data are integers.
    payload = json.dumps(
        run(), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    if args.output:
        args.output.write_bytes(gzip.compress(payload, mtime=0))
    elif gzip.decompress(OUTPUT.read_bytes()) != payload:
        raise ValueError("Deposited relaxed-charge record differs from the replay")
    else:
        print("All terminal records and decoded energies replayed.")


if __name__ == "__main__":
    main()
