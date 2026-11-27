"""Measure the first fugacity coefficient of signed static work from real-clock traces."""

import argparse
from functools import cache
import gzip
import json
from pathlib import Path
from math import sqrt, log

import numpy as np
from scipy.linalg import eigvalsh
from scipy.sparse import diags
from scipy.special import betainc
from model.executable_model import EMITTER_SEAM, federation
from ..free_bath import native
from . import measure, readback as read

OUTPUT = Path(__file__).with_name("virial.json.gz")
BETAS = (1 / 16, 1 / 8)
SHOTS = 2**42
FREQ = 2 * np.pi * np.arange(1, 193) / 256
TAIL = 0.01 / 2_000_000


def step(u):
    t = np.clip(u, 0.0, 1.0)
    return t**3 * (10 - 15 * t + 6 * t**2)


@cache
def weights(beta, band=None):
    n = 8192
    energy = 256 * np.arange(n) / n
    energy[energy >= 192] -= 256
    if band is None:
        lo, hi = -16, 112
        taper = step((energy + 32) / 16) * step((144 - energy) / 32)
    else:
        lo, hi = band
        taper = betainc(8, 8, np.clip((energy - lo + 24) / 24, 0, 1)) * betainc(
            8, 8, np.clip((hi + 24 - energy) / 24, 0, 1)
        )
    f = np.exp(-beta * energy) * taper
    c = np.fft.rfft(f)[:193] / n
    # On [-16,112] the target is the analytic exponential. Linear interpolation
    # between checks costs at most spacing² sup|f''-P''|/8.
    fine = 2**19
    coef = np.zeros(fine // 2 + 1, complex)
    coef[: len(c)] = c
    approx = np.fft.irfft(coef, n=fine) * fine
    q = np.arange(lo * fine // 256, hi * fine // 256 + 1)
    energies = 256 * q / fine
    abs_err = abs(approx[q] - np.exp(-beta * energies))
    curvature = 2 * np.abs(c[1:]) @ (FREQ**2) + beta**2 * np.exp(-lo * beta)
    return c, energies, abs_err, (256 / fine) ** 2 * curvature / 8 + 2e-12


def bounds(h):
    diag = h.diagonal()
    off = np.asarray(abs(h).sum(axis=1)).ravel() - abs(diag)
    lo = float(min(diag - off))
    hi = float(max(diag + off))
    if lo < -16 or hi > 112:
        raise ValueError("Energy outside enrolled thermal interval")
    return lo, hi


def trace_bank(lines, rng, shots):
    characteristic = np.exp(1j * np.outer(FREQ, lines)).mean(axis=1)
    proba = (1 + np.column_stack((characteristic.real, characteristic.imag))) / 2
    # Archive the actual probabilities on a dyadic grid to separate simulation from measurement.
    grid = 2**40
    numerators = np.rint(np.clip(proba, 0, 1) * grid).astype(np.int64)
    return dict(
        plus=rng.binomial(shots, numerators / grid).tolist(),
        probability_ticks=numerators.tolist(),
        shots=shots,
    )


def thermal(bank, beta, n, enclosure, *, band=None):
    c, x, err, interpolation = weights(beta, band)
    lo, hi = enclosure
    start = max(0, np.searchsorted(x, lo) - 1)
    stop = min(len(x), np.searchsorted(x, hi) + 1)
    bias = float(max(err[start:stop]) + interpolation)
    v = 2 * np.array(bank["plus"]) / bank["shots"] - 1
    value = c[0].real + 2 * np.real(c[1:] @ (v[:, 0] + 1j * v[:, 1]))
    # A weighted sum of independent two-quadrature Bernoulli draws.
    noise = 4 * np.linalg.norm(c[1:]) * sqrt(log(2 / TAIL) / (2 * bank["shots"]))
    rounding = 4 * np.abs(c[1:]).sum() / 2**40
    return float(n * value), float(n * (bias + noise + rounding))


def ramsey_line(energy, rng):
    z = np.exp(1j * read.TIMES * energy)
    ticks = np.rint((1 + np.column_stack([z.real, z.imag])) / 2 * 2**40)
    return rng.binomial(SHOTS, ticks / 2**40).tolist()


def blank_clock(k, rng, gain_rng):
    h = native.contacts(k)
    block = h[np.ix_([1, 3], [1, 3])]
    ordinary = h[4, 4]
    sources = [np.array([1, 1]) / sqrt(2), np.array([1, -1]) / sqrt(2)]
    energy = [float(v @ block @ v - ordinary) for v in sources]
    _, private = native.LAW.local_term_parts()
    # A missing degree has the private squared-load service with strength k.
    private_phase = -float(
        k * private[native.LAW.SRC.index((1, 0)), native.LAW.SRC.index((1, 0))]
    )
    energy.append(private_phase)
    banks = [ramsey_line(e, rng) for e in energy]
    # Subtract the ordinary and vacancy cost on the two native pin settings.
    e = h.diagonal().reshape(3, 3)
    gain_banks = {
        str(t): ramsey_line(e[0, t + 1] - e[1, t + 1] - e[0, 1] + e[1, 1], gain_rng)
        for t in (-1, 1)
    }
    observed = [read.clock_energy(b, shots=SHOTS) for b in banks]
    de = (
        np.arcsin(2 * sqrt(2) * (sqrt(log(2 / TAIL) / (2 * SHOTS)) + 1 / 2**40))
        / read.TIMES[-1]
    )
    return dict(
        plus=banks,
        gain_plus=gain_banks,
        gains={t: read.clock_energy(b, shots=SHOTS) for t, b in gain_banks.items()},
        diagonal=(observed[1] + observed[0]) / 2,
        hopping=(observed[0] - observed[1]) / 2,
        degree=observed[2],
        line_error=float(de + 1e-11),
    )


def predicted_generator(board, home, probe, k, blank, gain, s, t):
    free = np.array([v for v in range(board.n) if v not in (home, probe)])
    degree = board.degree[free]
    cost = degree * blank["diagonal"] + (6 - degree) * blank["degree"]
    for pin, sign in ((home, s), (probe, t)):
        cost += float(gain[str(sign)]) * np.isin(free, board.neighbours[pin])
    return (
        diags(cost) + blank["hopping"] * board.adjacency().tocsr()[free][:, free]
    ).tocsr()


def run():
    rng = np.random.Generator(np.random.PCG64(23042917).jumped(2))
    extra = np.random.Generator(np.random.PCG64(23042917).jumped(3))
    record = json.loads(gzip.decompress(measure.OUTPUT.read_bytes()))
    gain_rng = np.random.Generator(np.random.PCG64(23042917).jumped(5))
    calibration = {str(k): blank_clock(k, rng, gain_rng) for k in (8, 0)}
    panels = {}
    for level in (0, 1):
        board = federation(level)
        home = EMITTER_SEAM[level][0]
        dist = native.distances(board, home)
        for k in (8, 0):
            blank = calibration[str(k)]
            gain = blank["gains"]
            for src in record["boards"][f"{board.cells},{k}"]["shell_banks"]:
                shell, y = src["distance"], src["port"]
                rows = {}
                for s in (-1, 1):
                    for t in (-1, 1):
                        label = f"{s},{t}"
                        h, free = native.vacancy_operator(
                            board, k, {home: s + 1, y: t + 1}
                        )
                        n = len(free)
                        enc = bounds(h)
                        lines = eigvalsh(h.toarray(), check_finite=False)
                        h_fit = predicted_generator(
                            board, home, y, k, blank, gain, s, t
                        )
                        predicted = eigvalsh(h_fit.toarray(), check_finite=False)
                        exposures = dict(intact=trace_bank(lines, rng, SHOTS))
                        if shell in (3, dist.max()):
                            side = dist[free] <= native.distances(board, y)[free]
                            cut = np.concatenate(
                                [
                                    eigvalsh(
                                        h[mask][:, mask].toarray(), check_finite=False
                                    )
                                    for mask in (side, ~side)
                                ]
                            )
                            exposures["cut"] = trace_bank(cut, rng, SHOTS)
                        if k == 8 and shell == 3:
                            exposures["deep"] = trace_bank(lines, extra, 256 * SHOTS)
                        thermal_rows = {}
                        for beta in BETAS:
                            out = {
                                name: thermal(v, beta, n, enc)
                                for name, v in exposures.items()
                            }
                            # Six diagonal and six exchange contacts, one private term and two signed gains.
                            dh = 15 * blank["line_error"]
                            out.update(
                                forecast=float(np.exp(-beta * predicted).sum()),
                                no_hop=float(np.exp(-beta * h_fit.diagonal()).sum()),
                                forecast_error=float(
                                    n * beta * dh * np.exp(-beta * (enc[0] - 2 * dh))
                                ),
                            )
                            thermal_rows[str(beta)] = out
                        rows[label] = dict(
                            banks=exposures, gershgorin=enc, thermal=thermal_rows
                        )
                summaries = {}
                for beta in BETAS:
                    cols = {key: v["thermal"][str(beta)] for key, v in rows.items()}
                    v0 = cols["-1,-1"]
                    res = {}
                    for name in ("intact", "cut", "deep"):
                        if name in v0:
                            res[name] = (
                                -read.signed(
                                    {key: v[name][0] for key, v in cols.items()}
                                )
                                / beta
                            )
                            res[name + "_error"] = sum(
                                v[name][1] for v in cols.values()
                            ) / (4 * beta)
                    for name in ("forecast", "no_hop"):
                        res[name] = (
                            -read.signed({key: v[name] for key, v in cols.items()})
                            / beta
                        )
                    res["forecast_error"] = sum(
                        v["forecast_error"] for v in cols.values()
                    ) / (4 * beta)
                    summaries[str(beta)] = res
                panels[f"{board.cells},{k},{shell}"] = dict(
                    probe=y, signs=rows, readings=summaries
                )
    return dict(
        calibration=calibration,
        panels=panels,
        phase_probability_grid=2**40,
        shots=SHOTS,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.dumps(
        run(), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    if args.output:
        args.output.write_bytes(gzip.compress(data, mtime=0))
    elif gzip.decompress(OUTPUT.read_bytes()) != data:
        raise ValueError("Dilute charge receipt differs from replay")
    else:
        print("The full thermal trace counts and reconstructions replayed.")


if __name__ == "__main__":
    main()
