"""Independent sparse propagation and remote quadrature counts at all appointments."""

import argparse
import json

import numpy as np
from mpmath import iv
from scipy.sparse.linalg import expm_multiply

from . import measure


def wave_columns(board, route, q, odd):
    h = measure.action.native.generator(board, route, q, odd)
    ih = -0.5j * h
    x = np.zeros((ih.shape[0], 8))
    x[[4 * p + j for p in measure.PORTS[:2] for j in range(4)]] = np.eye(8)
    w = expm_multiply(ih, x, start=0.5, stop=6, num=12, traceA=ih.diagonal().sum())
    rows = [4 * p + j for p in measure.PORTS[2:] for j in range(4)]
    return w[:, rows].reshape(12, 2, 4, 2, 4)


def as_iv(v):
    return iv.mpc(
        measure.action.front.full.ri(v.as_real_imag()[0]),
        measure.action.front.full.ri(v.as_real_imag()[1]),
    )


def group_delay(fringes, times, weight, angles):
    bank = {}
    for name, samples in fringes.items():
        phases = [weight**n * x for n, x in enumerate(samples, 1)]
        bank[name] = [sum(t**k * a for t, a in zip(times, phases)) for k in (0, 1, 2)]
    f, g, g2 = bank["equilibrium"]
    mu = g / f
    omega = 2 * iv.atan2(weight.imag, weight.real)
    dispersion = -omega * (g2 / f - mu**2).imag
    out = dict(delay=str(mu.real), omega_delay_dispersion=str(dispersion))
    raw = {}
    slopes = {}
    for name, angle in angles.items():
        df, dg = (
            bank[name][j]
            - bank["equilibrium"][j]
            + angle * (bank["plus"][j] - bank["minus"][j]) / 2
            for j in (0, 1)
        )
        slopes[name] = ((dg - mu * df) / f).real
        out[name] = str(slopes[name])
        raw[name] = ((bank[name][1] - g - mu * (bank[name][0] - f)) / f).real
    h = raw["pairing"] - raw["ordinary"] * angles["pairing"] / angles["ordinary"]
    out["pairing_minus_ordinary_over_volume_ratio"] = str(h)
    return out, [iv.mpf(1), mu.real, dispersion, slopes["ordinary"], h]


SHOTS = 2**58
SEED = 344567


def detector_probabilities():
    native = measure.action.native
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    states = measure.action.preparations()
    output = {}
    for q in (-1, 1):
        waves = np.array([wave_columns(board, route, q, odd) for odd in (False, True)])
        laws = {}
        for n, chi in states.items():
            table = []
            for it in range(12):
                phases = []
                for z in (1, 1j, -1, -1j):
                    arm = (
                        waves[:, it, :, :, 0, :] + z * waves[:, it, :, :, 1, :]
                    ) / np.sqrt(2)
                    p = np.einsum("bryi,bij,bryj->r", arm, chi, arm.conj()).real
                    assert p.min() >= 0 and p.sum() < 1
                    phases.append([*p, 1 - p.sum()])
                table.append(phases)
            laws[n] = table
        output[str(q)] = laws
    return output


def acquire_counts(probabilities):
    rng = np.random.default_rng(SEED)
    raw = {}
    for charge, settings in probabilities.items():
        raw[charge] = {}
        for preparation, table in settings.items():
            raw[charge][preparation] = [
                [rng.multinomial(SHOTS, law).tolist() for law in phases] for phases in table
            ]
    return raw


def analyze_counts(raw):
    """Reconstruct the certified intervals from integer detector records alone."""
    iv.dps = 30
    source_rows = json.loads(
        (measure.HERE.parent.parent / "signed_transfers.json").read_text()
    )["tests"]
    angles = {
        ("ordinary" if x["source"] == "ordinary" else "pairing"): iv.mpf(
            x["angle_per_unit"]
        )
        for x in source_rows
    }
    # A binary marginal for each of two exclusive receiving ports.
    preparations = len(measure.action.preparations())
    L = iv.ln(2 * 2 * 2 * preparations * 12 * 4 / iv.mpf("0.01"))
    times = [iv.mpf(n) / 2 for n in range(1, 13)]
    result = {}
    frame_rows = {}
    for q in (-1, 1):
        counts = raw[str(q)]
        fringes = {0: {}, 1: {}}
        for n, table in counts.items():
            for j in (0, 1):
                samples = []
                for phases in table:
                    means = []
                    for phase in phases:
                        x = phase[j]
                        e = (iv.sqrt(2 * x * L) + 2 * L) / SHOTS
                        means.append(iv.mpf(x) / SHOTS + iv.mpf([-1, 1]) * e)
                    samples.append(
                        iv.mpc((means[0] - means[2]) / 2, (means[1] - means[3]) / 2)
                    )
                fringes[j][n] = samples
        read = {}
        for j, port in enumerate(measure.PORTS[2:]):
            read[str(port)] = {}
            for name, w in measure.PHASORS.items():
                read[str(port)][name], frame_rows[(str(q), str(port), name)] = (
                    group_delay(fringes[j], times, as_iv(measure.sp.sympify(w)), angles)
                )
        result[str(q)] = dict(counts=counts, nominal_intervals_99pct=read)
    return dict(
        nominal_affine_clock_witnesses_99pct=measure.frame_witnesses(frame_rows),
        shots_per_setting=SHOTS,
        seed=SEED,
        phase_order=["0", "pi/2", "pi", "3pi/2"],
        outcome_order=[*measure.PORTS[2:], "elsewhere"],
        charges=result,
    )


def run():
    return analyze_counts(acquire_counts(detector_probabilities()))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "counts.json").write_text(out)
    else:
        print(out, end="")
