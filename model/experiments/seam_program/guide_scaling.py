"""Low-momentum scattering and fixed-time two-load records on an addressed binary ring."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

from model.candidate_laws.interval_echo.candidate import KAPPA

HERE = Path(__file__).resolve().parent
ARMS = ("native", "exchange", "slow_detuning", "linear_detuning", "fast_detuning")


def contrast(arm, h):
    if arm == "native":
        return KAPPA
    if arm == "exchange":
        return 0.0
    exponent = {"slow_detuning": 0.5, "linear_detuning": 1.0, "fast_detuning": 2.0}[arm]
    return -(h**exponent) / 4


def exchange_amplitude(p1, p2, h, kappa):
    """Coefficient multiplying the exchanged Bethe plane wave on x<y."""
    delta = 1 + 2 * kappa
    common = np.expm1(1j * h * (p1 + p2)) - 4 * kappa
    numerator = common - 2 * delta * np.expm1(1j * h * p2)
    denominator = common - 2 * delta * np.expm1(1j * h * p1)
    if denominator == 0:
        raise ValueError("Degenerate scattering momenta")
    return -numerator / denominator


def ring(n, kappa):
    if n < 8 or n % 2:
        raise ValueError("Expected an even ring of at least eight ports")
    half = n // 2
    diag = np.full(half, 2.0)
    diag[0] = 1 - 2 * kappa
    off = np.full(half - 1, -1.0)
    off[-1] = -np.sqrt(2)
    energies, vectors = eigh_tridiagonal(diag, off)
    initial = np.ones(half)
    initial[-1] = 1 / np.sqrt(2)
    initial /= np.linalg.norm(initial)
    amplitudes = vectors.T @ initial
    rows = []
    for physical_time in (1 / 32, 1 / 8):
        state = vectors @ (amplitudes * np.exp(-1j * energies * n**2 * physical_time))
        law = abs(state) ** 2
        rows.append(
            dict(
                time=physical_time,
                separation_probabilities=law.tolist(),
                quarter_or_closer=float(law[: n // 4].sum()),
                uniform_return=float(abs(np.vdot(initial, state)) ** 2),
            )
        )
    below = energies < -1e-10
    return dict(
        lowest_native_energy=float(energies[0]),
        bound_weight=float(np.sum(abs(amplitudes[below]) ** 2)),
        readouts=rows,
    )


def scan():
    rows = []
    for n in (16, 32, 64, 128, 256):
        h = 1 / n
        for arm in ARMS:
            kappa = contrast(arm, h)
            data = []
            for relative in (0.5, 1.0, 2.0, 4.0):
                p1, p2 = (0.5 + relative) / 2, (0.5 - relative) / 2
                result = exchange_amplitude(p1, p2, h, kappa)
                limit = {
                    "native": -1.0,
                    "exchange": 1.0,
                    "slow_detuning": -1.0,
                    "linear_detuning": (relative - 1j) / (relative + 1j),
                    "fast_detuning": 1.0,
                }[arm]
                data.append(
                    dict(
                        relative_momentum=relative,
                        real=float(result.real),
                        imag=float(result.imag),
                        limit_error=float(abs(result - limit)),
                    )
                )
            rows.append(
                dict(n=n, arm=arm, kappa_eff=kappa, exchange=data, **ring(n, kappa))
            )
    return dict(
        protocol_sha256=sha256(
            (HERE / "guide_scaling_protocol.json").read_bytes()
        ).hexdigest(),
        spatial_dimension=1,
        dynamical_exponent=2,
        rows=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(
        json.dumps(scan(), indent=2, allow_nan=False) + "\n"
    )
