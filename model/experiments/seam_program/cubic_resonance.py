"""Find the cubic two-load threshold from native hopping and check it by finite-volume evolution."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.integrate import quad
from scipy.sparse.linalg import expm_multiply, spsolve
from scipy.special import ellipkm1, i0e

from model.candidate_laws.interval_echo.candidate import KAPPA
from .range_contact import native_check

HERE = Path(__file__).resolve().parent
ARMS = (
    "native",
    "exchange",
    "critical",
    "slow_detuning",
    "linear_detuning",
    "fast_detuning",
)


def green(mass):
    """Onsite resolvent of the infinite cubic relative Laplacian at E=-mass^2."""
    if not np.isfinite(mass) or mass < 0:
        raise ValueError("Expected a nonnegative inverse lattice length")

    # t=u^-2-1 regularizes the t^-3/2 large-time tail of the heat kernel.
    def integrand(u):
        t = 1 / (u * u) - 1
        return 2 * (i0e(2 * t) / u) ** 3 * np.exp(-mass * mass * t)

    breaks = [mass] if 0 < mass < 1 else None
    value, err = quad(
        integrand, 0, 1, points=breaks, epsabs=2e-13, epsrel=2e-13, limit=150
    )
    if err > 1e-11:
        raise ArithmeticError("Unresolved cubic return integral")
    return value


def coupling(arm, h, critical):
    if arm == "native":
        return KAPPA
    if arm == "exchange":
        return 0.0
    if arm == "critical":
        return critical
    power = {"slow_detuning": 0.5, "linear_detuning": 1.0, "fast_detuning": 2.0}[arm]
    return 1 / (1 / critical + 3 * h**power / np.pi)


def outgoing_green(energy):
    """Outgoing real-energy cubic Green function, integrating the third momentum exactly."""
    if not 0 < energy < 1:
        raise ValueError(
            "The enrolled real energies lie below the first lattice saddle"
        )

    def dos(s):
        return ellipkm1(((s - 4) / 4) ** 2) / (2 * np.pi**2)

    # At fixed two-dimensional energy s the remaining integral is
    # [(s-E-i0)(s+4-E-i0)]^(-1/2). Square-root substitutions remove its edge.
    def real(v):
        return 2 * dos(energy + v * v) / np.sqrt(4 + v * v)

    def imag(v):
        return 2 * dos(energy - v * v) / np.sqrt(4 - v * v)

    re, er = quad(
        real,
        0,
        np.sqrt(8 - energy),
        points=[np.sqrt(4 - energy)],
        epsabs=1e-12,
        epsrel=1e-12,
    )
    im, ei = quad(imag, 0, np.sqrt(energy), epsabs=1e-12, epsrel=1e-12)
    if er + ei > 1e-10:
        raise ArithmeticError("Unresolved outgoing lattice resolvent")
    return complex(re, im)


def scattering(h, kappa, energy, g0):
    delta = 1 + 2 * kappa
    # Keep the exchange-only zero and the near-threshold numerator resolved.
    numerator = 2 * kappa - delta * energy / 6
    denominator = delta / 6 - numerator * g0
    return h * numerator / (4 * np.pi * denominator)


def punctured_cube(n):
    indices = np.arange(n**3).reshape((n,) * 3)
    all_sites = indices.ravel()
    neighbors = np.stack(
        [
            np.roll(indices, sign, axis=axis).ravel()
            for axis in range(3)
            for sign in (-1, 1)
        ]
    )
    corners = np.minimum(np.indices((n,) * 3), n - np.indices((n,) * 3))
    distances = corners.max(axis=0).ravel()[1:] / n
    shell = np.any(neighbors[:, 1:] == 0, axis=0)
    # Relative jumps occur at rate one. Jumps onto the other load are absent.
    adjacency = sparse.coo_matrix(
        (np.ones(neighbors.size), (np.tile(all_sites, 6), neighbors.ravel())),
        shape=(n**3, n**3),
    ).tocsr()[1:, 1:]
    kinetic = sparse.diags(6.0 - shell.astype(float)) - adjacency
    return kinetic.tocsr(), shell, distances


def independent_shell_check(n, kappa, mass):
    kinetic, shell, _ = punctured_cube(n)
    hamiltonian = kinetic - sparse.diags(2 * kappa * shell)
    result = spsolve(
        hamiltonian + mass**2 * sparse.eye(len(shell), format="csr"),
        shell.astype(float),
    )
    angles = 2 * np.pi * np.arange(n) / n
    energies = 2 * (
        3
        - np.cos(angles[:, None, None])
        - np.cos(angles[None, :, None])
        - np.cos(angles[None, None, :])
    )
    g0 = np.mean(1 / (energies + mass**2))
    g1 = ((6 + mass**2) * g0 - 1) / 6
    expected = g1 / (g0 - (1 + 2 * kappa) * g1)
    return float(max(abs(result[shell] - expected)))


def packet(n, arm, critical):
    h = 1 / n
    kappa = coupling(arm, h, critical)
    kinetic, shell, distances = punctured_cube(n)
    H = kinetic - sparse.diags(2 * kappa * shell)
    initial = np.full(len(shell), 1 / np.sqrt(len(shell)))
    bins = np.minimum(9, np.floor(20 * distances).astype(int))
    times = (1 / 32, 1 / 8)
    rows = []
    for time in times:
        state = expm_multiply((-1j * time / h**2) * H, initial)
        law = np.bincount(bins, weights=abs(state) ** 2, minlength=10)
        if abs(law.sum() - 1) > 1e-9:
            raise ArithmeticError("Incomplete cubic packet read")
        rows.append(
            dict(
                time=time,
                separation_bins=law.tolist(),
                uniform_return=float(abs(np.vdot(initial, state)) ** 2),
            )
        )
    return dict(n=n, arm=arm, kappa_eff=kappa, readouts=rows)


def run():
    g0 = green(0.0)
    critical = 1 / (12 * (g0 - 1 / 6))
    rows = []
    for n in (16, 32, 64, 128, 256):
        h = 1 / n
        mus = (0.5, 1, 2, 4)
        kernels = [green(h * m) for m in mus]
        energies = [4 * np.sin(h * m / 2) ** 2 for m in mus]
        outgoing = [outgoing_green(e) for e in energies]
        for arm in ARMS:
            kappa = coupling(arm, h, critical)
            length = (
                None
                if arm == "critical"
                else (
                    0.0 if kappa == 0 else 3 * h / (np.pi * (1 / critical - 1 / kappa))
                )
            )
            rows.append(
                dict(
                    n=n,
                    arm=arm,
                    kappa_eff=kappa,
                    scattering_length=length,
                    threshold_pole=(arm == "critical"),
                    imaginary_momentum_response=[
                        scattering(h, kappa, -((h * m) ** 2), val)
                        for m, val in zip(mus, kernels, strict=True)
                    ],
                    real_momentum_response=[
                        [value.real, value.imag]
                        for value in (
                            scattering(h, kappa, e, val)
                            for e, val in zip(energies, outgoing, strict=True)
                        )
                    ],
                )
            )
    return dict(
        protocol_sha256=sha256(
            (HERE / "cubic_resonance_protocol.json").read_bytes()
        ).hexdigest(),
        green_at_zero=g0,
        kappa_star=critical,
        resonance_counterpulse=native_check(
            2 * (critical - KAPPA) / (64 * (2 * KAPPA + 1))
        ),
        kinetic_energy_ratio_at_double_momentum=[
            dict(n=n, ratio=4 * np.cos(0.5 / n) ** 2) for n in (16, 32, 64, 128, 256)
        ],
        physical_inverse_lengths=[0.5, 1, 2, 4],
        finite_shell_error=independent_shell_check(8, critical, 1),
        scattering=rows,
        packets=[
            packet(n, arm, critical)
            for n in (8, 16, 32)
            for arm in ("native", "exchange", "critical", "linear_detuning")
        ],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.parse_args().output.write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
