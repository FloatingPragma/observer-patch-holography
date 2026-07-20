#!/usr/bin/env python3
"""Exact W/Z resonance-convention and covariance map.

The public experimental parameters use a mass-dependent-width propagator

    D(s)^-1 = s - m_bar^2 + i s Gamma_bar / m_bar.

This module maps its pole to the energy-pole convention used by the package,

    s_p = (M - i Gamma/2)^2,  M > 0, Gamma >= 0.

It also retains the older ``sqrt(Re(s_p))`` mass so that the two inequivalent
uses of "pole mass" cannot be silently mixed.
"""

from __future__ import annotations

import cmath
import json
import math
from pathlib import Path
from typing import Iterable, Sequence


def running_width_to_energy_pole(mass: float, width: float) -> dict[str, float | list[float]]:
    """Map PDG/LEP running-width parameters to the exact complex pole.

    The returned ``mass`` and ``width`` obey
    ``s_pole == (mass - 0.5j * width)**2`` up to floating-point roundoff.
    """

    if not (math.isfinite(mass) and math.isfinite(width)):
        raise ValueError("mass and width must be finite")
    if mass <= 0.0 or width < 0.0:
        raise ValueError("require mass > 0 and width >= 0")

    ratio = width / mass
    s_pole = mass * mass / complex(1.0, ratio)
    energy_pole = cmath.sqrt(s_pole)
    if energy_pole.real < 0.0:
        energy_pole = -energy_pole
    pole_mass = energy_pole.real
    pole_width = -2.0 * energy_pole.imag

    old_mass = math.sqrt(s_pole.real)
    old_width = -s_pole.imag / old_mass
    return {
        "mass_GeV": pole_mass,
        "width_GeV": pole_width,
        "s_pole_GeV2": [s_pole.real, s_pole.imag],
        "sqrt_Re_s_pole_mass_GeV": old_mass,
        "legacy_constant_width_parameter_GeV": old_width,
    }


def convention_jacobian(mass: float, width: float) -> list[list[float]]:
    """Analytic Jacobian d(M, Gamma)/d(m_bar, Gamma_bar)."""

    if mass <= 0.0 or width <= 0.0:
        raise ValueError("Jacobian requires positive mass and width")
    r = width / mass
    u = 1.0 / math.sqrt(1.0 + r * r)
    f = math.sqrt((u + u * u) / 2.0)
    h = math.sqrt((u - u * u) / 2.0)
    du = -r * u**3
    df = du * (1.0 + 2.0 * u) / (4.0 * f)
    dh = du * (1.0 - 2.0 * u) / (4.0 * h)
    return [
        [f - r * df, df],
        [2.0 * (h - r * dh), 2.0 * dh],
    ]


def matmul(a: Sequence[Sequence[float]], b: Sequence[Sequence[float]]) -> list[list[float]]:
    if not a or not b or len(a[0]) != len(b):
        raise ValueError("incompatible matrix shapes")
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def transpose(a: Sequence[Sequence[float]]) -> list[list[float]]:
    return [list(row) for row in zip(*a)]


def propagate_covariance(jacobian: Sequence[Sequence[float]], covariance: Sequence[Sequence[float]]) -> list[list[float]]:
    """Return J C J^T after shape, symmetry, and PSD checks."""

    n = len(covariance)
    if n == 0 or any(len(row) != n for row in covariance):
        raise ValueError("covariance must be a non-empty square matrix")
    if any(len(row) != n for row in jacobian):
        raise ValueError("Jacobian input dimension does not match covariance")
    scale = max(1.0, max(abs(value) for row in covariance for value in row))
    for i in range(n):
        if covariance[i][i] < 0.0:
            raise ValueError("covariance has a negative diagonal")
        for j in range(n):
            if abs(covariance[i][j] - covariance[j][i]) > 1e-12 * scale:
                raise ValueError("covariance is not symmetric")
    # Cholesky-like PSD check with a tolerance that permits exact zero modes.
    lower = [[0.0] * n for _ in range(n)]
    tol = 1e-14 * scale
    for i in range(n):
        for j in range(i + 1):
            residual = covariance[i][j] - sum(lower[i][k] * lower[j][k] for k in range(j))
            if i == j:
                if residual < -tol:
                    raise ValueError("covariance is not positive semidefinite")
                lower[i][j] = math.sqrt(max(0.0, residual))
            elif lower[j][j] > math.sqrt(tol):
                lower[i][j] = residual / lower[j][j]
            elif abs(residual) > tol:
                raise ValueError("covariance is not positive semidefinite")
    return matmul(matmul(jacobian, covariance), transpose(jacobian))


def block_diagonal(blocks: Iterable[Sequence[Sequence[float]]]) -> list[list[float]]:
    blocks = list(blocks)
    size = sum(len(block) for block in blocks)
    out = [[0.0] * size for _ in range(size)]
    offset = 0
    for block in blocks:
        n = len(block)
        if any(len(row) != n for row in block):
            raise ValueError("each block must be square")
        for i in range(n):
            for j in range(n):
                out[offset + i][offset + j] = float(block[i][j])
        offset += n
    return out


def make_target_receipt(packet: dict) -> dict:
    """Build the frozen W/Z target receipt from a validated input packet."""

    target = packet["experimental_target"]
    order = target["input_order"]
    expected_order = ["W_mass_running_width", "W_width_running_width", "Z_mass_running_width", "Z_width_running_width"]
    if order != expected_order:
        raise ValueError(f"unexpected target order: {order!r}")
    values = [float(value) for value in target["central_GeV"]]
    covariance = [[float(value) for value in row] for row in target["covariance_GeV2"]]

    w = running_width_to_energy_pole(values[0], values[1])
    z = running_width_to_energy_pole(values[2], values[3])
    jw = convention_jacobian(values[0], values[1])
    jz = convention_jacobian(values[2], values[3])
    jacobian = block_diagonal([jw, jz])
    pole_covariance = propagate_covariance(jacobian, covariance)
    pole_order = ["W_mass_energy_pole", "W_width_energy_pole", "Z_mass_energy_pole", "Z_width_energy_pole"]
    pole_central = [w["mass_GeV"], w["width_GeV"], z["mass_GeV"], z["width_GeV"]]

    return {
        "schema": "physical_wz_target_receipt_v1",
        "release": target["release"],
        "input_convention": "mass_dependent_width_Breit_Wigner",
        "input_order": order,
        "input_central_GeV": values,
        "input_covariance_GeV2": covariance,
        "covariance_status": target["covariance_status"],
        "output_convention": "s_pole=(M-i*Gamma/2)^2",
        "output_order": pole_order,
        "output_central_GeV": pole_central,
        "output_covariance_GeV2": pole_covariance,
        "output_sigma_GeV": [math.sqrt(pole_covariance[i][i]) for i in range(4)],
        "jacobian": jacobian,
        "W": w,
        "Z": z,
        "joint_chi2_licensed": target["joint_chi2_licensed"],
        "sources": target["sources"],
    }


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = make_target_receipt(json.loads(args.packet.read_text(encoding="utf-8")))
    encoded = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
