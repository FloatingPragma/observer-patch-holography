"""One stationary lattice Weyl reference, one excitation, several regional readings."""

from dataclasses import dataclass
import json
from math import pi
from pathlib import Path

import numpy as np
from scipy import linalg
from scipy.special import xlogy


PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
DIPOLE = (0.6, -0.3, 0.2)


@dataclass
class Board:
    side: int
    spacing: float
    floor: float
    energy: float
    momentum_density: list
    vacuum_kernels: list
    source_kernels: list


def make_board(side, length, cutoff, floor, dipole=DIPOLE):
    """Anti-periodic hopping and two stationary occupation matrices in its energy basis."""
    if side % 2 or not 0 < floor < 0.5 or min(length, cutoff) <= 0:
        raise ValueError("use an even grid, positive lengths and a faithful occupation floor")
    if np.shape(dipole) != (3,) or np.linalg.norm(dipole) > 1:
        raise ValueError("the angular occupation requires a dipole of norm at most one")
    spacing = length / side
    momenta = 2 * pi * (side * np.fft.fftfreq(side) + 0.5) / length
    p = np.meshgrid(momenta, momenta, momenta, indexing="ij", sparse=True)
    d = [np.sin(spacing * component) / spacing for component in p]
    energy = np.sqrt(sum(component**2 for component in d))
    momentum_squared = sum(component**2 for component in p)
    occupation = 0.25 * np.maximum(0, 1 - momentum_squared / cutoff**2) ** 2
    difference = (1 - 2 * floor) * occupation
    odd_part = difference * sum(weight * value for weight, value in zip(dipole, p)) / cutoff
    vacuum_kernels = [np.fft.ifftn(-(0.5 - floor) * component / energy) for component in d]
    source_kernels = [np.fft.ifftn(difference * component / energy) for component in d]
    source_kernels.append(np.fft.ifftn(odd_part))
    density = 2 * np.sum(difference * energy) / length**3
    momentum_density = [float(2 * np.sum(odd_part * component) / length**3) for component in p]
    return Board(side, spacing, floor, float(density), momentum_density, vacuum_kernels, source_kernels)


def sites_in_ball(spacing, radius, center):
    center = np.asarray(center, dtype=float)
    corners = [np.arange(np.floor((x - radius) / spacing), np.ceil((x + radius) / spacing) + 1, dtype=int) for x in center]
    sites = np.stack(np.meshgrid(*corners, indexing="ij"), axis=-1).reshape(-1, 3)
    return sites[np.sum((spacing * sites - center) ** 2, axis=1) < radius**2]


def regional_covariances(board, sites):
    displacement = sites[:, None, :] - sites[None, :, :]
    index = tuple((displacement % board.side).transpose(2, 0, 1))
    twist = np.exp(1j * pi * displacement.sum(axis=-1) / board.side)

    def spin_matrix(kernels):
        terms = [np.kron(twist * kernel[index], sigma) for kernel, sigma in zip(kernels, [*PAULI, np.eye(2)])]
        return sum(terms)

    reference = 0.5 * np.eye(2 * len(sites)) + spin_matrix(board.vacuum_kernels)
    return reference, spin_matrix(board.source_kernels)


def entropy(eigenvalues):
    return -np.sum(xlogy(eigenvalues, eigenvalues) + xlogy(1 - eigenvalues, 1 - eigenvalues))


def read_ball(board, radius, center=(0, 0, 0)):
    sites = sites_in_ball(board.spacing, radius, center)
    reference, difference = regional_covariances(board, sites)
    occupations, basis = linalg.eigh(reference, check_finite=False, driver="evd")
    excited = linalg.eigvalsh(reference + difference, check_finite=False, driver="evd")
    if min(occupations[0], excited[0]) < board.floor - 1e-12 or max(occupations[-1], excited[-1]) > 1 - board.floor + 1e-12:
        raise ArithmeticError("regional eigenvalue outside the physical global occupation bounds")
    log_odds = np.log1p(-occupations) - np.log(occupations)
    modular = np.dot(log_odds, np.einsum("ij,ij->j", basis.conj(), difference @ basis).real)
    entropy_shift = entropy(excited) - entropy(occupations)
    product_spectrum = 0.5 + linalg.eigvalsh(difference, check_finite=False, driver="evd")
    product_shift = entropy(product_spectrum) - len(product_spectrum) * np.log(2)
    hs_squared = float(np.vdot(difference, difference).real)
    weights = (radius**2 - np.sum((board.spacing * sites - center) ** 2, axis=1)) / (2 * radius)
    measured_energy = 2 * pi * board.energy * board.spacing**3 * weights.sum()
    continuum_energy = 8 * pi**2 * board.energy * radius**4 / 15
    return {
        "radius": radius, "center": list(center), "sites": len(sites),
        "entropy_change": float(entropy_shift), "modular_change": float(modular),
        "product_reference_entropy_change": float(product_shift), "product_reference_modular_change": 0.0,
        "relative_entropy": float(modular - entropy_shift),
        "relative_entropy_upper_bound": hs_squared / (2 * board.floor * (1 - board.floor)),
        "weighted_hopping_energy": float(measured_energy), "ball_energy_from_density": float(continuum_energy),
        "entropy_over_ball_energy": float(entropy_shift / continuum_energy),
        "entropy_over_weighted_hopping": float(entropy_shift / measured_energy),
        "largest_onsite_source_change": float(max(abs(difference.reshape(len(sites), 2, len(sites), 2)[j, :, j]).max() for j in range(len(sites)))),
    }


def audit():
    # One physical box length, one fixed low-momentum occupation and one hopping clock.
    # The two balls at each resolution share the exact global reference and source.
    rows = []
    cutoff, length = 2.5, 8.0
    for side, radius in [(32, 0.75), (48, 0.60), (64, 0.50)]:
        floor = 0.001 * cutoff * radius
        board = make_board(side, length, cutoff, floor)
        measurements = [read_ball(board, radius), read_ball(board, radius * 0.8, (radius * 0.10, 0, 0))]
        rows.append({"grid_side": side, "spacing": board.spacing, "floor": floor, "energy_density": board.energy, "balls": measurements})
    return {"model": "stationary charge-symmetric Weyl comparison sector", "box_length": length, "momentum_cutoff": cutoff, "source_dipole": DIPOLE, "runs": rows}


if __name__ == "__main__":
    Path(__file__).with_name("results.json").write_text(json.dumps(audit(), indent=2) + "\n")
