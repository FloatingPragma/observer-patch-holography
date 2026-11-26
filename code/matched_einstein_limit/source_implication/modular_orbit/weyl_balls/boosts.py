"""Read tilted balls by unitary operations on one finite Cauchy record."""

import json
from math import pi
from pathlib import Path

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import expm_multiply

from .measure import DIPOLE, PAULI, entropy, make_board, sites_in_ball


def hopping(side, spacing):
    shift = sparse.diags([np.ones(side - 1), [-1]], [1, 1 - side], shape=(side, side), format="csr")
    derivative = (shift - shift.T) / (2j * spacing)
    identity = sparse.eye(side, format="csr")
    terms = []
    for axis, spin in enumerate(PAULI):
        factors = [identity] * 3
        factors[axis] = derivative
        lattice = sparse.kron(sparse.kron(*factors[:2]), factors[2])
        terms.append(sparse.kron(lattice, spin, format="csr"))
    return sum(terms).tocsr()


def smooth_step(value):
    """Rise smoothly from zero at zero to one at one."""
    value = np.asarray(value)
    result = np.zeros_like(value, dtype=float)
    result[value >= 1] = 1
    mask = (value > 0) & (value < 1)
    t = value[mask]
    left, right = np.exp(-1 / t), np.exp(-1 / (1 - t))
    result[mask] = left / (left + right)
    return result


def frame_generator(board, action, direction):
    side = board.side
    coordinates = board.spacing * (np.indices((side,) * 3).reshape(3, -1).T - side // 2)
    length = board.spacing * side
    taper = np.prod(smooth_step((0.5 * length - np.abs(coordinates)) / (0.25 * length)), axis=1)
    position = np.repeat((coordinates @ np.asarray(direction)) * taper, 2)
    return (action.multiply(position[:, None]) + action.multiply(position[None, :])).tocsr() / 2


def covariance_on_columns(board, columns):
    side = board.side
    position_sum = np.indices((side,) * 3).sum(axis=0)
    untwist = np.exp(-1j * pi * position_sum / side)
    twisted_columns = columns.reshape((side,) * 3 + (2, -1)) * untwist[..., None, None]
    fourier = np.fft.fftn(twisted_columns, axes=(0, 1, 2), norm="ortho")
    flat = fourier.reshape(-1, columns.shape[1])

    def compress(kernels):
        applied = np.zeros_like(fourier)
        for kernel, spin in zip(kernels, [*PAULI, np.eye(2)]):
            coefficient = np.fft.fftn(kernel)
            applied += coefficient[..., None, None] * np.einsum("ab,...bj->...aj", spin, fourier)
        result = flat.conj().T @ applied.reshape(flat.shape)
        return (result + result.conj().T) / 2

    return 0.5 * np.eye(columns.shape[1]) + compress(board.vacuum_kernels), compress(board.source_kernels)


def tilted_reading(board, action, radius, direction, rapidity):
    sites = sites_in_ball(board.spacing, radius, (0, 0, 0))
    offsets = np.concatenate((np.eye(3, dtype=int), -np.eye(3, dtype=int)))
    neighborhood = np.unique(np.concatenate([sites, (sites[:, None, :] + offsets).reshape(-1, 3)]), axis=0)
    site_ids = np.ravel_multi_index((neighborhood + board.side // 2).T, (board.side,) * 3)
    modes = (2 * site_ids[:, None] + [0, 1]).ravel()
    columns = np.zeros((2 * board.side**3, len(modes)), dtype=complex)
    columns[modes, np.arange(len(modes))] = 1
    if rapidity:
        generator = frame_generator(board, action, direction)
        columns = expm_multiply(-1j * rapidity * generator, columns, traceA=0)
    reference, difference = covariance_on_columns(board, columns)
    rr = np.sum((board.spacing * neighborhood) ** 2, axis=1)
    mask = np.repeat(rr < radius**2, 2)
    cov, delta = reference[np.ix_(mask, mask)], difference[np.ix_(mask, mask)]
    occupations, base = linalg.eigh(cov, check_finite=False)
    log_odds = np.log1p(-occupations) - np.log(occupations)
    modular = np.dot(log_odds, np.einsum("ij,ij->j", base.conj(), delta @ base).real)
    shift = entropy(linalg.eigvalsh(cov + delta, check_finite=False)) - entropy(occupations)
    weights = np.repeat(np.maximum(radius**2 - rr, 0) / (2 * radius), 2)
    local_action = action[modes][:, modes].toarray()
    energy = 2 * pi * np.einsum("i,ij,ji->", weights, local_action, difference).real
    radiation = (4 * np.cosh(rapidity) ** 2 - 1) * board.energy / 3
    radiation += np.sinh(2 * rapidity) * np.dot(direction, board.momentum_density)
    return {
        "direction": list(direction), "rapidity": rapidity, "sites": len(sites),
        "entropy_change": float(shift), "modular_change": float(modular),
        "relative_entropy": float(modular - shift), "weighted_hopping_energy": float(energy),
        "relative_entropy_upper_bound": float(np.vdot(delta, delta).real / (2 * board.floor * (1 - board.floor))),
        "continuum_radiation_prediction": float(8 * pi**2 * radius**4 * radiation / 15),
    }


def audit():
    side, length, cutoff, radius, rapidity = 12, 6.0, 2.5, 0.76, 0.32
    board = make_board(side, length, cutoff, 0.001 * cutoff * radius)
    action = hopping(side, board.spacing)
    axes = np.eye(3)
    directions = list(axes) + list(-axes) + [(axes[i] + axes[j]) / np.sqrt(2) for i, j in [(0, 1), (1, 2), (2, 0)]]
    rows = [tilted_reading(board, action, radius, (0., 0., 0.), 0.)]
    rows.extend(tilted_reading(board, action, radius, direction, rapidity) for direction in directions)
    return {"grid_side": side, "box_length": length, "momentum_cutoff": cutoff, "source_dipole": DIPOLE, "radius": radius, "reference_floor": board.floor, "rest_hopping_energy_density": board.energy, "rest_momentum_density": board.momentum_density, "frames": rows}


if __name__ == "__main__":
    Path(__file__).with_name("boosts.json").write_text(json.dumps(audit(), indent=2) + "\n")
