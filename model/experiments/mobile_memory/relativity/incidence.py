"""Local incidence transport in a declared Dirac extension of mobile memory."""

from collections import deque
import math

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import expm_multiply

from model.executable_model import carrier
from ..mobility import PORTS_PER_CELL, home_seam


SPEED = 1.0
CHARGE_MASS = 1.0
CAGE_MASS = 2 * math.sqrt(6)
TOLERANCE = 2e-10
WAITS = (0.0, 1.0, 2.0, 4.0)


def cell_incidence(car, weights=None):
    """Orient each actual external seam from its lower port to its higher port.

    A supplied mapping gives nonnegative link weights, with omitted links blocked.
    Each incidence column takes the square root of its link weight.
    """
    seams = tuple(tuple(sorted(seam)) for seam in car.external)
    if len(set(seams)) != len(seams):
        raise ValueError("the carrier has duplicate external seams")
    if weights is not None and not set(weights) <= set(seams):
        raise ValueError("a weight addresses a seam outside the carrier")
    left, right = np.asarray(seams, dtype=int).reshape(-1, 2).T // PORTS_PER_CELL
    if np.any(left == right) or np.any(left < 0) or np.any(right >= car.cells):
        raise ValueError("an external seam must connect two cells")
    values = np.ones(len(seams)) if weights is None else np.array(
        [weights.get(seam, 0.0) for seam in seams], dtype=float
    )
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("link weights must be finite and nonnegative")
    roots = np.sqrt(values)
    incidence = sparse.csr_matrix(
        (np.r_[roots, -roots], (np.r_[left, right], np.tile(np.arange(len(seams)), 2))),
        shape=(car.cells, len(seams)),
    )
    incidence.eliminate_zeros()
    return incidence, seams


def _check_parameters(mass, speed):
    if not math.isfinite(mass) or mass <= 0 or not math.isfinite(speed) or speed <= 0:
        raise ValueError("mass and speed must be finite and positive")


def dirac_operator(incidence, mass, speed=SPEED):
    """Build D on the actual cell and seam registers."""
    _check_parameters(mass, speed)
    cells, seams = incidence.shape
    return sparse.bmat(
        [[mass * sparse.eye(cells), speed * incidence],
         [speed * incidence.T, -mass * sparse.eye(seams)]],
        format="csr",
    )


def laplacian_modes(incidence):
    """Return the cell eigenvalues and corresponding orthonormal eigenvectors."""
    values, vectors = np.linalg.eigh((incidence @ incidence.T).toarray())
    if values.min(initial=0) < -TOLERANCE:
        raise ValueError("the incidence Laplacian must be positive semidefinite")
    return np.maximum(values, 0.0), vectors


def _shell(eigenvalues, mass, speed):
    _check_parameters(mass, speed)
    energies = np.hypot(mass, speed * np.sqrt(eigenvalues))
    return energies, np.sqrt((energies + mass) / (2 * energies))


def positive_mode(incidence, mass, eigenvalue, vertex_mode, speed=SPEED):
    """Dress a normalized cell eigenmode in the positive Dirac band."""
    vertex_mode = np.asarray(vertex_mode)
    if vertex_mode.shape != (incidence.shape[0],) or not math.isfinite(eigenvalue) or eigenvalue < 0:
        raise ValueError("the mode must give one finite Laplacian eigenpair")
    np.testing.assert_allclose(np.linalg.norm(vertex_mode), 1.0, atol=TOLERANCE, rtol=0)
    np.testing.assert_allclose(incidence @ (incidence.T @ vertex_mode),
                               eigenvalue * vertex_mode, atol=TOLERANCE, rtol=0)
    energy, amplitude = _shell(eigenvalue, mass, speed)
    cells = amplitude * vertex_mode
    seams = speed * (incidence.T @ cells) / (energy + mass)
    return float(energy), np.r_[cells, seams]


def positive_embedding(incidence, mass, speed=SPEED):
    """Isometrically dress an arbitrary vertex seed in the positive band."""
    values, vectors = laplacian_modes(incidence)
    energies, amplitudes = _shell(values, mass, speed)
    vertex_modes = vectors * amplitudes
    seam_modes = speed * (incidence.T @ vertex_modes) / (energies + mass)
    return np.vstack((vertex_modes, seam_modes)) @ vectors.T


def _largest_entry(matrix):
    return float(np.max(np.abs(matrix.data), initial=0.0))


def _spectrum_certificate(incidence, mass, speed):
    operator = dirac_operator(incidence, mass, speed)
    cells, seams = incidence.shape
    laplacian = incidence @ incidence.T
    squared = sparse.block_diag(
        (mass ** 2 * sparse.eye(cells) + speed ** 2 * laplacian,
         mass ** 2 * sparse.eye(seams) + speed ** 2 * (incidence.T @ incidence)),
        format="csr",
    )
    error = _largest_entry(operator @ operator - squared)
    if error > TOLERANCE:
        raise ValueError("the local Dirac square differs from the incidence shell")
    components = int(connected_components(laplacian, directed=False, return_labels=False))
    # Null cell vectors have energy +m. Null edge vectors have energy -m.
    return {
        "mass": mass,
        "speed_coefficient": speed,
        "components": components,
        "positive_band_size": cells,
        "unpaired_positive_cell_modes": components,
        "unpaired_negative_edge_modes": seams - cells + components,
        "square_residual": error,
    }


def addressed_records(car, mass, seam=None, speed=SPEED):
    """Compare the same positive-band preparation with one link written or blocked."""
    seam = home_seam(car) if seam is None else tuple(seam)
    written, _ = cell_incidence(car, {tuple(sorted(seam)): 1.0})
    blocked, _ = cell_incidence(car, {})
    sender, receiver = (port // PORTS_PER_CELL for port in seam)
    prep = positive_embedding(written, mass, speed)[:, sender]
    np.testing.assert_allclose(np.linalg.norm(prep), 1.0, atol=TOLERANCE, rtol=0)

    def occupation(incidence):
        operator = dirac_operator(incidence, mass, speed)
        result = []
        for wait in WAITS:
            step = -1j * wait * operator
            state = expm_multiply(step, prep, traceA=step.diagonal().sum())
            result.append(abs(state) ** 2)
        result = np.asarray(result)
        np.testing.assert_allclose(result.sum(axis=1), 1.0, atol=TOLERANCE, rtol=0)
        return result

    moving = occupation(written)
    parked = occupation(blocked)
    # On this link the symmetric mode has lambda=0 and the antisymmetric mode lambda=2.
    energy, vertex_weight = _shell(2.0, mass, speed)
    analytic = (1 + vertex_weight ** 2
                - 2 * vertex_weight * np.cos((energy - mass) * np.asarray(WAITS))) / 4
    np.testing.assert_allclose(moving[:, receiver], analytic, atol=TOLERANCE, rtol=0)
    np.testing.assert_allclose(parked, np.repeat(parked[:1], len(WAITS), axis=0),
                               atol=TOLERANCE, rtol=0)
    return {
        "mass": mass,
        "seam": list(seam),
        "sender_cell": sender,
        "receiver_cell": receiver,
        "waits": list(WAITS),
        "written": moving.tolist(),
        "blocked": parked.tolist(),
        "written_receiver_gain": float(moving[-1, receiver] - moving[0, receiver]),
    }


def carrier_ring(car):
    """Find a deterministic shortest loop through one actual removed seam."""
    neighbors = [set() for _ in range(car.cells)]
    for i, j in car.external:
        a, b = i // PORTS_PER_CELL, j // PORTS_PER_CELL
        neighbors[a].add(b)
        neighbors[b].add(a)
    start, stop = (port // PORTS_PER_CELL for port in home_seam(car))
    previous = {start: None}
    queue = deque([start])
    while queue and stop not in previous:
        a = queue.popleft()
        for b in sorted(neighbors[a]):
            if b not in previous and (a, b) not in ((start, stop), (stop, start)):
                previous[b] = a
                queue.append(b)
    if stop not in previous:
        raise ValueError("the addressed carrier seam is not part of a loop")
    path = [stop]
    while path[-1] != start:
        path.append(previous[path[-1]])
    return tuple(reversed(path))


def ring_incidence(car, spacing=1.0):
    """Use equal spacing on one supplied cycle of actual seams."""
    if not math.isfinite(spacing) or spacing <= 0:
        raise ValueError("the declared spacing must be finite and positive")
    ring = carrier_ring(car)
    pairs = {tuple(sorted(pair)) for pair in zip(ring, ring[1:] + ring[:1], strict=True)}
    weights = {tuple(sorted(seam)): 1 / spacing ** 2 for seam in car.external
               if tuple(sorted(port // PORTS_PER_CELL for port in seam)) in pairs}
    if len(weights) != len(ring):
        raise ValueError("the loop requires exactly one actual seam per ring link")
    incidence, seams = cell_incidence(car, weights)
    return ring, incidence, seams


def ring_certificate(car, spacing=1.0, speed=SPEED):
    """Check a Fourier shell and its finite-spacing clock-rate defect on a real loop."""
    ring, incidence, _ = ring_incidence(car, spacing)
    rows = []
    for mass in (CHARGE_MASS, CAGE_MASS):
        for mode in range(len(ring)):
            wave_number = math.tau * mode / len(ring)
            orbital = np.zeros(car.cells, complex)
            orbital[list(ring)] = np.exp(1j * wave_number * np.arange(len(ring))) / math.sqrt(len(ring))
            eigenvalue = 4 * math.sin(wave_number / 2) ** 2 / spacing ** 2
            energy, state = positive_mode(incidence, mass, eigenvalue, orbital, speed)
            operator = dirac_operator(incidence, mass, speed)
            np.testing.assert_allclose(operator @ state, energy * state, atol=TOLERANCE, rtol=0)
            # Physical velocity would also require a physical spacing. This is the derivative
            # dE/dp for the declared coordinate p=k/spacing on the periodic cut.
            group_speed = speed ** 2 * math.sin(wave_number) / (spacing * energy)
            rate = mass / energy
            defect = 4 * speed ** 2 * math.sin(wave_number / 2) ** 4 / (spacing * energy) ** 2
            np.testing.assert_allclose(rate ** 2 + (group_speed / speed) ** 2 + defect,
                                       1.0, atol=TOLERANCE, rtol=0)
            rows.append({
                "mass": mass,
                "mode": mode,
                "energy": energy,
                "group_speed": group_speed,
                "mass_derivative_rate": rate,
                "finite_spacing_defect": defect,
            })
    return {"cell_ring": list(ring), "declared_spacing": spacing, "records": rows}


def run():
    """Verify the incidence construction and retain all addressed detector outcomes."""
    carriers = [carrier.federation(level) for level in range(3)]
    spectra = []
    for car in carriers:
        incidence, _ = cell_incidence(car)
        spectra.append({
            "cells": car.cells,
            "seams": incidence.shape[1],
            "charge": _spectrum_certificate(incidence, CHARGE_MASS, SPEED),
            "cage": _spectrum_certificate(incidence, CAGE_MASS, SPEED),
        })
    return {
        "all_written": spectra,
        "addressed": {
            "charge": addressed_records(carriers[0], CHARGE_MASS),
            "cage": addressed_records(carriers[0], CAGE_MASS),
        },
        "ring": ring_certificate(carriers[0]),
    }
