"""Optional rotation-selected current and complete single-cell detector, using the shared native incidence."""
from __future__ import annotations

from itertools import combinations
import json

import networkx as nx
import numpy as np
from scipy.linalg import block_diag, expm

from model.executable_model import cell, rotations
from ..contrast_memory_spherical_v1 import port_response as port
from ..contrast_memory_spherical_v1.sectors import PairSector
from .candidate import KAPPA

TOL = 1e-9
STRENGTH = 0.41
SOURCE_PORT = 1


def _geometry(car):
    graph = nx.Graph(car.seams)
    planar, embedding = nx.check_planarity(graph)
    if not planar:
        raise ValueError('the service requires the native oriented cell')
    antipodes, turns = [], []
    for p in range(car.n):
        distance = nx.single_source_shortest_path_length(graph, p)
        far, = (v for v, d in distance.items() if d == max(distance.values()))
        antipodes.append(far)
        ring = list(embedding.neighbors_cw_order(p))
        j = np.zeros((car.n, car.n), complex)
        for a, b in zip(ring, ring[1:] + ring[:1]):
            j[b, a], j[a, b] = 1j, -1j
        turns.append(j)
    return antipodes, turns


def _permutation(permutation, sector):
    if sector is None:
        dest = permutation
    else:
        dest = [sector.index(permutation[i], permutation[j]) for i, j in zip(sector.aa, sector.bb)]
    return np.eye(len(dest), dtype=int)[np.argsort(dest)]


def _one_edge_turns(car):
    for perm in rotations():
        fixed = [p for p in range(car.n) if perm[p] == p]
        if len(fixed) == 2 and all(perm[q] in car.neighbours[q] for q in car.neighbours[fixed[0]]):
            yield perm


def _measured_bands(h, source, detector):
    energies, vec = np.linalg.eigh(h)
    levels = np.split(np.arange(len(h)), np.flatnonzero(np.diff(energies) > TOL) + 1)
    poles, retained = [], []
    for level in levels:
        v = vec[:, level]
        overlap = v.T @ source
        weights = (overlap**2).sum(axis=0)
        if weights.max() < TOL:
            continue
        if np.ptp(weights) > TOL or np.min(weights) < TOL:
            raise ArithmeticError('the rotated joint preparation lost its uniform source weight')
        agrees = ((v @ overlap)[detector]**2).sum(axis=0) / weights
        u, s, _ = np.linalg.svd(overlap, full_matrices=False)
        active = v @ u[:, s > TOL]
        keep = bool(agrees.min() > .5)
        poles.append(dict(energy=float(energies[level].mean()), source_weight=float(weights.mean()),
                          agreement=float(agrees.mean()), dimension=active.shape[1], retained=keep))
        if keep:
            retained.append(active)
    return poles, retained


def _service(car, pair, gain):
    anti, turns = _geometry(car)
    eye = np.eye(car.n)
    loads = [np.diag(e) for e in eye]
    sector = PairSector(car) if pair else None
    if pair:
        stars = np.zeros((sector.P, car.n))
        for p in range(car.n):
            for q in car.neighbours[p]:
                stars[sector.index(p, q), p] = 1 / np.sqrt(len(car.neighbours[p]))
        source = (stars - stars[:, anti]) / np.sqrt(2)
        h = sector.hamiltonian(KAPPA).toarray()
        detector = sector.adjacent
        turns = [port.second_quantized(sector, j) for j in turns]
        loads = [port.second_quantized(sector, n) for n in loads]
    else:
        source = (eye - eye[:, anti]) / np.sqrt(2)
        h = car.laplacian().toarray() / 2  # Subtract the constant one-load energy 6 kappa.
        detector = np.ones(car.n, bool)
    poles, bands = _measured_bands(h, source, detector)
    if [v.shape[1] for v in bands] != [3, 3]:
        raise ArithmeticError('selected source is outside the two-triplet response domain')
    group = tuple(_one_edge_turns(car))
    class_sum = sum(_permutation(p, sector) for p in group)
    if not np.array_equal(class_sum, class_sum.T):
        raise ArithmeticError('the incidence class sum is not inversion closed')
    eigenvalues, axes = np.linalg.eigh(class_sum)
    positive = axes[:, eigenvalues > TOL]
    projector = positive @ positive.T
    currents, bare = [], []
    for n, j in zip(loads, turns):
        service = gain * projector @ n @ projector + j
        currents.append(block_diag(*(v.T @ service @ v for v in bands)))
        bare.append(block_diag(*(v.T @ n @ v for v in bands)))
    return sector, source, h, poles, bands, group, class_sum, projector, currents, bare


def _response(car, pair, gain):
    sector, source, h, poles, bands, group, class_sum, projector, currents, bare = _service(car, pair, gain)
    algebra = [-1j * m for m in currents]
    basis = port.span_basis(algebra)
    residuals = []
    for a, b in combinations(algebra, 2):
        bracket = port.realvec(a @ b - b @ a)
        residuals.append(float(np.linalg.norm(bracket - basis.T @ (basis @ bracket))))
    cov = 0.
    for perm in rotations():
        rotation = _permutation(perm, sector)
        u = block_diag(*(v.T @ rotation @ v for v in bands))
        cov = max(cov, *(np.linalg.norm(u @ currents[p] @ u.T - currents[perm[p]])
                         for p in range(car.n)))
    mask_error = float(np.linalg.norm(projector @ bands[0] - bands[0]))
    leakage = float(np.linalg.norm(projector @ bands[1]))
    # Report a reproducible upper bound under the receipt's decimal convention.
    tolerance = 1e-11
    if max(mask_error, leakage, max(residuals), cov) >= tolerance:
        raise ArithmeticError('the bounded symmetry, commutator or band selection residual was exceeded')
    data = dict(basis_dimension=len(h), source_poles=poles, one_edge_rotations=len(group),
                class_eigenvalues_on_retained=[float(np.trace(v.T @ class_sum @ v) / v.shape[1]) for v in bands],
                guarded=port.describe(algebra), with_bare_port_zero=port.describe(algebra + [-1j * bare[0]]),
                brackets_checked=len(residuals), port_rotation_checks=car.n * len(rotations()),
                max_checked_mask_bracket_and_rotation_error_below=tolerance)
    if not pair:
        return data

    # The compiled ideal program acts as the identity outside the selected source bands.
    embedding = np.column_stack(bands)
    state = source[:, SOURCE_PORT]
    initial = embedding.T @ state
    # Preserve the full coherent input including the portion of the source outside the retained modes.
    untouched = state - embedding @ initial
    free = expm(-16j * h)
    rows = []
    for p, current in enumerate(currents):
        probs = []
        for sign in (+1, -1):
            out = free @ (untouched + embedding @ (expm(-1j * sign * STRENGTH * current) @ initial))
            weights = abs(out)**2
            reading = np.r_[weights[sector.seam_states], weights[~sector.adjacent].sum()]
            if abs(reading.sum() - 1) > 1.e-10:
                raise ArithmeticError('some local detector outcomes were lost')
            probs.append(reading)
        rows.append(dict(port=p, positive=probs[0], negative=probs[1],
                         total_variation=float(sum(abs(probs[0] - probs[1])) / 2)))
    data['detector'] = dict(source_port=SOURCE_PORT, pulse_magnitude=STRENGTH, free_duration=16,
                            seam_outcomes=[list(edge) for edge in car.seams], final_outcome='quiet',
                            rows=rows)
    return data


def _portable(data):
    """Eleven decimal places resolve the selection tolerance and remove eigenbasis-dependent roundoff."""
    if isinstance(data, dict):
        return {k: _portable(v) for k, v in data.items()}
    if isinstance(data, (list, tuple, np.ndarray)):
        return [_portable(v) for v in data]
    if isinstance(data, (float, np.floating)):
        return round(float(data), 11) + 0.0
    return data


def run():
    """Rebuild the optional one- and two-load current menus and the ideal patch readout from the native cell."""
    car = cell()
    return _portable(dict(kappa=KAPPA, gain=1, spectral_tolerance=TOL,
                          one_load=_response(car, False, 1.), two_load=_response(car, True, 1.)))


def pulse(port_number, strength, *, loads=2):
    """One ideal guarded service on its whole one-cell sector, with identity on the unused modes."""
    if loads not in (1, 2) or port_number not in range(12) or not np.isfinite(strength):
        raise ValueError('choose one or two loads, a native port and a finite pulse magnitude')
    _, _, _, _, bands, _, _, _, currents, _ = _service(cell(), loads == 2, 1.)
    embedding = np.column_stack(bands)
    change = expm(-1j * strength * currents[port_number]) - np.eye(embedding.shape[1])
    return np.eye(embedding.shape[0]) + embedding @ change @ embedding.T


if __name__ == '__main__':
    print(json.dumps(run(), allow_nan=False, indent=2))
