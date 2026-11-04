"""Independent finite source controls for the exact origin Schur interface."""
import numpy as np
from scipy.linalg import solve
from ..bounded_defect_index import finite as source
from . import bounds


def blocks(radius, spacing=1, velocity=1, mass=4, hopping_bound=2, sign=1):
    if type(radius) is not int or radius < 1:
        raise ValueError('integer symmetric cube radius at least one required')
    estimate = bounds.spectral_window(spacing, velocity, mass, hopping_bound)
    h = source.source_matrix(radius, spacing, velocity, mass, sign)
    onsite = source.source_matrix(radius, spacing, velocity, mass, sign, hopping_fraction=0)
    points = source.sites(radius)
    origin = points.index((0, 0, 0))
    origin_rows = np.array([2*origin, 2*origin+1])
    other_rows = np.delete(np.arange(h.shape[0]), origin_rows)
    return {'H': h, 'B': onsite, 'P_rows': origin_rows, 'Q_rows': other_rows,
            'D': h[other_rows][:, other_rows],
            'V': h[other_rows][:, origin_rows].toarray(),
            'B_Q': onsite[other_rows][:, other_rows],
            'estimate': estimate, 'points': points}


def kernel_frame(radius, spacing=1, velocity=1, mass=4, hopping_bound=2, sign=1):
    value = blocks(radius, spacing, velocity, mass, hopping_bound, sign)
    graph = -solve(value['D'].toarray(), value['V'], assume_a='her')
    raw = np.zeros((value['H'].shape[0], 2), complex)
    raw[value['P_rows']] = np.eye(2)
    raw[value['Q_rows']] = graph
    gram = raw.conj().T@raw
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    inverse_root = (eigenvectors/np.sqrt(eigenvalues))@eigenvectors.conj().T
    normalized = raw@inverse_root
    local = normalized.reshape(-1, 2, 2)
    density = np.einsum('nsk,nsk->n', local.conj(), local).real/2
    return {**value, 'graph': graph, 'raw_frame': raw, 'gram': gram,
            'frame': normalized, 'density': density}


def neumann_frame(radius, terms, spacing=1, velocity=1, mass=4, hopping_bound=2, sign=1):
    if type(terms) is not int or terms < 0:
        raise ValueError('nonnegative integer Neumann term count required')
    value = blocks(radius, spacing, velocity, mass, hopping_bound, sign)
    inverse_mass = value['B_Q']/float(value['estimate']['mass'])**2
    kinetic = value['D']-value['B_Q']
    step = -inverse_mass@kinetic
    term = -inverse_mass@value['V']
    graph = np.zeros_like(term)
    for _ in range(terms):
        graph += term
        term = step@term
    frame = np.zeros((value['H'].shape[0], 2), complex)
    frame[value['P_rows']] = np.eye(2)
    frame[value['Q_rows']] = graph
    return frame


def embed(radius_from, radius_to, values):
    if radius_from > radius_to:
        raise ValueError('target cube must contain the source cube')
    old, following = source.sites(radius_from), source.sites(radius_to)
    if values.shape[0] != 2*len(old):
        raise ValueError('every original spin/site row is required')
    addresses = {point: j for j, point in enumerate(following)}
    output = np.zeros((2*len(following), values.shape[1]), complex)
    for index, point in enumerate(old):
        target = addresses[point]
        output[2*target:2*target+2] = values[2*index:2*index+2]
    return output
