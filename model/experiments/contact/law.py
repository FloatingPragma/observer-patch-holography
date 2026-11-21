"""Local source map and finite instruments for the contrast-square contact."""

import math

import numpy as np
from scipy.linalg import eigh, expm

from model.candidate_laws.contrast_memory_spherical_v1 import collar, law as contrast
from model.candidate_laws.scalar_channel.carrier import require
from model.candidate_laws.scalar_channel.compiler import build, local_cycles, anti
from model.candidate_laws.scalar_channel.instruments import compile_matrix, decode_paulis
from model.experiments.cage import construction as cage, instruments as apparatus


def collar_states(sign):
    """Empty loads, with a blank memory or the sign eigenstate written by the collar fuel."""
    require(sign in (-1, 1), 'collar sign must be +1 or -1')
    arms = np.zeros((contrast.DIM, 2))
    arms[contrast.SRC.index((0, 0)), 0] = 1
    for d, amplitude in ((1, 1), (-1, sign)):
        label = contrast.LABELS.index((0, d))
        arms[contrast.N_SRC * (1 + label) + contrast.SRC.index((0, 0)), 1] = amplitude / math.sqrt(2)
    return arms


def arm_energies(kappa, coupling, filling, sign):
    exchange, square = contrast.local_term_parts()
    arms = collar_states(sign)
    action = (exchange + (kappa + coupling * filling) * square) @ arms
    reduced = arms.T @ action
    require(np.max(abs(action - arms @ reduced)) < 1e-12 and abs(reduced[0, 1]) < 1e-12,
            'the contrast law leaves the two prepared arms')
    return reduced.diagonal()


def interaction(graph, seam, coupling):
    """Decode the scalar factor n_s from its two parity words on the register graph."""
    site = graph['data']['index'][seam]
    word = graph['B'][site]
    require(all(not anti(word, loop['pauli']) for loop in local_cycles(graph)),
            'seam occupation changes the loop code')
    diagonal = np.zeros((graph['n'], graph['n']))
    diagonal[site, site] = 1
    words = compile_matrix(graph, diagonal)
    require(words == {(0, 0, 0): .5, word: -.5}, 'unexpected seam number words')
    decoded, residual = decode_paulis(graph, words)
    require(residual < 1e-13 and np.max(abs(decoded - diagonal)) < 1e-13,
            'seam occupation fails to decode')
    n = graph['nsource']
    return coupling * decoded[:n, :n].real


def propagate(h, t):
    values, modes = eigh(h)
    return (modes * np.exp(-1j * t * values)) @ modes.conj().T


def ramp(car, contact, steps):
    """Propagate every occupied wedge orbital, retaining the contact during every star ramp."""
    D = cage.law(car)
    seams, wedges = car['C'].shape
    frame = np.eye(len(D), dtype=complex)[:, seams:]
    base = np.diag(D.diagonal()) + contact
    centers = np.array([b for _, b, _ in car['wedges']])
    for center in range(car['ports']):
        sites = seams + np.flatnonzero(centers == center)
        star = np.zeros_like(base)
        star[:seams, sites] = D[:seams, sites]
        star[sites, :seams] = D[sites, :seams]
        for k in range(steps):
            s = (k + .5) / steps
            frame = propagate(base + s * s * (3 - 2 * s) * star, 16 / steps) @ frame
        base += star
    require(np.max(abs(base - D - contact)) < 1e-12, 'incomplete ramp source')
    require(np.linalg.norm(frame.conj().T @ frame - np.eye(wedges), 2) < 1e-10,
            'ramp lost occupied-frame normalization')
    return frame @ frame.conj().T


def memory_write(kappa):
    energies = arm_energies(kappa, 0, 0, 1)
    driven = np.diag(energies) + collar.G * np.array([[0, 1], [1, 0]])
    return abs(expm(-1j * math.pi / (2 * collar.G) * driven)[:, 0]) ** 2


def ramsey(kappa, coupling):
    """Full-fuel probabilities indexed by filling, collar sign, phase setting and query."""
    laws = np.zeros((2, 2, len(collar.SETTINGS), collar.QUERIES))
    for filling in (0, 1):
        for j, sign in enumerate((1, -1)):
            energy = arm_energies(kappa, coupling, filling, sign)
            driven = np.diag(energy) + collar.G * np.array([[0, 1], [1, 0]])
            pulse = expm(-1j * math.pi / (4 * collar.G) * driven)
            for k, theta in enumerate(collar.SETTINGS):
                tail = pulse[0] * np.array([np.exp(-1j * theta), 1])
                times = np.arange(collar.QUERIES) * collar.CYCLE
                laws[filling, j, k] = abs(np.exp(-1j * np.outer(times, energy)) @ (tail * pulse[:, 0])) ** 2
    return laws


def collision(protocol):
    """All joint histories of a retained write outcome and the two scalar pointers."""
    car = cage.cell()
    graph = build()
    ring = tuple(protocol['ring'])
    seam = tuple(protocol['marked_seam'])
    profile = cage.ring_profile(car, ring)
    pair, _, error = cage.modes(car, profile)
    require(error < 1e-13, 'unperturbed ring is not a cage')
    D = cage.law(car)
    unit = interaction(graph, seam, 1)
    p_write = memory_write(protocol['kappa'])
    laws = []
    conditional = {}
    for g in protocol['couplings_in_exchange_units']:
        term = g * unit
        prepared = ramp(car, term, protocol['ramp_midpoints_per_star'])
        full, _ = apparatus.persistence_law(D, pair[:, 0], pair, prepared, term=term)
        conditional[g] = apparatus.pointers(full)
    for g in protocol['couplings_in_exchange_units']:
        # Write fuel full means blank, fuel empty means written. The two terminal
        # bits copy the scalar fuel; 01 and 10 are retained zero-probability outcomes.
        laws.append(np.concatenate([p_write[1] * conditional[g], p_write[0] * conditional[0]], axis=-1))
    return np.array(laws)
