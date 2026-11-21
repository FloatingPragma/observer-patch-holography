"""Check the independent dense source for the co-active factorization."""

import itertools
from types import SimpleNamespace

import numpy as np

from model.executable_model.carrier import Carrier
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from .coactive import CURRENT_OUTCOMES, _joint_arms
from .contrast_binding import _all_links_record, _external_laplacian
from .current import one_charge_read


EDGES = ((0, 1), (1, 2), (2, 3))
LINK = (1, 2)
PAIRS = tuple(itertools.combinations(range(4), 2))
BASIS = tuple(itertools.product(range(4), PAIRS, range(2)))
KAPPA = 4
G = 0.75
WAIT = 0.29


def _full_source(record_positive):
    """Set matrix elements by physical moves in the entire tensor space."""
    index = {coordinate: i for i, coordinate in enumerate(BASIS)}
    source = np.zeros((len(BASIS), len(BASIS)))
    for column, (charge, loads, cage_cell) in enumerate(BASIS):
        source[column, column] = 1 - 2 * KAPPA * (loads in EDGES)
        for a, b in EDGES:
            if charge in (a, b):
                hopping = 1 / 8 if (a, b) == LINK else 1
                source[column, column] += hopping
                source[index[a + b - charge, loads, cage_cell], column] -= hopping
            if (a in loads) != (b in loads):
                moved = tuple(sorted(set(loads).symmetric_difference((a, b))))
                source[column, column] += 0.5
                source[index[charge, moved, cage_cell], column] -= 0.5
        weight = G * (((1 in loads) != (2 in loads)) + int(record_positive))
        source[column, column] += weight
        source[index[charge, loads, 1 - cage_cell], column] -= weight
    return source


def _propagator(hamiltonian):
    np.testing.assert_array_equal(hamiltonian, hamiltonian.conj().T)
    values, modes = np.linalg.eigh(hamiltonian)
    return (modes * np.exp(-1j * WAIT * values)) @ modes.conj().T


def test_factorized_coactive_records_match_dense_three_register_source_everywhere():
    native = PairSector(Carrier("two-cell fixture", 4, list(EDGES)))
    sparse_sector = SimpleNamespace(P=native.P, car=SimpleNamespace(cells=2), adjacent=native.adjacent)
    load_source = native.hamiltonian(KAPPA)
    activation = np.array([(1 in pair) != (2 in pair) for pair in PAIRS], dtype=int)
    # Use port 12 for the first address of the second cell in the production assembler.
    sparse_contacts = [
        _external_laplacian(len(PAIRS), 2, [((1, 12), activation + active)])
        for active in (False, True)
    ]
    dense_propagators = [_propagator(_full_source(active)) for active in (False, True)]
    charge_source = np.eye(4)
    for a, b in EDGES:
        difference = np.eye(4)[a] - np.eye(4)[b]
        charge_source += (1 / 8 if (a, b) == LINK else 1) * np.outer(difference, difference)
    charge_propagator = _propagator(charge_source)
    current = np.zeros((4, 4), complex)
    current[1, 2] = -1j
    current[2, 1] = 1j
    eigenvalues, vectors = np.linalg.eigh(current)
    seeds = np.vstack((np.eye(4), [1, -0.25j, 0.5 + 0.5j, 0.6j]))
    selectors = [(position // 2, int(loads in EDGES), cage_cell) for position, loads, cage_cell in BASIS]
    for pair_index, sender in itertools.product(range(len(PAIRS)), range(2)):
        load_start = np.eye(len(PAIRS))[pair_index]
        cage_start = np.eye(2)[sender]
        load_cage = np.array([
            _all_links_record(sparse_sector, load_source, contact, load_start, sender, G, WAIT)
            for contact in sparse_contacts
        ])
        for seed, sign in itertools.product(seeds, (-1, 1)):
            # Deliberately retain a discarded preparation probability of 0.03.
            prepared = np.sqrt(0.97) * seed / np.linalg.norm(seed)
            sparse_branches = one_charge_read(prepared, LINK, sign)
            charge_record = np.array([
                (abs(charge_propagator @ sparse_branches[r]) ** 2).reshape(2, 2).sum(axis=1)
                for r in CURRENT_OUTCOMES
            ])
            feedback, shuffled = _joint_arms(charge_record, load_cage)
            spaces = [vectors[:, np.isclose(eigenvalues, sign * record)] for record in CURRENT_OUTCOMES]
            dense_branches = [space @ (space.conj().T @ prepared) for space in spaces]
            probability_on = np.linalg.norm(dense_branches[-1]) ** 2 / 0.97
            dense_records = np.zeros((2, 3, 2, 2, 2))
            for j, branch in enumerate(dense_branches):
                product = np.kron(branch, np.kron(load_start, cage_start))
                for gate, propagator in enumerate(dense_propagators):
                    for selector, probability in zip(selectors, abs(propagator @ product) ** 2, strict=True):
                        dense_records[(gate, j, *selector)] += probability
            np.testing.assert_allclose(
                feedback, np.stack((dense_records[0, 0], dense_records[0, 1], dense_records[1, 2])),
                atol=1e-12, rtol=0,
            )
            np.testing.assert_allclose(
                shuffled, (1 - probability_on) * dense_records[0] + probability_on * dense_records[1],
                atol=1e-12, rtol=0,
            )
            np.testing.assert_allclose([feedback.sum(), shuffled.sum()], [0.97, 0.97], atol=1e-12, rtol=0)
