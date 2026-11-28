"""Classify a single stationary memory beside a fixed ternary boundary pin."""

from itertools import product

import numpy as np
from sympy import Matrix
from model.candidate_laws.interval_echo import candidate


def classify():
    law = candidate.native_law
    # Eight times the exact local term. The blank code reduces the term;
    # on written memories the two echo blocks are (2I-P P^T)/2,
    # where P is the integer write pattern, and the generator is G2/2.
    pattern = law._write_pattern()
    reverse, square = law._label_ops()
    swap, contrast = law._load_ops()
    k = int(candidate.KAPPA)
    assert k == candidate.KAPPA
    loads = np.eye(law.N_SRC, dtype=int)
    memories = np.eye(law.N_LAB, dtype=int)
    twice_echo = 2 * np.eye(len(pattern), dtype=int) - pattern @ pattern.T
    twice_gen = np.kron(memories - reverse + 2 * k * square, loads)
    h8 = np.zeros((law.DIM, law.DIM), dtype=int)
    h8[: law.N_SRC, : law.N_SRC] = 4 * (loads - swap) + 8 * k * contrast
    h8[law.N_SRC :, law.N_SRC :] = twice_echo @ twice_gen @ twice_echo
    exchange, cost = law.local_term_parts()
    np.testing.assert_allclose(h8, 8 * (exchange + k * cost), atol=1e-12)
    # The generator's contrast-square sectors 0,1,4 have symmetric and
    # antisymmetric reversal levels; at zero contrast only symmetry occurs.
    energies = (0, k, k + 1, 4 * k, 4 * k + 1)
    result = {}
    for reverse_word, pin in product((False, True), range(3)):
        grams = []
        for resident in (0, 1):
            word = (pin, resident) if reverse_word else (resident, pin)
            columns = np.eye(law.DIM, dtype=int)[:, law.SRC.index(word) :: law.N_SRC]
            acted = h8 @ columns
            defects = [acted - 8 * e * columns for e in energies]
            grams.append({e: a.T @ a for e, a in zip(energies, defects, strict=True)})
        pairs = []
        for e0, e1 in product(energies, repeat=2):
            # ker(M^T M + N^T N) is precisely ker(M) intersect ker(N).
            gram = grams[0][e0] + grams[1][e1]
            rank = Matrix(gram.tolist()).to_DM().rank()
            if rank < len(gram):
                pairs.append([e0, e1, len(gram) - rank])
        result[f"{'right' if reverse_word else 'left'}_resident_pin_{pin}"] = pairs
    return result
