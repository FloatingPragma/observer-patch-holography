"""Joint native control and a charge-selective receiver on a unit-filled cell."""

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from sympy import Matrix, diag

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell

SIGNS = (-1, 1)
HERE = Path(__file__).resolve().parent


def commutator(a, b):
    return a * b - b * a


def algebra_certificate():
    """Integer identities suffice to isolate X on either sign of every edge."""
    qminus, qplus = (-8, 24)
    d = diag(qminus, 0, qplus, 0)
    edge = Matrix([[17, -1], [-1, 17]])  # twice the active Hamiltonian
    h = diag(edge, edge)
    two = commutator(d, commutator(d, h))
    four = commutator(d, commutator(d, two))
    matrix = Matrix([[qminus**2, qplus**2], [qminus**4, qplus**4]])
    coeff = matrix.inv()
    isolated = [coeff[j, 0] * two + coeff[j, 1] * four for j in range(2)]
    x = Matrix([[0, -1], [-1, 0]])
    zero = Matrix.zeros(2)
    assert isolated == [diag(x, zero), diag(zero, x)]
    return dict(
        relative_reference_frequencies=[qminus, qplus],
        commutator_determinant=int(matrix.det()),
        isolation_coefficients=[[str(i) for i in row] for row in coeff.tolist()],
    )


def native_controls():
    """Check the endpoint and all 97 spectator rows of each pair block."""
    board = cell()
    i, j = 0, 1
    assert j in board.neighbours[i] and board.degree[i] == board.degree[j] == 5
    law = candidate.native_law
    checks = law.exact_checks()
    assert checks["checks"]["reference_pair_gives_n_squared_at_reference_zero"]
    h0, h1 = law.local_term_parts()
    full = h0 + candidate.KAPPA * h1
    h = []
    degree = []
    for sign in SIGNS:
        indices = [law.SRC.index(pair) for pair in ((1 + sign, 1), (1, 1 + sign))]
        embedded = np.eye(law.DIM)[:, indices]
        pair = embedded.T @ full @ embedded
        np.testing.assert_allclose(full @ embedded, embedded @ pair, atol=2e-13, rtol=0)
        np.testing.assert_allclose(pair, [[8.5, -0.5], [-0.5, 8.5]], atol=2e-13, rtol=0)
        np.testing.assert_allclose(full[:, law.SRC.index((1, 1))], 0, atol=2e-13)
        h.append(pair)
        degree.append(candidate.KAPPA * np.diag([(1 + sign) ** 2, 1]))
    return h, degree


def click(p, shots, rng):
    hits = int(rng.binomial(shots, p))
    return [hits, shots - hits]


def report():
    c = algebra_certificate()
    h, degree = native_controls()
    conditional = []
    conjugates = []
    for pair, degree_phase in zip(h, degree):
        splitter = expm(-0.5j * np.pi * pair)
        phase = expm(-1j * np.pi * degree_phase / 12)
        inverse = expm(-1.5j * np.pi * pair)
        np.testing.assert_allclose(inverse, splitter.conj().T, atol=4e-14)
        u = inverse @ phase @ splitter
        p = float(abs(u[1, 0]) ** 2)
        conditional.append(p)
        conjugates.append(u)
    np.testing.assert_allclose(conditional, [0.75, 0], atol=3e-14)
    shots = 2**18
    seed = 871031
    rng = np.random.default_rng(seed)
    # Native sign-load census on every copy, then the mobile-mark address.
    # Fixed public source tally with one in four vacancy starts.
    trials = [shots // 4, 3 * shots // 4]
    joint = [click(p, n, rng) for p, n in zip(conditional, trials)]
    for u in conjugates:
        np.testing.assert_allclose(u.conj().T @ u, np.eye(2), atol=3e-14)
    # Counterfactual with the same sign-conditional sharp laws but uncalibrated prior.
    t = Fraction(1, 4)
    w = t**2 / (t**2 + (1 - t) ** 2)
    # A different receiver violates disjoint-address invariance within one sign.
    # Before and after splitting only its two spectator rails, the first
    # address is untouched and has squared norm 1/2 in both shapes.
    shapes = [
        (Fraction(1, 2), Fraction(1, 2), Fraction(0)),
        (Fraction(1, 2), Fraction(1, 4), Fraction(1, 4)),
    ]
    spectator_read = [x[0] ** 2 / sum(v**2 for v in x) for x in shapes]
    return dict(
        algebra=c,
        cell_ports=cell().n,
        joint_lie_dimension=2 * cell().n ** 2,
        native_exposures=["pi/2", "pi/12", "3pi/2"],
        addressed_edge=[0, 1],
        source_flag_counts=trials,
        shots=shots,
        seed=seed,
        counts_order=["remote port", "initial port"],
        conditional_remote_probabilities=conditional,
        joint_counts=joint,
        tagged_remote_probability="3/16",
        unverified_prior_w_at_quarter=str(w),
        unverified_prior_remote_probability=str(w * Fraction(3, 4)),
        normalized_fourth_power_spectator_reads=[str(x) for x in spectator_read],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    text = json.dumps(report(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "results.json").write_text(text)
    else:
        print(text, end="")
