"""Informational completeness of an unmarked material return detector."""

from math import comb

import numpy as np
from scipy import sparse
from sympy import Matrix
from scipy.sparse.linalg import expm_multiply

from ..ranging.certificate import integer_flight
from ...mobile_collar_source import measure as formation


PRIME = 1009
AXES = (np.eye(2), *formation.PAULI)
PAULI = np.array([np.kron(p, q) for p in AXES for q in AXES])


def terminals(board, route):
    return tuple(
        min((set(board.neighbours[outer]) & set(board.neighbours[inner])) - set(route))
        for outer, inner in (route[:2], route[2:])
    )


def add_pivot(row, basis):
    for j in range(len(row)):
        x = int(row[j])
        if not x:
            continue
        if j in basis:
            row = (row - x * basis[j]) % PRIME
        else:
            basis[j] = row * pow(x, -1, PRIME) % PRIME
            return True
    return False


def derivative_rank(h, ports, degree):
    """Use the observable derivatives of sharp launches and port clicks modulo a prime."""
    addresses = [4 * p + j for p in ports for j in range(4)]
    launch = np.zeros((h.shape[0], len(addresses)), np.int64)
    launch[addresses] = np.eye(len(addresses), dtype=np.int64)
    symmetric, skew = np.triu_indices(4), np.triu_indices(4, 1)
    powers, certificate, basis = [], [], {}
    for n in range(degree + 1):
        blocks = (
            launch[addresses]
            .reshape(len(ports), 4, len(ports), 4)
            .transpose(0, 2, 1, 3)
        )
        powers.append(blocks)
        for q in range(len(ports)):
            for p in range(len(ports)):
                # The n-th derivative of U_qp^* U_qp is the sum below,
                # times (-i/2)^n. H here is twice the centered native generator.
                term = (
                    sum(
                        ((-1) ** j * comb(n, j) % PRIME)
                        * (powers[j][q, p].T @ powers[n - j][q, p])
                        for j in range(n + 1)
                    )
                    % PRIME
                )
                coords = np.zeros(16, np.int64)
                if n % 2:
                    coords[10:] = term[skew]
                else:
                    coords[:10] = term[symmetric]
                if add_pivot(coords, basis):
                    certificate.append([n, ports[p], ports[q]])
        if len(basis) == 16:
            break
        launch = h @ launch % PRIME
    return dict(rank=len(basis), highest_order=certificate[-1][0], pivots=certificate)


def integer_return_minor(h, port, orders):
    """Recompute the single-terminal witness with unbounded integer path sums."""
    matrix = h.tocoo()
    start = 4 * port
    waves = np.zeros((h.shape[0], 4), dtype=object)
    waves[start : start + 4] = np.eye(4, dtype=int)
    coeff = []
    observations = []
    for n in range(1 + max(orders)):
        coeff.append(waves[start : start + 4].copy())
        if n in orders:
            effect = sum(
                (-1) ** j * comb(n, j) * (coeff[j].T @ coeff[n - j])
                for j in range(n + 1)
            )
            np.testing.assert_array_equal(effect.T, (-1) ** n * effect)
            observations.append(
                ([0] * 10 + list(effect[np.triu_indices(4, 1)]))
                if n % 2
                else (list(effect[np.triu_indices(4)]) + [0] * 6)
            )
        following = np.zeros_like(waves)
        np.add.at(following, matrix.row, matrix.data[:, None] * waves[matrix.col])
        waves = following
    determinant = int(Matrix(observations).det())
    assert determinant % PRIME != 0
    return dict(
        orders=orders,
        determinant=str(determinant),
        residue_mod_1009=determinant % PRIME,
        residue_mod_1013=determinant % 1013,
    )


def port_records(board, route, states):
    ports = terminals(board, route)
    generator = formation.hamiltonian(board, route, 1, 8)
    active = integer_flight(board, route)
    ranks = [derivative_rank(active, ports[:count], 20) for count in (1, 2)]
    scaled_difference = active - 2 * generator
    scalar = scaled_difference.diagonal()[0]
    np.testing.assert_allclose(
        (scaled_difference - sparse.eye(active.shape[0]) * scalar).data, 0, atol=1e-12
    )
    orders = [row[0] for row in ranks[0]["pivots"]]
    proof = integer_return_minor(active, ports[0], orders)
    changed = formation.hamiltonian(board, route, 1, 0)
    position = changed[::4, ::4]
    assert (changed - sparse.kron(position, np.eye(4), format="csr")).nnz == 0

    ids = np.array([4 * p + j for p in ports for j in range(4)])
    start = np.zeros((generator.shape[0], 8), complex)
    start[ids] = np.eye(8)
    flight = expm_multiply(-1j * generator, start, start=0, stop=8, num=33, traceA=0)
    blocks = flight[1:, ids].reshape(32, 2, 4, 2, 4).transpose(0, 1, 3, 2, 4)
    held_remote = []
    rng = np.random.default_rng(7384)
    shots = 2**40
    width = 2 * np.sqrt(np.log(2 * 4 / 0.01) / shots)
    for t in (0.5, 1, 2, 4):
        block = blocks[int(4 * t) - 1, 0, 0]
        effect = block.conj().T @ block
        probs = [
            float(np.trace(effect @ states[key]).real)
            for key in ("formed", "independent")
        ]
        counts = [int(rng.binomial(shots, p)) for p in probs]
        delta = (counts[0] - counts[1]) / shots
        held_remote.append(
            dict(
                t=t,
                return_probabilities=probs,
                yes_counts=counts,
                population_difference=probs[0] - probs[1],
                confidence_interval=[delta - width, delta + width],
            )
        )
    out = []
    for count in (1, 2):
        packet = blocks[:, :count, :count].reshape(-1, 4, 4)
        effects = packet.conj().transpose(0, 2, 1) @ packet
        design = np.einsum("rij,kji->rk", effects, PAULI).real / 4
        singular = np.linalg.svd(design[:, 1:], compute_uv=False)
        errors = {}
        for key, rho in states.items():
            values, axes = np.linalg.eigh(rho)
            observed = (abs(packet @ axes) ** 2).sum(axis=1) @ values
            bloch = np.linalg.lstsq(design[:, 1:], observed - design[:, 0], rcond=None)[
                0
            ]
            inferred = (PAULI[0] + np.einsum("k,kij->ij", bloch, PAULI[1:])) / 4
            errors[key] = float(np.linalg.norm(inferred - rho, ord="nuc") / 2)
        out.append(
            dict(
                ports=list(ports[:count]),
                appointments=32,
                settings=len(design),
                singular_min=float(singular[-1]),
                singular_max=float(singular[0]),
                conditioned_on_known_trace=True,
                reconstruction_distance=errors,
            )
        )
    return dict(
        integer_prime=PRIME,
        exact_single_port_minor=proof,
        held_return=dict(
            source_order=["formed", "independent"],
            seed=7384,
            shots=shots,
            simultaneous_confidence=0.99,
            settings=held_remote,
        ),
        finite_rank_certificates=ranks,
        changed_readout_rank=1,
        finite_population_inversions=out,
    )
