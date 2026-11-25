"""Measure core-to-exterior conditional information and signed Gibbs response on one cell."""

import json
from fractions import Fraction
from itertools import combinations
from math import factorial

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import cell
from model.experiments.seam_program.collective_rate.measure import (
    clock_seam, reduced_gibbs_states, spectrum,
)
from model.experiments.seam_program.collective_rate.shell_response.measure import PARITY_PORTS
from .regional_moments import pair_transfers, embed


def entropy(probs):
    positive = probs[probs > 0]
    return float(-positive @ np.log(positive))


def marginal(probabilities, ports):
    words = np.arange(len(probabilities))
    local = sum(((words >> port) & 1) << j for j, port in enumerate(ports))
    return np.bincount(np.broadcast_to(local, words.shape), probabilities, minlength=2**len(ports))


def hidden_diagonal_response(region, modular_diagonal):
    """Largest signed slope at fixed complete bath and every one- and two-port state.

    The regional tangent is diagonal and has trace norm at most one. Normalized
    identity on exterior spectators extends it to a full admissible tangent.
    """
    words = np.arange(2**len(region))
    signs = 1 - 2*((words[None, :] >> np.arange(len(region))[:, None]) & 1)
    bath = words >> 2  # the first two ports form the core
    constraints = np.array([*(bath == i for i in np.unique(bath)), *signs,
                            *(x*y for x, y in combinations(signs, 2))], dtype=float)
    # Dual: minimize the maximum error of a predictor built from the retained readings.
    count = len(words)
    design = constraints.T
    inequalities = np.column_stack([np.vstack([design, -design]), -np.ones(2*count)])
    dual = linprog([*np.zeros(len(design.T)), 1.], A_ub=inequalities,
                   b_ub=np.r_[modular_diagonal, -modular_diagonal], bounds=(None, None),
                   options={"dual_feasibility_tolerance": 1e-10, "primal_feasibility_tolerance": 1e-10})
    if not dual.success:
        raise ValueError(dual.message)
    # The signed difference of the two inequality multipliers is a primal witness.
    u, v = np.split(dual.ineqlin.marginals, 2)
    tangent = u-v
    np.testing.assert_allclose(constraints @ tangent, 0, atol=1e-8)
    assert np.abs(tangent).sum() <= 1+1e-8
    slope = float(tangent @ modular_diagonal)
    np.testing.assert_allclose(slope, dual.fun, atol=1e-8)
    return slope


def native_front(board, a):
    """Leading powers of the native one-load generator, checked against shortest-path counts."""
    sector = BinarySector(board, 1)
    assert sector.states == [(port,) for port in range(board.n)]
    h = sector.hamiltonian(KAPPA, shifted=False)
    doubled = 2*h
    np.testing.assert_array_equal(doubled, doubled.astype(int))
    h2 = csr_matrix(doubled.astype(int))
    for i, neighbors in enumerate(board.neighbours):
        row = h2.getrow(i)
        off = {int(j): int(v) for j, v in zip(row.indices, row.data, strict=True) if j != i}
        assert off == {j: -1 for j in neighbors}
    counts = {a: 1}
    frontier = {a: 1}
    v = {a: 1}
    # At distance d only the previous frontier can contribute. Restricting to new
    # sites before each multiplication retains the leading term without integer overflow.
    distances = []
    while frontier:
        layer = len(distances)
        for i, count in frontier.items():
            assert v[i] == (-1)**layer * count
        distances.append({i: Fraction(n, 2**layer*factorial(layer)) for i, n in frontier.items()})
        following = {}
        for i, count in frontier.items():
            for j in board.neighbours[i]:
                if j not in counts:
                    following[j] = following.get(j, 0) + count
        # Use the entries of the actual one-load Hamiltonian, independently of neighbor counts.
        v = {j: sum(int(h2[j, i])*int(v[i]) for i in board.neighbours[j] if i in frontier)
             for j in following}
        counts.update(following)
        frontier = following
    return distances


def measure(beta=1/16, strength=KAPPA):
    board = cell()
    a, b = clock_seam()
    bath = tuple(p for p in range(board.n) if p not in (a, b))
    left, right = set(board.neighbours[a]), set(board.neighbours[b])
    buffers = [(), tuple(sorted(left & right)), tuple(sorted((left | right) - {a, b}))]
    regions = [(a, b, *buffer) for buffer in buffers]
    blocks = spectrum(strength)
    mu = strength * (6 - board.degree[a])
    weights = [np.exp(-beta * (energies - mu*n)) for n, (energies, _, _, _) in enumerate(blocks)]
    partition = sum(w.sum() for w in weights)
    probabilities = [w / partition for w in weights]
    full_entropy = entropy(np.concatenate(probabilities))
    reductions = reduced_gibbs_states(blocks, beta, [bath, *regions, *buffers[1:]], strength)
    bath_entropy = entropy(eigh(reductions[0], eigvals_only=True))
    buffer_entropies = [0., *(entropy(eigh(rho, eigvals_only=True)) for rho in reductions[4:])]

    # Both full tangents have zero bath marginal. Their squares are diagonal in native words.
    dim = 2**board.n
    floor = min(p.min() for p in probabilities)
    swap, _ = pair_transfers(strength)
    xor_trace = 0.
    for (_, vectors, masks, _), p in zip(blocks, probabilities, strict=True):
        xor = (((masks >> a) ^ (masks >> b)) & 1).astype(bool)
        xor_trace += (vectors[xor]**2).sum(axis=0) @ (1/p)
    q = dict(swap=float(xor_trace/(dim/4)**2),
             parity=float(sum((1/p).sum() for p in probabilities) / dim**2))
    c_bound = dict(swap=float(4/(dim * floor)), parity=float(1/(dim * floor)))
    # Commuting contrast-only control, identical carrier, all 30 seams, same chemical potential.
    words = np.arange(dim)
    unequal = sum(((words >> i) ^ (words >> j)) & 1 for i, j in board.seams)
    classical = np.exp(-beta * strength * unequal)
    classical /= classical.sum()
    classical_offset = entropy(marginal(classical, bath)) - entropy(classical)
    # The cell is edge-transitive; the uncorrelated endpoint marginals are both I/2.
    # At particle-hole symmetry the full grand energy is a sum of identical edge terms.
    ids = [native_law.SRC.index(w) for w in ((0, 0), (1, 0), (0, 1), (1, 1))]
    contact = native_law.seam_term(strength)[np.ix_(ids, ids)]
    pair = reductions[1]
    for subsystems in ((0, 2), (1, 3)):
        np.testing.assert_allclose(np.trace(pair.reshape(2, 2, 2, 2), axis1=subsystems[0],
                                            axis2=subsystems[1]), np.eye(2)/2, atol=1e-13)
    mean_contact = float(np.sum(contact.T*pair))
    grand_energy = sum(p @ (energies-mu*n) for n, ((energies, _, _, _), p)
                       in enumerate(zip(blocks, probabilities, strict=True)))
    np.testing.assert_allclose(grand_energy, len(board.seams)*mean_contact, atol=1e-12)
    crossings = sum((i in (a, b)) != (j in (a, b)) for i, j in board.seams)
    reset_global = float(beta*crossings*(np.trace(contact)/4-mean_contact))
    buffer_states = [np.ones((1, 1)), *reductions[4:]]
    rows = []
    for region, buffer, density, sb in zip(regions, buffers, reductions[1:4], buffer_entropies, strict=True):
        probs, vectors = eigh(density)
        modular = (vectors * -np.log(probs)) @ vectors.T
        info = entropy(probs) - sb + bath_entropy - full_entropy
        responses = dict(swap=float(np.sum(modular.T * embed(swap, (0, 1), len(region)))), parity=0.)
        if set(PARITY_PORTS) <= set(region):
            mask = sum(1 << region.index(p) for p in PARITY_PORTS)
            signs = np.array([(-1)**((w & mask).bit_count()) for w in range(len(probs))])
            responses['parity'] = float(signs @ modular.diagonal() / len(probs))
        bounds = {}
        for name, squared_norm in q.items():
            # This step minimizes I/t+t*q over the interval certified by the global floor.
            t = min(np.sqrt(max(info, 0.)/squared_norm), 1/c_bound[name])
            bounds[name] = float(info/t+t*squared_norm) if t else 0.
        diagonal_info = (entropy(marginal(classical, region)) - entropy(marginal(classical, buffer))
                         + classical_offset)
        diagonal_modular = -np.log(marginal(classical, region))
        reset_tangent = np.kron(buffer_states[len(rows)], pair)-density
        reset_regional = float(np.sum(modular.T*reset_tangent))
        rows.append(dict(buffer=buffer, conditional_information=info,
                         reset_signed_error=reset_regional-reset_global,
                         contrast_only_information=diagonal_info,
                         contrast_only_hidden_slope=hidden_diagonal_response(region, diagonal_modular),
                         max_hidden_diagonal_slope=hidden_diagonal_response(region, modular.diagonal()),
                         signed_error=dict(swap=responses['swap']-beta, parity=responses['parity']),
                         certified_absolute_error=bounds))
    return dict(beta=beta, contrast=strength, core=[a, b], relative_quadratic_cost=q,
                relative_norm_upper=c_bound, reset_global_response=reset_global, records=rows)


if __name__ == '__main__':
    print(json.dumps(measure(), indent=2))
