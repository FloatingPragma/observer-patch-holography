"""Compare two spatially uniform native energy variations in the same reduced Gibbs state."""

from itertools import combinations

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import cell
from model.experiments.seam_program.collective_rate.measure import clock_seam, reduced_gibbs_states, spectrum


def pair_transfers(strength):
    symmetric = np.array([0, 1, 1, 0])
    odd = np.array([0, 1, -1, 0])
    sym = np.outer(symmetric, symmetric) / 2
    swap = np.outer(odd, odd) / 2 - sym
    aligned = sym - np.diag([0.5, 0, 0, 0.5])
    ids = [native_law.SRC.index(word) for word in ((0, 0), (1, 0), (0, 1), (1, 1))]
    h = native_law.seam_term(strength)[np.ix_(ids, ids)]
    normalized = []
    for tangent in (swap, aligned):
        for pair in ((0, 2), (1, 3)):
            np.testing.assert_array_equal(np.trace(tangent.reshape(2, 2, 2, 2),
                                                    axis1=pair[0], axis2=pair[1]), 0)
        energy = np.trace(h @ tangent)
        if energy <= 0:
            raise ValueError("positive energy transfer required to normalize both launches")
        normalized.append(tangent / energy)
    return normalized


def embed(tangent, positions, count):
    """Tensor with normalized spectator identity on the declared low-bit basis."""
    a, b = positions
    masks = np.arange(2**count)
    fixed_mask = (1 << a) | (1 << b)
    bases = masks[(masks & fixed_mask) == 0]
    ids = bases[:, None] | np.array([0, 1 << a, 1 << b, fixed_mask])
    out = np.zeros((len(masks), len(masks)))
    out[ids[:, :, None], ids[:, None, :]] = tangent / len(bases)
    return out


def conserved_weights_witness():
    """A native two-load matrix entry enforcing proportionality of global conserved weights."""
    sector = BinarySector(cell(), 2)
    twice_exchange = 2 * sector.hamiltonian(0, shifted=False)
    contrast = sector.hamiltonian(1, shifted=False) - twice_exchange/2
    np.testing.assert_array_equal(twice_exchange, twice_exchange.astype(int))
    np.testing.assert_array_equal(contrast, contrast.astype(int))
    x, c = twice_exchange.astype(int), contrast.astype(int)
    i, j = sector.index[(0, 1)], sector.index[(0, 8)]
    return int((x@c - c@x)[i, j])


def measure(strength=KAPPA, beta=1/16):
    board = cell()
    a, b = clock_seam()
    region = tuple(sorted({a, b, *board.neighbours[a], *board.neighbours[b]}))
    positions = {port: j for j, port in enumerate(region)}
    blocks = spectrum(strength)
    density, = reduced_gibbs_states(blocks, beta, [region], strength)
    mu = strength * (6 - board.degree[0])
    all_weights = np.concatenate([np.exp(-beta * (energies - mu*n))
                                  for n, (energies, _, _, _) in enumerate(blocks)])
    global_floor = all_weights.min() / all_weights.sum()
    eigenvalues, vectors = eigh(density)
    modular = (vectors * -np.log(eigenvalues)) @ vectors.T
    transfers = pair_transfers(strength)
    max_norm = max(np.max(np.abs(np.linalg.eigvalsh(tangent))) for tangent in transfers)
    # Each placement tensored with the other ten maximally mixed ports contributes at most
    # max_norm / 2**10 to the full operator norm, including overlapping placements.
    global_step = global_floor * 2**(board.n - 2) / (len(board.seams) * max_norm)
    response = np.zeros((len(board.seams), 2))
    totals = np.zeros((2, *density.shape))
    for edge_index, edge in enumerate(board.seams):
        if not all(p in positions for p in edge):
            continue  # trace over either endpoint kills a fixed-marginal tangent
        for j, tangent in enumerate(transfers):
            delta = embed(tangent, tuple(positions[p] for p in edge), len(region))
            response[edge_index, j] = np.sum(modular * delta.T)
            totals[j] += delta
    outside = set(range(board.n)) - set(region)
    boundary_triangles = sum(len(set(board.neighbours[left]) & set(board.neighbours[right]) & outside)
                             for left, right in board.seams if left in positions and right in positions)
    direct = np.einsum('ij,kji->k', modular, totals)
    np.testing.assert_allclose(direct, response.sum(axis=0), rtol=0, atol=1e-12)
    return dict(region=region, internal_edges=sum(e in board.seams for e in combinations(region, 2)),
                high_temperature_coefficient=boundary_triangles*(1+2*strength)/4,
                min_probability=float(eigenvalues[0]), positive_step=float(global_step), central=response[board.seams.index((a,b))],
                slopes=direct, total_variation=np.abs(response).sum(axis=0),
                discrepancy=float(direct[1] - direct[0]), regional_tangents=totals, reference=density)


if __name__ == '__main__':
    print(f'Twice the exchange/contrast commutator at the selected words: {conserved_weights_witness()}')
    for strength in (KAPPA, KAPPA/2):
        report = measure(strength)
        print(f'contrast {strength:g}; region {report["region"]}; {report["internal_edges"]} internal seams')
        print(f'central seam slopes {report["central"]}; uniform energy slopes {report["slopes"]}')
        print(f'total absolute weights {report["total_variation"]}; difference {report["discrepancy"]:.9g}')
