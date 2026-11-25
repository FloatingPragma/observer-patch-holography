"""Inclusive clock quadratures and density transport in a stationary binary medium."""

import json
from itertools import combinations
from math import sqrt

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo.candidate import KAPPA
from model.executable_model import cell, federation


def clock_seam():
    host = federation(0)
    return next((a, b) for a, b in cell().seams if host.degree[a] == host.degree[b] == 5)


def spectrum(kappa):
    board = cell()
    blocks = []
    for count in range(board.n + 1):
        sector = BinarySector(board, count)
        energies, states = eigh(sector.hamiltonian(kappa, shifted=False), driver="evd")
        masks = np.array([sum(1 << v for v in word) for word in sector.states])
        blocks.append((energies, states, masks, dict(zip(masks, range(len(masks))))))
    return blocks


def clock_records(blocks, fractions, exposures):
    board = cell()
    a, b = clock_seam()
    bath = sorted(set(range(board.n)) - {a, b})
    records = []
    for t in exposures:
        propagators = [(states * np.exp(-1j * t * energies)) @ states.T
                       for energies, states, _, _ in blocks]
        clock_sums = []
        for count in range(len(bath) + 1):
            backgrounds = [sum(1 << p for p in word) for word in combinations(bath, count)]
            zero, one = blocks[count][3], blocks[count + 1][3]
            vac = [zero[x] for x in backgrounds]
            left = [one[x | (1 << a)] for x in backgrounds]
            right = [one[x | (1 << b)] for x in backgrounds]
            u0 = propagators[count][np.ix_(vac, vac)]
            u = propagators[count + 1]
            aa, ab = u[np.ix_(left, left)], u[np.ix_(left, right)]
            ba, bb = u[np.ix_(right, left)], u[np.ix_(right, right)]
            # Sum inclusive terminal quadratures over every bath word and launch.
            exchange = 2 * np.vdot(aa, ba).imag
            symmetric = -np.vdot(u0, (aa + ba + ab + bb) / 2).imag
            clock_sums.append((exchange, symmetric))

        for f in fractions:
            weights = [f**m * (1 - f)**(len(bath) - m) for m in range(len(bath) + 1)]
            exchange, symmetric = np.asarray(weights) @ clock_sums
            # Difference of two local occupancy readings at b with the launch at a
            # filled or vacant. Every sharp word on the other eleven ports is kept.
            density_transfer = 0.
            for count, (u, (_, _, masks, _)) in enumerate(zip(propagators, blocks)):
                at_a = (masks >> a) & 1
                at_b = (masks >> b) & 1
                bath_counts = count - at_a
                source = (2 * at_a - 1) * f**bath_counts * (1 - f)**(board.n - 1 - bath_counts)
                density_transfer += (np.abs(u[at_b.astype(bool)])**2).sum(axis=0) @ source
            records.append(dict(filling=f, exposure=t, exchange_quadrature=float(exchange),
                                symmetric_quadrature=float(symmetric),
                                adjacent_density_difference=float(density_transfer)))
    return records


def gibbs_read(blocks, beta, kappa=KAPPA):
    """Trace the grand Gibbs density onto one bond at particle-hole symmetry."""
    a, b = clock_seam()
    board = cell()
    mu = kappa * (6 - board.degree[a])
    density = np.zeros((4, 4))
    partition = 0.
    number = 0.
    for n, (energies, states, masks, ids) in enumerate(blocks):
        weights = np.exp(-beta * (energies - mu * n))
        partition += weights.sum()
        number += n * weights.sum()
        occupations = ((masks >> a) & 1) + 2 * ((masks >> b) & 1)
        populations = (states**2) @ weights
        density[np.diag_indices(4)] += np.bincount(occupations, populations, minlength=4)
        left = np.flatnonzero(occupations == 1)
        right = [ids[masks[i] ^ (1 << a) ^ (1 << b)] for i in left]
        density[1, 2] += np.sum((states[left] * states[right]) @ weights)
    density[2, 1] = density[1, 2]
    density /= partition
    probs, vectors = eigh(density)
    modular = (vectors * -np.log(probs)) @ vectors.T
    swap_slope = -2 * modular[1, 2]
    aligned_to_symmetric = ((modular[1, 1] + modular[2, 2] - modular[0, 0] - modular[3, 3]) / 2
                            + modular[1, 2])
    return dict(beta=beta, filling=float(number / (board.n * partition)),
                minimum_regional_probability=float(probs[0]),
                exchange_coherence=float(density[1, 2]), swap_slope=float(swap_slope),
                aligned_to_symmetric_slope=float(aligned_to_symmetric))


def reduced_gibbs_states(blocks, beta, regions, kappa=KAPPA):
    """Reduce the half-filled, regular-cell grand Gibbs state in the specified bit orders."""
    densities = [np.zeros((2**len(region),) * 2) for region in regions]
    mu = kappa * (6 - cell().degree[0])
    for n, (energies, states, masks, _) in enumerate(blocks):
        weights = np.exp(-beta * (energies - mu * n))
        block = (states * weights) @ states.T
        for region, density in zip(regions, densities, strict=True):
            outside = masks & ~sum(1 << site for site in region)
            inside = sum(((masks >> site) & 1) << j for j, site in enumerate(region))
            for background in np.unique(outside):
                indices = np.flatnonzero(outside == background)
                density[np.ix_(inside[indices], inside[indices])] += block[np.ix_(indices, indices)]

    return [density / np.trace(density) for density in densities]


def nested_regions(blocks, beta, kappa=KAPPA):
    """Trace the same Gibbs state onto three graph regions around the test seam."""
    from model.candidate_laws.interval_echo.candidate import native_law as law

    board = cell()
    a, b = clock_seam()
    near_a, near_b = set(board.neighbours[a]), set(board.neighbours[b])
    regions = [(a, b), (a, b, *sorted(near_a & near_b)),
               (a, b, *sorted((near_a | near_b) - {a, b}))]
    densities = reduced_gibbs_states(blocks, beta, regions, kappa)

    out = []
    binary = [law.SRC.index(word) for word in ((0, 0), (1, 0), (0, 1), (1, 1))]
    native_bond = law.seam_term(kappa)[np.ix_(binary, binary)]
    for region, density in zip(regions, densities, strict=True):
        probs, vectors = eigh(density)
        modular = (vectors * -np.log(probs)) @ vectors.T
        # A,B use the low bits. Trace the modular operator against the same
        # normalized spectator identity for each of the two launch tangents.
        spectators = len(density) // 4
        ka = np.trace(modular.reshape(spectators, 4, spectators, 4), axis1=0, axis2=2) / spectators
        diag, exchange_entry = ka.diagonal(), ka[1, 2]
        swap_slope = -2 * exchange_entry
        symmetric_slope = (diag[1] + diag[2] - diag[0] - diag[3]) / 2 + exchange_entry
        residual = ka - beta * native_bond
        two_sites = residual.reshape(2, 2, 2, 2)
        id_half = np.eye(2) / 2
        left_only = np.kron(np.trace(two_sites, axis1=1, axis2=3), id_half)
        right_only = np.kron(id_half, np.trace(two_sites, axis1=0, axis2=2))
        connected = residual - left_only - right_only + np.trace(residual) * np.eye(4) / 4
        out.append(dict(ports=region, minimum_regional_probability=float(probs[0]),
                        swap_per_energy=float(swap_slope), aligned_per_energy=float(symmetric_slope / kappa),
                        fixed_marginal_operator_norm=float(np.max(np.abs(np.linalg.eigvalsh(connected))))))
    return out


def modular_transfer(kappa=KAPPA):
    from model.candidate_laws.interval_echo.candidate import native_law as law

    exchange, contrast = law.local_term_parts()
    h = exchange + kappa * contrast
    v = np.eye(h.shape[0])
    left, right = (v[:, law.SRC.index(word)] for word in ((1, 0), (0, 1)))
    symmetric, odd = (left + right) / sqrt(2), (left - right) / sqrt(2)
    shift = np.outer(odd, odd) - np.outer(symmetric, symmetric)
    energy = np.trace(h @ shift)
    binary = [law.SRC.index(word) for word in ((0, 0), (0, 1), (1, 0), (1, 1))]
    diagonal = np.diag(shift)[binary]
    values = {}
    for f in (.25, .5, .75):
        probs = np.kron([1 - f, f], [1 - f, f])
        values[str(f)] = float(diagonal @ (-np.log(probs)))
    return dict(native_energy_slope=float(energy), modular_slopes=values)


def run():
    blocks = spectrum(KAPPA)
    return dict(seam=clock_seam(), active_kappa=KAPPA,
                clocks=clock_records(blocks, (.25, .5, .75), (1/512, 1/32)),
                exchange_only_control=clock_records(spectrum(0), (.25, .5, .75), (1/512,)),
                thermal_region=[gibbs_read(blocks, beta) for beta in (0., 1/32, 1/16, 1/8)],
                nested_gibbs_at_beta_one_sixteenth=nested_regions(blocks, 1/16),
                regional_transfer=modular_transfer())


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
