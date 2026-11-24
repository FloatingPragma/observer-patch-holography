"""Upper-energy three-wall read and pair-plus-belt breakup on the gated native path."""

import argparse
from collections import defaultdict
from itertools import product
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import candidate, currents
from ...row_13_walls.energy import is_wall, occupied_pairs
from . import belt, composite


def checked_ceiling(h, estimate):
    """A positive Gram factor checks the upper edge without assuming Lanczos found it."""
    bound = float(estimate) + 1e-7
    shifted = -np.array(h)
    shifted.flat[:: len(shifted) + 1] += bound
    factor = np.linalg.cholesky(shifted)
    error = np.linalg.norm(shifted - factor @ factor.T, ord=np.inf)
    assert error < 1e-8
    return float(bound + error + 1e-8)


def bins():
    groups = defaultdict(list)
    for word in range(1024):
        a = occupied_pairs(word, belt.ring(10).seams)
        groups[word.bit_count(), a].append(word)
    return {key: np.array(v, dtype=np.int64) for key, v in groups.items()}


def product_space(groups, count, bridges, max_pairs, cells):
    """All words meeting the cutoff, with the *full* diagonal on each word."""
    edges = [
        (10 * i + a, 10 * i + b) for i in range(cells) for a, b in belt.ring(10).seams
    ]
    edges.extend(bridges)
    words, pairs, degrees = [], [], []
    for keys in product(groups, repeat=cells):
        loads, costs = zip(*keys)
        intrinsic = sum(costs)
        if sum(loads) != count or intrinsic > max_pairs:
            continue
        factors = [groups[k] << (10 * i) for i, k in enumerate(keys)]
        w = factors[0]
        for factor in factors[1:]:
            w = np.bitwise_or.outer(w, factor).ravel()
        a = np.full(len(w), intrinsic)
        endpoints = np.zeros(len(w), dtype=np.int64)
        for i, j in bridges:
            left, right = (w >> i) & 1, (w >> j) & 1
            a += left * right
            endpoints += left + right
        ok = a <= max_pairs
        words.extend(w[ok])
        pairs.extend(a[ok])
        degrees.extend((4 * count + endpoints - 2 * a)[ok])
    order = np.argsort(words)
    return (
        np.asarray(words)[order],
        np.asarray(pairs)[order],
        np.asarray(degrees)[order],
        edges,
    )


def hopping(words, edges):
    rows, cols = [], []
    for a, b in edges:
        src = np.flatnonzero((((words >> a) & 1) == 1) & (((words >> b) & 1) == 0))
        dst = words[src] ^ ((1 << a) | (1 << b))
        pos = np.searchsorted(words, dst)
        ok = (pos < len(words)) & (words[pos.clip(max=len(words) - 1)] == dst)
        rows.append(src[ok])
        cols.append(pos[ok])
    x, y = np.concatenate(rows), np.concatenate(cols)
    half = sparse.csr_matrix((np.ones(len(x)), (x, y)), shape=(len(words), len(words)))
    return half + half.T


def top(h):
    v0 = np.random.default_rng(131212).normal(size=h.shape[0])
    _, vectors = eigsh(h, k=1, which="LA", v0=v0, tol=2e-12)
    vec = vectors[:, 0]
    energy = float(vec @ (h @ vec))
    return energy, vec, float(np.linalg.norm(h @ vec - energy * vec))


def fragment_cap(pair, sectors, spectra, right_spectra, alpha, cutoff):
    """Enclose the full pair edge with product spectra and a softened contact."""
    parts = [(i, pair - i) for i in range(11) if 0 <= pair - i <= 10]
    maximum = max(spectra[i][0][-1] + right_spectra[j][0][-1] for i, j in parts)
    bases = {}
    omitted = -np.inf
    for i, j in parts:
        energies = spectra[i][0][:, None] + right_spectra[j][0]
        selected = energies >= maximum - cutoff
        if (~selected).any():
            omitted = max(omitted, float(energies[~selected].max()))
        x, y = np.where(selected)
        if len(x):
            bases[i] = (x, y, energies[x, y])
    blocks = []
    for i, (x, y, energies) in bases.items():
        j = pair - i
        row = []
        for k, (xx, yy, _) in bases.items():
            right_loads = pair - k

            def at(a, b, indices, other, local):
                matrix = port_matrix(sectors, local, a, b, 0, None)
                return matrix[np.ix_(indices, other)]

            if i == k:
                left, right = at(i, i, x, xx, spectra), at(j, j, y, yy, right_spectra)
                block = (
                    0.5 * (left * (y[:, None] == yy) + right * (x[:, None] == xx))
                    - (1 + alpha) * left * right
                    + np.diag(energies)
                )
            elif abs(i - k) == 1:
                block = (
                    -0.5
                    * at(i, k, x, xx, spectra)
                    * at(j, right_loads, y, yy, right_spectra)
                )
            else:
                block = np.zeros((len(x), len(xx)))
            row.append(block)
        blocks.append(row)
    projected = np.block(blocks)
    upper_p = checked_ceiling(projected, np.linalg.eigvalsh(projected)[-1])
    upper_q = omitted + 1
    mixing = (1 + alpha) / 2
    cap = (upper_p + upper_q + np.hypot(upper_p - upper_q, 2 * mixing)) / 2
    return dict(
        alpha=alpha,
        cutoff=cutoff,
        dimension=sum(len(x) for x, _, _ in bases.values()),
        compressed_edge=upper_p,
        complement_ceiling=upper_q,
        cross_norm_ceiling=mixing,
        upper=float(cap),
    )


def fragment_edges(groups, boards, one, sectors, spectra, capped):
    difference = abs(boards[0][1] - boards[1][0])

    def with_caps(s, ports):
        bias = np.array([sum(p in w for p in ports) / 2 for w in s.states])
        return np.linalg.eigh(s.hamiltonian(candidate.KAPPA) + np.diag(bias))

    centre = [with_caps(s, [difference]) for s in sectors] if capped else spectra
    end = [with_caps(s, [0])[0][-1] for s in sectors] if capped else one
    middle = [with_caps(s, [0, difference])[0][-1] for s in sectors] if capped else one
    # At fixed pair load, both link orientations are isomorphic by belt rotations.
    cache = {}
    out = {}
    for n in (10, 11):
        bare = max(
            end[p[0]] + middle[p[1]] + end[p[2]]
            for p in product(range(11), repeat=3)
            if sum(p) == n
        )
        options = []
        for single in range(11):
            pair = n - single
            if not 0 <= pair <= 20:
                continue
            disjoint = max(
                one[i] + centre[pair - i][0][-1]
                for i in range(11)
                if 0 <= pair - i <= 10
            )
            # The shifted bridge (I-S)/2 - 2*kappa*n_i*n_j has maximum 1.
            ceiling = disjoint + 1 + end[single]
            if ceiling < bare - 1:
                options.append(dict(single=single, upper=ceiling))
                continue
            if pair not in cache:
                words, a, d, edges = product_space(groups, pair, boards[:1], 40, 2)
                assert len(words) == comb(20, pair)
                h = sparse.diags(d / 2 - 2 * candidate.KAPPA * a) - 0.5 * hopping(
                    words, edges
                )
                if capped:
                    h = h + sparse.diags(0.5 * ((words >> boards[1][0]) & 1))
                energy, _, residual = top(h)
                cap = fragment_cap(pair, sectors, spectra, centre, 0.0, 5.0)
                cache[pair] = dict(
                    energy=energy,
                    residual=residual,
                    dimension=len(words),
                    certificate=cap,
                )
            options.append(
                dict(
                    single=single,
                    upper=min(
                        ceiling, end[single] + cache[pair]["certificate"]["upper"]
                    ),
                    numerical_edge=end[single] + cache[pair]["energy"],
                )
            )
        out[n] = dict(
            separated_edge=bare,
            partitions=options,
            pair_plus_single_numerical_edge=max(
                x["numerical_edge"] for x in options if "numerical_edge" in x
            ),
            pair_plus_single_upper=max(x["upper"] for x in options),
        )
        assert (
            max(
                (x["upper"] for x in options if "numerical_edge" not in x),
                default=-np.inf,
            )
            < out[n]["pair_plus_single_numerical_edge"]
        )
    return out, cache


def ceiling_outside(words, a, d, edges, n, m, kappa):
    """Weighted absolute row sums on the complement of the minimum-pair words."""
    shell = words[a == m + 1]
    first_d = d[a == m + 1]
    shell_hops = hopping(shell, edges)
    same = np.asarray(shell_hops.sum(axis=0)).ravel()
    back = np.zeros(len(shell), dtype=int)
    tails = []
    for p, q in edges:
        src = np.flatnonzero(((shell >> p) ^ (shell >> q)) & 1)
        other = shell[src] ^ ((1 << p) | (1 << q))
        idx = np.searchsorted(words, other)
        present = (idx < len(words)) & (words[idx.clip(max=len(words) - 1)] == other)
        back[src[present]] += a[idx[present]] == m
        tails.append(other[~present])
    _, parent_count = np.unique(np.concatenate(tails), return_counts=True)
    max_parents = int(parent_count.max())
    # Tail row <= d-2*kappa*A, plus the increment for neighbours in the first shell.
    tail_plain = 4 * n + 2 * 2 - (2 + 2 * kappa) * (m + 2)

    def bounds(weight):
        near = np.max(
            first_d / 2
            - 2 * kappa * (m + 1)
            + 0.5 * (same + weight * (first_d - back - same))
        )
        far = tail_plain + 0.5 * (1 / weight - 1) * max_parents
        return max(near, far), near, far

    w = min(np.arange(1, 1001) / 1000, key=lambda x: bounds(x)[0])
    cap, near, far = bounds(w)
    return dict(
        first_shell=len(shell),
        tail_maximum_first_shell_neighbours=max_parents,
        tail_weight=float(w),
        first_shell_ceiling=float(near),
        tail_ceiling=float(far),
        complement_ceiling=float(cap),
    )


def port_matrix(sectors, spectra, a, b, port, width):
    """Matrix elements of the occupancy or a unit transfer in the upper belt bands."""
    src, dst = [], []
    for i, word in enumerate(sectors[b].states):
        occ = set(word)
        if (
            (0 == a - b and port in occ)
            or (a == b + 1 and port not in occ)
            or (a == b - 1 and port in occ)
        ):
            src.append(i)
            dst.append(sectors[a].index[tuple(sorted(occ if a == b else occ ^ {port}))])
    columns = slice(-width, None) if width is not None else slice(None)
    return spectra[a][1][dst, columns].T @ spectra[b][1][src, columns]


def joined_ceiling(n, sectors, spectra, bridges):
    """Bound the full joined operator above after dropping the negative bridge potentials."""
    orders = [p for p in product((3, 4), repeat=3) if sum(p) == n]
    peaks = [float(e[-1]) for e, _ in spectra]
    omitted = max(
        sum(peaks[i] for i in p)
        for p in product(range(11), repeat=3)
        if sum(p) == n and p not in orders
    )
    for p in orders:
        for i in range(3):
            omitted = max(
                omitted, sum(peaks[j] for j in p) - peaks[p[i]] + spectra[p[i]][0][-11]
            )
    identity = sparse.eye(10, format="csr")

    def tensor(factors):
        return sparse.kron(
            sparse.kron(sparse.csr_matrix(factors[0]), factors[1], format="csr"),
            factors[2],
            format="csr",
        )

    rows = []
    for left in orders:
        row = []
        for right in orders:
            block = sparse.csr_matrix((1000, 1000))
            if left == right:
                x, y, z = [spectra[i][0][-10:] for i in left]
                block = sparse.diags((x[:, None, None] + y[None, :, None] + z).ravel())
            for a, b in bridges:
                u, v = a // 10, b // 10
                other = 3 - u - v
                if left[other] != right[other]:
                    continue
                factors = [identity] * 3
                first = port_matrix(sectors, spectra, left[u], right[u], a % 10, 10)
                second = port_matrix(sectors, spectra, left[v], right[v], b % 10, 10)
                if left == right:
                    factors[u] = first
                    block = block + 0.5 * tensor(factors)
                    factors[u], factors[v] = identity, second
                    block = block + 0.5 * tensor(factors)
                    factors[u] = first
                    block = block - tensor(factors)
                elif abs(left[u] - right[u]) == 1:
                    factors[u], factors[v] = first, second
                    block = block - 0.5 * tensor(factors)
            row.append(block)
        rows.append(row)
    projected = sparse.bmat(rows, format="csr")
    compressed, _, residual = top(projected)
    compressed = checked_ceiling(projected.toarray(), compressed)
    q = omitted + 2
    estimate = (compressed + q + np.hypot(compressed - q, 2)) / 2
    return dict(
        product_dimension=3000,
        compressed_edge=compressed,
        residual=residual,
        complement_ceiling=q,
        off_block_norm_ceiling=1,
        upper=float(estimate),
    )


def trial_records(words, vector, edges, n):
    """On fresh copies read each seam's both-occupied and exchange-odd projectors."""
    rng = np.random.default_rng(131300 + n)
    shots = 1048576
    counts = []
    means = np.zeros(3)
    for p, q in edges:
        first, second = (words >> p) & 1, (words >> q) & 1
        pair = np.sum(vector[(first & second).astype(bool)] ** 2)
        src = np.flatnonzero((first == 1) & (second == 0))
        dst = words[src] ^ ((1 << p) | (1 << q))
        pos = np.searchsorted(words, dst)
        ok = (pos < len(words)) & (words[pos.clip(max=len(words) - 1)] == dst)
        odd = (
            0.5 * np.sum(vector[first != second] ** 2)
            - vector[src[ok]] @ vector[pos[ok]]
        )
        probabilities = np.array([pair, odd, 1 - pair - odd])
        counts.append(rng.multinomial(shots, probabilities).tolist())
        means += probabilities
    rates = np.sum(counts, axis=0) / shots
    return dict(
        shots_per_seam=shots,
        seed=131300 + n,
        both_odd_other=counts,
        energies={
            str(k): dict(
                mean=float(means[1] - 2 * k * means[0]),
                measured=float(rates[1] - 2 * k * rates[0]),
            )
            for k in (candidate.KAPPA, 0.0, -candidate.KAPPA)
        },
    )


def trial_state(groups, bridges, n):
    """Variational source restricted to the two lowest occupied-pair counts."""
    words, pairs, degree, edges = product_space(groups, n, bridges, n - 8, 3)
    h = sparse.diags(degree / 2 - 2 * candidate.KAPPA * pairs) - 0.5 * hopping(
        words, edges
    )
    return (words, pairs, degree, edges), top(h)


def joined(groups, bridges, n):
    kappa = candidate.KAPPA
    m = n - 9
    (words, a, d, edges), (energy, vec, residual) = trial_state(groups, bridges, n)
    wall = np.array([is_wall(i) for i in range(1024)])
    triple = wall[words & 1023] & wall[(words >> 10) & 1023] & wall[words >> 20]
    assert a.min() == m and triple[a == m].all()
    outside = ceiling_outside(words, a, d, edges, n, m, kappa)
    minimum = words[a == m]
    # Build the actual cut from minimum-pair words to every nonminimal word.
    dest, sources = [], []
    for p, q in edges:
        src = np.flatnonzero(((minimum >> p) ^ (minimum >> q)) & 1)
        other = minimum[src] ^ ((1 << p) | (1 << q))
        idx = np.searchsorted(minimum, other)
        ok = (idx == len(minimum)) | (minimum[idx.clip(max=len(minimum) - 1)] != other)
        dest.append(other[ok])
        sources.append(src[ok])
    dest, sources = np.concatenate(dest), np.concatenate(sources)
    hits = np.unique(dest, return_counts=True)[1]
    src_degree, dst_degree = int(np.bincount(sources).max()), int(hits.max())
    cut2 = src_degree * dst_degree / 4
    gap = energy - outside["complement_ceiling"]
    lower = 1 / (1 + cut2 / gap**2) if gap > 0 else 0
    return dict(
        loads=n,
        dimension=comb(30, n),
        minimum_pairs=m,
        minimum_words=len(minimum),
        trial_dimension=len(words),
        trial_energy=energy,
        trial_residual_inside=residual,
        trial_wall_probability=float(np.dot(vec[triple], vec[triple])),
        cut_degrees=[src_degree, dst_degree],
        **outside,
        top_eigenstate_wall_probability_lower=lower,
        contact_read=trial_records(words, vec, edges, n),
    )


def run():
    geometry, _ = composite.path()
    bridges = geometry["local_bridges"]
    groups = bins()
    sectors = [BinarySector(belt.ring(10), n) for n in range(11)]
    spectra = [np.linalg.eigh(s.hamiltonian(candidate.KAPPA)) for s in sectors]
    one = [float(e[-1]) for e, _ in spectra]
    comparison, pairs = fragment_edges(groups, bridges, one, sectors, spectra, False)
    capped, cpairs = fragment_edges(groups, bridges, one, sectors, spectra, True)
    exchange_only = [float(np.linalg.eigvalsh(s.K)[-1]) for s in sectors]
    records = []
    for n in (10, 11):
        result = joined(groups, bridges, n)
        upper = joined_ceiling(n, sectors, spectra, bridges)
        # Round the filter strictly inside the observed comparison gap.
        threshold = (
            comparison[n]["pair_plus_single_upper"] + result["trial_energy"]
        ) / 2
        threshold = round(threshold, 2)
        gap = threshold - result["complement_ceiling"]
        leak = np.prod(result["cut_degrees"]) / (4 * gap**2)
        assert (
            comparison[n]["pair_plus_single_upper"] < threshold < result["trial_energy"]
        )
        paint = dict(
            threshold=threshold,
            all_times_wall_probability_lower=1 - leak,
            coherent_trial_acceptance_lower=(result["trial_energy"] - threshold)
            / (upper["upper"] - threshold),
            uniform_minimum_word_acceptance_lower=result[
                "top_eigenstate_wall_probability_lower"
            ]
            / result["minimum_words"],
        )
        exchange_edge = max(
            sum(exchange_only[i] for i in assignment)
            for assignment in product(range(11), repeat=3)
            if sum(assignment) == n
        )
        records.append(
            dict(
                **result,
                **comparison[n],
                capped=capped[n],
                joined_ceiling=upper,
                high_energy_filter=paint,
                exchange_only_three_separated_edge=exchange_edge,
            )
        )
    neutral, positive = records
    differences = [
        positive["trial_energy"] - neutral["joined_ceiling"]["upper"],
        positive["joined_ceiling"]["upper"] - neutral["trial_energy"],
    ]
    unit = 6 * candidate.KAPPA
    return currents._portable(
        dict(
            geometry=geometry,
            single_belt_upper_edges=one,
            pair_sectors=pairs,
            capped_pair_sectors=cpairs,
            sectors=records,
            ordering=dict(
                unit_at_rest_on_separate_belt=unit,
                absolute_positive_minus_neutral=[unit + d for d in differences],
                positive_minus_neutral_with_rest_unit=differences,
            ),
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(result)
