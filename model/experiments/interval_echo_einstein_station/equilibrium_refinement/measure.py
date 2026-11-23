"""Native two-double-load regional propagation, monitoring and independent entropy readings."""

import argparse
from itertools import combinations, product
import json
from math import log
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import connected_components

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.clock import Registers
from model.executable_model import federation

SHOTS = 2**22
STEPS = 2**28
EXPOSURE = 0.5
ROOT_SEED = 783542


def seam():
    law = candidate.native_law
    h = law.seam_term(candidate.KAPPA)
    indices = [law.SRC.index(word) for word in product((0, 2), repeat=2)]
    others = np.setdiff1d(np.arange(len(h)), indices)
    if np.any(h[np.ix_(others, indices)]):
        raise ValueError("Sharp two-double-load words leak into other seam sectors")
    return h[np.ix_(indices, indices)]


def apparatus(car, cells):
    inside = frozenset(cells)
    ports = sorted(p for cell in inside for p in range(12*cell, 12*cell + 12))
    border = {}
    contacts = []
    for a, b in car.seams:
        a, b = int(a), int(b)
        if a//12 in inside and b//12 in inside:
            contacts.append((a, b))
        elif (a//12 in inside) != (b//12 in inside):
            entry, vacancy = (a, b) if a//12 in inside else (b, a)
            if entry in border:
                raise ValueError("Enrolled cross-cut endpoint lost its unique address")
            border[entry] = vacancy
    words = list(combinations(ports, 2))
    return ports, contacts, border, words


def capture_labels(car, cells, ports, border, words):
    """Read all sharp basis words through the actual cut and replay each undo."""
    regs = Registers(car, seed=0)
    anchor = min(border)
    # Endpoint transposition routes stay entirely in the enrolled cells.
    routes = {p: regs.route(p, anchor, cells) for p in ports if p not in border}
    if any(abs(rate-.5) > 1e-12 for rate in regs.rates.values()):
        raise ValueError("Native fair contact cannot supply the changed acceptance rule")
    groups, labels, occupancies = {}, [], []
    cut_labels, cut_words = [], {}
    reference_pair = set(words[0])
    for occupied in words:
        regs.loads.fill(0)
        regs.loads[list(occupied)] = 2
        cut = tuple(p for p in sorted(border)
                    if regs.flag(p, border[p], {p//12, border[p]//12}))
        interior = []
        for p, route in routes.items():
            for a, b in route:
                regs.swap(a, b, cells)
            if regs.flag(anchor, border[anchor], {anchor//12, border[anchor]//12}):
                interior.append(p)
            for a, b in reversed(route):
                regs.swap(a, b, cells)
        if tuple(sorted(interior + list(cut))) != occupied:
            raise ValueError("Echo words fail to identify the terminal sharp preparation")
        if tuple(np.flatnonzero(regs.loads)) != occupied:
            raise ValueError("Terminal reads disturb the sharp state")
        label = tuple(interior)
        labels.append(groups.setdefault(label, len(groups)))
        # Every reported word reaches a common sharp reference by a permutation
        # of native contacts. Composing two certified paths implements rejection.
        restoration = []
        for src, dest in zip(sorted(set(occupied)-reference_pair), sorted(reference_pair-set(occupied))):
            for a, b in regs.route(src, dest, cells):
                regs.swap(a, b, cells)
                restoration.append((a, b))
        if set(np.flatnonzero(regs.loads)) != reference_pair:
            raise ValueError("Controlled rejection cannot reprepare the reference pair")
        for a, b in reversed(restoration):
            regs.swap(a, b, cells)
        if tuple(np.flatnonzero(regs.loads)) != occupied:
            raise ValueError("Controlled rejection cannot reverse its native re-preparation")
        occupancies.append(len(cut))
        cut_labels.append(cut_words.setdefault(cut, len(cut_words)))
    k = np.array(occupancies)
    labels = np.array(labels)
    group_k = np.empty(len(groups), dtype=int)
    group_k[labels] = k
    if np.any(group_k[labels] != k):
        raise ValueError("Interior records do not determine boundary occupation")
    return labels, k, group_k, np.array(cut_labels), list(cut_words), {key: int(value) for key, value in regs.gates.items()}


def generator(car, words, contacts):
    local = seam()
    word_index = {w: i for i, w in enumerate(words)}
    h = np.zeros((len(words), len(words)))
    for i, occupied in enumerate(words):
        h[i, i] += 4*candidate.KAPPA*sum(6 - int(car.degree[v]) for v in occupied)
        for a, b in contacts:
            left, right = a in occupied, b in occupied
            local_state = 2*left + right
            h[i, i] += local[local_state, local_state]
            if left != right:
                moved = tuple(sorted(({a, b} ^ set(occupied))))
                h[word_index[moved], i] += local[2*right + left, local_state]
    np.testing.assert_array_equal(h, h.T)
    return h


def sampling_kernel(h):
    energy, modes = eigh(h)
    u = (modes * np.exp(-1j*EXPOSURE*energy)) @ modes.T
    p = abs(u)**2
    np.testing.assert_allclose(p.sum(axis=1), 1, rtol=0, atol=2e-12)
    np.testing.assert_allclose(p, p.T, rtol=0, atol=2e-13)
    if connected_components(p > 1e-12, directed=False)[0] != 1:
        raise ValueError("The monitored kernel is disconnected")
    unmonitored = u[:, 0]
    return p, dict(purity=float(abs(np.vdot(unmonitored, unmonitored))**2),
                   one_interval_sharp_entropy_nats=float(-np.sum(p[0]*np.log(p[0].clip(1e-300)))))


def repeated_endpoints(p, k, tilted):
    # The certified fair gate supplies one independent half-probability bit
    # for each outward increment.
    # P[i,j] is the next sharp word j from i. For the active control a rejected
    # proposal restores the previous sharp word by native endpoint swaps.
    weights = np.exp2(-k.astype(float)) if tilted else np.ones(len(p))
    target = weights/weights.sum()
    if tilted:
        p = p*np.minimum(1, weights[None, :]/weights[:, None])
        p[np.diag_indices_from(p)] += 1 - p.sum(axis=1)
    np.testing.assert_allclose(target @ p, target, rtol=0, atol=1e-12)
    root = np.sqrt(target)
    symmetric = root[:, None]*p/root[None, :]
    np.testing.assert_allclose(symmetric.T, symmetric, rtol=0, atol=1e-12)
    vals, vecs = eigh(symmetric)
    np.testing.assert_allclose(vals[-1], 1, rtol=0, atol=1e-12)
    if np.min(p.diagonal()) <= 0:
        raise ValueError("The sampled chain may be periodic")
    decay = max(abs(vals[0]), abs(vals[-2]))
    if decay >= 1:
        raise ValueError("No numerical spectral convergence")
    indices = [np.flatnonzero(k == count)[0] for count in (0, 2)]
    # The leading eigenvector is known exactly by detailed balance. Only
    # nonstationary components are powered, to avoid amplifying round-off at 1.
    powers = vals[:-1] ** STEPS
    rows = []
    for i in indices:
        row = target + ((vecs[i, :-1]*powers) @ vecs[:, :-1].T)*root/root[i]
        if np.min(row) < 0 or abs(row.sum() - 1) > 1e-12:
            raise ValueError("Finite-time endpoint probabilities are invalid")
        rows.append(row)
    log_tv = (log(.5) + .5*log(1/min(target)-1) + STEPS*log(decay))/log(10)
    return rows, dict(spectral_gap=float(1-decay), worst_start_log10_tv_upper=float(log_tv),
                      input_word_indices=list(map(int, indices)))


def matched_energy_and_load(energies, k, group_ids):
    sectors = np.array([.5, -1., .5])
    tangent = np.array([sectors[j]/np.count_nonzero(k == j) for j in k])
    interior = np.flatnonzero(k == 0)
    low = interior[np.argmin(energies[interior])]
    high = interior[np.argmax(energies[interior])]
    correction = -(tangent @ energies)/(energies[high]-energies[low])
    tangent[high] += correction
    tangent[low] -= correction
    moved = tangent[tangent > 0].sum()
    tangent /= moved
    if max(abs(tangent.sum()), abs(tangent @ energies), abs(tangent @ k)) > 1e-12:
        raise ValueError("Tangent changed one of the enrolled constraints")
    return dict(sector_transfers=(sectors/moved).tolist(), compensation=float(correction/moved),
                low_group=int(group_ids[low]), high_group=int(group_ids[high]),
                low_word=int(low), high_word=int(high))


def meters(counts_bulk, counts_edge, groups, group_k, energy_direction, cut_ids, cut_words, matched):
    nb = len(group_k)
    b = np.bincount(groups, weights=counts_bulk, minlength=nb)
    edge = np.bincount(groups, weights=counts_edge, minlength=nb)
    if np.min(counts_edge) <= 0 or np.min(b) <= 0:
        raise ValueError("Empirical support incomplete")
    information = -np.log(b/b.sum())
    edge_bits = np.bincount(groups, -counts_edge*np.log2(counts_edge/edge[groups]), minlength=nb)/edge
    # The pooled boundary observer has access only to the directly addressed
    # exterior-contact words. Interior labels enter the joint-access check above.
    cut_counts = np.bincount(cut_ids, weights=counts_edge)
    word_k = np.array([len(word) for word in cut_words])
    pooled_bits = []
    for number in range(3):
        distribution = cut_counts[word_k == number]
        distribution = distribution/distribution.sum()
        pooled_bits.append(float(-distribution @ np.log2(distribution)))
    slopes = {}
    for a, bclass in combinations(range(3), 2):
        before, after = group_k == a, group_k == bclass
        slopes[f"{a}_to_{bclass}"] = dict(
            bulk_nats=float(information[after].mean()-information[before].mean()),
            conditional_edge_bits=float(edge_bits[after].mean()-edge_bits[before].mean()),
            cut_only_bits=pooled_bits[bclass]-pooled_bits[a], occupied_cut_endpoints=bclass-a)
    first, last = np.flatnonzero(group_k == 0)[[0, -1]]
    slopes["localized_interior"] = dict(
        bulk_nats=float(information[last]-information[first]),
        conditional_edge_bits=float(edge_bits[last]-edge_bits[first]),
        cut_only_bits=0., occupied_cut_endpoints=0)
    slopes["matched_mean_energy"] = {
        measure: float(sum(energy_direction[j] * slopes[f"0_to_{j}"][measure] for j in (1, 2)))
        for measure in ("bulk_nats", "conditional_edge_bits", "cut_only_bits", "occupied_cut_endpoints")
    }
    t = matched["sector_transfers"]
    slopes["matched_mean_energy_and_cut_load"] = {
        measure: float(sum(t[j]*slopes[f"0_to_{j}"][measure] for j in (1, 2)))
        for measure in ("bulk_nats", "conditional_edge_bits", "cut_only_bits", "occupied_cut_endpoints")
    }
    slopes["matched_mean_energy_and_cut_load"]["bulk_nats"] += matched["compensation"]*(
        information[matched["high_group"]] - information[matched["low_group"]])
    return slopes, cut_counts.astype(int).tolist()


def run(levels=(0, 1, 2)):
    results = []
    for level in levels:
        car = federation(level)
        nearby = sorted({b//12 if a//12 == 0 else a//12
                         for a, b in car.external if 0 in (a//12, b//12)})
        for size in (1, 2, 3):
            cells = [0, *nearby[:size-1]]
            ports, contacts, border, words = apparatus(car, cells)
            labels, k, group_k, cut_ids, cut_words, gates = capture_labels(car, cells, ports, border, words)
            h = generator(car, words, contacts)
            energy_means = np.array([h.diagonal()[k == j].mean() for j in range(3)])
            direction = np.cross(np.ones(3), energy_means)
            direction *= np.sign(direction[np.flatnonzero(direction)[0]])
            direction /= direction[direction > 0].sum()
            np.testing.assert_allclose(direction @ energy_means, 0, atol=1e-12)
            matched = matched_energy_and_load(h.diagonal(), k, labels)
            p, pure = sampling_kernel(h)
            row = dict(board_cells=car.n//12, cells=cells, cut_ports=sorted(border),
                       word_order="lexicographic two-port combinations of sorted enrolled ports",
                       cut_word_order=cut_words,
                       enrolled_words=len(words), interior_words=len(group_k), exhaustive_read_gates=gates,
                       class_mean_native_energy=energy_means.tolist(),
                       mean_energy_matched_transfer=direction.tolist(), matched_energy_and_cut_load=matched,
                       unmonitored=pure, repair={})
            for mode in ("unbiased", "outward_rejection"):
                endpoints, diagnostics = repeated_endpoints(p, k, mode == "outward_rejection")
                acquisitions = {}
                for partition, salt in (("calibration", 31), ("holdout", 83)):
                    counts = [np.random.default_rng([ROOT_SEED, level, size, int(mode != "unbiased"), salt, i])
                              .multinomial(SHOTS, distribution) for i, distribution in enumerate(endpoints)]
                    slopes, cut_counts = meters(*counts, labels, group_k, direction, cut_ids, cut_words, matched)
                    acquisitions[partition] = dict(interior_start_bulk_counts=counts[0].tolist(),
                                                   cut_start_edge_counts=counts[1].tolist(),
                                                   cut_only_counts=cut_counts, slopes=slopes)
                row["repair"][mode] = dict(**diagnostics, acquisitions=acquisitions)
            results.append(row)
    train = results[0]["repair"]["unbiased"]["acquisitions"]["calibration"]["slopes"]["0_to_1"]
    eta = -train["bulk_nats"]/train["conditional_edge_bits"]
    eta_cut = -train["bulk_nats"]/train["cut_only_bits"]
    eta_count = -train["bulk_nats"]
    maxima = {mode: {} for mode in ("unbiased", "outward_rejection")}
    for row in results:
        for mode, experiment in row["repair"].items():
            slopes = experiment["acquisitions"]["holdout"]["slopes"]
            for slope in slopes.values():
                slope["capacity_residual_nats"] = slope["bulk_nats"] + eta*slope["conditional_edge_bits"]
                slope["cut_only_residual_nats"] = slope["bulk_nats"] + eta_cut*slope["cut_only_bits"]
                slope["occupation_residual_nats"] = slope["bulk_nats"] + eta_count*slope["occupied_cut_endpoints"]
                for meter in ("capacity", "cut_only", "occupation"):
                    value = abs(slope[f"{meter}_residual_nats"])
                    maxima[mode][meter] = max(maxima[mode].get(meter, 0), value)
    return dict(native_exposure=EXPOSURE, monitored_steps=STEPS, preparations_per_meter=SHOTS,
                seed=ROOT_SEED, fitted_nats_per_conditional_edge_bit=eta,
                fitted_nats_per_occupied_cut_endpoint=eta_count, fitted_nats_per_cut_only_bit=eta_cut,
                maximum_absolute_held_residual_nats=maxima, observations=results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), separators=(",", ":"), allow_nan=False)+"\n")
