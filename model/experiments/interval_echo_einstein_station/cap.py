"""Prepare a finite regional reference from overlapping local repair likelihoods."""

import argparse
import hashlib
import json
from itertools import combinations, product
from pathlib import Path

import mpmath as mp
import networkx as nx
import numpy as np
from scipy.linalg import eigh, expm, svd

from model.candidate_laws.contrast_memory_spherical_v1 import candidate, law
from model.executable_model import cell, federation
from model.experiments.interval_echo_consensus.modular_standard import standard_count
from model.experiments.interval_echo_consensus.offer_service import (
    exact_stationary_source,
)
from model.experiments.interval_echo_null_balance.experiment import (
    KAPPA,
    hardware,
    states,
)
from .pair import BLOCKS, marginal, rays

HERE = Path(__file__).resolve().parent
BINARY = [(0, 0), (0, 1), (1, 0), (1, 1)]
WORDS = [()] + [(p,) for p in range(12)] + list(combinations(range(12), 2))
WORD_INDEX = {w: i for i, w in enumerate(WORDS)}
SHOTS = 2**20


def cap_geometry(pole):
    car = cell()
    levels, axes = eigh(candidate.response(car))
    triplet = np.isclose(levels, levels[1])
    if sum(triplet) != 3:
        raise ValueError("The native first nonconstant eigenspace is not a triplet.")
    coordinates = axes[:, triplet]
    coordinates /= np.linalg.norm(coordinates, axis=1)[:, None]
    midpoints = np.array([(coordinates[i] + coordinates[j]) / 2 for i, j in car.seams])
    squared = np.sum((midpoints - coordinates[pole]) ** 2, axis=1)
    return (4 - squared) / 4, int(
        np.argmax(np.sum((coordinates - coordinates[pole]) ** 2, axis=1))
    )


def routes(links):
    """Bring each addressed pair to one reusable station in the neighboring cell, then undo."""
    sphere = federation(0)
    neighbour = links[0][1] // 12
    nodes = set(range(12)) | set(range(12 * neighbour, 12 * (neighbour + 1)))
    graph = nx.Graph([e for e in sphere.seams if set(e) <= nodes])
    inside = graph.subgraph(nodes - set(range(12)))
    entry = links[0][1]
    for a, b in sorted(inside.edges):
        for rail in sorted(set(inside[a]) & set(inside[b])):
            for other in sorted(set(inside[rail]) - {a, b}):
                if entry in (rail, other):
                    continue
                # Assign a private vacancy for an addressed single-port phase on each rail.
                vacancies = {
                    v: min(set(inside[v]) - {a, b, rail, other, entry}, default=-1)
                    for v in (a, b, rail)
                }
                if -1 in vacancies.values():
                    continue
                usable = graph.subgraph(nodes - {rail, other, *vacancies.values()})
                if not nx.is_connected(usable) or not nx.is_connected(
                    usable.subgraph(set(usable) - {a})
                ):
                    continue
                paths = []
                for i, j in cell().seams:
                    path_a = nx.shortest_path(usable, i, a)
                    j_after = path_a[path_a.index(j) - 1] if j in path_a else j
                    path_b = nx.shortest_path(
                        usable.subgraph(set(usable) - {a}), j_after, b
                    )
                    paths.append([path_a, path_b])
                for edge, (first, second) in zip(cell().seams, paths, strict=True):
                    positions = {port: port for port in nodes}
                    moves = list(zip(first, first[1:])) + list(zip(second, second[1:]))

                    def apply(moves):
                        for x, y in moves:
                            if not graph.has_edge(x, y):
                                raise ArithmeticError("Routing used a missing seam.")
                            positions[x], positions[y] = positions[y], positions[x]

                    apply(moves)
                    if (positions[a], positions[b]) != edge or any(
                        positions[v] != v for v in (rail, other, *vacancies.values())
                    ):
                        raise ArithmeticError(
                            "The gate would act on a different prepared register."
                        )
                    apply(reversed(moves))
                    if any(port != original for port, original in positions.items()):
                        raise ArithmeticError(
                            "The gate routing did not undo its carrier permutation."
                        )
                return dict(
                    cells=[0, neighbour],
                    data=[a, b],
                    ancilla=[rail, other],
                    phase_vacancies=vacancies,
                    paths=paths,
                )
    raise ValueError("No reachable two-cell station for a coherent local filter.")


class NativeFilter:
    """Binary blank-rail restrictions of native exchanges and addressed signed-memory phases."""

    def __init__(self):
        exchange, contrast = law.local_term_parts()
        term = exchange + KAPPA * contrast
        self.injection = np.eye(law.DIM)[:, [law.SRC.index(w) for w in BINARY]]
        self.binary = self.injection.T @ term @ self.injection
        np.testing.assert_allclose(
            term @ self.injection, self.injection @ self.binary, atol=1e-13
        )
        written = law.echo() @ self.injection
        labels = np.repeat([-9, *(d for _, d in law.LABELS)], law.N_SRC)

        def tag_phase(phase):
            gate = law.echo() @ (np.exp(1j * phase)[:, None] * written)
            reduced = self.injection.T @ gate
            np.testing.assert_allclose(gate, self.injection @ reduced, atol=1e-13)
            return reduced

        # With a vacancy at the second port, signed label +1 detects occupancy at the first.
        site = tag_phase((labels == 1) * (-np.pi / 2))[np.ix_([0, 2], [0, 2])]
        self.cz = np.kron(site, site) @ tag_phase((abs(labels) == 1) * (np.pi / 2))
        np.testing.assert_allclose(self.cz, np.diag([1, 1, 1, -1]), atol=1e-13)
        self.read_basis = expm(-0.5j * np.pi * self.binary) @ np.kron(
            site.conj().T, np.eye(2)
        )
        np.testing.assert_allclose(
            expm(-1j * np.pi * self.binary), np.eye(4)[:, [0, 2, 1, 3]], atol=2e-13
        )

    def success(self, odds, step):
        """Compile survival amplitudes using a coherent multiplexed dual-rail pulse."""
        even, odd = np.array(odds) ** step
        target = np.array([1.0, odd, even, 1.0])
        angles = 2 * np.arccos(target)
        signed = np.array(
            [[(-1) ** (a * x + b * y) for x, y in BINARY] for a, b in BINARY]
        )
        coefficients = signed @ angles / 4
        full = np.eye(8, dtype=complex)
        for (mask_a, mask_b), angle in zip(BINARY, coefficients, strict=True):
            phases = [
                self.cz[2 * a + 1, 2 * a + 1] ** mask_a
                * self.cz[2 * b + 1, 2 * b + 1] ** mask_b
                for a, b in BINARY
            ]
            control = np.diag(np.column_stack((np.ones(4), phases)).ravel())
            # The occupied ancilla has one load. Its left/right exchange has a common phase.
            pulse = expm(-1j * (angle % (2 * np.pi)) * self.binary[1:3, 1:3])
            full = control @ np.kron(np.eye(4), pulse) @ control @ full
        basis = np.kron(self.read_basis, np.eye(2))
        full = basis.conj().T @ full @ basis
        yes, no = full[::2, ::2], full[1::2, ::2]
        np.testing.assert_allclose(
            yes.conj().T @ yes + no.conj().T @ no, np.eye(4), atol=1e-12
        )
        yes = (
            yes / yes[0, 0]
        )  # unobservable phase; the equal-word acceptance is exactly one.
        expected = np.diag([1.0, 0.0, 0.0, 1.0])
        expected[1:3, 1:3] = [
            [(even + odd) / 2, (even - odd) / 2],
            [(even - odd) / 2, (even + odd) / 2],
        ]
        np.testing.assert_allclose(yes, expected, atol=2e-12, rtol=0)
        return yes.real


def seam_indices(edge):
    a, b = edge
    left = [i for i, w in enumerate(WORDS) if a in w and b not in w]
    right = [WORD_INDEX[tuple(sorted((set(WORDS[i]) - {a}) | {b}))] for i in left]
    return np.array(left), np.array(right)


def prepared_charge(filters, steps):
    # Accumulate one complete palindromic sweep before repeating it. No regional
    # energy or its eigensystem is used to prepare this reference.
    forward = np.eye(len(WORDS))
    for (left, right), local in filters:
        forward[np.r_[left, right]] = np.vstack(
            (
                local[2, 2] * forward[left] + local[2, 1] * forward[right],
                local[1, 2] * forward[left] + local[1, 1] * forward[right],
            )
        )
    amplitude = np.linalg.matrix_power(forward.T @ forward, steps)
    square_root_norm = np.linalg.norm(amplitude, "fro")
    state = amplitude @ amplitude.T / square_root_norm**2
    charge = np.zeros_like(state)
    for sector in BLOCKS:
        vectors, values, _ = svd(amplitude[sector, sector], lapack_driver="gesvd")
        charge[sector, sector] = (
            vectors * (-2 * np.log(values / square_root_norm))
        ) @ vectors.T
    # A scalar vanishes on every normalized difference and centering improves sampled variance.
    charge -= charge[0, 0] * np.eye(len(WORDS))
    weak_square = (forward.T @ forward) @ (forward.T @ forward)
    return (
        state,
        charge,
        square_root_norm**2 / len(WORDS),
        weak_square / weak_square.trace(),
        weak_square.trace() / len(WORDS),
    )


def native_energy(indices, weights):
    # This material meter reads mixed occupations and the antisymmetric click.
    mixed_word = np.zeros(len(WORDS))
    antisymmetric = []
    for (left, right), w in zip(indices, weights, strict=True):
        mixed_word[np.r_[left, right]] += KAPPA * w
        p = np.zeros((len(WORDS), len(WORDS)))
        p[left, left] = p[right, right] = 0.5
        p[left, right] = p[right, left] = -0.5
        antisymmetric.append(p)
    return (
        np.diag(mixed_word)
        + sum(w * p for w, p in zip(weights, antisymmetric, strict=True)),
        mixed_word,
        antisymmetric,
    )


def flights():
    sector, cells, links, _, internal, bonds = hardware()
    eigen, v = eigh((sum(internal) + sum(b.h for b in bonds)).toarray())
    prepared = {}
    for port, sign, wave in states(sector):
        coefficients = v.T @ wave
        for time in (1, 4, 16):
            state = marginal(sector, v @ (np.exp(-1j * time * eigen) * coefficients))
            read, _ = rays(state)
            np.testing.assert_allclose(read, state, atol=1e-13)
            prepared[port, sign, time] = read
    return cells, links, prepared


def read_charge(rho, observable, seed, seam_meter=None, shots=SHOTS):
    rng = np.random.default_rng(seed)
    diagonal = np.maximum(rho.diagonal().real, 0)
    diagonal /= diagonal.sum()
    word_weights = observable.diagonal().real if seam_meter is None else seam_meter[0]
    counts = rng.multinomial(shots, grid_distribution(diagonal))
    estimate = word_weights @ counts / shots
    variance = diagonal @ (word_weights - (word_weights @ diagonal)) ** 2

    if seam_meter is not None:
        _, weights, projectors = seam_meter
        rates = np.array([np.trace(p @ rho).real for p in projectors])
        clicks = rng.binomial(
            shots, grid_distribution(np.column_stack((rates, 1 - rates)))[:, 0]
        )
        return float(estimate + weights @ clicks / shots), float(
            (variance + (weights**2) @ (rates * (1 - rates))) / shots
        )

    for sector in BLOCKS[1:]:
        block = rho[sector, sector]
        i, j = np.triu_indices(len(block), 1)
        ray_base = (block.diagonal().real[i] + block.diagonal().real[j]) / 2
        entries = observable[sector, sector][i, j]
        # The real charge uses a direct word setting and the real quadratures.
        # Imaginary coefficients vanish for the palindromic, time-reversal-even reference.
        for offset, coef in (
            (block[i, j].real, entries.real),
            (-block[i, j].imag, -entries.imag),
        ):
            if not np.any(coef):
                continue
            probs = np.column_stack(
                [ray_base + offset, ray_base - offset, 1 - 2 * ray_base]
            )
            probs = np.maximum(probs, 0)
            counts = rng.multinomial(shots, grid_distribution(probs))
            estimate += coef @ (counts[:, 0] - counts[:, 1]) / shots
            variance += (coef**2) @ (2 * ray_base - 4 * offset**2)
    return float(estimate), float(variance / shots)


def grid_distribution(probs, scale=2**31):
    """Use a binary acquisition grid; keep all padding in the largest bin."""
    weights = np.rint(np.maximum(probs, 0) * scale).astype(np.int64)
    ids = np.argmax(weights, axis=-1)
    difference = scale - weights.sum(axis=-1)
    np.put_along_axis(
        weights,
        np.expand_dims(ids, -1),
        np.take_along_axis(weights, np.expand_dims(ids, -1), axis=-1)
        + np.expand_dims(difference, -1),
        axis=-1,
    )
    return weights / scale


def weak_tomography(reference, delta, sweeps, beta, acceptance):
    """Read the weak state and use its repeated positive success operator on a held tangent."""
    shots, batches, seed = 2**29, 16, 150721
    grids = [(None, grid_distribution(reference.diagonal()))]
    for sector in BLOCKS[1:]:
        block = reference[sector, sector]
        i, j = np.triu_indices(len(block), 1)
        base = (block.diagonal()[i] + block.diagonal()[j]) / 2
        for phase, offset in [(1.0, block[i, j]), (-1j, np.zeros(len(i)))]:
            probs = grid_distribution(
                np.column_stack([base + offset, base - offset, 1 - 2 * base])
            )
            grids.append(((sector, i, j, phase), probs))

    def extrapolated(rho):
        observable = np.zeros_like(rho, dtype=complex)
        minimum = 1.0
        for sector in BLOCKS:
            vals, vecs = eigh(rho[sector, sector])
            minimum = min(minimum, float(vals.min()))
            if vals.min() <= 0:
                raise ValueError(
                    "The weak-reference tomography has a nonpositive sector."
                )
            observable[sector, sector] = (
                vecs * (-(sweeps / beta) * np.log(vals))
            ) @ vecs.conj().T
        observable -= observable[0, 0] * np.eye(len(WORDS))
        charge = float(np.trace(observable @ delta).real)
        return charge, minimum, observable

    checksum = hashlib.sha256()
    estimates = []
    batch_states = []
    pooled = np.zeros_like(reference, dtype=complex)
    children = np.random.SeedSequence(seed).spawn(batches)
    for child in children:
        rng = np.random.Generator(np.random.PCG64(child))
        rho = np.zeros_like(pooled)
        for settings, probs in grids:
            counts = rng.multinomial(shots, probs)
            checksum.update(counts.astype("<u8").tobytes())
            if settings is None:
                np.fill_diagonal(rho, counts / shots)
                continue
            sector, i, j, phase = settings
            part = rho[sector, sector]
            part[i, j] += phase * (counts[:, 0] - counts[:, 1]) / (2 * shots)
            part[j, i] = part[i, j].conj()
        pooled += rho / batches
        batch_states.append(rho)
        estimates.append(extrapolated(rho)[:2])
    results = np.array(estimates)
    value, minimum, observable = extrapolated(pooled)
    ideal, weak_minimum, _ = extrapolated(reference)
    # Normalize each one-sweep reference on its successfully written branch.
    # The budget below counts all weak reference shots entering the tomography.
    return (
        dict(
            seed=seed,
            batches=batches,
            settings=4423,
            accepted_shots_per_batch_setting=shots,
            total_accepted_preparations=shots * batches * 4423,
            expected_input_preparations=shots * batches * 4423 / acceptance,
            one_sweep_acceptance=acceptance,
            ideal_minimum=weak_minimum,
            pooled_minimum=minimum,
            minimum_over_batches=float(results[:, 1].min()),
            exact_signed_energy=ideal,
            pooled_signed_energy=value,
            mean_raw_batch_signed_energy=float(results[:, 0].mean()),
            estimated_standard_error=float(
                np.std(results[:, 0], ddof=1) / np.sqrt(batches)
            ),
            batch_energies=results[:, 0].tolist(),
            packed_count_sha256=checksum.hexdigest(),
        ),
        observable,
        batch_states,
        pooled,
    )


def cut_entropy_test(reference, charge):
    """Measure an entropy variation at fixed total load and, for the interior exchange, fixed gateway state."""
    carrier = federation(0)
    gateway_ports = sorted(p for edge in carrier.external for p in edge if p < 12)
    first = cell().seams[0]
    protected = next(e for e in cell().seams if not set(e).intersection(gateway_ports))
    classes = {}
    for index, word in enumerate(WORDS):
        visible = tuple(p for p in word if p in gateway_ports)
        hidden = tuple(p for p in word if p not in gateway_ports)
        bits = sum(1 << gateway_ports.index(p) for p in visible)
        classes.setdefault(hidden, []).append((bits, index))

    def cut(state):
        rho = np.zeros((2 ** len(gateway_ports),) * 2, dtype=state.dtype)
        for group in classes.values():
            bits, indices = zip(*group, strict=True)
            rho[np.ix_(bits, bits)] += state[np.ix_(indices, indices)]
        return rho

    def entropy(state):
        vals = np.concatenate(
            [np.linalg.eigvalsh(state[part, part]) for part in BLOCKS]
        )
        if min(vals) <= 0:
            raise ValueError(
                "Entropy slope requires a faithful reference on the enrolled support."
            )
        return float(-np.dot(vals, np.log(vals)))

    base = entropy(reference)
    zero_cut = cut(reference)
    trials = []
    for edge in (first, protected, None):
        gate = np.eye(len(WORDS))
        if edge is None:
            gate = np.diag([(-1) ** len(w) for w in WORDS])
        else:
            left, right = seam_indices(edge)
            swaps = np.arange(len(WORDS))
            swaps[left], swaps[right] = right, left
            gate = gate[swaps]
        rotated = gate @ reference @ gate.T
        delta = rotated - reference
        slope = float(np.trace(charge @ delta).real)
        relative = -entropy(rotated)
        for part in BLOCKS:
            vals, axes = eigh(reference[part, part])
            relative -= np.dot(
                np.log(vals), np.diag(axes.T @ rotated[part, part] @ axes)
            ).real
        cut_distance = 0.5 * np.abs(np.linalg.eigvalsh(cut(rotated) - zero_cut)).sum()
        differences = [
            (entropy(reference + delta * step) - base) / step
            for step in (1 / 512, 1 / 4096)
        ]
        sector_difference = [float(np.trace(delta[part, part]).real) for part in BLOCKS]
        np.testing.assert_allclose(slope, relative, atol=4e-13, rtol=0)
        np.testing.assert_allclose(sector_difference, 0, atol=2e-13, rtol=0)
        if edge == protected and (slope <= 1e-6 or cut_distance > 1e-13):
            raise ArithmeticError("Interior gate failed the fixed-cut test.")
        trials.append(
            dict(
                edge=edge,
                tag="load_phase" if edge is None else "exchange",
                slope_nats=slope,
                relative_entropy_nats=float(relative),
                finite_entropy_slopes=differences,
                sector_changes=sector_difference,
                gateway_trace_distance=cut_distance,
            )
        )
    return dict(
        gateway_ports=gateway_ports,
        reference_entropy_nats=base,
        steps=[1 / 512, 1 / 4096],
        trials=trials,
    )


def collar_remainder(reference, support):
    """Recovery is on the full three-port algebra. Only the logarithmic reading is supported."""

    def spectral(state, function, allow_zero=False):
        vals, axes = eigh(state)
        if vals.min() < -2e-14 or (not allow_zero and vals.min() <= 0):
            raise ArithmeticError(
                "A nonpositive spectrum escaped the declared recovery support."
            )
        return (axes * function(np.maximum(vals, 0))) @ axes.conj().T

    ab = np.trace(reference.reshape(4, 2, 4, 2), axis1=1, axis2=3)
    bc = np.trace(reference.reshape(2, 4, 2, 4), axis1=0, axis2=2)
    b = np.trace(bc.reshape(2, 2, 2, 2), axis1=1, axis2=3)
    operator = spectral(bc, np.sqrt, True) @ np.kron(
        spectral(b, lambda x: 1 / np.sqrt(x)), np.eye(2)
    )
    # For every input on B, R(X)=W (X ⊗ I_C) W† is CP and TP.
    np.testing.assert_allclose(
        np.trace((operator.conj().T @ operator).reshape(2, 2, 2, 2), axis1=1, axis2=3),
        np.eye(2),
        atol=3e-12,
        rtol=0,
    )
    dressing = np.kron(np.eye(2), operator)
    recovered = dressing @ np.kron(ab, np.eye(2)) @ dressing.conj().T
    np.testing.assert_allclose(recovered.trace(), 1, atol=2e-12, rtol=0)
    outside = np.setdiff1d(np.arange(8), support)
    cross = np.linalg.norm(recovered[np.ix_(support, outside)])
    if cross > 2e-13:
        raise ArithmeticError(
            "The physical support no longer reduces the recovered state."
        )
    floor = min(
        np.linalg.eigvalsh(x[np.ix_(support, support)])[0]
        for x in (reference, recovered)
    )
    e = np.zeros_like(reference)
    e[np.ix_(support, support)] = spectral(
        recovered[np.ix_(support, support)], np.log
    ) - spectral(reference[np.ix_(support, support)], np.log)
    eta = np.abs(np.linalg.eigvalsh(recovered - reference)).sum() / 2
    return e, dict(
        eta=float(eta),
        outside_weight=float(recovered[outside, outside].real.sum()),
        support_cross_norm=float(cross),
        positive_floor=float(floor),
        operator_norm=float(np.linalg.norm(e, 2)),
        relative_entropy=float(-np.trace(e @ reference).real),
    )


def three_port_survey(state, shots=None, rng=None, checksum=None):
    """Read a complete three-port word or a balanced pair in all spectator fibres."""
    grid = 2**43

    def frequencies(probs):
        calibrated = grid_distribution(probs, grid)
        if shots is None:
            return calibrated
        counts = rng.multinomial(shots, calibrated)
        checksum.update(counts.astype("<u8").tobytes())
        return counts / shots

    word = frequencies(state.diagonal().real)
    fitted = np.diag(word).astype(complex)
    for sector in ((1, 2, 4), (3, 5, 6)):
        for i, j in combinations(sector, 2):
            base = (state[i, i] + state[j, j]).real / 2
            for offset, component in ((state[i, j].real, 1), (state[i, j].imag, 1j)):
                plus, minus, _quiet = frequencies(
                    np.array([base + offset, base - offset, 1 - 2 * base])
                )
                fitted[i, j] += component * ((plus - minus) / 2)
            fitted[j, i] = fitted[i, j].conj()
    return fitted


def direct_collar_survey(ref, physical, reductions, acceptance):
    from scipy.stats import t

    batches = 16
    matter_shots = 2**29
    axes = list(reductions)
    contrasts = [reductions[axis][0] - reductions[axis][1] for axis in axes]
    critical = float(t.ppf(1 - 0.05 / (2 * len(axes)), batches - 1))
    gridded, _ = collar_remainder(three_port_survey(ref), physical)
    exact, _ = collar_remainder(ref, physical)
    grid_bias = [
        float(
            np.trace(
                gridded
                @ (
                    three_port_survey(reductions[axis][0])
                    - three_port_survey(reductions[axis][1])
                )
                - exact @ contrast
            ).real
        )
        for axis, contrast in zip(axes, contrasts, strict=True)
    ]
    rows = []
    for dose, trials in enumerate((2**35, 2**47)):
        checksum = hashlib.sha256()
        ref_seeds = np.random.SeedSequence(150751 + dose).spawn(batches)
        read = [
            three_port_survey(ref, trials, np.random.default_rng(seed), checksum)
            for seed in ref_seeds
        ]
        pooled = np.mean(read, axis=0)
        floors = [np.linalg.eigvalsh(x[np.ix_(physical, physical)])[0] for x in read]
        if min(floors) <= 0:
            raise ArithmeticError(
                "Direct collar tomography lacks positive support at this dose."
            )
        instruments = [collar_remainder(r, physical)[0] for r in read]
        observable, diagnostics = collar_remainder(pooled, physical)
        projected = np.array(
            [[float(np.trace(k @ d).real) for d in contrasts] for k in instruments]
        )
        matter_seeds = iter(np.random.SeedSequence(150761 + dose).spawn(2 * len(axes)))
        estimates = []
        for index, axis in enumerate(axes):
            states = reductions[axis]
            sample = [
                three_port_survey(
                    r, matter_shots, np.random.default_rng(next(matter_seeds)), checksum
                )
                for r in states
            ]
            value = float(np.trace(observable @ (sample[0] - sample[1])).real)
            pooled_value = float(np.trace(observable @ contrasts[index]).real)
            # Direct words and twelve independently drawn ray settings have separate variances.
            variance = 0.0
            for raw in sample:
                diag = raw.diagonal().real
                coeff = observable.diagonal().real
                variance += diag @ (coeff - diag @ coeff) ** 2
                for sector in ((1, 2, 4), (3, 5, 6)):
                    for i, j in combinations(sector, 2):
                        mass = diag[i] + diag[j]
                        for weight, offset in (
                            (observable[i, j].real, raw[i, j].real),
                            (observable[i, j].imag, raw[i, j].imag),
                        ):
                            variance += weight**2 * (mass - 4 * offset**2)
            batch_readings = [
                float(np.trace(k @ (sample[0] - sample[1])).real) for k in instruments
            ]
            reference_se = float(np.std(batch_readings, ddof=1) / np.sqrt(batches))
            matter_se = float(np.sqrt(variance / matter_shots))
            half_width = float(critical * np.hypot(reference_se, matter_se))
            estimates.append(
                dict(
                    port=axis[0],
                    time=axis[1],
                    signed_nats=value,
                    pooled_reference_nats=pooled_value,
                    mean_batch_reference_nats=float(np.mean(projected[:, index])),
                    reference_standard_error=reference_se,
                    matter_standard_error=matter_se,
                    simultaneous_half_width=half_width,
                )
            )
        rows.append(
            dict(
                reference_trials_per_batch_setting=trials,
                reference_batches=batches,
                accepted_reference_preparations=trials * batches * 13,
                expected_input_preparations=trials * batches * 13 / acceptance,
                matter_trials_per_setting=matter_shots,
                reference_seed=150751 + dose,
                matter_seed=150761 + dose,
                packed_count_sha256=checksum.hexdigest(),
                minimum_batch_floor=float(min(floors)),
                recovery=diagnostics,
                estimates=estimates,
            )
        )
    return dict(
        settings=13,
        probability_grid=2**43,
        simultaneous_critical_t=critical,
        max_grid_charge_bias_nats=max(abs(x) for x in grid_bias),
        doses=rows,
    )


def collar_measurement(
    reference, matter, weak_batches, weak_pooled, sweeps, pole, acceptance
):
    gateways = {i for seam in federation(0).external for i in seam if i < 12}
    endports = sorted(set(cell().neighbours[pole]) - gateways)
    a, c = next(
        (a, c) for a, c in combinations(endports, 2) if c not in cell().neighbours[a]
    )
    path = (a, pole, c)
    groups = {}
    for i, word in enumerate(WORDS):
        exterior = tuple(v for v in word if v not in path)
        visible = sum(1 << (2 - j) for j, v in enumerate(path) if v in word)
        groups.setdefault(exterior, []).append((visible, i))
    pairs = [tuple(zip(*group, strict=True)) for group in groups.values()]

    def marginal3(rho):
        local = np.zeros((8, 8), dtype=rho.dtype)
        for visible, full in pairs:
            local[np.ix_(visible, visible)] += rho[np.ix_(full, full)]
        return local

    def readback(local):
        complete = np.zeros((len(WORDS), len(WORDS)), dtype=local.dtype)
        for visible, full in pairs:
            complete[np.ix_(full, full)] = local[np.ix_(visible, visible)]
        return complete

    ref = marginal3(reference)
    # Each local ray is a direct sum of disjoint, addressable regional two-word rays.
    # Coarsen the final word into plus/minus/quiet on every spectator fibre.
    inputs = [reference, *matter.values()]
    marginals = [marginal3(rho) for rho in inputs]
    n_fibres = {}
    for sector in ((1, 2, 4), (3, 5, 6)):
        for i, j in combinations(sector, 2):
            for phase in (1, -1j):
                lifted = []
                for sign in (1, -1):
                    ray = np.zeros(8, complex)
                    ray[[i, j]] = [1 / np.sqrt(2), sign * phase / np.sqrt(2)]
                    effect = readback(np.outer(ray, ray.conj()))
                    lifted.append(effect)
                    np.testing.assert_allclose(
                        effect @ effect, effect, atol=2e-15, rtol=0
                    )
                    np.testing.assert_allclose(
                        [np.trace(effect @ rho) for rho in inputs],
                        [ray.conj() @ rho @ ray for rho in marginals],
                        atol=2e-13,
                        rtol=0,
                    )
                np.testing.assert_allclose(lifted[0] @ lifted[1], 0, atol=2e-15, rtol=0)
                n_fibres[f"{i}-{j}"] = int(round(lifted[0].trace().real))
    physical = np.arange(7)
    remainder, summary = collar_remainder(ref, physical)
    classical, diag_summary = collar_remainder(np.diag(ref.diagonal()), physical)
    ab = np.trace(ref.reshape(4, 2, 4, 2), axis1=1, axis2=3)
    _, quiet = collar_remainder(
        np.kron(ab, np.diag([1.0, 0.0])), np.array([0, 2, 4, 6])
    )
    if quiet["eta"] > 2e-12 or quiet["operator_norm"] > 1e-9:
        raise ArithmeticError(
            "The independent vacancy recovery failed its switch-off control."
        )
    tangents = []
    reductions = {}
    for port in (0, 9):
        for time in (1, 4, 16):
            states = [marginal3(matter[port, sign, time]) for sign in (1, -1)]
            reductions[port, time] = states
            delta = states[0] - states[1]
            weight = np.abs(np.linalg.eigvalsh(delta)).sum()
            signed = float(np.trace(remainder @ delta).real)
            bound = weight * 2 * summary["eta"] / summary["positive_floor"]
            if abs(signed) > bound + 1e-13:
                raise ArithmeticError(
                    "The supported operator logarithm exceeded its trace-distance bound."
                )
            tangents.append(
                dict(
                    port=port,
                    time=time,
                    signed_nats=signed,
                    dephased_reference_nats=float(np.trace(classical @ delta).real),
                    tangent_trace_norm=float(weight),
                    eta_floor_bound_nats=float(bound),
                )
            )

    def weak_read(weak):
        fitted = np.linalg.matrix_power(weak, sweeps)
        fitted /= fitted.trace()
        e, data = collar_remainder(marginal3(fitted), physical)
        held = float(np.trace(e @ marginal3(matter[9, 1, 4] - matter[9, -1, 4])).real)
        return held, data["eta"], e

    sample = np.array([weak_read(w)[:2] for w in weak_batches])
    reference_estimate, fitted_eta, fitted = weak_read(weak_pooled)
    reference_se = float(np.std(sample[:, 0], ddof=1) / np.sqrt(len(sample)))
    held = [
        read_charge(matter[9, sign, 4], readback(fitted), 150741 + i, shots=2**24)
        for i, sign in enumerate((1, -1))
    ]
    sampled = held[0][0] - held[1][0]
    matter_se = float(np.sqrt(held[0][1] + held[1][1]))
    return dict(
        path=path,
        number_of_spectator_fibres=n_fibres,
        recovery=summary,
        classical_recovery=diag_summary,
        independent_vacancy=quiet,
        tangents=tangents,
        direct_survey=direct_collar_survey(ref, physical, reductions, acceptance),
        held_survey=dict(
            fitted_eta=fitted_eta,
            eta_batch_standard_error=float(
                np.std(sample[:, 1], ddof=1) / np.sqrt(len(sample))
            ),
            batch_reference_nats=sample[:, 0].tolist(),
            pooled_reference_nats=reference_estimate,
            reference_standard_error=reference_se,
            reference_and_matter_nats=sampled,
            matter_standard_error=matter_se,
            combined_standard_error=float(np.hypot(reference_se, matter_se)),
            shots_per_matter_setting=2**24,
        ),
    )


def wedge(beta):
    """Outward-rounded verification allowing an arbitrary constant in each load sector."""
    iv = mp.iv
    original = iv.dps
    try:
        iv.dps = 35
        b = iv.ln(2) / 2
        u = iv.mpf(17) / 2
        d = iv.sqrt(u * u / 4 + iv.mpf(1) / 2)
        scale = iv.exp(-3 * b * u / 2)
        plus, minus = iv.exp(b * d), iv.exp(-b * d)
        cosh, sinh = (plus + minus) / 2, (plus - minus) / 2
        leaf = (scale * (cosh + u * sinh / (2 * d)) + iv.exp(-b * u)) / 2
        middle = scale * (cosh - u * sinh / (2 * d))
        coherence = scale * sinh / (2 * d)
        ratio = (leaf + middle - 2 * coherence) / (leaf + middle + 2 * coherence)
        if not (
            ratio > iv.mpf("0.8089986126552") and ratio < iv.mpf("0.8089986126554")
        ):
            raise ArithmeticError("Wedge ratio escaped its certified interval.")
    finally:
        iv.dps = original
    native = NativeFilter().binary
    h = np.kron(native, np.eye(2)) + np.kron(np.eye(2), native)
    state = expm(-beta * h)
    loads = np.array([sum(w) for w in product((0, 1), repeat=3)])
    ratios = []
    for n in (1, 2):
        restricted = state * np.outer(loads == n, loads == n)
        local = np.trace(restricted.reshape(4, 2, 4, 2), axis1=1, axis2=3)
        diag = local[1, 1] + local[2, 2]
        interference = 2 * local[1, 2]
        ratios.append(float((diag - interference) / (diag + interference)))
    np.testing.assert_allclose(ratios, [float(ratio.mid)] * 2, atol=2e-13, rtol=0)
    return dict(
        sector1_and_2=ratios,
        certified_interval=[0.8089986126552, 0.8089986126554],
        isolated_ratio=float(np.exp(-beta)),
        empty_and_full_mixed_readings=[0, 0],
    )


def measure():
    cells, links, matter = flights()
    pole = links[0][0]
    weights, opposite = cap_geometry(pole)
    reversed_weights, _ = cap_geometry(opposite)
    indices = [seam_indices(edge) for edge in cell().seams]
    energy, word, projectors = native_energy(indices, weights)
    flat, _, _ = native_energy(indices, np.ones(30))
    reversed_energy, _, _ = native_energy(indices, reversed_weights)
    s = exact_stationary_source(1, standard_count())
    zero, one, odd, _, _ = s["mode_probabilities"]
    odds = np.array([one, odd]) / zero
    beta = np.log(one / odd)
    native = NativeFilter()
    rows = []
    for n in (1, 4, 16, 64):
        filters = [
            (pair, native.success(odds, w / (4 * n)))
            for pair, w in zip(indices, weights, strict=True)
        ]
        reference, k, acceptance, weak, weak_acceptance = prepared_charge(filters, n)
        rows.append(
            dict(
                sweeps=n,
                total_local_filters=60 * n,
                acceptance=acceptance,
                min_reference_eigenvalue=float(np.linalg.eigvalsh(reference).min()),
                max_operator_energy_error=float(np.linalg.norm(k / beta - energy, 2)),
            )
        )
    tangents = []

    def trace(operator, x):
        return float(np.trace(operator @ x).real)

    for p in (0, 9):
        for t in (1, 4, 16):
            delta = matter[p, 1, t] - matter[p, -1, t]
            tangents.append(
                dict(
                    port=p,
                    time=t,
                    read_energy=trace(energy, delta),
                    modular_energy=trace(k / beta, delta),
                    first_law_nats=trace(k, delta),
                    flat_meter=trace(flat, delta),
                    reversed_pole_meter=trace(reversed_energy, delta),
                    dephased_charge=trace(np.diag(np.diag(k)) / beta, delta),
                )
            )
    sample = {}
    for label, op, seed, seam_meter in [
        ("modular", k / beta, 150710, None),
        ("energy", energy, 150712, (word, weights, projectors)),
    ]:
        vals = [
            read_charge(matter[9, sign, 4], op, seed + i, seam_meter)
            for i, sign in enumerate([1, -1])
        ]
        sample[label] = dict(
            signed=vals[0][0] - vals[1][0], se=float(np.sqrt(vals[0][1] + vals[1][1]))
        )
    recovery, sampled_charge, weak_batches, weak_pooled = weak_tomography(
        weak, matter[9, 1, 4] - matter[9, -1, 4], 64, beta, weak_acceptance
    )
    np.testing.assert_allclose(
        recovery["exact_signed_energy"],
        trace(k / beta, matter[9, 1, 4] - matter[9, -1, 4]),
        atol=2e-9,
        rtol=0,
    )
    sampled = [
        read_charge(matter[9, sign, 4], sampled_charge, 150725 + i)
        for i, sign in enumerate((1, -1))
    ]
    recovery["sampled_reference_and_matter"] = sampled[0][0] - sampled[1][0]
    recovery["sampled_matter_se"] = float(np.sqrt(sampled[0][1] + sampled[1][1]))
    recovery["difference_from_sampled_material_meter"] = (
        recovery["sampled_reference_and_matter"] - sample["energy"]["signed"]
    )
    recovery["combined_se_including_material_meter"] = float(
        np.sqrt(
            recovery["estimated_standard_error"] ** 2
            + recovery["sampled_matter_se"] ** 2
            + sample["energy"]["se"] ** 2
        )
    )
    sample["residual"] = sample["modular"]["signed"] - sample["energy"]["signed"]
    sample["combined_se"] = float(
        np.hypot(sample["modular"]["se"], sample["energy"]["se"])
    )
    if not all(
        rows[i + 1]["max_operator_energy_error"] < rows[i]["max_operator_energy_error"]
        for i in range(3)
    ):
        raise ValueError("The cap charge failed the sweep convergence check.")
    return dict(
        pole=pole,
        opposite_pole=opposite,
        serviced_cells=cells,
        gateways=links,
        seam_weights=weights.tolist(),
        repair_odds=odds.tolist(),
        exchange_beta=beta,
        wedge=wedge(beta),
        three_port_recovery=collar_measurement(
            reference,
            matter,
            weak_batches,
            weak_pooled,
            64,
            pole,
            rows[-1]["acceptance"],
        ),
        fixed_cut_entropy=cut_entropy_test(reference, k),
        apparatus=routes(links),
        sweeps=rows,
        tangents=tangents,
        held_count_survey=sample,
        weak_reference_survey=recovery,
        shots_per_setting=SHOTS,
    )


def compare(actual, expected):
    if isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError("Changed cap evidence fields.")
        for key in expected:
            compare(actual[key], expected[key])
    elif isinstance(expected, list):
        for a, b in zip(actual, expected, strict=True):
            compare(a, b)
    elif isinstance(expected, float):
        np.testing.assert_allclose(actual, expected, rtol=2e-8, atol=2e-10)
    elif actual != expected:
        raise ValueError((actual, expected))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(measure(), allow_nan=False))
    path = HERE / "cap.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        compare(result, json.loads(path.read_text()))
    print("Conditional survival and cap energy errors:", result["sweeps"])
    print(
        "Held port-nine count difference and uncertainty:", result["held_count_survey"]
    )
