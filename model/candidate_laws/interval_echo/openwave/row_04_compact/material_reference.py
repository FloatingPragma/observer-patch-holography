"""Compare vacancy and excess transfer energies using a physical storage port and clock."""

import argparse
import json
from collections import defaultdict
from itertools import product
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import block_diag, expm
from scipy.optimize import minimize_scalar
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.experiments.material_inertia import neutral_defect
from . import anchored_vacancy as anchor

SHOTS = 2**24
WAITS = 0.10 + 0.013 * np.arange(64)
K = anchor.K
LAW = anchor.candidate.native_law


def layout(board, homes, launch):
    graph = board.adjacency()
    used = set()
    stores = []
    paths = []
    for home, other in (homes, homes[::-1]):
        local = np.arange(12 * home, 12 * home + 12)
        far = np.arange(12 * other, 12 * other + 12)
        source = 12 * home + launch
        dist, parent = shortest_path(
            graph,
            directed=False,
            indices=source,
            unweighted=True,
            return_predecessors=True,
        )
        gap = shortest_path(graph, directed=False, indices=local, unweighted=True).min(
            axis=0
        )
        other_gap = shortest_path(
            graph, directed=False, indices=far, unweighted=True
        ).min(axis=0)
        pool = sorted(range(board.n), key=lambda p: (dist[p], p))
        choices, routes = {}, {}
        for degree in (5, 6):
            for p in pool:
                if board.degree[p] != degree or min(gap[p], other_gap[p]) < 2:
                    continue
                route = [p]
                while route[-1] != source:
                    route.append(int(parent[route[-1]]))
                if not used.intersection(route) and min(other_gap[route]) >= 2:
                    choices[degree], routes[degree] = p, route[::-1]
                    break
            else:
                raise ValueError("No disjoint storage path")
        for route in routes.values():
            used.update(route)
        used.update(local)
        stores.append(choices)
        paths.append(routes)
    near = used | {p for a in used for p in board.neighbours[a]}
    clock = next((a, b) for a, b in board.seams if a not in near and b not in near)
    return stores, paths, clock


def optics(sign):
    labels = [(1 + sign, 1), (1, 1 + sign)]
    idx = [LAW.SRC.index(pair) for pair in labels]
    seam = anchor.free_seam(K)
    # At integer kappa=8 a full addressed exchange is a pure permutation on all nine columns.
    swaps = [LAW.SRC.index((b, a)) for a, b in LAW.SRC]
    np.testing.assert_allclose(expm(-1j * np.pi * seam), np.eye(9)[swaps], atol=4e-14)
    splitter = expm(-0.5j * np.pi * seam)[np.ix_(idx, idx)]
    combiner = expm(-1.5j * np.pi * seam)[np.ix_(idx, idx)]
    np.testing.assert_allclose(combiner @ splitter, np.eye(2), atol=3e-14)
    blank = np.eye(LAW.DIM)[:, [*idx, LAW.SRC.index((1, 1))]]
    echo = LAW.echo()
    for delta in (-1, 0, 1):
        mask = np.repeat([False, *(d == delta for _, d in LAW.LABELS)], LAW.N_SRC)
        np.testing.assert_allclose(
            echo.T @ (mask[:, None] * (echo @ blank)),
            blank * (np.array([sign, -sign, 0]) == delta),
            atol=3e-14,
        )
    return splitter[:, 0], combiner


def sample_generator(board, sign, ref, cuts):
    h = anchor.defect_generator(board, sign, K, ref)
    # Remove each gated seam once. Its diagonal difference relative to a filled seam is K+1/2.
    a, b = np.asarray([e for e in board.seams if cuts.intersection(e)]).T
    lap = sparse.coo_matrix(
        (np.full(2 * len(a), -0.5), (np.r_[a, b], np.r_[b, a])), shape=h.shape
    )
    diag = np.bincount(np.r_[a, b], minlength=board.n) * (K + 0.5)
    return (h - lap - sparse.diags(diag)).tocsr()


def energy(h, v):
    hv = h @ v
    result = float(np.vdot(v, hv).real)
    np.testing.assert_allclose(hv, result * v, atol=2e-12, rtol=0)
    return result


def single_paths(h, sample, success, launch, store, sign, ref):
    # basis states are the sampled port, storage port and inaccessible spatial complement.
    # 1,024 clear tests act on the sample. A stored defect sees the filled source on every test.
    split, combine = optics(sign)
    tip = sample[launch]
    kept = np.zeros((3, 2), complex)
    kept[0, 0] = split[0] * np.sqrt(success) * tip
    kept[2, 0] = split[0] * np.sqrt(success * (1 - abs(tip) ** 2))
    kept[1, 1] = split[1]
    analyzer = np.eye(3, dtype=complex)
    analyzer[:2, :2] = combine
    store_state = np.zeros_like(sample)
    store_state[store] = 1
    levels = [energy(h, sample), energy(h, store_state)]
    # Implement the Y setting by a physical degree reference on the sharp degree-five launch port.
    y = analyzer @ np.diag([phase_quadrature(sign, ref), 1, 1])
    return kept, (analyzer, y), levels


def dimer_clock(rng):
    ids = [LAW.SRC.index(pair) for pair in ((0, 1), (1, 0))]
    h = anchor.free_seam(K)[np.ix_(ids, ids)]
    counts = []
    for t in WAITS:
        swap = abs(expm(-1j * t * h)[1, 0]) ** 2
        counts.append(int(rng.binomial(SHOTS, round(swap, 13))))
    return counts


def fit(readings, clock):
    # The fit has only the timer's local binary counts. The acquisition schedule is not an input.
    tau = 2 * np.arcsin(np.sqrt(np.asarray(clock) / SHOTS))

    def residual(freq):
        design = np.column_stack(
            [np.ones_like(tau), np.cos(freq * tau), np.sin(freq * tau)]
        )
        coeff = np.linalg.lstsq(design, readings, rcond=None)[0]
        return float(np.linalg.norm(design @ coeff - readings) ** 2)

    grid = np.arange(0, 160.125, 0.125)
    center = grid[min(range(len(grid)), key=lambda i: residual(grid[i]))]
    found = minimize_scalar(
        residual, bounds=(max(0.0, center - 0.125), center + 0.125), method="bounded"
    )
    f = found.x
    design = np.column_stack([np.ones_like(tau), np.cos(f * tau), np.sin(f * tau)])
    coef = np.linalg.lstsq(design, readings, rcond=None)[0]
    slope = tau * (-coef[1] * design[:, 2] + coef[2] * design[:, 1])
    jac = np.column_stack([design, slope])
    residuals = readings - design @ coef
    # Heteroskedasticity-consistent first-order standard error across independent appointments.
    error = np.linalg.norm(np.linalg.pinv(jac)[-1] * residuals) * np.sqrt(
        len(tau) / (len(tau) - 4)
    )
    return {
        "frequency_per_clock_gap": float(found.x),
        "frequency_standard_error": float(error),
        "rms_probability_residual": float(np.sqrt(found.fun / len(clock))),
    }


def joint_columns(board, samples, stores, sign, ref, cuts, pair_potential=None):
    """Apply actual native two-defect seam columns to each of the four coherent paths."""
    choices = []
    for sample, store in zip(samples, stores, strict=True):
        # Include every represented amplitude, even numerical residue on spectator ports.
        choices.append(
            [{int(p): sample[p] for p in np.flatnonzero(sample)}, {store: 1.0}]
        )
    seam = anchor.free_seam(K)
    energies = []
    for left, right in product(*choices):

        def key(a, b):
            return (min(a, b), max(a, b)) if sign == -1 else (a, b)

        state = {key(a, b): v * w for a, v in left.items() for b, w in right.items()}
        action = defaultdict(complex)
        for (a, b), amp in state.items():
            if a == b:
                raise ValueError("Overlapping sharp loads")
            shifts = {a: -1, b: sign}
            edges = {
                tuple(sorted((p, int(q))))
                for p in (a, b)
                for q in board.neighbours[p]
                if p not in cuts and q not in cuts
            }
            diag = sum(
                K * (6 - board.degree[p]) * ((1 + d - ref) ** 2 - (1 - ref) ** 2)
                for p, d in shifts.items()
            )
            if pair_potential is not None:
                diag += pair_potential[a // 12, b // 12]
            for u, v in edges:
                x, y = 1 + shifts.get(u, 0), 1 + shifts.get(v, 0)
                col = LAW.SRC.index((x, y))
                diag += seam[col, col]
                if x != y:
                    swap = {u: v, v: u}
                    action[key(swap.get(a, a), swap.get(b, b))] += (
                        seam[LAW.SRC.index((y, x)), col] * amp
                    )
            action[a, b] += diag * amp
        value = sum(np.conj(v) * action[k] for k, v in state.items())
        err = (
            sum(abs(v - value * state.get(k, 0)) ** 2 for k, v in action.items()) ** 0.5
        )
        if abs(value.imag) > 1e-12 or err > 1e-11:
            raise ValueError("Two-defect path is not a stationary joint state")
        energies.append(float(value.real))
    return energies


def record(arms, rng, joint=None, total_shift=0.0):
    amps, (matrices, y), levels = arms[0]
    ops = [matrices]
    if len(arms) == 2:
        second, (bx, by), _ = arms[1]
        amps = np.kron(amps, second)
        ops = [np.kron(matrices, bx), np.kron(y, by)]
        if joint is None:
            raise ValueError("Missing two-defect action")
        levels = joint
    largest_change = 0.0
    clock = dimer_clock(rng)
    counts = []
    signs = np.array([1, -1, 0])
    weights = signs if len(arms) == 1 else np.kron(signs, signs)
    for analyzer in ops:
        read = []
        for t in WAITS:
            probs = abs(analyzer @ amps @ np.exp(-1j * t * np.array(levels))) ** 2
            shifted = (
                abs(
                    analyzer @ amps @ np.exp(-1j * t * (np.array(levels) + total_shift))
                )
                ** 2
            )
            largest_change = max(largest_change, float(max(abs(shifted - probs))))
            p = np.maximum(np.round(np.r_[probs, 1 - probs.sum()], 12), 0)
            read.append(rng.multinomial(SHOTS, p / p.sum()).tolist())
        counts.append(read)
    parity = np.asarray(counts)[:, :, :-1] @ weights / SHOTS
    signal = parity[0] if len(arms) == 1 else parity[0] - parity[1]
    return {
        "timer_swaps": clock,
        "detector_counts": counts,
        "common_sector_shift": float(total_shift),
        "max_probability_change_under_sector_shift": largest_change,
        **fit(signal, clock),
    }


def phase_quadrature(sign, ref):
    cost = K * ((1 + sign - ref) ** 2 - (1 - ref) ** 2)
    pulse = (np.pi / 2 if cost < 0 else 3 * np.pi / 2) / abs(cost)
    return np.exp(-1j * pulse * cost)


def storage_calibration(h, stores, sign, ref, rng):
    a, b = stores[5], stores[6]
    delta = float(h[a, a] - h[b, b])
    split, read = optics(sign)
    arm = (
        np.vstack([np.diag(split), np.zeros(2)]),
        (block_diag(read, 1), None),
        [delta, 0],
    )
    data = record([arm], rng)
    y = read @ np.diag([phase_quadrature(sign, ref), 1])
    counts = []
    for t in WAITS:
        probs = abs(y @ (split * np.array([np.exp(-1j * t * delta), 1]))) ** 2
        counts.append(int(rng.binomial(SHOTS, float(np.clip(probs[0], 0, 1)))))
    tau = 2 * np.arcsin(np.sqrt(np.asarray(data["timer_swaps"]) / SHOTS))
    freq = data["frequency_per_clock_gap"]
    orientation = np.sign(
        np.dot(2 * np.asarray(counts) / SHOTS - 1, np.sin(freq * tau))
    )
    data.update(y_return_counts=counts, signed_frequency=float(orientation * freq))
    return data


def phase_flip(board, sample, success, launch, sign, ref, rng):
    h = anchor.defect_generator(board, sign, K, ref)
    state = sample.copy()
    cost = K * ((1 + sign - ref) ** 2 - (1 - ref) ** 2)
    state[launch] *= np.exp(-1j * np.pi * cost / abs(cost))
    delta = float((np.vdot(state, h @ state) - np.vdot(sample, h @ sample)).real)
    # Select the quadrature that reads 2 Re(psi_a* psi_b) with its sign.
    split, _ = optics(sign)
    splitter = np.array([[split[0], split[1]], [split[1], split[0]]])
    analyzer = splitter @ np.diag([phase_quadrature(sign, ref), 1])
    z = np.diag([1.0, -1.0])
    read_matrix = analyzer.conj().T @ z @ analyzer
    orientation = float(np.sign(read_matrix[0, 1].real))
    np.testing.assert_allclose(
        read_matrix, orientation * np.array([[0, 1], [1, 0]]), atol=1e-12
    )
    reads = []
    estimate = 0.0
    for other in board.neighbours[launch]:
        pair_counts = []
        for weight, preparation in ((0.5, sample), (-0.5, state)):
            xy = abs(analyzer @ preparation[[launch, other]]) ** 2
            probs = np.r_[success * xy, success * (1 - xy.sum()), 1 - success]
            counts = rng.multinomial(SHOTS, np.maximum(probs, 0)).tolist()
            pair_counts.append(counts)
            estimate += weight * orientation * (counts[0] - counts[1]) / sum(counts[:3])
        reads.append([int(other), pair_counts])
    difference = np.linalg.norm(h @ state - np.vdot(state, h @ state) * state)
    return {
        "bonds_and_counts": reads,
        "mean_energy_change_from_counts": float(estimate),
        "native_mean_change": delta,
        "post_flip_energy_standard_deviation": float(difference),
    }


def panel(board, ref, rng):
    samples, success, prep = anchor.preparation(board)
    homes = prep["cells"]
    tips = [home * 12 + prep["launch_local_port"] for home in homes]
    stores, paths, timer = layout(board, homes, prep["launch_local_port"])
    out = {
        "home_cells": homes,
        "source_success_if_released": float(success),
        "launch_ports": tips,
        "storage_paths": [{str(k): v for k, v in p.items()} for p in paths],
        "clock_seam": list(timer),
        "reference_only": {},
        "panels": {},
    }
    # Branches have disjoint support and no shared active seam, including the storage alternatives.
    support = [
        set(np.flatnonzero(abs(v) > 1e-12)) | set(st.values())
        for st, v in zip(stores, samples)
    ]
    if support[0] & support[1] or any(
        support[0].intersection(board.neighbours[p]) for p in support[1]
    ):
        raise ArithmeticError("Two release branches share a contact")
    for degree in (5, 6):
        cuts = {s[degree] for s in stores} | set(timer)
        h = {sign: sample_generator(board, sign, ref, cuts) for sign in (-1, 1)}
        arms = {
            (i, sign): single_paths(
                h[sign], samples[i], success, tips[i], stores[i][degree], sign, ref
            )
            for i, sign in product(range(2), (-1, 1))
        }
        if degree == 5:
            for sign in (-1, 1):
                both = cuts | {s[6] for s in stores}
                out["reference_only"][str(sign)] = storage_calibration(
                    sample_generator(board, sign, ref, both), stores[0], sign, ref, rng
                )
        settings = {
            "vacancy": [arms[0, -1]],
            "excess": [arms[0, 1]],
            "two_vacancies": [arms[0, -1], arms[1, -1]],
            "opposite_signs": [arms[0, -1], arms[1, 1]],
        }
        joint = {
            name: joint_columns(
                board, samples, [s[degree] for s in stores], sign, ref, cuts
            )
            for name, sign in (("two_vacancies", -1), ("opposite_signs", 1))
        }
        out["panels"][str(degree)] = {
            name: record(
                arms_for_read,
                rng,
                joint.get(name),
                {
                    "vacancy": 16,
                    "excess": -16,
                    "two_vacancies": 32,
                    "opposite_signs": 0,
                }[name],
            )
            for name, arms_for_read in settings.items()
        }
        if board.cells == 20 and ref == 0 and degree == 6:
            a, b = np.array(board.external).T // 12
            adjacency = sparse.coo_matrix(
                (np.ones(len(a)), (a, b)), shape=(board.cells, board.cells)
            ).tocsr()
            distances = shortest_path(adjacency, directed=False, unweighted=True)
            modified = joint_columns(
                board,
                samples,
                [s[degree] for s in stores],
                -1,
                ref,
                cuts,
                -2 / (1 + distances),
            )
            out["separated_interaction_control"] = record(
                settings["two_vacancies"], rng, modified
            )
            out["separated_interaction_control"]["changed_joint_levels"] = modified
    out["local_phase_flip"] = {
        str(sign): phase_flip(board, samples[0], success, tips[0], sign, ref, rng)
        for sign in (-1, 1)
    }
    positive = anchor.defect_generator(board, 1, K, ref)
    off_radius = abs(positive).sum(axis=1).A1 - abs(positive.diagonal())
    out["all_single_excess_formation_energy_bounds"] = [
        float(np.min(positive.diagonal() - off_radius)),
        float(np.max(positive.diagonal() + off_radius)),
    ]
    if board.cells == 20 and ref == 0:
        # Independently check both sample-sample energies on the interacting connected sectors.
        pair = PairSector(board)
        a, b = pair.aa, pair.bb
        va = K * (2 * board.degree - 6)
        full = pair.exchange() + sparse.diags(va[a] + va[b] - 2 * K * pair.adjacent)
        left, right = samples
        x, y, _, _, mixed, _ = neutral_defect.generators(board)
        out["full_two_sample_energies"] = [
            energy(full, left[a] * right[b] + right[a] * left[b]),
            energy(mixed + 12 * K * sparse.eye(len(x)), right[x] * left[y]),
        ]
    return out


def run():
    rng = np.random.default_rng(178031)
    return {
        "shots": SHOTS,
        "seed": 178031,
        "native": [panel(anchor.federation(level), 0, rng) for level in (0, 2)],
        "unit_centered_control": panel(anchor.federation(0), 1, rng),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, allow_nan=False, separators=(",", ":"))
            stream.write("\n")
    else:
        from ..electron_field.replay import compare

        compare(
            result,
            json.loads(Path(__file__).with_suffix(".json").read_text()),
            "material_reference",
        )
        print(
            "Reproduced the storage calibration, clock records and both two-defect spectra."
        )
