"""Joint energy and complete two-port propagation reads on the registered seam law."""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label

HERE = Path(__file__).resolve().parent
TIMES = (0.5, 1.0, 2.0, 4.0)
N = 2**28
SEED = 390177
EPS = float(np.sqrt(np.log(2 * 100000 / 0.01) / (2 * N)))
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
UNIT = np.eye(2)
TRAIN = np.vstack([np.eye(3), -np.eye(3)])
HELD = np.array(
    [(a, b, a * b) for a, b in itertools.product((-1, 1), repeat=2)]
) / np.sqrt(3)
INPUTS = np.vstack([TRAIN, HELD, np.zeros(3)])
ARMS = {
    "low": (1.0, 1, 1, 8.0),
    "high": (-1.0, 1, 1, 8.0),
    "quarter": (0.5, 1, 1, 8.0),
    "collar": (1.0, -1, 1, 8.0),
    "disconnected_low": (1.0, 1, -1, 8.0),
    "disconnected_high": (-1.0, 1, -1, 8.0),
    "changed_low": (1.0, 1, 1, 0.0),
    "changed_high": (-1.0, 1, 1, 0.0),
}


def layouts(board):
    """Ordered six-port instruments; the second dimer has equal native degree."""
    rows = []
    for a, b in sorted(e for e in board.external if e[0] // 12 == 0):
        c, d, u, v = min(
            (c, d, u, v)
            for c in board.neighbours[b]
            if c // 12 == b // 12
            for d in board.neighbours[c]
            if d // 12 == c // 12 and d != b and board.degree[c] == board.degree[d]
            for u in board.neighbours[a]
            if u // 12 == a // 12
            for v in board.neighbours[d]
            if len({a, b, c, d, u, v}) == 6
        )
        rows.append((a, b, c, d, u, v))
    if len(rows) != 3:
        raise ValueError("Missing bridge")
    return rows


def reduced_operators(strength, reversal, coin, q):
    """Verify invariance in all 99 directions, including the inactive zero guards."""
    law = candidate.native_law
    kinetic, contrast = law.local_term_parts()
    h = kinetic + strength * contrast
    eye = np.eye(law.DIM)
    rail = eye[:, [law.SRC.index((1 + q, 1)), law.SRC.index((1, 1 + q))]]
    rail_h = rail.T @ h @ rail
    np.testing.assert_allclose(h @ rail, rail @ rail_h, atol=1e-13)
    modes = np.kron(label(2, reversal, coin)[:, None], np.eye(law.N_SRC))
    contact = modes.T @ h @ modes
    np.testing.assert_allclose(h @ modes, modes @ contact, atol=1e-13)
    np.testing.assert_allclose(contact, np.diag(np.diag(contact)), atol=1e-13)
    guards = np.kron(label(0, 1, 1)[:, None], np.eye(law.N_SRC))
    np.testing.assert_allclose(h @ guards, 0, atol=1e-13)
    return rail_h, dict(zip(law.SRC, contact.diagonal(), strict=True))


def source_meter(board, sites, s, arm, rng):
    """Read the source alone; count parity and the total stored contrast on its marks."""
    mu, parity, coin, k = ARMS[arm]
    c, d = sites[2:4]
    x = int(rng.binomial(N, (1 - mu) / 2))  # odd exchange on the source rail
    law = candidate.native_law
    operators = law.generator_parts()
    memory = label(2, parity, coin)
    # The actual inverse-echo read is diagonal in the source position.
    # Sum all three mark energies on each source position; verify a sharp sum.
    records = []
    for where in (c, d):
        rs, sq = 0.0, 0.0
        for left, right in ((sites[1], c), (sites[0], sites[4]), (d, sites[5])):
            word = (1 + s * int(left == where), 1 + s * int(right == where))
            basis = np.kron(memory, np.eye(9)[:, law.SRC.index(word)])
            lifted = law.echo() @ basis
            for idx, operator in enumerate(operators):
                mean = float(lifted @ operator @ lifted)
                np.testing.assert_allclose(operator @ lifted, mean * lifted, atol=1e-13)
                if idx == 0:
                    rs += mean
                else:
                    sq += mean
        records.append([round(rs), round(sq)])
    if records[0] != records[1]:
        raise ValueError("Source mark sum is uncertain")
    references = int(6 - board.degree[c]) * ((1 + s) ** 2 - 1)
    intercept = k * (1 + records[0][1] + references) + records[0][0]
    return dict(
        odd_exchange_count=x,
        mark_reversal=records[0][0],
        mark_contrast_squared=records[0][1],
        material_reference_increment=references,
        mean_energy_above_uniform=intercept + x / N,
    )


def generator(board, sites, s, q, arm):
    mu, parity, coin, k = ARMS[arm]
    rail_p, contact = reduced_operators(k, parity, coin, q)
    rail_s, _ = reduced_operators(k, parity, coin, s)
    a, b, c, d, u, v = sites
    h = np.kron(rail_p, UNIT) + np.kron(UNIT, rail_s)
    for i, (p, src) in enumerate(itertools.product((a, b), (c, d))):

        def occupation(port):
            return 1 + q * int(p == port) + s * int(src == port)

        h[i, i] += sum(
            contact[occupation(x), occupation(y)] for x, y in ((b, c), (a, u), (d, v))
        )
        h[i, i] += k * sum(
            (6 - board.degree[z]) * (occupation(z) ** 2 - 1) for z in (a, b, c, d)
        )
    h -= h.trace() / 4 * np.eye(4)
    return h, (UNIT + mu * PAULI[0]) / 2


def probabilities(h, rho, time):
    u = expm(-1j * time * h)
    effects = (UNIT[None] + PAULI) / 2
    effects = np.array([u.conj().T @ np.kron(e, UNIT) @ u for e in effects])
    inputs = (UNIT + np.einsum("ka,amn->kmn", INPUTS, PAULI)) / 2
    launches = np.array([np.kron(p, rho) for p in inputs])
    prob = np.einsum("aij,kji->ka", effects, launches).real
    if prob.min() < -1e-12 or prob.max() > 1 + 1e-12:
        raise ValueError("Invalid pointer probability")
    return np.clip(prob, 0, 1)


def map_predict(counts, vectors):
    p = np.asarray(counts[:6]) / N
    offset = p.mean(axis=0)
    return offset + np.asarray(vectors) @ ((p[:3] - p[3:]) / 2)


def acquire():
    rng = np.random.default_rng(SEED)
    shots, source_rows, clock = [], [], []
    for level in range(3):
        board = federation(level)
        for bridge, sites in enumerate(layouts(board)):
            # A separate, unmarked blank dimer clock starts at its first endpoint.
            excluded_owners = {site // 12 for site in sites}
            clock_ports = next(
                (i, j)
                for i, adjacent in enumerate(board.neighbours)
                if i // 12 not in excluded_owners
                for j in sorted(adjacent)
                if j > i and j // 12 == i // 12 and board.degree[j] == board.degree[i]
            )
            h, _ = reduced_operators(8, 1, 1, -1)
            readings = []
            for t in TIMES:
                psi = expm(-1j * t * h) @ np.array([1, 0])
                yz = [float((psi.conj() @ op @ psi).real) for op in PAULI[1:]]
                counts = rng.binomial(N, np.clip((np.array(yz) + 1) / 2, 0, 1)).tolist()
                readings.append(counts)
            clock.append(
                dict(
                    cells=board.cells,
                    bridge=bridge,
                    ports=clock_ports,
                    quadrature_plus=readings,
                )
            )
            for s in (-1, 1):
                for name in ARMS:
                    source_rows.append(
                        dict(
                            cells=board.cells,
                            bridge=bridge,
                            s=s,
                            arm=name,
                            **source_meter(board, sites, s, name, rng),
                        )
                    )
                    for q in (-1, 1):
                        h, rho = generator(board, sites, s, q, name)
                        counts = [
                            rng.binomial(N, probabilities(h, rho, t)).tolist()
                            for t in TIMES
                        ]
                        shots.append(
                            dict(
                                cells=board.cells,
                                bridge=bridge,
                                sites=list(sites),
                                q=q,
                                s=s,
                                arm=name,
                                plus_counts=counts,
                            )
                        )
    return shots, source_rows, clock


def summarize(shots, source, clock):
    by_arm = {
        (r["cells"], r["bridge"], r["s"], r["q"], r["arm"]): np.array(r["plus_counts"])
        for r in shots
    }
    meter = {
        (r["cells"], r["bridge"], r["s"], r["arm"]): r["mean_energy_above_uniform"]
        for r in source
    }
    calibration = {arm: by_arm[20, 0, -1, -1, arm] for arm in ("low", "high")}
    held = {}
    for arm in ("low", "high", "quarter", "collar"):
        errors, detected = [], []
        for r in shots:
            if r["arm"] != arm:
                continue
            c, b, s = r["cells"], r["bridge"], r["s"]
            if (c, b, s, r["q"]) == (20, 0, -1, -1) and arm in calibration:
                continue
            alpha = (meter[c, b, s, arm] - meter[c, b, s, "low"]) / (
                meter[c, b, s, "high"] - meter[c, b, s, "low"]
            )
            for ix, actual in enumerate(r["plus_counts"]):
                forecast = map_predict(
                    (1 - alpha) * calibration["low"][ix]
                    + alpha * calibration["high"][ix],
                    INPUTS,
                )
                errs = np.abs(forecast - np.array(actual) / N)
                errors.append(float(errs.max()))
                # Sum of absolute calibration weights <= 1 + sqrt(3), plus held read;
                # alpha within EPS for the quarter-mixture; exact on all other arms.
                detected.append(int((errs > (3 + np.sqrt(3)) * EPS).sum()))
        held[arm] = dict(
            largest_frequency_error=max(errors),
            rejected_pointer_components=sum(detected),
            compared_pointer_components=len(errors) * len(INPUTS) * 3,
        )
    comparisons = {}
    for arm in ARMS:
        results = []
        for row in shots:
            if row["arm"] == arm:
                mixed = np.array(row["plus_counts"])[:, -1, :] / N - 0.5
                results.append(np.max(abs(mixed), axis=1))
        observed = np.array(results)
        comparisons[arm] = dict(
            minimum_mixed_state_distance_by_time=observed.min(axis=0).tolist(),
            maximum_mixed_state_component_by_time=observed.max(axis=0).tolist(),
        )
    composition = max(
        np.max(np.abs(cnt - by_arm[c, b, s, 1, arm])) / N
        for (c, b, s, q, arm), cnt in by_arm.items()
        if q == -1
    )
    clock_error = max(
        np.max(
            np.abs(
                (2 * np.array(r["quadrature_plus"]) / N - 1)
                - np.array([[np.sin(t), np.cos(t)] for t in TIMES])
            )
        )
        for r in clock
    )
    readings = [2 * np.asarray(row["quadrature_plus"]) / N - 1 for row in clock]
    phase_error = max(
        np.max(np.abs(np.unwrap(np.arctan2(yz[:, 0], yz[:, 1])) - TIMES))
        for yz in readings
    )
    return dict(
        held_energy_calibration=held,
        calibration_and_held_radius=(3 + np.sqrt(3)) * EPS,
        max_two_charge_probability_difference=float(composition),
        mixed_probe=comparisons,
        clock_max_quadrature_error=float(clock_error),
        clock_max_exposure_error=float(phase_error),
    )


def run():
    shots, source, clock = acquire()
    total = (
        len(shots) * len(TIMES) * len(INPUTS) * 3
        + len(source)
        + len(clock) * len(TIMES) * 2
    )
    if total > 100000:
        raise ValueError("Simultaneous budget exceeded")
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        shots_per_binary_cohort=N,
        single_frequency_radius=EPS,
        times=TIMES,
        probe_axes_xyz=INPUTS.tolist(),
        number_of_binary_cohorts=total,
        source=source,
        clocks=clock,
        trials=shots,
        summary=summarize(shots, source, clock),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    result = run()
    with p.parse_args().output.open("x") as out:
        json.dump(result, out, separators=(",", ":"), allow_nan=False)
        out.write("\n")
    print(json.dumps(result["summary"], indent=2))
