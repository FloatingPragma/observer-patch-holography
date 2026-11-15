#!/usr/bin/env python3
"""Independent finite review calculations; imports no producer result code."""
import argparse
import ast
from fractions import Fraction
import gzip
import itertools
import json
import math
from pathlib import Path
import random
import tempfile
from datetime import datetime, timezone

import numpy as np

from reproduce import HERE, FREE, extract, save, sha, verify_manifest


def require(condition, description):
    if not condition:
        raise ValueError(description)


def endpoint(position, edge, coin):
    """Floor/ceiling averaging, independently of native.py's difference formula."""
    a, b = edge[:2]
    left, right = int(position == a), int(position == b)
    total = left + right
    new_left = total // 2 if coin == 0 else (total + 1) // 2
    current = left - new_left
    return (b if current == 1 else a if current == -1 else position), current


def geometry():
    golden = (1 + math.sqrt(5)) / 2
    vertices = np.array([
        [-1, golden, 0], [1, golden, 0], [-1, -golden, 0], [1, -golden, 0],
        [0, -1, golden], [0, 1, golden], [0, -1, -golden], [0, 1, -golden],
        [golden, 0, -1], [golden, 0, 1], [-golden, 0, -1], [-golden, 0, 1]])
    edges = [(a, b) for a in range(12) for b in range(a + 1, 12)
             if abs(float(np.sum((vertices[a] - vertices[b]) ** 2)) - 4) < 1e-12]
    lap = np.zeros((12, 12))
    for a, b in edges:
        lap[a, a] += 1
        lap[b, b] += 1
        lap[a, b] -= 1
        lap[b, a] -= 1
    require(len(edges) == 30 and np.all(np.diag(lap) == 5), "icosahedral geometry")
    return edges, lap


def curl_from_counts(counts, cells, width, weight):
    """Explicit adjacent-cell comparisons, rather than the producer's np.roll."""
    count = np.asarray(counts)
    y = math.sqrt(width / (4 * weight * 32768)) * (count - 32768 / width)
    potentials = y.reshape(4, len(cells), 12).sum(axis=2) / math.sqrt(12)
    lookup = {tuple(c): i for i, c in enumerate(cells)}
    fields = []
    for mu in range(4):
        for nu in range(mu + 1, 4):
            values = []
            for i, cell in enumerate(cells):
                cm, cn = list(cell), list(cell)
                dm, dn = (3, 0, 1, 2)[mu], (3, 0, 1, 2)[nu]
                cm[dm] = (cm[dm] + 1) % 3
                cn[dn] = (cn[dn] + 1) % 3
                # This helper is used for the free 0..2 coordinate chart only.
                values.append(potentials[nu, lookup[tuple(cm)]] - potentials[nu, i]
                              - potentials[mu, lookup[tuple(cn)]] + potentials[mu, i])
            fields.append(values)
    return np.array(fields), y


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--output", type=Path, default=HERE / "results/INDEPENDENT.json")
    args = parser.parse_args()
    if args.bundle is None:
        # Reuse only the reviewer's Git-byte extraction, never producer computations.
        import os
        scratch = Path(os.environ.get("LOCALAPPDATA", tempfile.gettempdir())) / "Temp/opencode"
        scratch.mkdir(parents=True, exist_ok=True)
        commit, _, bundle = extract(Path(tempfile.mkdtemp(prefix="oph_wave_math_", dir=scratch)))
    else:
        bundle, commit = args.bundle, "caller-supplied; bound by manifest below"
    pins = {}

    def load(relative, compressed=False):
        path = bundle / relative
        pins[relative] = sha(path)
        data = path.read_bytes()
        return json.loads(gzip.decompress(data) if compressed else data)

    report = dict(subject_commit=commit, started_utc=datetime.now(timezone.utc).isoformat(),
                  script_sha256=sha(Path(__file__)), numpy=np.__version__,
                  manifest_sha256=sha(bundle / "MANIFEST.sha256"),
                  manifest_entries=verify_manifest(bundle, bundle / "MANIFEST.sha256"),
                  scope="Independent finite arithmetic and tape review, not a formal proof or G6.")
    records = load(FREE + "/RECORDS.json.gz", True)
    receipt = load(FREE + "/RECEIPT.json")
    internal, lap = geometry()
    spectrum = np.linalg.eigvalsh(lap)
    bands = [(0., 1), (5 - math.sqrt(5), 3), (6., 5), (5 + math.sqrt(5), 3)]
    target = np.concatenate([np.full(n, value) for value, n in bands])
    spectrum_error = float(np.max(abs(spectrum - target)))
    require(spectrum_error < 1e-12, "native internal spectrum")
    cells = records["graph"]["cells"]
    require(cells == [list(c) for c in itertools.product(range(3), repeat=4)], "free cell chart")
    lookup = {tuple(c): i for i, c in enumerate(cells)}
    edges = []
    for c, cell in enumerate(cells):
        edges += [[12 * c + a, 12 * c + b, -1] for a, b in internal]
        for axis in range(4):
            neighbor = list(cell)
            neighbor[axis] = (neighbor[axis] + 1) % 3
            edges += [[12 * c + p, 12 * lookup[tuple(neighbor)] + p, axis] for p in range(12)]
    require(edges == records["graph"]["edges"], "full 972-vertex product routing")
    require(len(edges) == 6318, "seam count")
    for edge_id, slot, coin, current, target_slot in records["branches"]:
        predicted, j = endpoint(slot, (0, 1), coin)
        require(0 <= edge_id < len(edges) and predicted == target_slot and j == current,
                "native branch census")
    # Inspect just the pinned primitive as an isolated function, not its result code.
    source = bundle / FREE / "inputs/native.py"
    pins[FREE + "/inputs/native.py"] = sha(source)
    node = next(n for n in ast.parse(source.read_text(encoding="utf-8")).body
                if isinstance(n, ast.FunctionDef) and n.name == "native")
    namespace = {}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(source), "exec"), namespace)
    for position, coin in itertools.product(range(3), (0, 1)):
        result, rec = namespace["native"](tuple(int(i == position) for i in range(3)), (0, 1), coin)
        predicted, current = endpoint(position, (0, 1), coin)
        require(result.index(1) == predicted and rec == (coin, current), "isolated native primitive")
    rng = random.Random(2026111432)
    all_counts, all_clocks = [], []
    for start, tape, expected_counts, expected_clock in zip(
            records["channel_starts"], records["tapes"], records["occupation_counts"],
            records["same_tape_clock_records"]):
        require(start == rng.randrange(972), "declared uniform starting draw")
        position, counts, clock, ticks = start, [0] * 972, [], 0
        origin = cells[start // 12][3]
        for step, (edge_id, coin, current) in enumerate(tape, 1):
            require(edge_id == rng.randrange(6318) and coin == rng.randrange(2), "free RNG schedule")
            position, j = endpoint(position, edges[edge_id], coin)
            require(j == current, "current record agrees with independent endpoint")
            counts[position] += 1
            if cells[position // 12][3] == (origin + ticks + 1) % 3:
                ticks += 1
                clock.append([step, ticks, position])
        require(counts == expected_counts and clock == expected_clock, "full occupation/clock readout")
        all_counts.append(counts)
        all_clocks.append(clock)
    field, y = curl_from_counts(all_counts, cells, 972, 6318)
    field_error = float(np.max(abs(field - np.array(records["field_strength_readouts"]))))
    band_error = 0.
    for (value, _), recorded in zip(bands, records["internal_band_readouts"]):
        projection = np.eye(12)
        for other, _ in bands:
            if other != value:
                projection = projection @ (lap - other * np.eye(12)) / (value - other)
        error = float(np.max(abs(y.reshape(4, 81, 12) @ projection - recorded)))
        band_error = max(error, band_error)
    require(field_error < 1e-12 and band_error < 1e-12, "local field and spectral readouts")
    report["free_tapes"] = dict(branches=len(records["branches"]), isolated_primitive_cases=6,
                              actual_current_records=sum(map(len, records["tapes"])),
                              clock_ticks=list(map(len, all_clocks)),
                              graph_vertices=972, graph_seams=6318, field_max_error=field_error,
                              band_readout_max_error=band_error, spectrum_max_error=spectrum_error)

    max_covariance_error = max_hit_error = 0.
    for entry in receipt["population"]:
        lam = min((x for x, _ in bands), key=lambda x: abs(x - entry["internal_lambda"]))
        A = lam + 3 * sum(x != 0 for x in entry["spatial"])
        values = [A, A + 3, A + 3]
        finite = []
        n, W = 100_000_000, 6318
        for value in values:
            z = 1 - value / (2 * W)
            # Closed sum of all ordered time pairs; no eigenvector Rayleigh quotient.
            total = n * (1 + z) / (1 - z) - 2 * z * (1 - z ** n) / (1 - z) ** 2
            finite.append(total / (4 * W * n))
        cov = [(finite[0] + 2 * finite[1]) / 3, (finite[0] - finite[1]) / 3]
        max_covariance_error = max(max_covariance_error,
                                   max(abs(a - b) for a, b in zip(cov, entry["finite_block_covariance"])))
        max_hit_error = max(max_hit_error, abs(1 / (A + 1) - entry["first_hit_transform"]))
    require(max_covariance_error < 1e-10 and max_hit_error < 1e-10, "15 analytic populations")
    report["population"] = dict(modes=len(receipt["population"]), block_length=100_000_000,
                                actual_tape=False, max_covariance_difference=max_covariance_error,
                                max_hitting_difference=max_hit_error)

    # Check occupation CLT ingredients with a different, explicit moment recursion.
    P = np.eye(12) - lap / 60
    eigenvalues, vectors = np.linalg.eigh(lap)
    green = sum(np.outer(vectors[:, i], vectors[:, i]) / eigenvalues[i] for i in range(1, 12))
    projection = np.eye(12) - np.ones((12, 12)) / 12
    g = np.zeros(12)
    g[0], g[1] = 1, -1
    variance_target = float(g @ (green - projection / 120) @ g)
    moments = np.zeros((5, 12))
    moments[0] = 1 / 12
    moment_results = []
    for n in range(1, 4097):
        propagated = moments @ P
        moments = np.array([sum(math.comb(k, j) * g ** (k - j) * propagated[j]
                                for j in range(k + 1)) for k in range(5)])
        if n in (64, 256, 1024, 4096):
            scale = math.sqrt(12 / (120 * n))
            m2, m4 = float(moments[2].sum() * scale ** 2), float(moments[4].sum() * scale ** 4)
            moment_results.append(dict(n=n, variance=m2, fourth_moment=m4,
                                       variance_target=variance_target,
                                       fourth_target=3 * variance_target ** 2))
    require(abs(moment_results[-1]["variance"] - variance_target) <
            abs(moment_results[0]["variance"] - variance_target), "variance convergence diagnostic")
    require(abs(moment_results[-1]["fourth_moment"] - 3 * variance_target ** 2) <
            abs(moment_results[0]["fourth_moment"] - 3 * variance_target ** 2),
            "fourth moment convergence diagnostic")
    # Direct finite sum verifies the contact correction and finite-window formula.
    covariance = projection.copy()
    power = np.eye(12)
    window = 64
    for k in range(1, window):
        power = power @ P
        covariance += 2 * (1 - k / window) * (projection @ power @ projection)
    covariance /= 120
    exact_variance = float(g @ covariance @ g)
    require(abs(exact_variance - moment_results[0]["variance"]) < 1e-12, "finite moment/covariance join")
    report["moment_diagnostics"] = moment_results

    dispersion_tests = 0
    max_group_speed = 0.
    for spacing, (lam, _), p in itertools.product(
            (1., .5, .25, .0625), bands, itertools.product((-2., -.7, 0., .9, 2.), repeat=3)):
        p = np.array(p)
        q = float(np.sum(4 * np.sin(spacing * p / 2) ** 2) / spacing ** 2)
        energy = 2 * math.asinh(spacing * math.sqrt(lam + q) / 2) / spacing
        error = lam + float(p @ p) - energy ** 2
        bound = spacing ** 2 * (float(np.sum(p ** 4)) + (lam + float(p @ p)) ** 2) / 12
        require(-1e-12 <= error <= bound + 1e-12, "uniform KG bound")
        if lam + q > 0:
            speed = float(np.linalg.norm(np.sin(spacing * p)) /
                          (spacing * math.sqrt((lam + q) * (1 + spacing ** 2 * (lam + q) / 4))))
            require(speed <= 1 + 1e-12, "finite group speed")
            max_group_speed = max(max_group_speed, speed)
        dispersion_tests += 1
    report["dispersion"] = dict(sampled_cases=dispersion_tests, max_group_speed=max_group_speed,
                                note="Finite corroboration; analytic inequalities reviewed separately.")

    clock_tests, clock_error, mean_error = 0, 0., 0.
    for period, A in itertools.product((3, 5, 7, 17), (0., 5 - math.sqrt(5), 3., 6., 9.)):
        d = period - 1
        killed = (A + 2) * np.eye(d) - np.eye(d, k=1) - np.eye(d, k=-1)
        boundary = np.zeros(d)
        boundary[0] += 1
        boundary[-1] += 1
        hit = float(np.linalg.solve(killed, boundary)[0])
        omega = 2 * math.asinh(math.sqrt(A) / 2)
        predicted = math.cosh(omega * (period / 2 - 1)) / math.cosh(omega * period / 2)
        clock_error = max(clock_error, abs(hit - predicted))
        if A:
            require(-1e-12 <= hit / math.exp(-omega) - 1 <= math.exp(-omega * (period - 2)) + 1e-12,
                    "cycle clock convergence bound")
        W = 78 * 27 * period
        mean = float(np.linalg.solve(2 * np.eye(d) - np.eye(d, k=1) - np.eye(d, k=-1),
                                    np.full(d, 2 * W))[0])
        mean_error = max(mean_error, abs(mean - W * (period - 1)))
        clock_tests += 1
    require(clock_error < 1e-12 and mean_error < 1e-8, "killed-chain clock and cost")
    report["clock"] = dict(cases=clock_tests, max_multiplier_error=clock_error,
                           max_mean_attempt_error=mean_error)

    reflected_cases, factor_error, least_negative = 0, 0., 0.
    for p in ((1., 0., 0.), (.4, -.7, 1.1), (2., 3., -1.)):
        k = np.array(p)
        Q = float(np.linalg.norm(k))
        transverse = np.eye(3) - np.outer(k, k) / Q ** 2
        x, y, z = k
        cross = 1j * np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
        gram = np.block([[Q ** 2 * transverse, -Q * cross], [-Q * cross, Q ** 2 * transverse]])
        factor = np.concatenate([Q * transverse, -cross], axis=1)
        factor_error = max(factor_error, float(np.max(abs(gram - factor.conj().T @ factor))))
        spectrum = np.linalg.eigvalsh(gram)
        require(spectrum[0] > -1e-12 and sum(spectrum > 1e-10) == 2, "Maxwell positive rank-two quotient")
        wrong = gram.copy()
        wrong[:3, :3] *= -1
        least_negative = min(least_negative, float(np.linalg.eigvalsh(wrong)[0]))
        require(np.linalg.eigvalsh(wrong)[0] < -.1, "wrong reflection sign is detected")
        reflected_cases += 1
    report["maxwell_quotient"] = dict(cases=reflected_cases, rank=2,
                                     max_factorization_error=factor_error,
                                     wrong_sign_control_min_eigenvalue=least_negative)

    unified = load("evidence/unified/RECEIPT.json")
    # For the two registered directions A=3 and B²=3 or 0.
    expected = [Fraction(3, 2 * (9 + 3)), Fraction(3, 2 * 9)]
    contrast = 6 * (expected[1] - expected[0])
    require(contrast == Fraction(1, 4), "exact directed-coin contrast")
    for row, value in zip(unified["changed_law_same_F12_population"], expected):
        require(abs(row["biased_covariance"] - (float(value) - 1 / (4 * 6318))) < 1e-10,
                "directed population response")
    report["changed_law"] = dict(exact_F12_contrast=str(contrast),
                                ratios=["3/4", "1"], measured_precision_claim=False)
    report["source_inputs_sha256"] = pins
    report["all_review_checks_pass"] = True
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    save(args.output, report)
    print(json.dumps(dict(report=str(args.output), all_review_checks_pass=True,
                          native_tape_records=report["free_tapes"]["actual_current_records"],
                          mathematical_dispersion_checks=dispersion_tests), indent=2))


if __name__ == "__main__":
    main()
