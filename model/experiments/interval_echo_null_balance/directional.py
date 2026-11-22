"""Exact and finite-record tests of additive, time-odd single-seam matter currents."""

import argparse
from collections import Counter
from itertools import combinations, combinations_with_replacement
import json
from math import log, sqrt

import numpy as np
import sympy as sp
from scipy import sparse
from scipy.linalg import eigvalsh, expm
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from model.candidate_laws.interval_echo import candidate
from model.evaluation.patch_spatial_readout import transitions
from model.executable_model import cell, federation, rotations
from model.experiments.interval_echo_angular.__main__ import reconstruct
from .experiment import HERE, portable
from .regional import verify

SHOTS = 2**20
TIMES = (0, 1 / 8, 1 / 4)
SEED = 253119


def local_map(gateways):
    """Columns for B real skew, rows for the symmetric [L,B] on actual cell paths."""
    car = cell()
    g = len(gateways)
    n = 12 + g
    laplacian = np.zeros((n, n), dtype=np.int64)
    laplacian[:12, :12] = car.laplacian().toarray().astype(int)
    mixed = []
    bridges = []
    for k, gate in enumerate(gateways, 12):
        bridges.append((gate, k))
        laplacian[gate, gate] += 1
        laplacian[k, k] = 6
        laplacian[gate, k] = laplacian[k, gate] = -1
        mixed.extend((v, k) for v in car.neighbours[gate])
    rows = list(combinations_with_replacement(range(12), 2)) + mixed
    columns = []
    for a, b in car.seams + bridges:
        skew = np.zeros_like(laplacian)
        skew[a, b], skew[b, a] = 1, -1
        residual = laplacian @ skew - skew @ laplacian
        columns.append([residual[i, j] for i, j in rows])
    # Each diagonal output occurs once. Every other symmetric output occurs twice.
    # Frobenius input counts an internal skew entry twice and each bridge once per owner.
    output = np.array([1 if i == j else 2 for i, j in rows])
    input_mass = np.diag([2] * len(car.seams) + [1] * g)
    return np.array(columns, dtype=np.int64).T, output, input_mass


def canon(gateways):
    return min(tuple(sorted(rotation[g] for g in gateways)) for rotation in rotations())


def exact_ward_bound():
    types = Counter(canon(g) for g in combinations(range(12), 3))
    by_class = []
    for gates, size in sorted(types.items()):
        operator, weight, mass = local_map(gates)
        gram = operator.T @ (weight[:, None] * operator)
        # For H=6*kappa I+L/2, requiring a > 1/8 norm ratio is
        # equivalent to positivity of 16*gram - mass. Strengthen this to
        # 9/40 by proving positivity of 400*gram - 81*mass.
        certificate = sp.Matrix(400 * gram - 81 * mass)
        _, diagonal = certificate.LDLdecomposition(hermitian=False)
        if any(d <= 0 for d in diagonal.diagonal()):
            raise ValueError("The registered uniform bound failed")
        numerical_min = sqrt(eigvalsh(gram, mass, subset_by_index=[0, 0])[0]) / 2
        by_class.append(
            dict(
                gateways=list(gates),
                triples=size,
                positive_exact_pivots=len(gram),
                approximate_local_minimum=numerical_min,
            )
        )
    isolated, _, _ = local_map(())
    kernel = isolated.shape[1] - sp.Matrix(isolated).rank()
    if kernel != 6 or sum(types.values()) != 220:
        raise ValueError("The isolated-cell control or gateway inventory changed")
    return dict(
        closed_cell_conserved_current_dimension=kernel,
        registered_frobenius_lower_bound=1 / 8,
        strict_frobenius_lower_bound=9 / 40,
        gateway_triples=220,
        rotation_orbits=by_class,
    )


def incidence(board):
    owners = {(a // 12, b // 12) for a, b in board.external}
    used = {p for e in board.external for p in e}
    if len(owners) != len(board.external) or len(used) != 2 * len(board.external):
        raise ValueError("The mixed path may have another contribution")
    gates = [
        tuple(np.flatnonzero(board.degree[12 * i : 12 * (i + 1)] == 6))
        for i in range(board.cells)
    ]
    if any(len(g) != 3 for g in gates):
        raise ValueError("A cell requires exactly three distinct gateway ports")
    return Counter(canon(g) for g in gates)


def check_global_decomposition(board, level):
    """Compare the entire full-board [L,B] norm with separately assembled cell tiles."""
    rng = np.random.default_rng([72843, level])
    edge_values = rng.integers(-7, 8, len(board.seams))
    ports = np.array(board.seams)
    a, b = ports.T
    skew = sparse.csr_matrix(
        (np.r_[edge_values, -edge_values], (np.r_[a, b], np.r_[b, a])),
        shape=(board.n, board.n),
    )
    lap = board.laplacian().astype(np.int64)
    residual = lap @ skew - skew @ lap
    actual_norms = [
        int(skew.multiply(skew).sum()),
        int(residual.multiply(residual).sum()),
    ]
    local_norms = np.zeros(2, dtype=np.int64)
    by_edge = dict(zip(board.seams, edge_values, strict=True))
    maps = {}
    for owner in range(board.cells):
        offset = 12 * owner
        gates = tuple(np.flatnonzero(board.degree[offset : offset + 12] == 6))
        if gates not in maps:
            maps[gates] = local_map(gates)
        operator, weight, mass = maps[gates]
        internal = [by_edge[offset + i, offset + j] for i, j in cell().seams]
        external = []
        for g in gates:
            a = offset + g
            b = next(p for p in board.neighbours[a] if p // 12 != owner)
            external.append(by_edge[min(a, b), max(a, b)] * (1 if a < b else -1))
        column = np.array(internal + external)
        output = operator @ column
        local_norms += [column @ mass @ column, output @ (weight * output)]
    if list(local_norms) != actual_norms:
        raise ValueError("Full-board and per-cell current norm decompositions differ")
    return dict(
        input_frobenius_squared=actual_norms[0],
        commutator_laplacian_squared=actual_norms[1],
    )


def frame_and_meter(gates):
    """Infer the slow frame from local transition counts, with a chart-fixed orientation."""
    car = cell()
    probabilities = transitions(candidate.generator(car).toarray(), 0.25)
    rng = np.random.default_rng([74291, 65536, 0])
    labels = rng.permutation(car.n)
    renamed = probabilities[np.ix_(labels, labels)]
    counts = np.column_stack([rng.multinomial(65536, column) for column in renamed.T])
    rebuilt = reconstruct(counts)
    if rebuilt is None:
        raise ValueError("The earlier short-exposure angular instrument abstained")
    adjacency, measured = rebuilt
    order = np.argsort(labels)
    if not np.array_equal(adjacency[np.ix_(order, order)], car.adjacency().toarray()):
        raise ValueError("Reconstructed angular instrument has different seams")
    directions = measured[order]
    orientation = next(
        np.sign(value)
        for triple in combinations(range(12), 3)
        if abs(value := np.linalg.det(directions[list(triple)])) > 1e-8
    )
    edges = np.array(car.seams)
    first = next(e for e in car.seams if (e[0] in gates) != (e[1] in gates))
    crosses = orientation * np.cross(directions[edges[:, 0]], directions[edges[:, 1]])
    strengths = crosses[car.seams.index(first)] @ directions.T
    axis = int(np.argmax(np.round(abs(strengths), 10)))
    coeff = np.round(crosses @ directions[axis], 13)
    generator = np.zeros((12, 12), dtype=complex)
    for (a, b), c in zip(edges, coeff, strict=True):
        generator[a, b], generator[b, a] = 1j * c, -1j * c
    isolated = car.laplacian().toarray() / 2
    if np.max(abs(isolated @ generator - generator @ isolated)) > 1e-11:
        raise ValueError("Closed-cell circulation control failed")
    return first, axis, coeff, generator


def native_analyzer():
    """Sharp left/right after the native half exchange measure the signed current."""
    law = candidate.native_law
    source = [law.SRC.index((1, 0)), law.SRC.index((0, 1))]
    unitary = expm(-0.5j * np.pi * law.seam_term(candidate.KAPPA))
    cols = unitary[:, source]
    ideal = np.array([[1, 1j], [1j, 1]]) / sqrt(2)
    phase = cols[source[0], 0] / ideal[0, 0]
    np.testing.assert_allclose(cols[source, :], phase * ideal, atol=3e-13, rtol=0)
    np.testing.assert_allclose(np.delete(cols, source, axis=0), 0, atol=3e-13, rtol=0)
    difference = ideal.conj().T @ np.diag([1, -1]) @ ideal
    np.testing.assert_allclose(difference, [[0, 1j], [-1j, 0]], atol=1e-14, rtol=0)
    return ideal


def prepare(board, launch, relative_phase):
    """Prepare equal endpoint weights by half exchange and a native degree-reference phase."""
    a, b = launch
    law = candidate.native_law
    ports = [law.SRC.index((1, 0)), law.SRC.index((0, 1))]
    half = expm(-0.5j * np.pi * law.seam_term(candidate.KAPPA))
    prepared = half[ports, ports[0]]
    # The off-gateway rail has a unit degree-reference coefficient. The
    # pulse phase on that rail sets the relative phase at positive exposure.
    lower = 0 if board.degree[a] == 5 else 1
    if board.degree[launch[1 - lower]] != 6:
        raise ValueError("Source seam requires one native degree-reference rail")
    change = relative_phase - np.angle(prepared[1] / prepared[0])
    angle = change if lower == 1 else -change
    exposure = float(np.remainder(-angle, 2 * np.pi) / candidate.KAPPA)
    prepared[lower] *= np.exp(-1j * candidate.KAPPA * exposure)
    target = np.array([1.0, np.exp(1j * relative_phase)]) / sqrt(2)
    np.testing.assert_allclose(
        prepared / prepared[0], target / target[0], atol=5e-13, rtol=0
    )
    source = np.zeros(board.n, dtype=complex)
    source[list(launch)] = target
    return source, exposure


def counts_for(rng, probs):
    # Sampling on a fixed 53-bit uniform grid avoids platform-sensitive choices
    # between complementary binomial algorithms at equiprobable outcomes.
    p_left, p_right = np.maximum(probs, 0)
    draws = rng.random(SHOTS)
    left = int(np.count_nonzero(draws < p_left))
    right = int(np.count_nonzero(draws < p_left + p_right)) - left
    return [left, right, SHOTS - left - right]


def patch_read(wave, coeff, meter, analyzer, seed):
    local = wave[:12]
    rng = np.random.default_rng(seed)
    histograms = [
        counts_for(rng, abs(analyzer @ local[[a, b]]) ** 2) for a, b in cell().seams
    ]
    counts = np.array(histograms)
    response = float(coeff @ (counts[:, 0] - counts[:, 1]) / SHOTS)
    expectation = float(np.vdot(local, meter @ local).real)
    outside = float(np.sum(abs(wave[12:]) ** 2))
    observed_outside = int(np.count_nonzero(rng.random(SHOTS) < outside))
    return dict(
        estimate=response,
        exact_patch=expectation,
        outside_probability=outside,
        outside_count=observed_outside,
        histograms=histograms,
    )


def gateway_torques():
    board = federation(0)
    h = board.laplacian() / 2
    first, axis, coeff, meter = frame_and_meter(
        tuple(np.flatnonzero(board.degree[:12] == 6))
    )
    source, pulse = prepare(board, first, np.pi / 3)
    time = 3 / 16
    wave = expm_multiply(-1j * time * h, source)
    analyzer = native_analyzer()
    law = candidate.native_law
    lived = [law.SRC.index((1, 0)), law.SRC.index((0, 1))]
    responses = []
    home_gates = [edge for edge in board.external if edge[0] < 12]
    for e, (a, b) in enumerate(home_gates):
        readings = []
        for j, angle in enumerate((np.pi / 2, 3 * np.pi / 2)):
            full = expm(-1j * angle * law.seam_term(candidate.KAPPA))
            unitary = full[np.ix_(lived, lived)]
            swap = np.array([[0, 1], [1, 0]])
            even, odd = (np.eye(2) + swap) / 2, (np.eye(2) - swap) / 2
            np.testing.assert_allclose(
                unitary, even + np.exp(-1j * angle) * odd, atol=5e-13, rtol=0
            )
            exposed = wave.copy()
            exposed[[a, b]] = unitary @ wave[[a, b]]
            pair = []
            for k, owner in enumerate((a // 12, b // 12)):
                start = 12 * owner
                local = exposed[start : start + 12]
                read = patch_read(local, coeff, meter, analyzer, [SEED + 1, e, j, k])
                pair.append(
                    {
                        key: read[key]
                        for key in ("histograms", "exact_patch", "estimate")
                    }
                )
            readings.append(pair)
        exact = [
            (readings[0][k]["exact_patch"] - readings[1][k]["exact_patch"]) / 2
            for k in range(2)
        ]
        sampled = [
            (readings[0][k]["estimate"] - readings[1][k]["estimate"]) / 2
            for k in range(2)
        ]
        # Four independent ensembles with coefficient magnitude |c|/2.
        # Three common-owner sums and one total-home measurement share alpha.
        radius = sqrt(2 * float(coeff @ coeff) * log(2 * 4 / 0.01) / SHOTS)
        responses.append(
            dict(
                ports=[a, b],
                exact_owner_torques=exact,
                estimated_owner_torques=sampled,
                reciprocal_sum_interval=[sum(sampled) - radius, sum(sampled) + radius],
                half_exchange_reads=readings,
            )
        )
    slopes = []
    for step in (1 / 256, 1 / 512):
        forward = expm_multiply(-1j * step * h, wave)[:12]
        back = expm_multiply(1j * step * h, wave)[:12]
        slopes.append(
            float(
                (np.vdot(forward, meter @ forward) - np.vdot(back, meter @ back)).real
                / (2 * step)
            )
        )
    total = sum(r["exact_owner_torques"][0] for r in responses)
    observed = sum(r["estimated_owner_torques"][0] for r in responses)
    home_radius = sqrt(3 * float(coeff @ coeff) * log(2 * 4 / 0.01) / SHOTS)
    if abs(slopes[-1] - total) > 1e-5:
        raise ValueError(
            "The finite independent free slope disagrees with the exact contact sum"
        )
    # At release the other owner is entirely vacant. Its three angular reads
    # all annihilate the preparation, independently of its choice of chart.
    a, b = home_gates[0]
    launch_response = np.zeros(board.n, dtype=complex)
    launch_response[[a, b]] = (source[a] - source[b]) * np.array([1.0, -1.0]) / 2
    initial_torque = float(-2 * np.vdot(launch_response[:12], meter @ source[:12]).imag)
    if np.any(source[b // 12 * 12 : (b // 12 + 1) * 12]) or abs(initial_torque) < 0.1:
        raise ValueError("The vacant-neighbour frame control failed")
    return dict(
        cells=board.cells,
        source_ports=list(first),
        axis=axis,
        exposure=time,
        relative_phase_radians=np.pi / 3,
        reference_pulse_exposure=pulse,
        gateway_response=responses,
        initial_gateway_frame_control=dict(
            ports=[a, b],
            home_axis_torque=initial_torque,
            neighbouring_angular_vector=[0, 0, 0],
        ),
        sources_seed=SEED + 1,
        uninterrupted_regional_slopes=slopes,
        summed_native_home_torque=total,
        sampled_home_torque=observed,
        sampled_home_interval=[observed - home_radius, observed + home_radius],
    )


def experiment():
    analyzer = native_analyzer()
    rows = []
    # Twelve signed drift claims. Each uses a two-sided bounded-mean bound and
    # an independent one-sided outside-mass bound for the final appointment.
    alpha = 0.01 / (12 * 2)
    for level in range(3):
        board = federation(level)
        classes = incidence(board)
        gates = tuple(np.flatnonzero(board.degree[:12] == 6))
        launch, axis, coeff, meter = frame_and_meter(gates)
        h = board.laplacian() / 2
        native_difference = candidate.generator(board) - h
        native_difference.setdiag(native_difference.diagonal() - 6 * candidate.KAPPA)
        if native_difference.count_nonzero():
            raise ValueError(
                "Free-flight sector departed from the native one-load dynamics"
            )
        norm = float(max(abs(np.linalg.eigvalsh(meter))))
        radius = sqrt(4 * float(coeff @ coeff) * log(2 / alpha) / SHOTS)
        arms = []
        for sign in (1, -1):
            source, pulse = prepare(board, launch, 0 if sign == 1 else np.pi)
            reads = []
            for j, t in enumerate(TIMES):
                wave = expm_multiply(-1j * t * h, source)
                read = patch_read(
                    wave, coeff, meter, analyzer, [SEED, level, int(sign < 0), j]
                )
                if j:
                    leaked = read["outside_count"]
                    upper = (
                        beta.ppf(1 - alpha, leaked + 1, SHOTS - leaked)
                        if leaked < SHOTS
                        else 1
                    )
                    delta = read["estimate"] - reads[0]["estimate"]
                    half_width = float(radius + norm * upper)
                    full = sum(
                        float(np.vdot(v, meter @ v).real) for v in wave.reshape(-1, 12)
                    )
                    read.update(
                        estimated_change=delta,
                        simultaneous_interval=[delta - half_width, delta + half_width],
                        exact_whole_sphere_change=full,
                        max_omitted_expected_current=norm * read["outside_probability"],
                    )
                    if (
                        abs(full - read["exact_patch"])
                        > read["max_omitted_expected_current"] + 1e-11
                    ):
                        raise ValueError("Omitted-cell bound failed")
                read["exposure"] = t
                reads.append(read)
            arms.append(
                dict(relative_sign=sign, reference_pulse_exposure=pulse, readings=reads)
            )
        rows.append(
            dict(
                cells=board.cells,
                native_seams=len(board.seams),
                gateway_classes={
                    ",".join(map(str, k)): v for k, v in sorted(classes.items())
                },
                whole_board_integer_check=check_global_decomposition(board, level),
                source_ports=list(launch),
                axis=axis,
                weights=coeff.tolist(),
                complement_norm=norm,
                sampling_radius=radius,
                arms=arms,
            )
        )
    return portable(
        dict(
            exact_operator=exact_ward_bound(),
            seed=SEED,
            shots_per_setting=SHOTS,
            family_coverage=0.99,
            screens=rows,
            gateway_torques=gateway_torques(),
        )
    )


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = experiment()
    path = HERE / "directional.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        verify(json.loads(path.read_text()), result)
    print("Gateway triples:", result["exact_operator"]["gateway_triples"])
    for row in result["screens"]:
        print(
            row["cells"],
            [
                (arm["relative_sign"], arm["readings"][1]["simultaneous_interval"])
                for arm in row["arms"]
            ],
        )


if __name__ == "__main__":
    main()
