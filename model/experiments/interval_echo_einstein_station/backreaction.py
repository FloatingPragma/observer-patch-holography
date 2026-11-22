"""Distinguish the stored energy of a prepared collar from its action on material probes."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

law = candidate.native_law
KAPPA = candidate.KAPPA
RECEIPT = Path(__file__).with_name("backreaction.json")


def coin_basis(magnitude, parity):
    """Two unnormalized integer eigenvectors on the occupied memory, one for each coin."""
    vectors = np.zeros((law.N_LAB, 2), dtype=int)
    for coin in (0, 1):
        vectors[law.LABELS.index((coin, magnitude)), coin] = 1
        if magnitude:
            vectors[law.LABELS.index((coin, -magnitude)), coin] = parity
    return vectors


def collar_vector(magnitude, parity, coin_phase):
    lab = coin_basis(magnitude, parity) @ np.array([1, coin_phase])
    return np.r_[0, lab] / np.linalg.norm(lab)


def classify():
    """Compute the complete scalar-action subspaces on occupied collars at the pin."""
    import sympy as sp

    W = law._write_pattern()  # sqrt(2) times the native write, exactly integral
    reflection, contrast = law._label_ops()
    twice_record_h = (
        np.eye(law.N_LAB, dtype=int) - reflection + int(2 * KAPPA) * contrast
    )
    rows = []
    exchange, agreement = law.local_term_parts()
    native = exchange + KAPPA * agreement
    for magnitude, parity in [(0, 1), (1, 1), (1, -1), (2, 1), (2, -1)]:
        basis = coin_basis(magnitude, parity)
        energy = int(KAPPA * magnitude**2 + (1 - parity) / 2)
        if np.any(twice_record_h @ basis != 2 * energy * basis):
            raise ValueError("Record eigenbasis failed")
        # For an eigenvector with positive E, native scalar action on every
        # load is equivalent to W* (record tensor load) = 0 for every load.
        restrictions = np.column_stack(
            [
                (W.T @ np.kron(column[:, None], np.eye(law.N_SRC, dtype=int))).ravel()
                for column in basis.T
            ]
        )
        gram = restrictions.T @ restrictions
        constraints = sp.Matrix(energy * gram)
        nullspace = [list(map(int, v)) for v in constraints.nullspace()]
        for coins in nullspace:
            record = np.r_[0, basis @ coins]
            columns = np.kron(record[:, None], np.eye(law.N_SRC))
            np.testing.assert_allclose(
                native @ columns, energy * columns, atol=2e-13, rtol=0
            )
        rows.append(
            dict(
                energy=energy,
                write_overlap_gram=gram.tolist(),
                scalar_action_coin_vectors=nullspace,
            )
        )
    if [len(r["scalar_action_coin_vectors"]) for r in rows] != [2, 0, 0, 1, 1]:
        raise ValueError("Unexpected scalar-action spectrum at the pinned coupling")
    return rows


def source_reads():
    exchange, contrast = law.local_term_parts()
    h = exchange + KAPPA * contrast
    vacant = np.eye(law.N_SRC)[:, law.SRC.index((0, 0))]
    zero = np.kron(collar_vector(0, 1, -1), vacant)
    records = []
    for magnitude, parity in [(1, 1), (1, -1), (2, 1), (2, -1)]:
        record = np.kron(collar_vector(magnitude, parity, -1), vacant)
        preparation = (zero + record) / np.sqrt(2)
        arms = {}
        for name, generator, kappa in (
            ("native", h, KAPPA),
            ("exchange_only", exchange, 0),
        ):
            readings = []
            for duration in (1 / 32, 3 / 64):
                state = expm_multiply(-1j * duration * generator, preparation)
                x = abs(np.vdot(preparation, state)) ** 2
                y = abs(np.vdot((zero - 1j * record) / np.sqrt(2), state)) ** 2
                gap = np.arctan2(2 * y - 1, 2 * x - 1) / duration
                np.testing.assert_allclose(
                    gap, kappa * magnitude**2 + (1 - parity) / 2, atol=2e-12
                )
                readings.append(
                    dict(exposure=duration, x=float(x), y=float(y), gap=float(gap))
                )
            arms[name] = readings
        records.append(dict(magnitude=magnitude, sign_parity=parity, arms=arms))
    return records


def lift_one_load(board, edge):
    """Keep all eleven states on the marked seam during a full-sphere unit-load flight."""
    a, b = edge
    memories = 1 + law.N_LAB
    free = candidate.generator(board).tocsr()
    unit = [law.SRC.index(pair) for pair in [(1, 0), (0, 1)]]
    edge_h = law.seam_term(KAPPA)[np.ix_(unit, unit)]
    marked_ports = np.array(edge)
    seam = sparse.coo_matrix(
        (edge_h.ravel(), (np.repeat(marked_ports, 2), np.tile(marked_ports, 2))),
        shape=free.shape,
    )
    external = sparse.kron(free - seam, sparse.eye(memories), format="csr")
    exchange, contrast = law.local_term_parts()
    h = exchange + KAPPA * contrast
    offsets = np.arange(memories) * law.N_SRC
    zero = h[np.ix_(offsets, offsets)]
    # The local vacuum block acts also when the probe is elsewhere.
    total = external + sparse.kron(sparse.eye(board.n), zero, format="csr")
    local = np.concatenate([offsets + i for i in unit])
    correction = h[np.ix_(local, local)] - np.kron(np.eye(2), zero)
    indices = np.r_[
        a * memories + np.arange(memories), b * memories + np.arange(memories)
    ]
    total += sparse.coo_matrix(
        (
            correction.ravel(),
            (np.repeat(indices, len(indices)), np.tile(indices, len(indices))),
        ),
        shape=total.shape,
    ).tocsr()
    return total


def flight(level):
    board = federation(level)
    edge = min(board.external)
    receiver = list(range((edge[1] // 12) * 12, (edge[1] // 12 + 1) * 12))
    memories = 1 + law.N_LAB
    sources = {
        "blank": np.eye(memories)[:, 0],
        "zero": collar_vector(0, 1, -1),
        "even_double": collar_vector(2, 1, -1),
        "odd_double": collar_vector(2, -1, -1),
        "ordinary_even_double": collar_vector(2, 1, 1),
    }
    states = np.zeros((board.n, memories, len(sources)), dtype=complex)
    states[edge[0]] = np.column_stack(list(sources.values()))
    states = states.reshape(board.n * memories, len(sources))
    total = lift_one_load(board, edge)
    # Remove a common scalar to reduce the cost of exponentiation.
    total = total - sparse.eye(total.shape[0]) * (6 * KAPPA)
    appointments = []
    for duration in (3 / 4, 3 / 2):
        evolved = expm_multiply(-1j * duration * total, states).reshape(
            board.n, memories, len(sources)
        )
        loads = (abs(evolved) ** 2).sum(axis=1)
        np.testing.assert_allclose(loads.sum(axis=0), 1, atol=3e-12, rtol=0)
        patch = loads[receiver].T
        outcomes = np.column_stack((patch, 1 - patch.sum(axis=1)))
        np.testing.assert_allclose(
            outcomes[1:],
            np.broadcast_to(outcomes[1], outcomes[1:].shape),
            rtol=0,
            atol=3e-12,
        )
        distance = 0.5 * np.abs(outcomes[0] - outcomes[1]).sum()
        if distance < 0.005:
            raise ValueError("The blank-memory sensitivity control lost its signal")
        appointments.append(
            dict(
                exposure=duration,
                blank_vs_marked_tv=float(distance),
                patch_plus_quiet={
                    key: list(row) for key, row in zip(sources, outcomes)
                },
            )
        )
    return dict(
        cells=board.cells,
        launch=edge[0],
        seam=list(edge),
        receiver_ports=receiver,
        appointments=appointments,
    )


def ordered_collar_control(
    board, a, b, c, native, energies, label, splitter, inverse, shots
):
    """A directed change with every unordered physical source word held fixed."""
    from scipy.linalg import expm

    if board.degree[a] != 5 or board.degree[c] != 5:
        raise ValueError("Both outer ports need a native degree reference.")
    gaps = np.array([0.0, energies[(1, 1)]["gap"], energies[(-1, 1)]["gap"]])
    records = [label(0, 1), label(2, 1), label(2, -1)]
    # Orient each seam with its outer, degree-five port first. With the shared
    # port parked, sign-parity reversal of an oriented odd mode is a global phase.
    both = np.eye(law.N_SRC)[:, law.SRC.index((1, 1))]
    checker = np.tile([u**2 for u, _ in law.SRC], 1 + law.N_LAB)
    read_phase = np.exp(-1j * np.pi * checker / 8)
    reader = expm(-1j * np.pi / 2 * native) @ (read_phase[:, None] * law.echo())
    loads11 = np.column_stack([np.kron(record, both) for record in records])
    words = [law.SRC.index(v) for v in ((1, 1), (0, 2), (2, 0))]
    outcome = abs(reader @ loads11) ** 2
    np.testing.assert_allclose(outcome[words], np.eye(3), atol=2e-13, rtol=0)
    np.testing.assert_allclose(outcome.sum(axis=0), 1, atol=2e-13, rtol=0)
    # Preparation reverses this detector with positive-duration native gates.
    undo_phase = np.exp(-1j * (2 * np.pi - np.pi / 8) * checker)
    prep = law.echo() @ (undo_phase[:, None] * expm(-1j * 3 * np.pi / 2 * native))
    np.testing.assert_allclose(
        abs(loads11.T @ prep[:, words]) ** 2, np.eye(3), atol=2e-13
    )
    spare = []
    for v in (a, b, c):
        spare.append(
            min(
                u
                for u in board.neighbours[v]
                if u // 12 == a // 12 and u not in (a, b, c) and u not in spare
            )
        )
    # Feed or clear a unit on each endpoint by swapping with its private spare;
    # address one collar at a time and park the other, including their shared port.
    swap = expm(-1j * np.pi * native)
    single = [law.SRC.index(v) for v in ((1, 0), (0, 1))]
    np.testing.assert_allclose(
        abs(swap[np.ix_(single, single)]) ** 2, [[0, 1], [1, 0]], atol=2e-13
    )
    read_time = np.pi / (3 * (gaps[2] - gaps[1]))
    contact = expm(-1j * read_time * native)
    parked_one = np.eye(law.N_SRC)[:, law.SRC.index((0, 1))]
    single_controls = np.eye(law.N_SRC)[
        :, [law.SRC.index(w) for w in ((0, 0), (1, 0), (0, 1))]
    ]
    contact_phases = []
    for record in records:
        states = np.column_stack(
            [np.kron(record, loads) for loads in (both, parked_one)]
        )
        phases = np.diag(states.T @ contact @ states)
        np.testing.assert_allclose(
            contact @ states, states * phases, atol=3e-13, rtol=0
        )
        control = np.kron(record[:, None], single_controls)
        np.testing.assert_allclose(
            contact @ control, control * phases[1], atol=3e-13, rtol=0
        )
        contact_phases.append(phases)
    contact_phases = np.array(contact_phases)
    output = {}
    one_load = {}
    for setting, ref_time in (("cosine", 0.0), ("quadrature", np.pi / (2 * KAPPA))):
        bias = np.exp(-1j * KAPPA * (6 - board.degree[a]) * ref_time * np.array([1, 0]))
        one_load[setting] = (abs(inverse @ (bias * splitter[:, 0])) ** 2).tolist()
        response = np.empty((3, 3, 2))
        for i, first in enumerate(contact_phases):
            for j, second in enumerate(contact_phases):
                final = inverse @ (bias * first * second[::-1] * splitter[:, 0])
                response[i, j] = abs(final) ** 2
        np.testing.assert_allclose(response.sum(axis=2), 1.0, atol=3e-13, rtol=0)
        output[setting] = response
    np.testing.assert_allclose(
        output["cosine"], output["cosine"].transpose(1, 0, 2), atol=3e-13
    )

    forward = np.roll(np.eye(3, dtype=int), 1, axis=1)
    reference_units = 4 * np.eye(3, dtype=int) + 6 * forward + 2 * forward.T
    total_energy = np.rint(gaps[:, None] + gaps[None, :]).astype(int)
    np.testing.assert_allclose(total_energy, gaps[:, None] + gaps[None, :], atol=3e-12)
    hist = np.bincount(total_energy.ravel(), weights=reference_units.ravel())
    reference = reference_units / 36
    score = np.log((1 / 9) / reference)
    rng = np.random.default_rng(2026112209)

    def draw(mix):
        return rng.multinomial(shots, mix.ravel()).reshape(3, 3)

    calibration_counts = draw(reference)
    frequency = calibration_counts / shots
    calibrated_score = np.log(
        np.outer(frequency.sum(axis=1), frequency.sum(axis=0)) / frequency
    )

    # Conditional categorical sources compile to two successive binary native
    # dimer interrogations with supplied reset and classical feed-forward.
    def check_draw(prob):
        p0, p1, p2 = prob
        for p in (p0, p1 / (p1 + p2)):
            step = expm(-2j * np.arcsin(np.sqrt(p)) * native)
            np.testing.assert_allclose(
                abs(step[single, single[0]]) ** 2, [1 - p, p], atol=2e-13
            )

    check_draw(np.ones(3) / 3)
    arms = []
    for direction in (-1, 0, 1):
        units = reference_units + direction * (forward.T - forward)
        np.testing.assert_array_equal(units.sum(axis=0), [12, 12, 12])
        np.testing.assert_array_equal(units.sum(axis=1), [12, 12, 12])
        np.testing.assert_array_equal(
            units + units.T, reference_units + reference_units.T
        )
        np.testing.assert_array_equal(
            np.bincount(total_energy.ravel(), weights=units.ravel()), hist
        )
        mix = units / 36
        for row in mix:
            check_draw(row / row.sum())
        charge = float(np.sum((mix - reference) * score))
        np.testing.assert_allclose(charge, direction * np.log(3) / 12, atol=2e-14)
        counts = draw(mix)
        readings = {}
        for key, response in output.items():
            probs = np.einsum("ij,ijk->k", mix, response)
            readings[key] = dict(
                probabilities=probs.tolist(),
                counts=rng.multinomial(shots, probs).tolist(),
            )
        arms.append(
            dict(
                direction=direction,
                source_units_of_1_over_36=units.tolist(),
                charge_change=charge,
                sampled_charge_change=float(
                    np.sum((counts / shots - frequency) * calibrated_score)
                ),
                source_word_counts=counts.tolist(),
                readers=readings,
            )
        )
    signed_return_change = (
        arms[2]["readers"]["quadrature"]["probabilities"][0]
        - arms[0]["readers"]["quadrature"]["probabilities"][0]
    )
    np.testing.assert_allclose(
        signed_return_change, np.sqrt(3) / 18, rtol=0, atol=3e-13
    )
    return dict(
        source_energies=gaps.tolist(),
        total_energy_units={str(k): int(hist[k]) for k in np.flatnonzero(hist)},
        detector_word_order=["1,1", "0,2", "2,0"],
        spare_ports=spare,
        contact_exposure=read_time,
        quadrature_reference_port=a,
        quadrature_reference_exposure=np.pi / (2 * KAPPA),
        phase_source_read_exposure=np.pi / (8 * KAPPA),
        outcomes=["return on AB", "other on BC"],
        one_load_outcomes_a_or_c=one_load,
        shots_per_setting=shots,
        seed=2026112209,
        calibration_counts=calibration_counts.tolist(),
        calibrated_scores=calibrated_score.tolist(),
        arms=arms,
    )


def correlated_contacts():
    """Change a two-collar history at fixed marginals and compare three local probes."""
    from scipy.linalg import expm
    from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
    from model.experiments.material_inertia.dark_collar import (
        common_pair_generator,
        label,
        native_readouts,
    )
    from model.experiments.material_inertia.federation import propagate

    board = federation(0)
    b, opposite = min(board.external)
    a = min(
        i for i in board.neighbours[b] if i // 12 == b // 12 and board.degree[i] == 5
    )
    c = min(
        i
        for i in set(board.neighbours[a]) & set(board.neighbours[b])
        if i // 12 == b // 12
    )
    edges = [(a, b), (b, c)]
    exchange, contrast = law.local_term_parts()
    native = exchange + KAPPA * contrast
    binary = [law.SRC.index((u, v)) for u, v in ((0, 0), (0, 1), (1, 0), (1, 1))]
    loads = np.eye(law.N_SRC)[:, binary]
    energies, _, _ = native_readouts()
    energy = energies[(1, 1)]["gap"]
    read_time = np.pi / energy
    for parity in (1, -1):
        for bit in (0, 1):
            lift = np.kron(label(2 * bit, 1, parity)[:, None], loads)
            active = np.diag([0, 0, 0, int(parity == 1)])
            action = energy * bit * (np.eye(4) - active)
            np.testing.assert_allclose(native @ lift, lift @ action, atol=1e-13)
    # The first and third triangle ports remain connected by a blank native seam.
    word_pair = [law.SRC.index(w) for w in ((1, 0), (0, 1))]
    splitter = expm(-1j * np.pi / 2 * native)[np.ix_(word_pair, word_pair)]
    inverse = expm(-1j * 3 * np.pi / 2 * native)[np.ix_(word_pair, word_pair)]
    start = splitter[:, 0]
    sectors = PairSector(board)
    marks = [sectors.index(*edge) for edge in edges]
    free = sectors.hamiltonian(KAPPA)
    base = sum((common_pair_generator(sectors, e) - free for e in edges), free)
    launch = np.zeros(sectors.P, complex)
    launch[marks] = start
    remote = (sectors.aa // 12 == opposite // 12) | (sectors.bb // 12 == opposite // 12)
    if np.any(remote[marks]):
        raise ValueError("Remote and marked outcomes must be disjoint")
    words = np.array([(0, 0), (0, 1), (1, 0), (1, 1)])
    # Fixed comparator log(P_product/P_ref) with a blank boundary record before flight.
    reference = np.array([3, 1, 1, 3]) / 8
    score = np.log(0.25 / reference)
    appointments = (0.25, 1.0, 2.0)
    outputs = []
    fringes = []
    for z, old in words:
        diagonal = np.zeros(sectors.P)
        diagonal[marks] = -energy * np.array([z, old])
        h = base + sparse.diags(diagonal)
        state = launch
        previous = 0
        tab = []
        for t in appointments:
            state = propagate(h, state, t - previous)
            prob = abs(state) ** 2
            np.testing.assert_allclose(prob.sum(), 1, rtol=0, atol=3e-10)
            loc = [*prob[marks], prob[remote].sum()]
            tab.append([*loc, 1 - sum(loc)])
            previous = t
        outputs.append(tab)
        phases = np.exp(1j * energy * read_time * np.array([z, old]))
        fringes.append(abs((inverse @ (phases * start))[0]) ** 2)
    np.testing.assert_allclose(fringes, [1, 0, 0, 1], atol=2e-13)
    # The antisymmetric memory control has the base Hamiltonian in every joint branch.
    one = candidate.generator(board).tocsr()
    one_local = law.seam_term(KAPPA)[np.ix_(word_pair, word_pair)]
    for u, v in edges:
        one -= sparse.coo_matrix(
            (one_local.ravel(), ([u, u, v, v], [u, v, u, v])), shape=one.shape
        ).tocsr()
    single = np.zeros(board.n, complex)
    single[[a, c]] = start
    one_remote = []
    t_prev = 0
    for t in appointments:
        single = propagate(one, single, t - t_prev)
        np.testing.assert_allclose(np.vdot(single, single), 1, atol=1e-12)
        one_remote.append(
            float((abs(single.reshape(-1, 12)[opposite // 12]) ** 2).sum())
        )
        t_prev = t
    arms = []
    rng = np.random.default_rng(2026112203)
    shots = 131072
    calibration_counts = rng.multinomial(shots, reference)
    calibrated = calibration_counts.reshape(2, 2) / shots
    independent = np.outer(calibrated.sum(axis=1), calibrated.sum(axis=0))
    calibrated_score = np.log(independent / calibrated).ravel()
    for p in (1 / 8, 1 / 4, 3 / 8):
        mix = np.array([1 - p, p, p, 1 - p]) / 2
        held_reads = rng.multinomial(shots, mix)
        np.testing.assert_allclose(mix @ words, [0.5, 0.5], atol=0)
        patch = np.einsum("w,wto->to", mix, outputs)
        fringe_return = float(mix @ fringes)
        # Fresh independent histories for each material interrogation.
        arms.append(
            dict(
                flip_probability=p,
                signed_charge_change=float((mix - reference) @ score),
                source_word_counts=held_reads.tolist(),
                measured_charge_change=float(
                    (held_reads / shots - calibration_counts / shots) @ calibrated_score
                ),
                isolated_interference_return=fringe_return,
                isolated_return_count=int(rng.binomial(shots, fringe_return)),
                free_patch_probabilities=patch.tolist(),
                free_patch_counts=[
                    rng.multinomial(shots, row).tolist() for row in patch
                ],
            )
        )
    return dict(
        ordered_three_levels=ordered_collar_control(
            board, a, b, c, native, energies, label, splitter, inverse, shots
        ),
        ports=[a, b, c],
        remote_gateway=opposite,
        states=sectors.P,
        native_energy_per_mark=energy,
        mean_total_stored_energy=energy,
        exposures=list(appointments),
        shots_per_setting=shots,
        seed=2026112203,
        joint_reference=reference.tolist(),
        log_reference_score=score.tolist(),
        source_calibration_counts=calibration_counts.tolist(),
        source_calibrated_scores=calibrated_score.tolist(),
        interference_exposure=read_time,
        one_load_remote=one_remote,
        branch_order=words.tolist(),
        free_outcomes=["both on Z", "both on C", "any at remote cell", "other"],
        conditional_free_probabilities=np.asarray(outputs).tolist(),
        arms=arms,
    )


def run():
    return dict(
        correlated_contacts=correlated_contacts(),
        scalar_sectors=classify(),
        collar_ramsey=source_reads(),
        flights=[flight(level) for level in (0, 1, 2)],
    )


def compare(actual, saved):
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            raise ValueError("Receipt fields differ")
        for key in actual:
            compare(actual[key], saved[key])
    elif isinstance(actual, list):
        for x, y in zip(actual, saved, strict=True):
            compare(x, y)
    elif isinstance(actual, float):
        np.testing.assert_allclose(actual, saved, atol=1e-13, rtol=2e-10)
    elif actual != saved:
        raise ValueError((actual, saved))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    observed = run()
    if args.write:
        RECEIPT.write_text(json.dumps(observed, indent=2, allow_nan=False) + "\n")
    else:
        compare(observed, json.loads(RECEIPT.read_text()))
    print(
        "Scalar-action dimensions:",
        [len(r["scalar_action_coin_vectors"]) for r in observed["scalar_sectors"]],
    )
    print(
        "Patch sensitivity:",
        [
            [p["blank_vs_marked_tv"] for p in f["appointments"]]
            for f in observed["flights"]
        ],
    )
