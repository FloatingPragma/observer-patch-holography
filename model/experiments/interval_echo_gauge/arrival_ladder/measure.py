"""A charge-selective, reference-compensated ladder on arrived ternary material."""

import argparse
from functools import cache
from itertools import permutations, product
import json
from pathlib import Path
from types import SimpleNamespace

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model.carrier import federation
from model.experiments.interval_echo_consensus.signed_effects.collar_storage import (
    MEMORIES,
    conserved_columns,
)
from model.experiments.interval_echo_encounter.conjugation import (
    defect_generator,
    free_seam,
)
from model.experiments.interval_echo_lepton_register.measure import count_bank, SHOTS

HERE = Path(__file__).resolve().parent
LAW = candidate.native_law
K = candidate.KAPPA
WORDS = list(product(range(3), repeat=3))
SOURCES = {
    "negative": (0, 1, 2),
    "neutral": (0, 2, 1),
    "positive": (2, 1, 2),
    "filled": (1, 1, 2),
    "vacancy_at_other_port": (1, 0, 2),
    "donor_depleted": (0, 1, 1),
}
T = 0.25
DELAYS = (0.25, 1, 4, 16)
RELEASE_SHOTS = 2**34


def layout():
    board = federation(0)
    exits = {u: v for u in range(12) for v in board.neighbours[u] if v >= 12}
    for incoming, x in permutations(sorted(exits), 2):
        if exits[incoming] // 12 == exits[x] // 12:
            continue
        for a in sorted(
            j for j in board.neighbours[x] if j < 12 and board.degree[j] == 5
        ):
            tips = sorted(set(board.neighbours[a]) & set(board.neighbours[x]))
            for s, middle in permutations(tips):
                if s >= 12 or board.degree[s] != 5:
                    continue
                for b in sorted(j for j in board.neighbours[a] if j < 12):
                    selected = {incoming, x, a, s, b, middle}
                    if len(selected) < 6:
                        continue
                    guards = sorted(set(range(12)) - selected)
                    for compensators in permutations(guards, 2):
                        offline = {x, a, s, b, *compensators}
                        local = nx.Graph(
                            (u, v)
                            for u, v in board.seams
                            if max(u, v) < 12 and {u, v}.isdisjoint(offline - {x})
                        )
                        if not nx.has_path(local, incoming, x):
                            continue
                        route = nx.shortest_path(local, incoming, x)
                        # Both deposited collars must remain parked during the relay.
                        if set(route) & {a, s}:
                            continue
                        dest = exits[x]
                        v = min(
                            j for j in board.neighbours[dest] if j // 12 == dest // 12
                        )
                        return dict(
                            source=exits[incoming],
                            entry=incoming,
                            x=x,
                            a=a,
                            guard=s,
                            donor=b,
                            compensators=list(compensators),
                            relay=route,
                            offline=sorted(offline),
                            neutral_path=[a, middle, x, dest, v],
                        )
    raise ArithmeticError("No guarded workstation on the carrier")


def admission(site, source):
    """Assemble the measured laboratory word after a capture or a local calibration."""
    board = federation(0)
    w = np.ones(board.n, dtype=int)
    x, a, b = site["x"], site["a"], site["donor"]
    w[site["guard"]] = w[b] = 2
    w[site["compensators"]] = 0
    if source == "neutral":
        w[x], w[a], w[b] = 0, 2, 1
    elif source != "filled":
        w[site["entry"]] = 2 if source == "positive" else 0
        if source == "donor_depleted":
            w[b] = w[site["compensators"][0]] = 1
        paths = [site["relay"]]
        if source == "vacancy_at_other_port":
            paths.append([x, site["neutral_path"][1], a])
        written = {frozenset((x, a)), frozenset((a, site["guard"]))}
        for route in paths:
            edges = list(zip(route, route[1:]))
            for u, v in edges + edges[:-1][::-1]:
                if v not in board.neighbours[u] or frozenset((u, v)) in written:
                    raise ArithmeticError(
                        "An admission relay uses an unavailable blank seam"
                    )
                w[u], w[v] = w[v], w[u]
    read = tuple(int(w[p]) for p in (x, a, b))
    if read != SOURCES[source]:
        raise ArithmeticError("Heralded material does not supply the calibrated word")
    return read


@cache
def stored_energies(kappa):
    exchange, contrast = LAW.local_term_parts()
    h = exchange + kappa * contrast
    embedded = np.kron(MEMORIES[:, [0]], np.eye(9))
    energies = np.diag(embedded.conj().T @ h @ embedded).real
    np.testing.assert_allclose(h @ embedded, embedded * energies, atol=8e-14)
    np.testing.assert_allclose(
        energies, [4 * kappa * (w != (1, 1)) for w in LAW.SRC], atol=8e-14
    )
    return energies.reshape(3, 3)


def selective_phase(phi, kappa):
    """Echo the physical (x,2) input, dwell with both stored collars, and undo it."""
    echo = LAW.echo()
    energy = stored_energies(kappa)
    active = np.array(LAW.SRC)
    phase = np.zeros((3, 3), complex)
    for x, a in product(range(3), repeat=2):
        incoming = LAW.SRC.index((x, 2))
        half = phi / (8 * K)
        inverse_half = (-phi / 2 % (2 * np.pi)) / (4 * K)
        dwell = np.exp(-1j * half * (energy[active[:, 0], a] + energy[active[:, 1], a]))
        temporary = echo[:, incoming] * np.tile(dwell, LAW.DIM // 9)
        final = (echo @ temporary) * np.exp(-1j * inverse_half * energy[x, a])
        phase[x, a] = final[incoming]
        np.testing.assert_allclose(np.delete(final, incoming), 0, atol=3e-14)
    return phase


def interaction(phi, kappa):
    opening = expm(-0.5j * np.pi * free_seam(kappa))
    closing = expm(-1.5j * np.pi * free_seam(kappa))
    phases = np.repeat(selective_phase(phi, kappa).ravel(), 3)
    return np.kron(np.eye(3), closing) @ np.diag(phases) @ np.kron(np.eye(3), opening)


def aperture_map(ports, angle=np.pi, kappa=K):
    """Return amplitudes on single-defect rails before any sharp position read.

    A blank work port starts and ends filled. Two work rails start (1,2)
    for the charged input or (2,1) for the neutral input. The source rails
    carry one vacancy; the neutral output includes the common target double.
    """
    site = layout()
    board = federation(0)
    ports = tuple(int(p) for p in ports)
    parked = set(site["offline"]) - {site["x"]}
    if (
        not ports
        or len(set(ports)) != len(ports)
        or any(p in parked or p == site["x"] or not 0 <= p < board.n for p in ports)
    ):
        raise ValueError("Use distinct material rails disjoint from the workstation")
    graph = nx.Graph(e for e in board.seams if set(e).isdisjoint(parked))
    # The full contact, tested on all blank loads, is the plain transposition
    # including its relative phase. The unguarded path is used only for these
    # addressed relay pulses; each internal transit register is restored.
    native_swap = expm(-1j * np.pi * free_seam(kappa))
    permutation = [LAW.SRC.index((b, a)) for a, b in LAW.SRC]
    np.testing.assert_allclose(native_swap, np.eye(9)[:, permutation], atol=1e-13)
    routes = []
    for p in ports:
        path = nx.shortest_path(graph, p, site["x"])
        edges = list(zip(path, path[1:]))
        word = edges + edges[-2::-1] if len(edges) > 1 else edges
        labels = list(range(board.n))
        for i, j in word:
            labels[i], labels[j] = labels[j], labels[i]
        expected = list(range(board.n))
        expected[p], expected[site["x"]] = expected[site["x"]], expected[p]
        if labels != expected:
            raise ArithmeticError("Relay left a record in its transit rails")
        routes.append(path)
    phase = selective_phase(angle, kappa)
    opening = expm(-0.5j * np.pi * free_seam(kappa))
    closing = expm(-1.5j * np.pi * free_seam(kappa))
    pair = [LAW.SRC.index((1, 2)), LAW.SRC.index((2, 1))]
    keys = (
        "charged_to_charged",
        "neutral_to_charged",
        "charged_to_neutral",
        "neutral_to_neutral",
    )
    result = {key: np.zeros((len(ports), len(ports)), complex) for key in keys}
    singlet = []
    for j in range(len(ports)):
        for value in (0, 2):
            dwells = np.ones(9, complex)
            for i in range(len(ports)):
                source = value if i == j else 1
                dwells *= np.repeat(phase[source], 3)
            gate = closing @ np.diag(dwells) @ opening
            code = gate[np.ix_(pair, pair)]
            np.testing.assert_allclose(
                gate[:, pair], np.eye(9)[:, pair] @ code, atol=1e-12
            )
            if value == 0:
                for key, amplitude in zip(keys, code.ravel(), strict=True):
                    result[key][j, j] = amplitude
            else:
                singlet.append(code[0, 0])
    result.update(routes=routes, singlet_survival=np.diag(singlet))
    return result


def aperture_records(rng):
    board = federation(0)
    owner = layout()["source"] // 12
    ports = list(range(12 * owner, 12 * owner + 12))
    result = aperture_map(ports)
    convert = result["charged_to_neutral"]
    # Prepare and reread a relative phase on an actual adjacent pair of rails.
    left, right = next(
        (p, q) for p in ports for q in board.neighbours[p] if q in ports and p < q
    )
    pair = [LAW.SRC.index((0, 1)), LAW.SRC.index((1, 0))]
    h = free_seam(K)[np.ix_(pair, pair)]
    launch = expm(-0.5j * np.pi * h)[:, 0]
    receiver = expm(-1.5j * np.pi * h)
    ix = [ports.index(left), ports.index(right)]
    joint = convert[np.ix_(ix, ix)] @ launch
    coherent = abs(receiver @ joint) ** 2
    recorded = abs(receiver) ** 2 @ abs(joint) ** 2

    def sample(p):
        return count_bank(rng, np.r_[p, max(0, 1 - sum(p))]).tolist()

    return dict(
        ports=ports,
        routes=result["routes"],
        interferometer=[left, right],
        amplitudes={
            key: [[float(v.real), float(v.imag)] for v in np.diag(matrix)]
            for key, matrix in result.items()
            if key != "routes"
        },
        coherent_counts=sample(coherent),
        address_recorded_counts=sample(recorded),
        min_conversion_probability=float(
            min(np.linalg.eigvalsh(convert.conj().T @ convert))
        ),
        max_conversion_probability=float(
            max(np.linalg.eigvalsh(convert.conj().T @ convert))
        ),
    )


def read_loads():
    """One blank reference and the signed contrast meter resolve all three loads."""
    echo = LAW.echo()
    outcomes = []
    for load in range(3):
        c = LAW.SRC.index((load, 0))
        probs = []
        for d in range(-2, 3):
            bright = np.repeat([False] + [v == d for _, v in LAW.LABELS], 9)
            resolved = echo[:, c] * bright
            probs.append(float(np.vdot(resolved, resolved).real))
            if d == load:
                np.testing.assert_allclose(
                    echo @ resolved, np.eye(LAW.DIM)[:, c], atol=3e-14
                )
        np.testing.assert_allclose(probs, np.eye(5)[load + 2], atol=3e-14)
        outcomes.append(int(np.argmax(probs)) - 2)
    return outcomes


def arrival(site, sign, kappa):
    board = federation(0)
    edges = [e for e in board.seams if set(e).isdisjoint(site["offline"])]
    active = SimpleNamespace(n=board.n, seams=edges, degree=board.degree)
    h = defect_generator(active, sign, kappa, 0)
    launch = np.zeros(board.n)
    launch[site["source"]] = 1
    state = expm_multiply(-1j * T * h, launch)
    return float(abs(state[site["entry"]]) ** 2)


def path_generator(word, vertices, kappa):
    """Reconstruct the five-port open path from actual blank seam columns."""
    board = federation(0)
    degrees = board.degree[vertices]
    basis = sorted(set(permutations(word)))
    where = {w: i for i, w in enumerate(basis)}
    h = np.zeros((len(basis), len(basis)))
    seam = free_seam(kappa)
    for j, w in enumerate(basis):
        h[j, j] += kappa * np.dot(6 - degrees, np.square(w) - 1)
        for k in range(len(vertices) - 1):
            if vertices[k + 1] not in board.neighbours[vertices[k]]:
                raise ArithmeticError("Five-port receiver path leaves the carrier")
            pair = w[k : k + 2]
            i = LAW.SRC.index(pair)
            h[j, j] += seam[i, i]
            if pair[0] != pair[1]:
                h[where[w[:k] + pair[::-1] + w[k + 2 :]], j] += seam[
                    LAW.SRC.index(pair[::-1]), i
                ]
    return h, basis


def release(site, input_word, kappa, t):
    # The rest of the workstation is parked; only these five ports and the four path contacts run.
    x, a, _donor = input_word
    initial = (a, 1, x, 1, 1)
    h, basis = path_generator(initial, site["neutral_path"], kappa)
    amps = expm(-1j * t * h)[:, basis.index(initial)]
    terminal = np.zeros(9)
    for w, amp in zip(basis, amps, strict=True):
        terminal[LAW.SRC.index(w[-2:])] += abs(amp) ** 2
    return terminal


def measure():
    # The source-preparation check explicitly uses the native 99-state collar encoder.
    conserved_columns()
    site = layout()
    rng = np.random.default_rng(5403)
    banks = []
    flights = []
    for kappa in (K, K / 2):
        rates = [arrival(site, sign, kappa) for sign in (-1, 1)]
        flights.append(dict(kappa=kappa, negative=rates[0], positive=rates[1]))
        for angle in (0, np.pi / 2, np.pi):
            unitary = interaction(angle, kappa)
            np.testing.assert_allclose(
                unitary.conj().T @ unitary, np.eye(27), atol=7e-13
            )
            for source in SOURCES:
                w = admission(site, source)
                # The neutral calibration and filled sham are prepared at the reader.
                efficiency = (
                    (rates[1] if source == "positive" else rates[0])
                    if source
                    in (
                        "positive",
                        "negative",
                        "donor_depleted",
                        "vacancy_at_other_port",
                    )
                    else 1.0
                )
                probabilities = abs(unitary[:, WORDS.index(w)]) ** 2
                complete = np.r_[1 - efficiency, efficiency * probabilities]
                counts = count_bank(rng, complete)
                banks.append(
                    dict(
                        kappa=kappa,
                        angle_over_pi=float(angle / np.pi),
                        source=source,
                        per_herald_probabilities=probabilities.tolist(),
                        missed_and_all_27_counts=counts.tolist(),
                    )
                )
    terminal = []
    # Preserve the donor outcome in every branch. Reading it independently removes
    # interference between the charged and converted alternatives in the outgoing patch.
    for source, phi in (
        ("negative", np.pi),
        ("negative", 0),
        ("positive", np.pi),
        ("neutral", 0),
        ("filled", np.pi),
        ("donor_depleted", np.pi),
    ):
        w = admission(site, source)
        outcomes = abs(interaction(phi, K)[:, WORDS.index(w)]) ** 2
        for t in DELAYS:
            prob = sum(
                p * release(site, out, K, t)
                for out, p in zip(WORDS, outcomes, strict=True)
                if p > 1e-22
            )
            efficiency = (
                flights[0]["positive"]
                if source == "positive"
                else flights[0]["negative"]
                if source in ("negative", "donor_depleted")
                else 1
            )
            inclusive = np.r_[1 - efficiency, efficiency * prob]
            weights = np.round(inclusive, 14)
            counts = rng.multinomial(RELEASE_SHOTS, weights / weights.sum())
            pair_ix = [LAW.SRC.index((0, 2)), LAW.SRC.index((2, 0))]
            terminal.append(
                dict(
                    source=source,
                    angle_over_pi=float(phi / np.pi),
                    time=t,
                    per_herald_probabilities=prob.tolist(),
                    missed_and_all_9_counts=counts.tolist(),
                    neutral_count=int(sum(counts[1:][pair_ix])),
                    per_herald_neutral=float(sum(prob[pair_ix])),
                )
            )
    # Rank-one phase reflections at pi exchange the two codes; a round trip recovers the source.
    ladder = interaction(np.pi, K)
    ids = [WORDS.index(SOURCES[s]) for s in ("negative", "neutral", "positive")]
    block = ladder[np.ix_(ids, ids)]
    return dict(
        schema=1,
        site=site,
        duration_in=T,
        shots_per_bank=SHOTS,
        read_loads=read_loads(),
        flight_probabilities=flights,
        code_words=WORDS,
        admitted_words={s: admission(site, s) for s in SOURCES},
        joint_charge_block_real=block.real.tolist(),
        joint_charge_block_imag=block.imag.tolist(),
        two_pulse_negative_return=float(abs((ladder @ ladder)[ids[0], ids[0]]) ** 2),
        banks=banks,
        independent_receiver=terminal,
        release_shots_per_bank=RELEASE_SHOTS,
        coherent_aperture=aperture_records(rng),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = measure()
    out = json.dumps(result, indent=2) + "\n"
    if args.write:
        (HERE / "results.json").write_text(out)
    else:
        print(out, end="")


if __name__ == "__main__":
    main()
