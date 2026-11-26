"""Material parity-phase programs and the two stations' interference readings."""

from fractions import Fraction
from itertools import combinations, permutations
from math import log, pi, sqrt

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_consensus.signed_effects import (
    record_transfer,
    sequential,
)

from .. import coherent, material

SHOTS = 4096
SEED = 550056
WORDS = np.arange(32)


def charges(word):
    colour, weak = (word & 7).bit_count(), (word >> 3).bit_count()
    b3 = (0, 1, -1, 0)[colour]
    ell = int(colour == 3 and weak < 2) - int(colour == 0 and weak > 0)
    return b3, ell


def pair_gates():
    """Derive CZ on two marked ports; spectator rails factor at the end of the word."""
    _, _, h = material.dimer_gates()
    law = candidate.native_law
    ex, con = law.local_term_parts()
    full = ex + candidate.KAPPA * con
    physical = np.eye(law.DIM)[
        :, [law.SRC.index((a, b)) for a in (1, 0) for b in (1, 0)]
    ]
    contact = physical.T @ full @ physical
    np.testing.assert_allclose(full @ physical, physical @ contact, atol=2e-13)
    reference = candidate.KAPPA * np.diag([1, 0])
    refs = [np.kron(reference, np.eye(2)), np.kron(np.eye(2), reference)]
    unitary = np.eye(4)
    for kind, a, _, duration in record_transfer.controlled_phase(0, 1):
        generator = contact if kind == "contact" else refs[a]
        unitary = expm(-1j * pi * float(duration) * generator) @ unitary
    unitary /= unitary[0, 0]
    np.testing.assert_allclose(unitary, np.diag([1, 1, 1, -1]), atol=2e-13)
    hadamards = [coherent.at(h, j) for j in range(5)]
    cz, cx = {}, {}
    for i, j in combinations(range(5), 2):
        bits = 2 * ((WORDS >> i) & 1) + ((WORDS >> j) & 1)
        gate = np.diag(unitary.diagonal()[bits])
        cz[i, j] = cz[j, i] = gate
        for a, b in ((i, j), (j, i)):
            cx[a, b] = hadamards[b].conj().T @ gate @ hadamards[b]
            expected = np.eye(32)[:, WORDS ^ (((WORDS >> a) & 1) << b)]
            np.testing.assert_allclose(cx[a, b], expected, atol=2e-13)
    return hadamards, cz, cx


def connectors(geometry):
    board = federation(0)
    graph = nx.Graph(board.seams)
    records = []
    for side, owner in geometry["owners"].items():
        local = graph.subgraph(range(12 * owner, 12 * (owner + 1)))
        ports = min(
            (min(a, b), max(a, b))
            for a, b in local.edges()
            if board.degree[a] == board.degree[b] == 5
        )
        for i, j in combinations(range(5), 2):
            tags = list(range(board.n))
            paths = []
            marked = [geometry["dimers"][side][m][1] for m in (i, j)]
            for label, port in zip(marked, ports, strict=True):
                path = material.shared_report.path(local, tags.index(label), port)
                material.shared_report.transpose(tags, path)
                paths.append(path)
            assert [tags[p] for p in ports] == marked
            for path in reversed(paths):
                material.shared_report.transpose(tags, path)
            assert tags == list(range(board.n))
            records.append([side, i, j, list(ports), paths])
    return records


def relative_phase(j, turns):
    # A filled marked port has reference energy 8 and a vacancy has zero.
    dwell = 2 * pi * float(turns % 1) / candidate.KAPPA
    local = np.diag(np.exp(-1j * dwell * np.array([candidate.KAPPA, 0])))
    return coherent.at(local / local[0, 0], j)


def charge_program(values, cx):
    operator = np.eye(32)
    pulses = []
    for mask in range(1, 32):
        sign_sum = sum(
            v * (-1) ** ((w & mask).bit_count()) for w, v in enumerate(values)
        )
        coefficient = Fraction(-sign_sum, 16)
        if coefficient == 0:
            continue
        parity = [j for j in range(5) if mask & (1 << j)]
        target = parity[-1]
        compute = np.eye(32)
        for j in parity[:-1]:
            compute = cx[j, target] @ compute
        operator = (
            compute.conj().T
            @ relative_phase(target, coefficient / 4)
            @ compute
            @ operator
        )
        pulses.append([mask, str(coefficient)])
    np.testing.assert_allclose(operator, np.diag((1j) ** np.array(values)), atol=4e-12)
    return operator, pulses


def phase_banks(raw, gates, delivery):
    h, _, cx = gates
    rng = np.random.Generator(np.random.PCG64(SEED))
    b3, ell = zip(*(charges(w) for w in range(32)), strict=True)
    programs = {
        tag: charge_program(values, cx)
        for tag, values in (("3B", b3), ("L", ell), ("opposite_L", [-v for v in ell]))
    }
    records = []
    for word in range(1, 32):
        modes = [j for j in range(5) if word & (1 << j)]
        pivot = modes[0]
        fanout = np.eye(32)
        for j in modes[1:]:
            fanout = cx[pivot, j] @ fanout
        prepared = fanout @ h[pivot][:, 0]
        for tag, (operator, _) in programs.items():
            state = operator @ prepared
            for side, path in (("sender", np.eye(32)), ("receiver", delivery)):
                arms = fanout.conj().T @ path @ state
                for axis in ("X", "Y"):
                    measured = (
                        arms
                        if axis == "X"
                        else relative_phase(pivot, Fraction(-1, 4)) @ arms
                    )
                    probabilities = abs(h[pivot].conj().T @ measured) ** 2
                    np.testing.assert_allclose(probabilities.sum(), 1, atol=5e-12)
                    counts = sequential.grid_counts(
                        rng, SHOTS, probabilities / probabilities.sum()
                    )
                    raw_counts = [
                        [raw[k][side], int(counts[k])] for k in np.flatnonzero(counts)
                    ]
                    records.append([word, tag, side, axis, raw_counts])
    return records, {tag: p[1] for tag, p in programs.items()}


def phase_read(records):
    radius = sqrt(log(2 * len(records) / 0.01) / (2 * SHOTS)) + 3e-12
    by_setting = {}
    for word, tag, side, axis, raw_counts in records:
        assert sum(n for _, n in raw_counts) == SHOTS
        bright = sum(n for raw, n in raw_counts if material.decode(raw) == 0)
        by_setting.setdefault((word, tag, side), {})[axis] = bright / SHOTS
    values = {}
    for key, seen in by_setting.items():
        passing = [
            q
            for q in (-1, 0, 1)
            if all(
                abs(seen[axis] - p) <= radius
                for axis, p in (
                    ("X", (1 + (1j**q).real) / 2),
                    ("Y", (1 + (1j**q).imag) / 2),
                )
            )
        ]
        assert len(passing) == 1, (key, seen)
        values[key] = passing[0]
    return values, radius


def sorting_trials(raw, gates, delivery):
    h, cz, cx = gates
    trials = []
    for order in permutations(range(4)):
        for signed in (True, False):
            state = h[4][:, 0].copy()
            for j in order:
                # Both interferometer arms pass the local absence query.
                np.testing.assert_allclose(
                    state[(WORDS & (1 << j)) != 0], 0, atol=2e-13
                )
                if signed:
                    for preceding in range(j):
                        state = cz[4, preceding] @ state
                state = cx[4, j] @ state
            for j in range(4):
                state = cx[4, j] @ state
            out = coherent.endpoint_for_unit_branch(h[4].conj().T @ delivery @ state)
            assert out in (0, 16)
            trials.append([list(order), signed, raw[out]["receiver"]])
    duplicates = []
    for j in range(4):
        word = coherent.endpoint_for_unit_branch(cx[4, j][:, 16])
        occupied = raw[word]["sender"][4 * j : 4 * j + 4]
        assert occupied == "0001"
        duplicates.append([j, occupied])
    return trials, duplicates
