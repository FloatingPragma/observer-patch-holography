"""Exact contact support, reachable source coherences and a one-arm vacancy read."""

import argparse
import json
from functools import cache
from itertools import product
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import eigh, expm
from mpmath.ctx_iv import MPIntervalContext

from model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure import (
    LAW,
    blocks,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.material_interferometer.apparatus import (
    PHASE_PULSES,
    Terminal,
)
from model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder.release.geometry import (
    layout,
)
from ...eraser import TIMES, gram_by_owner, launch


@cache
def wiring():
    lab = Terminal("tip")
    marks = layout()[1]
    arms = [next((a, c, g) for a, c, g in lab.arms if (a, g) == mark) for mark in marks]
    ports = (lab.hub, *(p for arm in arms for p in arm))
    assert len(set(ports)) == len(ports) == 7
    return ports, arms


@cache
def contact(a, b, bit, strength):
    """Integer 8H in the full seven-position/two-direction-word input chart."""
    ports, _ = wiring()
    mat = blocks(strength)["blank" if bit is None else (1, 1)]
    term = np.zeros((28, 28), dtype=int)
    left, right = sorted((a, b))
    for pos, head in product(ports, range(4)):
        col = 4 * ports.index(pos) + head
        word = (int(pos != left), int(pos != right))
        coin = 0 if bit is None else ((head >> bit) & 1) ^ int(a > b)
        source = 9 * coin + LAW.SRC.index(word)
        background = 9 * coin + LAW.SRC.index((1, 1))
        term[col, col] = mat[source, source] - mat[background, background]
        for dest in np.flatnonzero(mat[:, source]):
            if dest == source:
                continue
            new_word = LAW.SRC[dest % 9]
            assert 0 in word and new_word == word[::-1]
            assert dest // 9 == (coin ^ int(bit is not None))
            new_pos = (left, right)[new_word.index(0)]
            new_head = head if bit is None else head ^ (1 << bit)
            term[4 * ports.index(new_pos) + new_head, col] = mat[dest, source]
    np.testing.assert_array_equal(term, term.T)
    return term


def partitions(open_detours, strength):
    ports, arms = wiring()
    graph = nx.Graph()
    graph.add_nodes_from(range(28))
    for bit, (a, c, g) in enumerate(arms):
        edges = [(ports[0], a, None), (a, c, None), (a, g, bit)]
        if open_detours & (1 << bit):
            edges.append((c, g, None))
        for u, v, mark in edges:
            graph.add_edges_from(
                zip(*np.nonzero(np.triu(contact(u, v, mark, strength), 1)), strict=True)
            )
    cohorts = sorted(
        sorted(nodes & set(range(4)))
        for nodes in nx.connected_components(graph)
        if nodes & set(range(4))
    )
    expected = [
        [v for v in range(4) if (v & ~open_detours) == r]
        for r in sorted({v & ~open_detours for v in range(4)})
    ]
    assert cohorts == sorted(expected)
    shortest = dict(nx.single_source_shortest_path_length(graph, 1))
    return dict(
        hub_components=cohorts,
        distinct_one_step_exchanges=graph.number_of_edges(),
        shortest_from_01_to_10=shortest.get(2),
    )


@cache
def sharp_read(bit):
    """Borrow the vacancy, read a sharp mark at a physical port, then return the probe."""
    ports, arms = wiring()
    a, _, g = arms[bit]
    u = expm((-1j * (np.pi / 8) / 8) * contact(a, g, bit, 8)) @ expm(
        (-1j * np.pi / 8) * contact(ports[0], a, None, 8)
    )
    start = np.eye(28, 4, dtype=complex)
    operators = []
    for value, pos in enumerate((g, a)):
        ended = u @ start
        ended[np.repeat(np.arange(7) != ports.index(pos), 4)] = 0
        restored = u.conj().T @ ended
        ideal = np.diag(((np.arange(4) >> bit) & 1) == value)
        np.testing.assert_allclose(restored, start @ ideal, atol=3e-13)
        operators.append(ideal)
    return operators


def phase_gate_error():
    """Outward-rounded action of the blank three-port pulse on the occupied corner."""
    ctx = MPIntervalContext()
    ctx.dps = 45
    ports, arms = wiring()
    first, hub, spare = arms[0][0], ports[0], arms[1][0]
    for u, v in ((first, hub), (hub, spare)):
        indexes = [4 * ports.index(u), 4 * ports.index(v)]
        term = contact(u, v, None, 8)
        np.testing.assert_array_equal(
            term[np.ix_(indexes, indexes)], [[68, -4], [-4, 68]]
        )
    vector = [ctx.mpc(1), ctx.mpc(0), ctx.mpc(0)]
    unit_i = ctx.mpc(0, 1)
    for j, duration in enumerate(PHASE_PULSES):
        num, den = duration.as_integer_ratio()
        half = ctx.mpf(num) / (2 * den)
        i = j % 2
        a, b = vector[i : i + 2]
        co, si = ctx.cos(half), unit_i * ctx.sin(half)
        common = ctx.exp(-17 * unit_i * half)
        vector[i : i + 2] = [common * (co * a + si * b), common * (si * a + co * b)]
    vector[0] -= ctx.exp(2 * unit_i * ctx.pi / 3)
    error = sum(abs(x) for x in vector)
    assert error < ctx.mpf(1) / 10**12
    return 1e-12


@cache
def collection(arm, strength, panel):
    ports, arms = wiring()
    a, c, g = arms[arm]
    gates = [(ports[0], a, None, np.pi), (a, g, arm, np.pi / 8)]
    if panel != "cut":
        gates.append((g, c, None, np.pi))
    if panel == "rotated":
        spare = arms[1 - arm][0]
        phase_edges = [(a, ports[0]), (ports[0], spare)]
        seq = [(*phase_edges[i % 2], None, t) for i, t in enumerate(PHASE_PULSES)]
        rot = np.eye(28, dtype=complex)
        for u, v, mark, time in seq:
            rot = expm((-1j * time / 8) * contact(u, v, mark, strength)) @ rot
        probe_phase = np.exp(2j * np.pi / 3)
        for place, phase in ((a, probe_phase), (c, 1)):
            for head in range(4):
                index = 4 * ports.index(place) + head
                np.testing.assert_allclose(
                    rot[:, index], phase * np.eye(28)[:, index], atol=1e-13
                )
        gates.extend(seq)
    if panel != "park":
        gates.append((c, a, None, (3 if panel == "inverse" else 1) * np.pi / 2))
    amplitudes = np.eye(28, 4, dtype=complex)
    for u, v, bit, time in gates:
        amplitudes = expm((-1j * time / 8) * contact(u, v, bit, strength)) @ amplitudes
    np.testing.assert_allclose(amplitudes.conj().T @ amplitudes, np.eye(4), atol=3e-13)
    return amplitudes


def record(rho, arm, strength, panel):
    a = collection(arm, strength, panel)
    filtered = rho
    if panel == "dephased":
        sector = (np.arange(4) >> arm) & 1
        filtered = rho * (sector[:, None] == sector)
    freq = np.einsum("is,st,it->i", a, filtered, a.conj()).real.reshape(7, 4)
    spectator = (np.arange(4) >> (1 - arm)) & 1
    out = np.stack([freq[:, spectator == b].sum(axis=1) for b in (0, 1)], axis=-1)
    assert out.min() >= -1e-13 and abs(out.sum() - 1) < 2e-12
    out = np.maximum(out, 0)
    return out / out.sum()


def tv(left, right):
    return float(np.abs(left - right).sum() / 2)


def bounds(rhos, setting):
    mask = (np.arange(4)[:, None] ^ np.arange(4)) & ~setting
    diff = (rhos[0] - rhos[1]) * (mask == 0)
    np.testing.assert_allclose(diff, diff.conj().T, atol=1e-12)
    eig, vec = eigh(diff)
    positive = vec[:, eig > 0]
    effect = positive @ positive.conj().T
    upper = float(abs(eig).sum() / 2)
    # A declared sign-of-difference measurement attains the binary decision bound.
    np.testing.assert_allclose(np.trace(effect @ diff).real, upper, atol=1e-12)
    return upper


def calculations():
    ports, arms = wiring()
    for bit in (0, 1):
        sharp_read(bit)
    structure = {k: {s: partitions(s, k) for s in range(4)} for k in (4, 8)}
    source = launch()
    rows = []
    for owner, gram in gram_by_owner().items():
        for time, overlap in zip(TIMES, gram, strict=True):
            received = source[1:] * overlap
            mass = np.trace(received, axis1=1, axis2=2).real
            np.testing.assert_allclose(mass[0], mass[1], rtol=1e-10, atol=0)
            rho = received / mass[:, None, None]
            ceiling = [bounds(rho, s) for s in range(4)]
            records = []
            for arm, strength, panel in product(
                range(2),
                (8, 4),
                ("open", "inverse", "park", "cut", "dephased", "rotated"),
            ):
                if panel == "rotated" and strength != 8:
                    continue
                r, t = (record(d, arm, strength, panel) for d in rho)
                joint_tv = tv(r, t)
                assert joint_tv <= ceiling[1 << arm] + 1e-11
                records.append(
                    dict(
                        arm=arm,
                        analyzer_strength=strength,
                        panel=panel,
                        plus=r.tolist(),
                        minus=t.tolist(),
                        tv_joint=joint_tv,
                        tv_vacancy=tv(r.sum(axis=1), t.sum(axis=1)),
                    )
                )
            rows.append(
                dict(
                    owner=owner,
                    time=time,
                    original_arrival=mass.mean(),
                    conditional_tv_ceiling=ceiling,
                    original_guess_advantage_ceiling=(
                        mass.mean() * np.array(ceiling) / 2
                    ).tolist(),
                    conditional_source_coherence_difference=[
                        [i, j, z.real, z.imag]
                        for i in range(4)
                        for j in range(i + 1, 4)
                        for z in [rho[0, i, j] - rho[1, i, j]]
                    ],
                    records=records,
                )
            )
    standards = []
    for arm in range(2):
        for label, coef in (
            ("zero", (1, 0)),
            ("one", (0, 1)),
            ("+i", (1, 1j)),
            ("-i", (1, -1j)),
        ):
            psi = np.zeros(4, complex)
            psi[[0, 1 << arm]] = coef
            psi /= np.linalg.norm(psi)
            h = np.outer(psi, psi.conj())
            standards.append(
                dict(
                    arm=arm,
                    label=label,
                    native=record(h, arm, 8, "open").sum(axis=1).tolist(),
                    rotated=record(h, arm, 8, "rotated").sum(axis=1).tolist(),
                )
            )
    return dict(
        ports=ports,
        arms=arms,
        source_support=structure,
        three_port_phase_vector_error_upper=phase_gate_error(),
        standards=standards,
        panels=rows,
        conditional_probabilities="port (listed order), then sharp spectator 0/1",
        prior="two opposite charged windings with equal prior, all misses pooled",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    out = calculations()
    if args.write:
        args.write.write_text(json.dumps(out, separators=(",", ":")) + "\n")
    for row in out["panels"]:
        print(
            row["owner"],
            row["time"],
            row["conditional_tv_ceiling"],
            [
                (r["tv_vacancy"], r["tv_joint"])
                for r in row["records"]
                if r["panel"] == "open" and r["analyzer_strength"] == 8
            ],
        )
