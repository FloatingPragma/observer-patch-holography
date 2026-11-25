"""Matched-host formation spectroscopy on mixed blank and written physical rings."""

import json
from fractions import Fraction
from math import floor, prod
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import eigh_tridiagonal

from model.candidate_laws.interval_echo.openwave.row_11_spectrum.marked_bands import (
    material,
)
from . import measure


HERE = Path(__file__).resolve().parent
ELECTRON_MEV = Fraction("0.51099895069")
TAU_RATIO = Fraction("1776.93") / ELECTRON_MEV
LIMIT = 198 / TAU_RATIO
MU, TAU = 206.7682827, float(TAU_RATIO)  # CODATA and PDG, used for comparisons.


def contacts():
    tables = measure.blocks()
    out = {}
    for mark, name in [(None, "blank"), (1, "even"), (-1, "odd")]:
        term = tables["blank" if mark is None else (1, mark)]
        for charge in (-1, 1):
            # Conjugation swaps the marked-direction coin along with the signed load.
            c = int(charge > 0 and mark is not None)
            pair = (1 + charge, 1)
            i = 9 * c + measure.LAW.SRC.index(pair)
            j = 9 * (1 - c if mark is not None else 0) + measure.LAW.SRC.index(
                pair[::-1]
            )
            f = 9 * c + measure.LAW.SRC.index((1, 1))
            col = term[:, i].copy()
            assert set(np.flatnonzero(col)) <= {i, j}
            assert set(np.flatnonzero(term[:, f])) <= {f}
            if mark is not None:
                blocked = 9 * (1 - c) + measure.LAW.SRC.index(pair)
                assert set(np.flatnonzero(term[:, blocked])) == {blocked}
                assert term[blocked, blocked] == term[f, f]
            block = (int(col[i] - term[f, f]), int(col[j]))
            if name in out:
                assert out[name] == block
            out[name] = block
    assert out == {"blank": (68, -4), "even": (-32, -32), "odd": (-36, -36)}
    return out


def lifted(pattern):
    """All reachable positions from source zero, with unit direction zero at each mark."""
    n = len(pattern)
    if all(x == "0" for x in pattern):
        raise ValueError("All blank has a cycle, with no winding memory")
    lo, hi = 0, 0
    while pattern[(lo - 1) % n] == "0":
        lo -= 1
    while pattern[hi % n] == "0":
        hi += 1
    return "".join(pattern[t % n] for t in range(lo, n + hi)), lo


def coefficients(word, parity, cyclic=False):
    value = -32 if parity == "even" else -36
    site = [-36 if mark == "2" else value if mark == "1" else 68 for mark in word]
    hop = [-36 if mark == "2" else value if mark == "1" else -4 for mark in word]
    diag = [a + b for a, b in zip([0, *site], [*site, 0], strict=True)]
    if cyclic:
        diag[0] += diag.pop()
    return diag, hop


def inertia(diag, hop, energy, cyclic=False):
    """Sturm sign variations, including an exact bordered final cyclic determinant."""
    den, num = energy.denominator, energy.numerator
    d = [den * a - 8 * num for a in diag]
    e = [den * a for a in hop]

    def minors(ds, es):
        vals = [1, ds[0]]
        for a, b in zip(ds[1:], es, strict=True):
            vals.append(a * vals[-1] - b * b * vals[-2])
        return vals

    m = minors(d, e[: len(d) - 1])
    if cyclic:
        m[-1] += -(e[-1] ** 2) * minors(d[1:-1], e[1:-2])[-1] + 2 * (-1) ** (
            len(d) - 1
        ) * prod(e)
    signs = [x > 0 for x in m if x]
    return sum(a != b for a, b in zip(signs, signs[1:])), int(m[-1] == 0)


def bracket(diag, hop, minimum, cyclic=False, index=0):
    scale = 10**9
    below = floor(minimum * scale)
    for lower in (below, below - 1):
        lo, hi = Fraction(lower, scale), Fraction(lower + 2, scale)
        if (
            inertia(diag, hop, lo, cyclic) == (index, 0)
            and inertia(diag, hop, hi, cyclic)[0] > index
        ):
            return [float(lo), float(hi)]
    raise ArithmeticError("Root bracket failed")


def sweep():
    n = 15
    classes = set()
    for mask in range(1, 2**n):
        path, _ = lifted(format(mask, "015b"))
        classes.add(min(path, path[::-1]))
    result = {}
    for parity in ("even", "odd"):
        hits, best = 0, (float("inf"), None)
        for path in sorted(classes):
            d, e = coefficients(path, parity)
            neg, zero = inertia(d, e, Fraction(0))
            end, on_cut = inertia(d, e, LIMIT)
            assert not on_cut
            if end > neg + zero:
                hits += 1
                k = neg + zero
                level = float(
                    eigh_tridiagonal(
                        np.array(d) / 8,
                        np.array(e) / 8,
                        select="i",
                        select_range=(k, k),
                        eigvals_only=True,
                    )[0]
                )
                if level < best[0]:
                    best = level, (path, k)
        result[parity] = dict(
            canonical_open_words=len(classes),
            with_line_below_ceiling_over_tau_ratio=hits,
        )
        if best[1]:
            path, k = best[1]
            d, e = coefficients(path, parity)
            result[parity].update(
                minimum_positive=bracket(d, e, best[0], index=k),
                attaining_lifted_word=path,
            )
    return result


def mixed_parities():
    # All 128 independently assigned parity words on the selected seven marks.
    hits = 0
    for i in range(128):
        parities = format(i, "07b")
        motif = "".join("0" + ("1" if p == "0" else "2") for p in parities) + "0"
        path, _ = lifted(motif)
        d, e = coefficients(path, "even")
        lower, zero = inertia(d, e, Fraction(0))
        upper, bound = inertia(d, e, LIMIT)
        assert not bound
        hits += int(upper > lower + zero)
    return dict(assignments=128, with_line_below_cutoff=hits)


def periodic(motif, periods):
    d, e = coefficients(motif, "even", cyclic=True)
    h = np.diag(d) + np.diag(e[:-1], 1) + np.diag(e[:-1], -1)
    h[-1, 0] = h[0, -1] = e[-1]
    zero = np.linalg.eigvalsh(h / 8)[0]
    assert inertia(d, e, Fraction(0), True) == (0, 0)
    band = bracket(d, e, zero, cyclic=True)
    h[-1, 0] *= -1
    h[0, -1] *= -1
    zone_edge = float(np.linalg.eigvalsh(h / 8)[0])
    word, offset = lifted(motif * periods)
    d, e = coefficients(word, "even")
    levels, v = eigh_tridiagonal(np.array(d) / 8, np.array(e) / 8)
    line = float(levels[0])
    certified = bracket(d, e, line)
    n = len(motif) * periods
    ports = (np.arange(len(d)) + offset) % n
    # Select a native undispersed sharp source by routing the charge to the
    # middle of one well before release, and select the same port one period away.
    well = len(motif) // 2
    j = int(np.nonzero(np.arange(len(d)) + offset == well)[0][0])
    peak = v[j, 0] ** 2
    detector = (well + len(motif)) % n
    remote = np.flatnonzero(ports == detector)
    if periods == 1:
        detector = (well + n // 2) % n
        remote = np.flatnonzero(ports == detector)
    # Energy pole residue in the propagator, plus the unfiltered arrival of a
    # sharp material vacancy. Time set by the independently computed bandwidth.
    dwell = 1 / (zone_edge - zero)
    amp = v[remote] @ (np.exp(-1j * (levels - line) * dwell) * v[j])
    matched = {}
    for name, goal in [("middle", MU), ("heavy", TAU)]:
        i = int(np.argmin(abs(levels / line - goal)))
        matched[name] = dict(
            ratio=float(levels[i] / line),
            index=i,
            relative_error=float(levels[i] / (line * goal) - 1),
        )
    return dict(
        period=len(motif),
        repeats=periods,
        phase_zero_edge=band,
        zone_edge=zone_edge,
        bandwidth=zone_edge - zero,
        finite_line=certified,
        source_offset=well,
        detector_offset=detector,
        source_line_weight=float(peak),
        source_to_detector_pole_norm=float(np.linalg.norm(v[remote, 0]) * abs(v[j, 0])),
        free_time=dwell,
        free_arrival=float(np.sum(abs(amp) ** 2)),
        free_record=measure.counts(
            float(np.sum(abs(amp) ** 2)), np.random.default_rng(2689 + len(motif))
        ),
        nearest_anchored_levels=matched,
    )


def stock_routes(ports, motifs):
    board = material.federation(2)
    graph = nx.Graph(board.seams)
    ring = set(ports)
    distance = nx.multi_source_dijkstra_path_length(graph, ring)
    donor = max(graph, key=lambda p: (distance[p], -p))
    assert distance[donor] > 2
    for motif in motifs:
        marked = [
            (a, b)
            for i, (a, b) in enumerate(zip(ports, ports[1:] + ports[:1], strict=True))
            if motif[i % len(motif)] == "1"
        ]
        usable = graph.copy()
        usable.remove_node(donor)
        usable.remove_edges_from(marked)
        assert nx.is_connected(usable)
        # No removal of a single held endpoint can cut the borrowed phase path.
        assert not set(nx.articulation_points(usable))
    reference = next(
        v
        for v in sorted(graph)
        if v not in ring and v != donor and board.degree[v] == 5
    )
    return dict(
        ports=len(ports),
        owners=len({v // 12 for v in ports}),
        private_reference=reference,
        erased_private=sum(board.degree[p] == 5 for p in ports),
        compensator=donor,
        all_motif_blank_bypasses_biconnected=True,
    )


def run():
    terms = contacts()
    ring, _ = measure.apparatus()
    long, _ = material.geometry()
    assert len(ring) == 15 and len(long) == 144
    # The cyclic control with strict even alternation has a constant zero mode.
    control = {}
    for n in (15, 144):
        pat = ("01" * n)[:n]
        word, _ = lifted(pat)
        d, e = coefficients(word, "even")
        z = float(
            eigh_tridiagonal(np.array(d) / 8, np.array(e) / 8, eigvals_only=True)[0]
        )
        control[str(n)] = bracket(d, e, z)
    rows = [periodic("01" * 7 + "0", 1)]
    motifs = []
    for p in [12, 16, 18, 24, 36, 48, 72]:
        motif = "01" * ((p - 2) // 2) + "00"
        motifs.append(motif)
        rows.append(periodic(motif, 144 // p))
    return dict(
        native_eighth_units=terms,
        alternating_prefix_finite=control,
        laboratory_ratios=dict(middle=MU, heavy=TAU),
        mixed_parities_selected_motif=mixed_parities(),
        physical_route_lengths=[len(ring), len(long)],
        long_preparation=stock_routes(long, motifs),
        sweep_cutoff=str(LIMIT),
        homogeneous_parity_sweep=sweep(),
        periodic_rows=rows,
    )


if __name__ == "__main__":
    import sys

    output = json.dumps(run(), indent=2, sort_keys=True, default=int) + "\n"
    if "--write" in sys.argv:
        (HERE / "mixed_masses.json").write_text(output)
    else:
        print(output, end="")
