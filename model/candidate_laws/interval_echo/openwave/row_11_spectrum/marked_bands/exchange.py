"""Exchange of identical, interacting vacancies on a calibrated pair of native rails."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import eigh, expm

from model.candidate_laws.interval_echo.openwave.row_11_spectrum.conjugate_flight.measure import (
    LAW,
    edge,
    seam_data,
    binomial_interval,
)
from model.experiments.interval_echo_gauge.flag_measure.measure import meter
from . import paired
from .material import programmed_mark, synthesize

HERE = Path(__file__).resolve().parent
N = paired.N
SIZE = 2 * N
ROUTE = (*paired.rail.RAILS[0], *reversed(paired.rail.RAILS[1]))
PAIRS = tuple(combinations(range(SIZE), 2))
INDEX = {e: j for j, e in enumerate(PAIRS)}
FIRST, SECOND = np.array(PAIRS).T
SHOTS = 2**23
DELAY = 2**27


def hardware():
    board = paired.federation(0)
    graph = nx.Graph(board.seams)
    assert all(b in board.neighbours[a] for a, b in zip(ROUTE, ROUTE[1:] + ROUTE[:1]))
    endpoints = [edge(ROUTE[-1], ROUTE[0]), edge(ROUTE[N - 1], ROUTE[N])]
    assert all(e in paired.rail.preparation()["cuts"] for e in endpoints)
    blank = set(ROUTE)
    stores = []
    for _ in range(2):
        dist = nx.multi_source_dijkstra_path_length(graph, blank)
        a, b = next((a, b) for a, b in board.seams if min(dist[a], dist[b]) > 3)
        stores.append([a, b])
        cuts = [e for e in board.seams if (a in e or b in e) and e != (a, b)]
        assert all((b if a in e else a) not in e for e in cuts)
        blank.update((a, b))
    return dict(
        route=ROUTE, reopened_filled_shutters=endpoints, stock_excess_dimers=stores
    )


def matrices(closed):
    one = paired.guide(-1)
    if closed:
        pair_block = paired.rail.pair(8, -1)
        for a, b in [(0, N), (N - 1, SIZE - 1)]:
            ids = [a, b]
            one[np.ix_(ids, ids)] += pair_block
    order = np.r_[np.arange(N), np.arange(SIZE - 1, N - 1, -1)]
    one = one[np.ix_(order, order)]
    adjacent = lambda a, b: a != b and bool(one[a, b])
    # The scalar edge correction compares two true holes with one hole in each independent copy.
    local = seam_data(8)[0]
    cost = lambda a, b: local[LAW.SRC.index((a, b)), LAW.SRC.index((a, b))]
    correction = cost(0, 0) - cost(0, 1) - cost(1, 0) + cost(1, 1)
    assert abs(correction + 17) < 1e-13 and abs(programmed_mark(8) - 32) < 1e-13
    two = np.zeros((len(PAIRS), len(PAIRS)))
    fermion = np.zeros_like(two)
    for j, (a, b) in enumerate(PAIRS):
        two[j, j] = fermion[j, j] = one[a, a] + one[b, b] + correction * adjacent(a, b)
        for old, spectator in [(a, b), (b, a)]:
            for dest in np.flatnonzero(one[:, old]):
                if dest == old or dest == spectator:
                    continue
                i = INDEX[tuple(sorted((dest, spectator)))]
                two[i, j] = one[dest, old]
                # Exterior product in increasing order. Crossing a spectator reverses orientation.
                crossed = min(old, dest) < spectator < max(old, dest)
                fermion[i, j] = (-1 if crossed else 1) * one[dest, old]
    for h in (two, fermion):
        np.testing.assert_array_equal(h, h.T)
    return one, two, fermion


def source(band, side="left", rail=0):
    h = paired.guide(-1)[:N, :N]
    _, orbitals = eigh(h[np.ix_(paired.LEFT, paired.LEFT)])
    v = synthesize(orbitals[:, band].astype(complex))
    sites = paired.LEFT if side == "left" else paired.RIGHT
    if side == "right":
        v = v[::-1]
    s = np.zeros(SIZE, complex)
    s[(SIZE - 1 - sites) if rail else sites] = v
    return s


def propagate(spectral, s, t):
    e, v = spectral
    # Cancel a common phase before multiplying out long free waits.
    return v @ (np.exp(-1j * t * (e - e[0])) * (v.conj().T @ s))


def positions(v):
    return np.abs(v) ** 2 / np.vdot(v, v).real


def bank(rng, p):
    p = np.maximum(0.0, p)
    return rng.multinomial(SHOTS, p / sum(p))


def coarse(counts):
    both = SECOND < N
    split = (FIRST < N) & (SECOND >= N)
    left = np.isin(FIRST, paired.LEFT)
    right = np.isin(SECOND, paired.RIGHT)
    return dict(
        both_rail_a=int(sum(counts[both])),
        one_each_rail=int(sum(counts[split])),
        both_rail_b=int(sum(counts[FIRST >= N])),
        one_each_original_well=int(sum(counts[left & right])),
        both_original_left_well=int(sum(counts[left & np.isin(SECOND, paired.LEFT)])),
        both_original_right_well=int(sum(counts[right & np.isin(FIRST, paired.RIGHT)])),
    )


def prep_two(left, right):
    # Independently synthesized on disjoint physical wells, hence no order/sign convention in the source.
    assert not np.any((abs(left) > 1e-12) & (abs(right) > 1e-12))
    plus = left[FIRST] * right[SECOND] + left[SECOND] * right[FIRST]
    minus = left[FIRST] * right[SECOND] - left[SECOND] * right[FIRST]
    # If the order of the supports is reversed an unobservable global sign can result.
    assert abs(abs(np.vdot(plus, minus)) - 1) < 1e-12
    return plus, minus


def split_controls(one_e, left, right, t):
    a = propagate(one_e, left, t)
    b = propagate(one_e, right, t)
    det = a[FIRST] * b[SECOND] - a[SECOND] * b[FIRST]
    # Soft boson foil permits two quanta at one address; report these separately.
    sym = a[FIRST] * b[SECOND] + a[SECOND] * b[FIRST]
    diagonal = 2 * abs(a * b) ** 2
    ordinary = float(
        np.sum(
            abs(det[np.isin(FIRST, paired.LEFT) & np.isin(SECOND, paired.RIGHT)]) ** 2
        )
    )
    matched = (
        abs(np.vdot(left, a) * np.vdot(right, b) - np.vdot(left, b) * np.vdot(right, a))
        ** 2
    )
    soft = (
        abs(np.vdot(left, a) * np.vdot(right, b) + np.vdot(left, b) * np.vdot(right, a))
        ** 2
    )
    assert abs(sum(abs(sym) ** 2) + sum(diagonal) - 1) < 2e-11
    return dict(
        independent_fermion_one_each_well=ordinary,
        free_fermion_one_each_matched_mode=float(matched),
        soft_boson_one_each_matched_mode=float(soft),
        soft_same_site=float(sum(diagonal)),
    )


def phase_witness(rng):
    # Independent, one-address bare contact calibrations for every loop edge.
    U = expm(-1j * np.pi * seam_data(8)[0])
    matrix_element = lambda a, b: U[LAW.SRC.index((b, a)), LAW.SRC.index((a, b))]

    def phase_run(n):
        starts = (2,) if n == 1 else (2, 8)
        words = np.ones(SIZE, dtype=int)
        words[list(starts)] = 0
        labels = np.arange(n)
        tracks = dict(zip(starts, labels.tolist()))
        if n == 1:
            currents = [
                (a, (a + 1) % SIZE) for a in [(2 + i) % SIZE for i in range(SIZE)]
            ]
        else:
            currents = [(a, (a + 1) % SIZE) for a in [*range(8, SIZE), 0]]
            currents += [(a, a + 1) for a in range(2, 8)] + [(1, 2)]
        amplitude = 1.0 + 0j
        for x, y in currents:
            assert words[x] == 0 and words[y] == 1
            amplitude *= matrix_element(words[x], words[y])
            words[x], words[y] = words[y], words[x]
            tracks[y] = tracks.pop(x)
        assert set(tracks) == set(starts) and len(currents) == SIZE
        assert [tracks[x] for x in starts] == list(reversed(labels.tolist()))
        counts = [
            int(rng.binomial(SHOTS, np.clip((1 + v) / 2, 0, 1)))
            for v in (amplitude.real, amplitude.imag)
        ]
        return dict(
            start=starts,
            hops=currents,
            final_traced_port_labels=[tracks[s] for s in starts],
            xy_bright=counts,
            predicted_complex=[amplitude.real, amplitude.imag],
        )

    one = phase_run(1)
    two = phase_run(2)
    _, modes = eigh(paired.guide(-1)[np.ix_(paired.LEFT, paired.LEFT)])
    labels = []
    for b in range(3):
        waves = []
        for sites, cols, src in [
            (paired.LEFT, modes, source(b)),
            (paired.RIGHT, modes[::-1, :], source(b, "right")),
        ]:
            wave = np.zeros(SIZE, complex)
            wave[sites] = paired.orbital_analyzer(cols, -1) @ src[sites]
            for j in range(b, 0, -1):
                a, c = sites[j - 1 : j + 1]
                wave[[a, c]] = wave[[c, a]]
            assert abs(wave[sites[0]]) ** 2 > 1 - 2e-12
            waves.append(wave)
        input_state, _ = prep_two(*waves)
        # The π contact is the ordinary positive load permutation on all nine pairs.
        # Apply its full word even to the residual amplitudes of the source compiler.
        final = np.arange(SIZE)
        for x, y in two["hops"]:
            final[[x, y]] = final[[y, x]]
        addresses = np.array(
            [INDEX[tuple(sorted((final[i], final[j])))] for i, j in PAIRS]
        )
        amplitude = np.vdot(input_state, input_state[addresses])
        measured = [
            int(rng.binomial(SHOTS, np.clip((1 + t) / 2, 0, 1)))
            for t in (amplitude.real, amplitude.imag, -amplitude.real)
        ]
        labels.append(dict(band=b, real_imag_and_pi_shifted_bright=measured))
    return dict(
        one_around_loop=one, two_around_each_other=two, same_orbital_records=labels
    )


def run():
    hw = hardware()
    assert meter() == [0, 1, 2]
    rng = np.random.default_rng(401513)
    opened = matrices(False)
    ring = matrices(True)
    tree = [eigh(h) for h in opened]
    cycle = [eigh(h) for h in ring]
    h_pi = ring[0].copy()
    h_pi[0, -1] *= -1
    h_pi[-1, 0] *= -1
    flux_pi = eigh(h_pi)
    samples = []
    pilot_times = np.arange(257) * 32768
    for band in range(3):
        left = source(band)
        right = source(band, "right")
        msrc, fsrc = prep_two(left, right)
        # Frozen from the one-vacancy bank alone, before the two-vacancy beam splitter is read.
        e, v = tree[0]
        amp = np.exp(-1j * np.outer(pilot_times, e - e[0])) @ (
            (right.conj() @ v) * (v.T @ left)
        )
        pilot = rng.binomial(2**16, np.clip(abs(amp) ** 2, 0, 1))
        j = np.argmin(abs(pilot / 2**16 - 0.5))
        t = int(pilot_times[j])
        controls = split_controls(tree[0], left, right, t)
        open_counts = bank(rng, positions(propagate(tree[1], msrc, t)))
        other, _ = prep_two(left, source(band, "right", rail=1))
        disjoint = bank(rng, positions(propagate(tree[1], other, t)))
        different, _ = prep_two(left, source((band + 1) % 3, "right"))
        mismatch = bank(rng, positions(propagate(tree[1], different, t)))
        single = bank(rng, positions(propagate(cycle[0], left, DELAY)))
        stock = bank(rng, positions(propagate(cycle[1], msrc, DELAY)))
        minus = positions(propagate(cycle[2], fsrc, DELAY))
        alternate = positions(propagate(flux_pi, left, DELAY))
        oneprob = positions(propagate(cycle[0], left, DELAY))
        roundtrip = positions(propagate(cycle[1], msrc, DELAY))
        # Two independently prepared rails on the opened guide: project either rail into its L/R positions.
        physical_left = np.isin(
            np.arange(SIZE), np.r_[paired.LEFT, SIZE - 1 - paired.LEFT]
        )
        physical_right = np.isin(
            np.arange(SIZE), np.r_[paired.RIGHT, SIZE - 1 - paired.RIGHT]
        )
        lr = (physical_left[FIRST] & physical_right[SECOND]) | (
            physical_right[FIRST] & physical_left[SECOND]
        )
        samples.append(
            dict(
                band=band,
                split_exposure=t,
                pilot_counts=pilot.tolist(),
                calibrated_single_transfer=float(abs(amp[j]) ** 2),
                open_path_joint=open_counts.tolist(),
                open_path_regions=coarse(open_counts),
                independent_rail_regions=coarse(disjoint),
                independent_rail_joint=disjoint.tolist(),
                independent_rail_one_each_well=int(sum(disjoint[lr])),
                next_label_joint=mismatch.tolist(),
                next_label_regions=coarse(mismatch),
                **controls,
                one_loop_by_port=single.tolist(),
                two_loop_by_port_pair=stock.tolist(),
                loop_regions=coarse(stock),
                native_two_a=float(sum(roundtrip[SECOND < N])),
                fermion_same_interactions_two_a=float(sum(minus[SECOND < N])),
                fermion_same_interactions_all=minus.tolist(),
                single_native_a=float(sum(oneprob[:N])),
                single_pi_flux_a=float(sum(alternate[:N])),
                two_a_interval=binomial_interval(
                    int(sum(stock[SECOND < N])), SHOTS, 0.01 / 6
                ).tolist(),
                one_a_interval=binomial_interval(
                    int(sum(single[:N])), SHOTS, 0.01 / 6
                ).tolist(),
            )
        )
    return dict(
        shots_per_setting=SHOTS,
        pilot_shots=2**16,
        pilot_times=pilot_times.tolist(),
        loop_exposure=DELAY,
        hardware=hw,
        native_phase=phase_witness(rng),
        cases=samples,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    dest = p.parse_args().output
    result = run()
    if dest:
        with dest.open("x") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        for row in result["cases"]:
            print(
                row["band"], row["loop_regions"], row["fermion_same_interactions_two_a"]
            )
