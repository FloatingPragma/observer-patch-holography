"""Test memory in angularly equilibrated native radial preparations."""

import argparse
from collections import deque
from functools import lru_cache
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.special import jv

from model.executable_model import cell, federation
from model.executable_model.refinement import load_refinement
from ...native_tidal_response.measure import build_station, local_calibration


HERE = Path(__file__).resolve().parent
SHOTS = 2**29
SEED = 631080
PHASES = (0.0, math.pi / 2, math.pi, -math.pi / 2, math.pi / 3)
PHASE_RADIUS = math.sqrt(math.log(2 * 80 / 0.01) / (2 * SHOTS))


class ActiveRoutes:
    """Matched-ray preparation with one closed seam at each mobile port."""

    def __init__(self, car):
        self.external = np.full(12 * car.cells, -1, dtype=np.int32)
        u, v = car.external.T
        self.external[u], self.external[v] = v, u
        self.local = cell().neighbours

    @lru_cache(maxsize=256)
    def pattern(self, gates):
        ports = set(range(12))
        resident = min(
            ports - set(gates), key=lambda p: (-len(set(gates) & set(self.local[p])), p)
        )
        close = {tuple(sorted((resident, v))) for v in self.local[resident]}
        other = ports - {resident} - set(self.local[resident])
        allowed = sorted(
            (u, v) for u in other for v in self.local[u] if u < v and v in other
        )
        matching = next(
            es for es in combinations(allowed, 3) if {v for e in es for v in e} == other
        )
        close.update(matching)
        return resident, tuple(
            tuple(v for v in ns if tuple(sorted((p, v))) not in close)
            for p, ns in enumerate(self.local)
        )

    @lru_cache(maxsize=100000)
    def cell_pattern(self, c):
        gates = tuple(np.flatnonzero(self.external[12 * c : 12 * c + 12] >= 0).tolist())
        return self.pattern(gates)

    def anchor(self):
        p, _ = self.cell_pattern(0)
        return min(v for v in self.local[p] if self.external[v] >= 0)

    def neighbours(self, site):
        owner, p = divmod(site, 12)
        _, adjacent = self.cell_pattern(owner)
        out = tuple(12 * owner + v for v in adjacent[p])
        across = self.external[site]
        return out if across < 0 else (*out, int(across))

    def distances(self, start, stop):
        q = deque([start])
        distance = {start: 0}
        while q:
            site = q.popleft()
            step = distance[site] + 1
            if step <= stop:
                for dest in self.neighbours(site):
                    if dest not in distance:
                        distance[dest] = step
                        q.append(dest)
        return distance


def chebyshev(scales):
    """Uniform operator-norm tail bound from the real-argument Bessel integral."""
    z = max(scales)
    order = math.ceil(z)
    while True:
        first = (z / 2) ** (order + 1) / math.factorial(order + 1)
        bound = 2 * first / (1 - z / (2 * (order + 2)))
        if bound < 1e-13:
            break
        order += 1
    indices = np.arange(order + 1)
    # The physically irrelevant common kinetic shift dmax/2 is excluded.
    series = np.array(
        [jv(indices, s) * np.where(indices, 2, 1) * (-1j) ** indices for s in scales]
    )
    return series, bound


def stencil(routes, width, window, order):
    last = window * width + order
    # Extra fringe completes the count of the outermost reachable readout bin.
    distances = routes.distances(routes.anchor(), ((last + width - 1) // width) * width)
    sites = np.array(sorted(k for k, d in distances.items() if d <= last))
    ids = {int(v): i for i, v in enumerate(sites)}
    bins = np.array([(distances[v] + width - 1) // width for v in sites], dtype=int)
    populations = np.bincount([(v + width - 1) // width for v in distances.values()])
    rows, cols, diagonal = [], [], []
    for i, v in enumerate(sites):
        near = routes.neighbours(int(v))
        if not 4 <= len(near) <= 5:
            raise ValueError("Mobile port lost an incorrect number of seams")
        diagonal.append(len(near) / 5 - 1)
        js = [ids[k] for k in near if k in ids]
        rows.extend([i] * len(js))
        cols.extend(js)
    scaled = sparse.diags(diagonal) - sparse.csr_matrix(
        (np.full(len(rows), 1 / 5), (rows, cols)), shape=(len(sites), len(sites))
    )
    return scaled.tocsr(), bins, populations


def columns(operator, source, coefficients, block=24):
    """One real Chebyshev recurrence for all time points and address columns."""
    for start in range(0, len(source), block):
        launches = source[start : start + block]
        prev = np.zeros((operator.shape[0], len(launches)))
        prev[launches, np.arange(len(launches))] = 1
        cur = operator @ prev
        waves = [c[0] * prev + c[1] * cur for c in coefficients]
        for m in range(2, coefficients.shape[1]):
            prev, cur = cur, 2 * (operator @ cur) - prev
            for wave, c in zip(waves, coefficients, strict=True):
                if c[m].real:
                    wave.real += c[m].real * cur
                elif c[m].imag:
                    wave.imag += c[m].imag * cur
        yield launches, waves


def phase_receipts(probabilities, level, window):
    receipts = []
    for charge in (1, 2):
        frequencies = []
        counts = []
        for setting, probability in enumerate(probabilities):
            rng = np.random.default_rng([SEED, level, window, 20 + setting, charge])
            returned = int(rng.binomial(SHOTS, probability))
            counts.append([returned, SHOTS - returned])
            frequencies.append(returned / SHOTS)
        c = (frequencies[0] + frequencies[2]) / 2
        a = (frequencies[0] - frequencies[2]) / 2
        b = frequencies[1] - c
        holdouts = []
        for i in (3, 4):
            co, si = math.cos(PHASES[i]), math.sin(PHASES[i])
            predicted = c + a * co + b * si
            weights = ((1 + co - si) / 2, si, (1 - co - si) / 2)
            bound = PHASE_RADIUS * (1 + sum(abs(v) for v in weights))
            holdouts.append(
                dict(
                    setting=i,
                    predicted=predicted,
                    observed=frequencies[i],
                    absolute_error=abs(predicted - frequencies[i]),
                    simultaneous_bound=bound,
                )
            )
        receipts.append(dict(load=charge, counts=counts, holdouts=holdouts))
    return receipts


def read_window(routes, width, window, level):
    time = width / 2
    c, bound = chebyshev([2.5 * time, 2.5 * 2 * time])
    h, bins, populations = stencil(routes, width, window, c.shape[1] - 1)
    source = np.flatnonzero(bins == window)
    # Store one-leg columns on the same finite propagation stencil as the direct control.
    # The two-leg cohort is computed directly on the larger reachable subgraph.
    leg = np.empty((h.shape[0], len(source)), dtype=complex)
    direct = 0.0
    for ids, (first, full) in columns(h, source, c):
        leg[:, np.searchsorted(source, ids)] = first
        direct += np.square(np.abs(full[source, :])).sum()
    direct /= len(source)
    probabilities = np.bincount(
        bins,
        weights=np.square(abs(leg)).sum(axis=1) / len(source),
        minlength=len(populations),
    )
    np.testing.assert_allclose(probabilities.sum(), 1.0, atol=4e-12, rtol=0)
    probabilities /= probabilities.sum()
    refresh = float(len(source) * np.sum(probabilities**2 / populations))
    # Symmetric U: restricted two-leg amplitudes are U_BA^T U_BA, not U_BA^† U_BA.
    # Distinct midpoint outcomes erase interference between these matrices.
    monitored = 0.0
    amplitude = np.zeros((len(source), len(source)), dtype=complex)
    for b in np.flatnonzero(probabilities > 1e-25):
        plane = leg[bins == b]
        through = plane.T @ plane
        if b == window:
            inside = through
        amplitude += through
        monitored += np.square(abs(through)).sum() / len(source)
    coherent = np.square(abs(amplitude)).sum() / len(source)
    np.testing.assert_allclose(coherent, direct, atol=5e-12, rtol=0)
    phase_probs = [
        float(
            np.square(abs(amplitude + (np.exp(1j * theta) - 1) * inside)).sum()
            / len(source)
        )
        for theta in PHASES
    ]
    outcomes = {}
    for arm, p in enumerate([direct, monitored, refresh]):
        rng = np.random.default_rng([SEED, level, window, arm])
        n = int(rng.binomial(SHOTS, p))
        double_rng = np.random.default_rng([SEED, level, window, arm, 2])
        double_returned = int(double_rng.binomial(SHOTS, p))
        outcomes[("uninterrupted", "window_observed", "angular_refresh")[arm]] = dict(
            return_probability=float(p),
            counts=[n, SHOTS - n],
            double_load_counts=[double_returned, SHOTS - double_returned],
        )
    rng = np.random.default_rng([SEED, level, window, 3])
    histogram = rng.multinomial(SHOTS, probabilities)
    estimated = histogram / SHOTS
    return dict(
        window=window,
        width=width,
        time=time,
        marked_gateway=routes.anchor(),
        launched_bearings=len(source),
        pinned_stencil_ports=h.shape[0],
        chebyshev_degree=c.shape[1] - 1,
        single_flight_norm_tail=bound,
        populations=populations.tolist(),
        one_leg_probabilities=probabilities.tolist(),
        one_leg_counts=histogram.tolist(),
        double_load_one_leg_counts=np.random.default_rng([SEED, level, window, 3, 2])
        .multinomial(SHOTS, probabilities)
        .tolist(),
        phase_probabilities=phase_probs,
        phase_acquisitions=phase_receipts(phase_probs, level, window),
        refreshed_prediction_from_one_leg_counts=float(
            len(source) * sum(estimated**2 / populations)
        ),
        coherent_composition=float(coherent),
        outcomes=outcomes,
        absolute_refresh_error=abs(float(direct - refresh)),
        error_fraction_of_departure=abs(float(direct - refresh)) / (1 - direct),
    )


def complete_control(routes, rows, blanks):
    # Reconstruct the whole native sector independently of the truncated stencil.
    free, _, enabled, _, _ = build_station(federation(1), blanks)
    distance = routes.distances(routes.anchor(), len(enabled))
    if set(distance) != set(enabled):
        raise ValueError(
            "Deposited mobile addresses differ from the whole native station"
        )
    bins = np.array([distance[int(s)] for s in enabled])
    values, vectors = np.linalg.eigh(free.toarray())
    flight = (vectors * np.exp(-0.5j * values)) @ vectors.T
    # Read a complete kernel directly in both launch directions.
    masks = [np.flatnonzero(bins == b) for b in range(int(bins.max()) + 1)]
    transitions = np.square(abs(flight))
    kernel = np.array(
        [
            [transitions[np.ix_(dest, src)].sum() / len(src) for dest in masks]
            for src in masks
        ]
    )
    checks = []
    for row in rows:
        window = row["window"]
        src = masks[window]
        reset = float(kernel[window, :] @ kernel[:, window])
        monitored = 0.0
        for mid in masks:
            wave = flight[np.ix_(src, mid)] @ flight[np.ix_(mid, src)]
            monitored += np.square(abs(wave)).sum() / len(src)
        phases = []
        for theta in PHASES:
            gate = np.where(bins == window, np.exp(1j * theta), 1.0)
            final = flight[src, :] @ (gate[:, None] * flight[:, src])
            phases.append(float(np.square(abs(final)).sum() / len(src)))
        expected = row["outcomes"]
        np.testing.assert_allclose(
            [reset, monitored],
            [
                expected["angular_refresh"]["return_probability"],
                expected["window_observed"]["return_probability"],
            ],
            rtol=0,
            atol=8e-12,
        )
        np.testing.assert_allclose(
            phases, row["phase_probabilities"], rtol=0, atol=8e-12
        )
        checks.append(
            dict(
                window=window,
                reset=reset,
                monitored=float(monitored),
                phase_probabilities=phases,
            )
        )
    reversed_probability = float(kernel[6, 3])
    reverse_returned = int(
        np.random.default_rng([SEED, 1, 3, 10]).binomial(SHOTS, reversed_probability)
    )
    predicted = rows[0]["one_leg_counts"][6] / SHOTS * len(masks[3]) / len(masks[6])
    return dict(
        complete_mobile_ports=len(enabled),
        windows=checks,
        reverse_window=6,
        read_window=3,
        measured_reverse_probability=reversed_probability,
        observed_counts=[reverse_returned, SHOTS - reverse_returned],
        predicted_from_forward_counts=predicted,
    )


def run(directory):
    blanks, cal = local_calibration()
    assert cal["coin_odd"]["source_energy"] == 32
    np.testing.assert_allclose(
        cal["coin_odd"]["measured_potential"], [0, 0], atol=1e-12
    )
    panels = []
    for j, level in enumerate((1, 3, 5, 7)):
        car = load_refinement(directory / f"FEDERATION_LEVEL{level}.npz", level=level)
        routes = ActiveRoutes(car)
        rows = [read_window(routes, 2**j, k, level) for k in (3, 6)]
        if level == 1:
            complete = complete_control(routes, rows, blanks)
        panels.append(
            dict(level=level, cells=car.cells, sha256=car.sha256, windows=rows)
        )
        routes.cell_pattern.cache_clear()
        routes.pattern.cache_clear()
        print(
            car.cells,
            [
                (r["absolute_refresh_error"], r["error_fraction_of_departure"])
                for r in rows
            ],
            flush=True,
        )
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        phase_protocol_sha256=hashlib.sha256(
            (HERE / "PHASE_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        phase_settings_radians=PHASES,
        phase_single_frequency_simultaneous_radius=PHASE_RADIUS,
        complete_board_control=complete,
        marker_calibration=cal["coin_odd"],
        shots=SHOTS,
        seed=SEED,
        panels=panels,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    packet = run(args.directory)
    with args.output.open("x") as stream:
        json.dump(packet, stream, indent=2, allow_nan=False)
        stream.write("\n")
