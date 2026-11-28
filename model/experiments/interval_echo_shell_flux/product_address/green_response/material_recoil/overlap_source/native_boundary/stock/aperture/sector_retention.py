"""Track modular and ordinary calorimeter records at successive native returns."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np
from scipy.sparse import eye
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu

from . import peer_acquisition as acquisition
from . import landmark_arrival as arrival
from . import peer_population
from . import sector_cadence
from . import sector_sampling
from . import sector_thermal as thermal

HERE = acquisition.HERE


class ReturnClock:
    """Adjoint return operator with checked dyadic harmonic-extension error."""

    den = 128**3
    grid = 2**32

    def __init__(self, sweep, subset):
        outside = np.ones(sweep.shape[0], bool)
        outside[subset] = False
        self.tt = sweep[outside][:, outside].tocsr()
        self.ts = sweep[outside][:, subset].tocsr()
        self.ss = sweep[subset][:, subset].tocsr()
        self.st = sweep[subset][:, outside].tocsr()
        self.solve = splu((eye(self.tt.shape[0]) - self.tt / self.den).tocsc())
        mean = self.solve.solve(np.ones(self.tt.shape[0]))
        cap = np.ceil(2 * mean).astype(np.int64)
        assert np.all(self.den * cap - self.tt @ cap >= self.den) and min(cap) >= 1
        self.exit_bound = int(cap.max())

    def advance(self, ticks):
        # Each row has exactly den tickets. The largest intermediate sum,
        # including a residual subtraction, must fit int64.
        assert (abs(ticks).max() + self.grid) * self.den < 2**61
        rhs = self.ts @ ticks
        harmonic = np.rint(self.solve.solve(rhs.astype(float) / self.den)).astype(
            np.int64
        )
        assert (abs(harmonic).max() + self.grid) * self.den < 2**61
        residual = abs(self.den * harmonic - self.tt @ harmonic - rhs).max(0)
        raw = self.ss @ ticks + self.st @ harmonic
        updated = (raw + self.den // 2) // self.den
        error = (1 + self.exit_bound * residual / self.den) / self.grid
        return updated, error


def run():
    receipt = json.loads((HERE / "sector_calorimetry.json").read_bytes())
    w, _, groups, blocks, ops = thermal.apparatus()
    ordinary = peer_population.fibers.source.records()[1]["original"]
    _, arms = peer_population.histogram(16)
    d = np.array(arms["agreed"] - ordinary, float)
    full, initial, _, service, _ = acquisition.returns.retention.setup()
    lookup = {tuple(row): i for i, row in enumerate(full)}
    subset = np.array([lookup[tuple(row)] for row in w])
    assert set(subset) == set(np.flatnonzero(initial))
    M = [mat for mat, _ in service]
    sweep = (M[2] @ M[1] @ M[0]).T.tocsr()
    den = ReturnClock.den
    assert np.all(sweep.sum(1) == den) and np.all(sweep.sum(0) == den)
    assert connected_components(sweep, directed=True, connection="strong")[
        0
    ] == 1 and np.all(sweep.diagonal() > 0)
    clock = ReturnClock(sweep, subset)
    e = ops[0].diagonal(axis1=1, axis2=2).sum(0) / 4
    gmean = np.array([e[inds].mean() for inds in blocks])
    prices = np.array(receipt["sector_log_enclosures"])
    centered = prices.mean(1)[groups]
    sector = centered - thermal.BETA * gmean[groups]
    # Feed the very same returned words to two independent material clocks.
    board, marks, probes, _ = arrival.pulse_inputs()
    adj, slots = arrival.passage.contacts(board, marks)
    pins = [25, 36, *marks[-1]]
    F, denominator = arrival.blank_exits(adj, pins)
    stop = np.array(pins)
    exterior = np.ones(board.n, bool)
    exterior[stop] = False
    bearings = []
    stopped_chart = None
    checked_tickets = 0
    for load in (1, 2):
        table = probes[load]
        ports = (table % board.n).reshape(11, board.n, 24)
        assert np.all(ports[:, exterior] == ports[0, exterior])
        assert np.array_equal(
            np.sort(ports[0, exterior], axis=1),
            np.sort(np.repeat(slots[exterior], 4, axis=1), axis=1),
        )
        # Observers stop at either receiving beacon or first source contact.
        # Agreement for each ticket also covers complete stopped transcripts
        # when both arms share an ordinary launch and sampling rule.
        stopped = ports.copy()
        stopped[:, stop, :] = stop[None, :, None]
        if stopped_chart is None:
            stopped_chart = stopped[0]
        assert np.all(stopped == stopped_chart)
        checked_tickets += int(ports[:, exterior].size)
        field = arrival.eventual(table, pins, F, denominator)
        bearings.append(np.array(field[w[:, -1], 104], float))
    scores = np.column_stack(
        (e, sector, thermal.BETA * e + sector, centered, *bearings)
    )
    names = (
        "contact_energy",
        "conserved_sector_term",
        "retained_modular_nats",
        "centered_spectral_price",
        "single_arrival",
        "double_arrival",
    )
    ticks = np.rint(scores * clock.grid).astype(np.int64)
    err = np.full(len(names), 1 / clock.grid + np.ptp(prices, axis=1).max() + 1e-12)
    deadlines = (0, 1, 2, 4, 8, 16, 32, 64, 128, 256)
    bounds = {}
    end = 2048
    series = []
    sector_series = []
    flights = []
    for k in range(end + 1):
        if k in deadlines:
            vals = ticks / clock.grid
            change = d @ vals
            width = abs(d).sum() * err + 1e-11
            bounds[str(k)] = {
                name: [v - z, v + z]
                for name, v, z in zip(names, change, width, strict=True)
            }
        q = ticks[:, 2] / clock.grid
        value = d @ q
        span = float(np.ptp(q) + 2 * err[2])
        width = abs(d).sum() * err[2] + 1e-11
        series.append([value - width, value + width])
        sv = ticks[:, 1] / clock.grid
        sector_val = d @ sv
        sector_span = float(np.ptp(sv) + 2 * err[1])
        sector_width = abs(d).sum() * err[1] + 1e-11
        sector_series.append([sector_val - sector_width, sector_val + sector_width])
        read = ticks[:, 4:] / clock.grid
        flight = d @ read
        fend = abs(d).sum() * err[4:] + 1e-11
        flight_tail = abs(d).sum() * (np.ptp(read, axis=0) + 2 * err[4:]) / 2 + 1e-11
        flights.append(np.stack((flight - fend, flight + fend), axis=1))
        ticks, step = clock.advance(ticks)
        err += step
    # A stochastic return kernel cannot increase the oscillation of a
    # bounded readout. The zero-sum contrast bounds every later deadline.
    tail = span * abs(d).sum() / 2 + 1e-11
    rate_scan = {}
    sector_rate_scan = {}
    flight_rates = {}
    flights = np.array(flights)
    sector_tail = sector_span * abs(d).sum() / 2 + 1e-11
    sector_series = np.array(sector_series)
    series = np.array(series)
    age = np.arange(1, end + 1)
    for j in range(13):
        p = 2.0**-j
        measure = p * (1 - p) ** (age - 1)
        rem = (1 - p) ** end
        assert abs(1 - measure.sum() - rem) < 1e-12
        edges = measure @ series[1:] + np.array([-tail, tail]) * rem
        key = str(Q(1, 2**j))
        rate_scan[key] = edges.tolist()
        sector_rate_scan[key] = (
            measure @ sector_series[1:] + np.array([-sector_tail, sector_tail]) * rem
        ).tolist()
        measured = np.einsum("k,kib->ib", measure, flights[1:]) + rem * np.stack(
            (-flight_tail, flight_tail), axis=1
        )
        flight_rates[key] = dict(
            single=measured[0].tolist(), double=measured[1].tolist()
        )
    envelope = [min(series[:, 0].min(), -tail), max(series[:, 1].max(), tail)]
    # Independent pulse trajectories, stratified by first word, retaining
    # the same carrier throughout all subsequent accepted returns.
    tables = acquisition.returns.ticket_tables(M)
    rng = np.random.default_rng(325677)
    copies = 128
    label = np.repeat(np.arange(len(w)), copies)
    sites = subset[label]
    native = np.full(len(full), -1)
    native[subset] = np.arange(len(w))
    held = {}
    waits = 0
    for k in range(1, 17):
        pending = np.arange(len(sites))
        for sweep_number in range(1, 4097):
            for table in tables:
                sites[pending] = table[
                    sites[pending], rng.integers(128, size=len(pending))
                ]
            pending = pending[native[sites[pending]] < 0]
            if not len(pending):
                break
        assert not len(pending)
        waits = max(waits, sweep_number)
        if k in (1, 4, 16):
            values = scores[native[sites], 2].reshape(len(w), copies)
            mean = float(d @ values.mean(1))
            width = float(
                np.ptp(scores[:, 2])
                * math.sqrt((d @ d) * math.log(6 / 0.01) / (2 * copies))
            )
            lo, hi = bounds[str(k)][names[2]]
            assert mean - width <= lo <= hi <= mean + width
            held[str(k)] = dict(
                mean=mean,
                bounds=[mean - width, mean + width],
                terminal_word_counts=np.bincount(
                    native[sites], minlength=len(w)
                ).tolist(),
            )
    calibrations = np.array(
        [row["centered_log_band"] for row in receipt["ordinary_clock_calibration"]]
    )
    after_subtraction = (calibrations - thermal.BETA * gmean[:, None])[groups]
    detected = sector_sampling.acquire(tables, subset, native, d, after_subtraction)
    for key, row in detected.items():
        lo, hi = sector_rate_scan[key]
        assert row["bounds"][0] <= lo <= hi <= row["bounds"][1]
    contraction = sector_cadence.mixing(clock, len(w))
    cadence_total = sector_cadence.all_rates(series, tail, contraction)
    cadence_sector = sector_cadence.all_rates(sector_series, sector_tail, contraction)
    return dict(
        sha256={
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in ("sector_calorimetry.json", "landmark_returns.json")
        },
        phase_order=[0, 1, 2],
        ordinary_reference="uniform native first return",
        uniform_limit_change="0",
        held_geometric_records=detected,
        held_geometric_and_clock_coverage_at_least=0.993,
        maximum_mean_outside_sweeps_bound=clock.exit_bound,
        dyadic_grid=clock.grid,
        contrasts=bounds,
        last_age=end,
        all_later_modular_contrast_bound=tail,
        geometric_sampling=rate_scan,
        sector_after_ordinary_heat_subtraction=sector_rate_scan,
        material_arrival_cadence=flight_rates,
        before_source_contact_common_contrast_all_ages="0",
        stopped_position_transcripts=dict(
            exterior_ticket_entries_checked=checked_tickets,
            source_contacts=list(map(int, pins[2:])),
            absorbing_beacons=pins[:2],
            loads=[1, 2],
            projection="position history through first arrival at a beacon or source contact",
            same_stopped_update_for_every_ticket=True,
        ),
        all_age_sector_excess_nats=[
            min(sector_series[:, 0].min(), -sector_tail),
            max(sector_series[:, 1].max(), sector_tail),
        ],
        continuous_cadence=cadence_total,
        continuous_sector_excess=cadence_sector,
        all_word_independent_age_mixtures_nats=envelope,
        largest_upper_at_computed_age=int(np.argmax(series[:, 1])),
        smallest_lower_at_computed_age=int(np.argmin(series[:, 0])),
        held_native_copies_per_word=copies,
        held_contrasts=held,
        unreturned=0,
        longest_observed_excursion=waits,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "sector_retention.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        from .closed_balls import radar

        radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_text()), out
        )
    print(
        {
            k: out["contrasts"][k]["retained_modular_nats"]
            for k in ("0", "1", "16", "256")
        }
    )
