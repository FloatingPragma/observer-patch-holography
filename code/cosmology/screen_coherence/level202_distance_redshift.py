"""Observer-zero pulse/radar relation and a comoving de Sitter clock comparison."""

import argparse
from collections import Counter
from decimal import Decimal as D, ROUND_CEILING, localcontext
from pathlib import Path
import json

from distance_redshift import PERIOD, TIMEOUT, inward, timed_sources
from observer_zero_radar import local_graph
from strict_history.wide_replay import WideReplay, SCHEME

LEVEL = 202
SEED = 20261126
LATE = 96
PI = D("3.1415926535897932384626433832795028841971693993751058209749445923078164")
PIXEL = D("1.6309682414")


def clock_comparison(distance, h, span):
    """Flat expanding de Sitter, comoving clocks. Inputs in receiver ticks, c=1."""
    with localcontext() as ctx:
        # h is ~1e-61. The centered finite-pair logarithms need ~125 guard places.
        ctx.prec = 210
        d = D(distance)
        length = h * d
        dx = (2 * length).exp()
        u = (dx - 1) / (dx + 1)
        half = D(PERIOD) * h / 2
        z = ((1 - u * (-half).exp()) / (1 - u * half.exp())).ln() / (2 * half)
        growth = (h * span).exp()
        d2 = -(1 - growth * (1 - (-2 * length).exp())).ln() / (2 * h)
        excess = ((1 - u) / (1 - u * growth)).ln() / h
        return dict(
            z_centered=float(z),
            next_radar_distance=float(d2),
            extra_received_ticks_for_long_source_period=float(excess),
        )


def trace(model):
    ports, depths, edges = local_graph(model.links, 0, 2 * TIMEOUT)
    sources = [
        p for p, depth in zip(ports, depths, strict=True) if p // 12 and depth <= 12
    ]
    lookup = {p: i for i, p in enumerate(ports)}
    blank = [{v: 0 for v in row} for row in edges]
    panels = {}
    late = []
    for epoch in (0, 4, 16, LATE, "transparent"):
        rows = timed_sources(
            ports,
            blank if epoch == "transparent" else edges,
            sources,
            0 if epoch == "transparent" else epoch,
        )
        for row in rows:
            del row["flash"]
            row["source"][0] = str(row["source"][0])
        panels[str(epoch)] = rows
        if epoch == LATE:
            late = rows
    assert all(
        r["z_pulse"] == 0 and r["radar_distance"] == depths[lookup[p]]
        for p, r in zip(sources, panels["transparent"], strict=True)
    )
    last = max(t for connections in edges for t in connections.values())
    if last > LATE - PERIOD // 2:
        raise ValueError("Contact changes could affect a comparison at the late epoch")
    with localcontext() as ctx:
        ctx.prec = 140
        h = (PI / (5 * PIXEL * 4**LEVEL)).sqrt()
        span = int((D(2).ln() / h).to_integral_value(rounding=ROUND_CEILING))
    after = timed_sources(ports, edges, sources, LATE + span)
    # Store differences rather than huge absolute timestamps for numeric interoperability.
    continuations = []
    for local, remote, p in zip(late, after, sources, strict=True):
        for key in ("reflection_at", "echo_at"):
            t1, t2 = local[key], remote[key]
            assert t2 is None if t1 is None else t2 - t1 == span
        assert (remote["z_pulse"], remote["radar_distance"]) == (
            local["z_pulse"],
            local["radar_distance"],
        )
        if local["z_pulse"] is None or local["radar_distance"] is None:
            continue
        m = local["reflection_at"]
        # Independently triggered source pulses at m, m+span in matched copies.
        first = inward(edges, lookup[p], m, TIMEOUT)
        second = inward(edges, lookup[p], m + span, TIMEOUT)
        assert first is not None and second is not None
        continuations.append(
            dict(
                source=local["source"],
                initial_radar_distance=local["radar_distance"],
                later_roundtrip_ticks=remote["echo_at"] - remote["query_at"],
                first_source_pulse=m,
                first_received=first,
                late_source_pulse_offset=0,
                late_receipt_offset=second - first - span,
                predicted=clock_comparison(local["radar_distance"], h, span),
            )
        )
    return dict(
        level=LEVEL,
        seed=SEED,
        key_scheme=SCHEME,
        receiver=[0, 0],
        query_epochs=[0, 4, 16, LATE],
        four_tick_period=PERIOD,
        initial_aperture_relay_units=12,
        maximum_round_trip_ticks=TIMEOUT,
        enrolled=len(sources),
        complete_dependency_ball_ports=len(ports),
        last_closure_in_ball=last,
        panels=panels,
        paired_z_counts={
            k: dict(
                Counter(
                    str(r["z_pulse"])
                    for r in rows
                    if r["z_pulse"] is not None and r["radar_distance"] is not None
                )
            )
            for k, rows in panels.items()
        },
        vacuum_comparison=dict(
            h_per_tick_if_tick_is_planck_time=str(h),
            duration_ticks=str(span),
            input_pixel_in_planck_areas=str(PIXEL),
            assumptions="Physical de Sitter horizon entropy, expanding flat vacuum chart, comoving source clocks and observer, common null radar ruler, tau=t_P.",
            late_source_pulse_offset_definition="Late source pulse time minus first source pulse time minus duration. Receipt offset is the analogous local receiver difference.",
            continuation_basis="Stationary finite dependency ball. Two widely separated local services queried directly. Clock and beacon identities assumed to persist between services.",
            matched_long_baseline=continuations,
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=110)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(
        json.dumps(
            trace(WideReplay(LEVEL, args.digits, SEED)), indent=1, allow_nan=False
        )
        + "\n"
    )
