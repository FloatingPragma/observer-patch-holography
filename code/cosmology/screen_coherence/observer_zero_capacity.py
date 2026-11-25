"""Measure the settled material-relay boundary and range every surviving remote port."""

import argparse
import json
import math
from collections import deque
from pathlib import Path

from observer_zero_radar import MANIFEST_SHA, ArchiveRadar, echo_table, local_graph


def component(links, origin):
    distance = {origin: 0}
    queue = deque([origin])
    while queue:
        u = queue.popleft()
        for v, tick in links(u):
            if tick == 0 and v not in distance:
                distance[v] = distance[u] + 1
                queue.append(v)
    return distance


def measure(archive):
    source = ArchiveRadar(archive)
    distances = component(source.links, 0)
    cut = sorted(
        (u, v, t) for u in distances for v, t in source.links(u) if v not in distances
    )
    if any(not t or (u, t) not in source.links(v) for u, v, t in cut):
        raise ValueError("Boundary write missing from peer history")
    if len({p for u, v, t in cut for p in (u, v)}) != 2 * len(cut):
        raise ValueError("Boundary writes do not have disjoint endpoints")
    radius = max(distances.values())
    # Every finally reachable reflector always has its static round trip available.
    # Any faster route must also stay inside the intact ball of that radius.
    ports, _, edges = local_graph(source.links, 0, 2 * radius)
    targets = [i for i, p in enumerate(ports) if p // 12 != 0 and p in distances]
    if len(targets) != sum(p // 12 != 0 for p in distances):
        raise ValueError("A reachable remote reflector was omitted")
    last = source.manifest["settlement_round"]
    ticks = list(range(last + 17))
    echoes = echo_table(edges, targets, ticks, 2 * radius, last)
    expected = [2 * distances[ports[i]] for i in targets]
    if any(None in row for row in echoes) or any(
        row != expected for row in echoes[last:]
    ):
        raise ValueError("Static return-distance check failed")
    totals = [sum(row) for row in echoes]
    mean = [n / (2 * len(targets)) for n in totals]

    def derivatives(center, step):
        left, mid, right = (mean[center + d] for d in (-step, 0, step))
        return dict(
            center=center,
            half_window=step,
            mean_range_speed=(right - left) / (2 * step),
            log_range_rate=math.log(right / left) / (2 * step),
            mean_range_acceleration=(right - 2 * mid + left) / step**2,
            relative_range_acceleration=(right - 2 * mid + left) / (step**2 * mid),
        )

    # Displayed gauge-width branch in the discrete spherical capacity paper.
    pixel = 1.6309682414
    n_screen = source.manifest["observers"] * pixel / 4
    lambda_lp2 = 3 * math.pi / n_screen
    return dict(
        manifest_sha256=MANIFEST_SHA,
        conditional_whole_screen=dict(
            pixel_area_in_planck_areas=pixel,
            assumptions="Each observer represents one cell of that physical area and the whole-screen area entropy is the correctable de Sitter horizon entropy.",
            source="paper/tex_fragments/DISCRETE_SPHERICAL_CAPACITY.tex",
            capacity_nats=n_screen,
            lambda_times_planck_length_squared=lambda_lp2,
            flat_vacuum_H_times_planck_time=math.sqrt(lambda_lp2 / 3),
            vacuum_relative_acceleration_times_planck_time_squared=lambda_lp2 / 3,
        ),
        origin=[0, 0],
        component_owners=len({p // 12 for p in distances}),
        component_ports=len(distances),
        max_one_way_relay_length=radius,
        cut_links=len(cut),
        cut=[[list(divmod(u, 12)), list(divmod(v, 12)), t] for u, v, t in cut],
        ideal_native_fair_cut_coin_entropy_nats=len(cut) * math.log(2),
        whole_sphere_binary_preparation_log_cardinality_nats=source.manifest["ports"]
        * math.log(2),
        physical_horizon_capacity_nats=None,
        cosmological_constant=None,
        all_terminal_remote_ports=len(targets),
        origin_emission_ticks=ticks,
        total_round_trip_ticks=totals,
        mean_radar_range=mean,
        last_echo_change=max(
            (t for t in ticks[1:] if echoes[t] != echoes[t - 1]), default=0
        ),
        estimates=[derivatives(4, 4), derivatives(last + 8, 4)],
        reflectors=[list(divmod(ports[i], 12)) for i in targets],
        echo_waits=echoes,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    print(json.dumps(measure(parser.parse_args().archive), indent=1, allow_nan=False))
