"""All labelled replies to observer zero, with no enrollment radius or timeout."""

import argparse
from collections import deque
from functools import cache
import json
from pathlib import Path

from observer_zero_capacity import component
from observer_zero_reach import echo_waits, return_delays
from strict_history.wide_replay import SCHEME, WideReplay


def antipodal_challenge(model, labels):
    faces = model.base_faces.tolist()
    vertices = model.base_vertices.tolist()
    opposite = {vertices.index([-x for x in vertices[i]]) for i in faces[0]}
    target = next(i for i, face in enumerate(faces) if set(face) == opposite)
    height = {v: 0 for v in faces[0]}
    queue = deque(height)
    while queue:
        v = queue.popleft()
        for face in faces:
            if v in face:
                for w in face:
                    if w not in height:
                        height[w] = height[v] + 1
                        queue.append(w)
    total = sum(height[v] for v in faces[target])

    # The three-vertex height sum doubles each refinement on the two
    # constant-height end faces. It changes by at most two across a shared edge.
    def bound(level):
        crossings = (total * 2**level + 1) // 2
        return 2 * (2 * crossings - 1)

    controls = []
    for level in (0, 1, 2):
        small = WideReplay(level)

        def transparent(p):
            return [(q, 0) for q, _ in small.links(p)]

        distance = component(transparent, 0)
        first = 2 * min(
            t for p, t in distance.items() if p // (12 * 4**level) == target
        )
        if first < bound(level):
            raise ValueError("Transparent echo precedes the architecture bound")
        controls.append(
            dict(level=level, lower_bound_ticks=bound(level), first_echo_wait=first)
        )
    return dict(
        source_base_face=0,
        target_base_face=target,
        target_vertex_heights=[height[v] for v in faces[target]],
        minimum_roundtrip_ticks=str(bound(model.level)),
        target_echoes_ever=sum(p // (12 * 4**model.level) == target for p in labels),
        controls_all_contacts_transmitting=controls,
    )


def measure(digits=110):
    model = WideReplay(202, digits, 20261126)
    links = cache(model.links)
    home = component(links, 0)
    # After every incident closure the home component is stationary. Outside
    # repairs can no longer deliver a signal through its closed boundary.
    last = max(t for u in home for _, t in links(u))
    boundary = [(u, v, t) for u in home for v, t in links(u) if v not in home]
    if any(t == 0 or (u, t) not in links(v) for u, v, t in boundary):
        raise ValueError("Endpoints disagree about the sealed boundary")
    delays = return_delays(links, last)
    rows = [echo_waits(links, delays, s) for s in range(last + 1)]
    expected = {p: 2 * d for p, d in home.items() if p // 12}
    if rows[-1] != expected:
        raise ValueError("Forward echoes disagree with the stationary return clock")
    for early, late in zip(rows, rows[1:]):
        if late.keys() - early.keys() or any(
            early[p] > wait for p, wait in late.items()
        ):
            raise ValueError("A closed path reopened or shortened a fixed pulse flight")
    labels = sorted(rows[0])
    antipode = antipodal_challenge(model, labels)
    epochs = []
    for send, row in enumerate(rows):
        if send and row == rows[send - 1]:
            epochs[-1]["last_send_tick"] = send
        else:
            epochs.append(
                dict(
                    first_send_tick=send,
                    last_send_tick=send,
                    responding_remote_observers=len({p // 12 for p in row}),
                    responding_remote_ports=len(row),
                    minimum_echo_wait=min(row.values()),
                    maximum_echo_wait=max(row.values()),
                    echo_wait_by_label=[row.get(p) for p in labels],
                )
            )
    return dict(
        level=202,
        seed=model.seed,
        key_scheme=SCHEME,
        receiver=["0", 0],
        apparatus="Individually labelled immediate reflections on matched copies. One tick per material transfer. Observer zero records send and return times.",
        source_labels=[[str(p // 12), p % 12] for p in labels],
        antipodal_fresh_challenge=antipode,
        surveys=epochs,
        stationary_from_send=epochs[-1]["first_send_tick"],
        stationary_continuation="For every later send with the same material clock and relay, the final row repeats.",
        completeness_check=dict(
            final_component_ports=len(home),
            closed_external_boundary_links=len(boundary),
            last_incident_closure_tick=last,
            proof="At and after the last incident closure, only the transmitting home component can reach the receiver. The backward arrival recurrence exhausts all earlier departures. Forward queries are followed only if a return remains possible. No time or distance cutoff.",
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--digits", type=int, default=110)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.write_text(json.dumps(measure(args.digits), indent=1) + "\n")
