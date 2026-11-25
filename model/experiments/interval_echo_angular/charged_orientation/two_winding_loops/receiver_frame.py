"""Geometric frame controls for a fixed gateway and two local triangle head words."""

from itertools import combinations
import json

from model.executable_model.symmetry import rotations
from ..live_trimer.measure import geometry


def compose(outer, inner):
    return tuple(outer[j] for j in inner)


def moved(frame, rotation):
    gateway, arrows = frame
    return rotation[gateway], tuple(
        sorted((rotation[a], rotation[b]) for a, b in arrows)
    )


def inspect():
    board, ports = geometry()
    target = ports[3:]
    assert all(p // 12 == target[0] // 12 for p in target)
    assert all(
        b in board.neighbours[a] for a, b in zip(target, target[1:] + target[:1])
    )
    g, i, j = (int(p % 12) for p in target)
    arms = [
        (g, tuple(sorted(arrows)))
        for arrows in [[(g, i), (i, j), (j, g)], [(i, g), (j, i), (g, j)]]
    ]
    group = rotations()
    ident = tuple(range(12))
    halves = [r for r in group if r != ident and compose(r, r) == ident]
    orthogonal = [
        (r, s) for r, s in combinations(halves, 2) if compose(r, s) == compose(s, r)
    ]
    orbits = [{moved(arm, r) for r in group} for arm in arms]
    assert orbits[0].isdisjoint(orbits[1])
    frames = orbits[0] | orbits[1]
    for first, second in orthogonal:
        for frame in frames:
            assert moved(moved(frame, first), second) == moved(
                moved(frame, second), first
            )
    fixed_port = [r for r in group if r[g] == g]
    return dict(
        global_ports=list(map(int, target)),
        gateway=g,
        turns=len(group),
        proper_half_turns=len(halves),
        perpendicular_pairs=len(orthogonal),
        fixed_port_turns=len(fixed_port),
        fixed_port_and_face_turns=sum({r[i], r[j]} == {i, j} for r in fixed_port),
        fixed_port_half_turns=sum(r[g] == g for r in halves),
        frame_orbit_sizes=[len(orbit) for orbit in orbits],
        compared_frame_orderings=len(frames) * len(orthogonal),
    )


if __name__ == "__main__":
    print(json.dumps(inspect(), indent=2))
