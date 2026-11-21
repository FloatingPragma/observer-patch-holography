"""An independent spherical ruler and a scheduled massless incidence radar."""

import itertools
import math

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.executable_model import carrier
from ..cage_transfer import PORTS_PER_CELL
from ..relativity.clock_read import two_read_records


TOLERANCE = 2e-10
DETUNING = 1.05


def _edge(a, b):
    return (a, b) if a < b else (b, a)


def _unit(vector):
    return vector / np.linalg.norm(vector)


def _angle(left, right):
    return math.atan2(float(np.linalg.norm(np.cross(left, right))), float(left @ right))


def _dual(faces):
    graph = nx.Graph()
    graph.add_nodes_from(range(len(faces)))
    owners = {}
    for index, face in enumerate(faces):
        for a, b in itertools.combinations(face, 2):
            owners.setdefault(_edge(a, b), []).append(index)
    if any(len(cells) != 2 for cells in owners.values()):
        raise ValueError("the spherical atlas has an unpaired triangle side")
    graph.add_edges_from(owners.values())
    return graph


def _subdivide(vertices, faces):
    """Construct geodesic midpoints and keep each parent’s center child fourth."""
    vertices = list(vertices)
    midpoints = {}

    def midpoint(a, b):
        edge = _edge(a, b)
        if edge not in midpoints:
            midpoints[edge] = len(vertices)
            vertices.append(_unit(vertices[a] + vertices[b]))
        return midpoints[edge]

    children = []
    for a, b, c in faces:
        ab, bc, ac = midpoint(a, b), midpoint(b, c), midpoint(a, c)
        children.extend(((a, ab, ac), (ab, b, bc), (ac, bc, c), (ab, bc, ac)))
    return vertices, children


def _child_map(parent_map, native, atlas):
    result = {}
    for parent, atlas_parent in parent_map.items():
        first, atlas_first = 4 * parent, 4 * atlas_parent
        result[first + 3] = atlas_first + 3
        # A corner is uniquely specified by its two sides crossing parent boundaries.
        atlas_corners = {
            frozenset(j // 4 for j in atlas[child] if j // 4 != atlas_parent): child
            for child in range(atlas_first, atlas_first + 3)
        }
        if len(atlas_corners) != 3 or any(len(sides) != 2 for sides in atlas_corners):
            raise ValueError("the atlas corners lack distinct parent-side signatures")
        for child in range(first, first + 3):
            sides = frozenset(parent_map[j // 4] for j in native[child] if j // 4 != parent)
            if sides not in atlas_corners:
                raise ValueError("a stored corner lacks its atlas parent-side signature")
            result[child] = atlas_corners[sides]
    return result


def _atlas(carriers):
    points, seams = carrier.icosahedron_from_coordinates()
    vertices = [_unit(point) for point in points]
    faces = [
        (a, b, c)
        for a, b, c in itertools.combinations(range(len(vertices)), 3)
        if {(a, b), (b, c), (a, c)}.issubset(seams)
    ]
    native_graphs, centers = [], []
    cell_map = None
    for level, native in enumerate(carriers):
        if level:
            vertices, faces = _subdivide(vertices, faces)
        atlas = _dual(faces)
        graph = nx.Graph()
        graph.add_nodes_from(range(native.cells))
        graph.add_edges_from((a // PORTS_PER_CELL, b // PORTS_PER_CELL)
                             for a, b in native.external)
        if cell_map is None:
            cell_map = nx.vf2pp_isomorphism(graph, atlas)
        else:
            cell_map = _child_map(cell_map, graph, atlas)
        if (cell_map is None or len(cell_map) != len(atlas)
                or len(set(cell_map.values())) != len(atlas)
                or len(graph.edges) != len(native.external)
                or {_edge(cell_map[a], cell_map[b]) for a, b in graph.edges}
                != {_edge(a, b) for a, b in atlas.edges}):
            raise ValueError("the parent-consistent atlas does not realize every native seam")
        centers.append(np.array([
            _unit(np.sum([vertices[v] for v in faces[cell_map[cell]]], axis=0))
            for cell in range(native.cells)
        ]))
        native_graphs.append(graph)

    for level, level_centers in enumerate(centers[1:], start=1):
        # Nested central triangles of each regular base face keep its spherical centroid.
        np.testing.assert_allclose(level_centers[4**level - 1::4**level], centers[0],
                                   atol=TOLERANCE, rtol=0)
    return native_graphs, centers


def _hop_obstruction(coarse, fine):
    """Test all distance-one and distance-two pairs as exact integer statements."""
    checked = {1: 0, 2: 0}
    for source in coarse:
        coarse_distances = nx.single_source_shortest_path_length(coarse, source, cutoff=2)
        fine_distances = nx.single_source_shortest_path_length(fine, 4 * source + 3, cutoff=4)
        for target, distance in coarse_distances.items():
            if source < target:
                # First and last central seams plus one or two parent crossings.
                if fine_distances.get(4 * target + 3) != distance + 2:
                    raise ValueError("a central-child hop distance differs from the coarse distance plus two")
                checked[distance] += 1
    if not all(checked.values()):
        raise ValueError("both distinct hop refinements require a native witness")
    return [{"coarse_hops": d, "fine_hops": d + 2, "checked_pairs": count}
            for d, count in checked.items()]


def _rest_clock(duration, frequency):
    # At most pi/2 between quadrature reads allows the elapsed phase to be unwrapped.
    intervals = max(1, math.ceil(duration * frequency / (math.pi / 2)))
    phases = []
    for index in range(intervals + 1):
        phase = frequency * duration * index / intervals
        record = two_read_records(np.array([0.5]), np.array([0.5 * np.exp(-1j * phase)]), 0)
        phases.append(record["phase"])
    steps = [math.atan2(math.sin(b - a), math.cos(b - a))
             for a, b in itertools.pairwise(phases)]
    return {"elapsed": math.fsum(steps) / frequency,
            "read_intervals": intervals, "terminal_pointer_record": record}


def _pulse(length, speed, wait):
    """Freeze a single real incidence seam and exponentiate the three-mode operator."""
    rate = speed / length
    generator = rate * np.array([[0, 0, 1], [0, 0, -1], [1, -1, 0]])
    return expm(-1j * wait * generator)


def _c_from_read(length, speed, frequency, wait):
    # Select the first rising half-period when enrolling the physical reference.
    if not 0 < math.sqrt(2) * speed * wait / length < math.pi:
        raise ValueError("the reference pulse must be on the first rising half-period")
    clock = _rest_clock(wait, frequency)
    occupation = np.abs(_pulse(length, speed, wait)[:, 0]) ** 2
    probability = float(occupation[1])
    measured_omega = math.acos(1 - 2 * math.sqrt(probability)) / clock["elapsed"]
    measured_speed = measured_omega * length / math.sqrt(2)
    return measured_speed, {
        "clock": clock,
        "receiver_probability": probability,
        "seam_probability": float(occupation[2]),
        "sender_probability": float(occupation[0]),
    }


def _route_read(route, graph, cell_count, true_speed, calibrated_speed, frequency):
    # Propagate on the whole cell and seam space so missed swaps remain in the state.
    edges = [_edge(a, b) for a, b in itertools.pairwise(route)]
    state = np.zeros(cell_count + len(edges), dtype=complex)
    state[route[0]] = 1.0
    waits = [math.pi * graph[a][b]["length"] / (math.sqrt(2) * calibrated_speed)
             for a, b in edges]
    schedules = [(a, b, cell_count + index, wait, graph[a][b]["length"])
                 for index, ((a, b), wait) in enumerate(zip(edges, waits, strict=True))]
    for a, b, seam, wait, length in schedules:
        state[[a, b, seam]] = _pulse(length, true_speed, wait) @ state[[a, b, seam]]
    received = float(abs(state[route[-1]]) ** 2)
    for a, b, seam, wait, length in reversed(schedules):
        state[[a, b, seam]] = _pulse(length, true_speed, wait) @ state[[a, b, seam]]
    returned = float(abs(state[route[0]]) ** 2)
    np.testing.assert_allclose([received, returned, float(np.vdot(state, state).real)],
                               1.0, atol=TOLERANCE, rtol=0)
    clock = _rest_clock(2 * math.fsum(waits), frequency)
    relay_distance = clock["elapsed"] * math.sqrt(2) * calibrated_speed / (2 * math.pi)
    atlas_distance = math.fsum(graph[a][b]["length"] for a, b in edges)
    np.testing.assert_allclose(relay_distance / atlas_distance, 1, atol=TOLERANCE, rtol=0)
    return {"cells": route, "receiver_probability": received, "return_probability": returned,
            "clock": clock, "relay_distance": relay_distance, "atlas_route_distance": atlas_distance}


def run(carriers=None, *, reference_length=1.0, speed=1.0, reference_frequency=1.0,
        calibration_wait=0.5):
    """Assign native seam lengths independently of the simulated packet evolution.

    The supplied base arc is the length standard, and a separate oscillator is
    the time standard. The pulse read determines the common incidence coefficient
    in these units. The scheduled relay rate is sqrt(2)/pi times that coefficient.
    """
    if not all(math.isfinite(value) and value > 0 for value in
               (reference_length, speed, reference_frequency, calibration_wait)):
        raise ValueError("the ruler, speed, reference frequency and wait must be positive and finite")
    if carriers is None:
        carriers = tuple(carrier.federation(level) for level in range(3))
    carriers = tuple(carriers)
    if tuple(native.cells for native in carriers) != (20, 80, 320):
        raise ValueError("use the three stored native refinement levels")
    graphs, centers = _atlas(carriers)
    reference_ports = tuple(carriers[0].external[0])
    reference_cells = tuple(port // PORTS_PER_CELL for port in reference_ports)
    radius = reference_length / _angle(*(centers[0][cell] for cell in reference_cells))
    calibrated_speed, calibration = _c_from_read(
        reference_length, speed, reference_frequency, calibration_wait
    )
    calibration.update({"reference_ports": list(reference_ports),
                        "reference_cells": list(reference_cells), "reference_length": reference_length})
    # Two fixed base stations persist as central descendants, including their coordinate positions.
    station_a = min(graphs[0])
    station_b = max(nx.single_source_shortest_path_length(graphs[0], station_a).items(),
                    key=lambda item: (item[1], -item[0]))[0]
    beacons = (station_a, station_b)
    geodesic = radius * _angle(*(centers[0][cell] for cell in beacons))

    levels = []
    for level, (native, graph, level_centers) in enumerate(zip(carriers, graphs, centers, strict=True)):
        seam_lengths = []
        max_swap_error = 0.0
        count = 0
        for a, b in native.external:
            u, v = a // PORTS_PER_CELL, b // PORTS_PER_CELL
            length = radius * _angle(level_centers[u], level_centers[v])
            if not math.isfinite(length) or length <= 0:
                raise ValueError("a native seam has no finite positive atlas arc")
            graph[u][v]["length"] = length
            seam_lengths.append([a, b, length])
            if level == 0 and (a, b) == reference_ports:
                continue
            wait = math.pi * length / (math.sqrt(2) * calibrated_speed)
            result = _pulse(length, speed, wait)[:, 0]
            max_swap_error = max(max_swap_error, abs(1 - float(abs(result[1]) ** 2)))
            count += 1

        # Keep the one reference seam out of the base-level held-out relay.
        relay_graph = graph.copy()
        if level == 0:
            relay_graph.remove_edge(*reference_cells)
        source, target = [4**level * cell + 4**level - 1 for cell in beacons]
        route = nx.shortest_path(relay_graph, source, target, weight="length")
        relay = _route_read(route, relay_graph, native.cells, speed, calibrated_speed, reference_frequency)
        if relay["atlas_route_distance"] / geodesic < 1 - TOLERANCE:
            raise ValueError("the native relay route is shorter than the direct great-circle arc")
        relay["great_circle_distance"] = geodesic

        control_length = relay_graph[route[0]][route[1]]["length"]
        control_wait = math.pi * control_length / (math.sqrt(2) * calibrated_speed)
        detuned = _pulse(control_length, speed * DETUNING, control_wait)[:, 0]
        detuned_probability = float(abs(detuned[1]) ** 2)
        if max_swap_error > TOLERANCE or detuned_probability >= 0.99:
            raise ValueError("the held-out seam swaps or the detuned control failed")
        levels.append({"level": level, "cell_centers": level_centers.tolist(),
                       "seam_lengths": seam_lengths, "held_out_seams": count,
                       "maximum_swap_probability_error": max_swap_error, "relay": relay,
                       "detuned_receiver_probability": detuned_probability})

    return {"radius": radius, "calibrated_speed": calibrated_speed,
            "reference_frequency": reference_frequency, "calibration": calibration,
            "relay_speed": math.sqrt(2) * calibrated_speed / math.pi,
            "base_beacons": list(beacons), "levels": levels,
            "uniform_hop_obstructions": [_hop_obstruction(a, b) for a, b in itertools.pairwise(graphs)]}
