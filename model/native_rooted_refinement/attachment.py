"""Bind finite rooted windows, canonical states and the complete record bound."""
from fractions import Fraction as F
import hashlib
import json
import numpy as np

from model.native_program import interpreter as native
from model.native_local_history_limit import limit
from . import refinement as r, routing


def source_prefix(level, radius=3):
    r.integer(level, upper=3); r.integer(radius)
    local = routing.LocalSource()
    mesh = local.ico.build_geodesic_icosahedral_tower(level).levels[level]
    raw = routing.source._source().build_federation(level, 'port_pair')
    cells = local.ball(level, radius)
    source_routes = {c: {} for c in cells}
    for a, p, b, q in raw.inter_pairs:
        if int(a) in source_routes:
            source_routes[int(a)][int(b)] = int(p)
        if int(b) in source_routes:
            source_routes[int(b)][int(a)] = int(q)
    for c in cells:
        points = np.asarray([local.point(level, v) for v in local.cell(level, c)])
        if not np.array_equal(points, mesh.vertices[mesh.faces[c]]):
            raise ArithmeticError('local midpoint expressions differ from literal source arrays')
        if dict(local.route(level, c)[0]) != source_routes[c]:
            raise ArithmeticError('certified local routing differs from source output')
        for k in range(3):
            other, side = local.neighbour(level, c, k)
            if local.neighbour(level, other, side) != (c, k):
                raise ArithmeticError('source neighbour reciprocity failed')
    return {'level': level, 'cell_radius': radius, 'cells_checked': len(cells),
            'midpoint_coordinates_bit_identical': True,
            'complete_routed_neighbour_rows_identical': True,
            'global_geometry_evaluated': True,
            'expanded_native_geometry_supported': level <= native.MAX_EVALUATED_LEVEL}


def occupation(level, density):
    r.integer(level)
    if type(density) is not F or not 0 <= density <= 1:
        raise ValueError('an exact declared density in [0,1] required')
    population = 240*4**level
    return population, (population*density.numerator)//density.denominator


def history_budget(window, density, partition):
    """Conditional full local records on this authenticated finite geometry.

    The initial state is the declared fair canonical number sector with an
    independent reference. Observe records only on seams with both endpoints
    in the twelve-port root, and the final state on those ports. Other window
    records are traced. The target is the same finite Poisson observation.
    """
    routing.authenticate_window(window)
    level = window['level']; population, occupied = occupation(level, density)
    seams = 630*4**level
    budget = limit.local_budget(population, occupied, seams,
        len(window['canonical_ports']), len(window['internal_seams']), 12, 6,
        window['cell_radius'], partition, density)
    conditional = budget.pop('native_to_matching_infinite_local_process_upper')
    budget['infinite_process_upper_if_all_radius_routing_is_supplied'] = conditional
    # An arbitrary port path of length <=radius crosses at most radius
    # cells, so the whole-cell ball contains the required port ball.
    attempts = limit.aligned_clock(seams, partition)[-1]
    program = native.program(federation={'level': level, 'total': occupied},
        biased=False, commands=[{'op': 'repair', 'selector': {'kind': 'uniform'},
                               'law': 'mixture', 'repeat': attempts}],
        capacity=1, overflow='export_oldest')
    native.admit(program)
    root = set(window['root_ports'])
    observed = tuple(e for e in window['internal_seams'] if set(e) <= root)
    return {'level': level, 'population': population, 'occupied': occupied,
            'declared_density': density, 'seams': seams,
            'clock_partition': partition, 'attempts': attempts,
            'program_digest': native.digest(program), 'budget': budget,
            'observed_root_internal_seams': observed,
            'final_observed_root_ports': window['root_ports'],
            'other_window_record_values_traced': True,
            'arbitrary_correlated_exterior_reference_asserted': False,
            'expanded_complete_history_execution_claimed': False}


def host_and_arithmetic_controls():
    ico, _ = routing.modules()
    value = native.program(federation={'level': 200, 'total': 0}, biased=False,
                           commands=[{'op': 'hold', 'repeat': 1}])
    sentinel = {'unevaluated': True}
    history = native.execute(value, sentinel)
    if history.status != 'RESOURCE_STOP_UNEVALUATED_CARRIER' or history.attempted_steps != 0:
        raise ArithmeticError('unsupported geometry must stop before its constructor')
    # This tiny call exercises the literal hash operation, without allocating
    # a high-level graph or pretending that the host reached this operation.
    overflow = False
    try:
        ico._geometry_hash(2**32, np.zeros((1, 3)),
                           np.zeros((0, 3), dtype=np.int64),
                           np.zeros((0, 2), dtype=np.int64), None)
    except OverflowError:
        overflow = True
    if not overflow:
        raise ArithmeticError('the four-byte source level tag changed')
    # Literal child0 arithmetic, with inherited vertices re-cleaned at each
    # level. Vertex identities are kept even when coordinates coincide.
    base, faces = ico._base_icosahedron()
    triangle = ico._clean_coordinates(base)[faces[0]]
    collapse = None
    for level in range(1, 65):
        a, b, c = triangle
        ab = a+b; ab /= np.linalg.norm(ab)
        ca = c+a; ca /= np.linalg.norm(ca)
        triangle = ico._clean_coordinates(ico._clean_coordinates(np.asarray([a, ab, ca])))
        orient = ico._orient_faces_outward(triangle, np.asarray([[0, 1, 2]]))[0]
        triangle = triangle[orient]
        area = ico._spherical_face_areas(triangle, np.asarray([[0, 1, 2]]))[0]
        if area == 0 and collapse is None:
            collapse = level
    if collapse is None:
        raise ArithmeticError('the bounded source corner-collapse control changed')
    return {'admitted_level200_host_status': history.status,
            'attempts_before_stop': history.attempted_steps,
            'geometry_constructor_called_for_host_control': False,
            'four_byte_level_tag_overflow': overflow,
            'corner_zero_spherical_area_within64_steps': True,
            'corner_arithmetic_steps_checked': 64,
            'astronomic_global_graph_execution_claimed': False,
            'exact_arithmetic_continuation_selected_by_this_package': False}


def window_digest(window):
    routing.authenticate_window(window)
    keys = ('cell_radius', 'canonical_ports', 'internal_seams', 'root_ports', 'routing')
    raw = json.dumps(routing.typed({k: window[k] for k in keys}), sort_keys=True).encode()
    return hashlib.sha256(raw).hexdigest()
