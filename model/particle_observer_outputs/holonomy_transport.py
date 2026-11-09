"""Gauge-invariant loop-defect creation, transport and complete records.

The link bits and local work register are declared additions to the native
icosahedral carrier. Existing Z6 topology and transport calculations supply
the incidence and the order-two subsector. This code does not identify its
curvature excitations with Standard Model particles or native memory charges.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import importlib
import json
from pathlib import Path
import sys

from model.common_source import slow_geometry
from model.observer import native, native_loop_readout


ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/common_source/slow_geometry.py',
    'model/observer/native_loop_readout.py',
    'code/a5_closure/native_defect_stability_gate.py',
    'code/a5_closure/flux_defect_criterion_certificate.py',
    'code/a5_closure/routed_seam_grammar_certificate.py',
    'model/common_source/native_gauge_object.py',
)


@lru_cache(maxsize=1)
def _upstream():
    location = ROOT/'code/a5_closure'
    sys.path.insert(0, str(location))
    try:
        gate = importlib.import_module('native_defect_stability_gate')
        flux = importlib.import_module('flux_defect_criterion_certificate')
    finally:
        sys.path.remove(str(location))
    for module in (gate, flux):
        if Path(module.__file__).resolve().parent != location.resolve():
            raise ValueError('a foreign finite holonomy module is loaded')
    return gate, flux


def _mask(value, width, label):
    if type(value) is not int or not 0 <= value < 1 << width:
        raise ValueError(label+' must be a complete binary address of the stated width')
    return value


def _indices(mask):
    return [i for i in range(mask.bit_length()) if (mask >> i) & 1]


@lru_cache(maxsize=1)
def complex_data():
    """Reuse the native cell complex and its established sphere homology."""
    gate, _ = _upstream()
    support, punctures = gate.support_and_punctures()
    edges = tuple(map(tuple, support['edges']))
    faces = tuple(map(tuple, support['face_rows']))
    face_masks = tuple(sum(1 << e for e in range(30) if support['boundary'][e][f])
                       for f in range(20))
    vertex_masks = tuple(sum(1 << e for e, ends in enumerate(edges) if v in ends)
                         for v in range(12))
    edge_faces = tuple(tuple(f for f, mask in enumerate(face_masks) if (mask >> e) & 1)
                       for e in range(30))
    if any(len(pair) != 2 for pair in edge_faces):
        raise ArithmeticError('the native edge has lost its two incident faces')
    if any((f & v).bit_count() % 2 for f in face_masks for v in vertex_masks):
        raise ArithmeticError('a vertex gauge operation changes face curvature')
    ranks = gate.homology_ranks(support, punctures)['F2']
    if (ranks['boundary_one_rank'], ranks['full_boundary_two_rank'],
            ranks['full_H1_dimension']) != (11, 19, 0):
        raise ArithmeticError('the native closed-sphere gauge quotient differs')
    return {'support': support, 'edges': edges, 'faces': faces,
            'face_masks': face_masks, 'vertex_masks': vertex_masks,
            'edge_faces': edge_faces, 'punctures': punctures, 'ranks': ranks}


def curvature(link_mask):
    """Twenty face holonomy bits; bit1 means Wilson product minus1."""
    a = _mask(link_mask, 30, 'link field')
    return sum(((a & row).bit_count() % 2) << f
               for f, row in enumerate(complex_data()['face_masks']))


def gauge_update(vertex_mask):
    result = 0
    for v in _indices(_mask(vertex_mask, 12, 'vertex gauge')):
        result ^= complex_data()['vertex_masks'][v]
    return result


def region_boundary(face_region):
    result = 0
    for f in _indices(_mask(face_region, 20, 'face region')):
        result ^= complex_data()['face_masks'][f]
    return result


def region_parity(link_mask, face_region):
    links = _mask(link_mask, 30, 'link field')
    region = _mask(face_region, 20, 'face region')
    inside = (curvature(links) & region).bit_count() % 2
    boundary = (links & region_boundary(region)).bit_count() % 2
    if inside != boundary:
        raise ArithmeticError('finite regional Stokes parity failed')
    return inside


def preserves_region(update_mask, face_region):
    return (_mask(update_mask, 30, 'link update') & region_boundary(face_region)).bit_count() % 2 == 0


def flat_gauge_potential(link_mask):
    """Recover the vertex gauge using the existing exact Z6 theorem."""
    links = _mask(link_mask, 30, 'link field')
    if curvature(links):
        raise ValueError('nonzero curvature is not a flat vertex gauge')
    gate, _ = _upstream()
    potential = gate.flat_gauge_potential(complex_data()['support'],
                                         {e: 3 for e in _indices(links)})
    if any(v not in (0, 3) for v in potential):
        raise ArithmeticError('the order-two subsector left its coefficient subgroup')
    result = sum((v == 3) << i for i, v in enumerate(potential))
    if gauge_update(result) != links:
        raise ArithmeticError('the recovered gauge does not reconstruct the flat link field')
    return result


def _state(face_mask, work):
    face_mask = _mask(face_mask, 20, 'curvature record')
    work = _mask(work, 1, 'work register')
    if face_mask.bit_count() % 2:
        raise ValueError('the closed sphere has no odd total-curvature sector')
    return face_mask, work


def edge_syndrome(edge):
    if type(edge) is not int or not 0 <= edge < 30:
        raise ValueError('an actual carrier edge index is required')
    return sum(1 << f for f in complex_data()['edge_faces'][edge])


def work_exchange(face_mask, work, edge):
    """Complete reversible pair creation/annihilation; hold every other input.

    The work quantum costs 2J and each face defect costs J. This permutation
    preserves their total energy on every admitted curvature/work state.
    """
    faces, work = _state(face_mask, work)
    pair = edge_syndrome(edge)
    occupied = (faces & pair).bit_count()
    if (occupied, work) in ((0, 1), (2, 0)):
        return faces ^ pair, 1-work
    return faces, work


def move_or_hold(face_mask, work, edge):
    """Move exactly one adjacent defect; keep all zero/two-defect inputs."""
    faces, work = _state(face_mask, work)
    pair = edge_syndrome(edge)
    return (faces ^ pair if (faces & pair).bit_count() == 1 else faces), work


def energy_units(face_mask, work):
    """Total energy divided by the declared positive defect cost J."""
    faces, work = _state(face_mask, work)
    return faces.bit_count()+2*work


def creation_pulse(face_mask, work, edge, cosine=F(3, 5), sine=F(4, 5)):
    """Exact Gaussian-rational columns of a full energy-preserving unitary.

    On an active work-exchange pair, U=cI-isW; it is identity on holds.
    The returned amplitudes retain every nonzero output. All other labels
    in the specified full Hilbert space have exactly zero amplitudes.
    """
    if type(cosine) not in (int, F) or type(sine) not in (int, F):
        raise ValueError('exact rational pulse coordinates are required')
    c, a = F(cosine), F(sine)
    if c*c+a*a != 1:
        raise ValueError('the pulse coordinates must lie on the unit circle')
    initial = _state(face_mask, work)
    final = work_exchange(*initial, edge)
    if final == initial:
        return {initial: (F(1), F(0))}
    result = {initial: (c, F(0)), final: (F(0), -a)}
    return {k: v for k, v in result.items() if v != (0, 0)}


def _snapshot(label, links, work, edge=None, region=1):
    faces = curvature(links)
    return {'operation': label, 'flipped_edge': edge,
            'negative_link_edges': _indices(links),
            'link_string_is_a_gauge_representative': True,
            'defect_faces': _indices(faces),
            'all_face_holonomy_signs': [1-2*((faces >> f) & 1) for f in range(20)],
            'work_register': work, 'defect_energy_in_J': faces.bit_count(),
            'work_energy_in_J': 2*work, 'total_energy_in_J': energy_units(faces, work),
            'protected_region_faces': _indices(region),
            'boundary_holonomy_sign': 1-2*region_parity(links, region),
            'update_preserves_region': None if edge is None else preserves_region(1 << edge, region)}


def graphical_path():
    """Literal source-labelled creation, separation, return and annihilation."""
    data = complex_data(); gate, flux = _upstream()
    start, end = data['punctures']
    path = flux.dual_path(data['support'], start, end)
    gate.verify_local_link_path(data['support'], (start, end), path, 3)
    first_edge = path[0][1]; links, work = 0, 1
    region = 1 << start
    rows = [_snapshot('flat, work available', links, work, region=region)]
    expected = work_exchange(0, work, first_edge)
    links ^= 1 << first_edge; work = expected[1]
    if curvature(links) != expected[0]:
        raise ArithmeticError('creation did not match the actual face incidence')
    rows.append(_snapshot('create adjacent pair', links, work, first_edge, region))
    for _, edge in path[1:]:
        expected = move_or_hold(curvature(links), work, edge)
        links ^= 1 << edge
        if curvature(links) != expected[0]:
            raise ArithmeticError('transport failed to move just one endpoint')
        rows.append(_snapshot('move one endpoint', links, work, edge, region))
    for _, edge in reversed(path[1:]):
        expected = move_or_hold(curvature(links), work, edge)
        links ^= 1 << edge
        if curvature(links) != expected[0]:
            raise ArithmeticError('return failed to retain its other endpoint')
        rows.append(_snapshot('return one endpoint', links, work, edge, region))
    expected = work_exchange(curvature(links), work, first_edge)
    links ^= 1 << first_edge; work = expected[1]
    if links != 0 or expected != (0, 1):
        raise ArithmeticError('annihilation failed to return the work quantum')
    rows.append(_snapshot('annihilate pair and restore work', links, work, first_edge, region))
    return {'start_face': start, 'separated_face': end,
            'dual_path': [{'face': f, 'edge': e} for f, e in path], 'stages': rows}


def complete_history():
    """One coherent birth pulse, full readout, record-controlled transport.

    Each read measures the entire even-curvature and work address. Absent
    labels have zero CP maps on the specified input, not deleted outcomes.
    The input's two nonzero histories are retained without postselection.
    """
    illustration = graphical_path()
    path = illustration['dual_path']; edge = path[0]['edge']
    pulse = creation_pulse(0, 1, edge)
    histories = []
    for state, amplitude in sorted(pulse.items()):
        born = state == (edge_syndrome(edge), 0)
        records = [state]
        for step in path[1:]:
            state = move_or_hold(*state, step['edge']) if born else state
            records.append(state)
        for step in reversed(path[1:]):
            state = move_or_hold(*state, step['edge']) if born else state
            records.append(state)
        state = work_exchange(*state, edge) if born else state
        records.append(state)
        probability = amplitude[0]**2+amplitude[1]**2
        if state != (0, 1) or any(energy_units(*r) != 2 for r in records):
            raise ArithmeticError('a complete history lost its field/work accounting')
        histories.append({'birth_record': born, 'probability': str(probability),
                          'records': [{'defect_faces': _indices(f), 'work_register': b}
                                      for f, b in records]})
    if sum((F(h['probability']) for h in histories), F(0)) != 1:
        raise ArithmeticError('the complete experiment lost an outcome')
    return {'initial_state': 'flat gauge-orbit state, work=1',
            'pulse_cosine': '3/5', 'pulse_sine': '4/5',
            'face_and_work_readout_alphabet_size': 1 << 21,
            'odd_face_parity_pointer_labels_with_zero_projector': 1 << 20,
            'additional_outside_Gauss_outcome': True,
            'full_readout_alphabet_including_outside': (1 << 21)+1,
            'other_first_records_including_outside_select_hold': True,
            'readout_count': len(histories[0]['records']),
            'nonzero_complete_histories': histories,
            'all_other_histories_have_exact_zero_weight': True,
            'outside_Gauss_sector_outcome_retained_for_ambient_inputs': True,
            'postselection_used': False, 'arbitrary_untouched_references_admitted': True,
            'record_control_is_declared_apparatus': True}


def certificate():
    data = complex_data()
    coordinates = slow_geometry.raw_frame()
    # This independently reruns the actual native current readout on its
    # full finite load domain. Its contrast is a vertex gradient, not a link bit.
    native_checks = 0
    for loads in native.states():
        contrasts = native_loop_readout.read_contrasts(loads, data['edges'])
        if any(contrasts[a, b] != loads[a]-loads[b] for a, b in data['edges']):
            raise ArithmeticError('the native current/gradient identity changed')
        if any(native_loop_readout.integer_holonomy(contrasts, face+(face[0],))
               for face in data['faces']):
            raise ArithmeticError('native current records produced face holonomy')
        native_checks += 60
    return {'schema': 'oph.holonomy_transport.v1', 'source_inputs': list(SOURCE_INPUTS),
        'coefficient_choice': 'Z2 order-two subgroup {0,3} of the retained Z6 link construction',
        'carrier': {'vertices': [{'id': i,
                     'coordinates_exact_Qsqrt5': [[str(a), str(b)] for a, b in v],
                     'coordinates_display': [float(a)+float(b)*5**0.5 for a, b in v]}
                    for i, v in enumerate(coordinates)],
                    'edges': [list(e) for e in data['edges']],
                    'faces': [list(f) for f in data['faces']],
                    'edge_incident_faces': [list(fs) for fs in data['edge_faces']],
                    'Euler_characteristic': 2, 'homology_F2': data['ranks']},
        'complete_spaces': {'link_qubits': 30, 'ambient_link_dimension': 1 << 30,
                            'independent_vertex_Gauss_constraints': 11,
                            'gauge_invariant_dimension': 1 << 19,
                            'work_dimension': 2, 'face_curvature_rank': 19,
                            'global_face_holonomy_product': 1,
                            'nontrivial_flat_sphere_holonomy': False},
        'protection': {'exact_condition': 'every allowed update u has even intersection with boundary(R)',
                       'preserved_quantity': 'product of face holonomies in R equals boundary Wilson product',
                       'odd_region_cannot_be_repaired_to_all_flat_under_this_condition': True,
                       'arbitrary_local_link_flips_permanently_protect_defects': False,
                       'gauge_invariance_alone_protects_defects': False},
        'native_boundary': {'actual_current_branches_checked': native_checks,
                            'actual_native_current_face_holonomy': 0,
                            'link_bits_and_work_register_added_to_this_program': True,
                            'unextended_C_E_creates_these_spatial_defects': False,
                            'protected_memory_charge_identified_with_curvature': False},
        'graphical_path': graphical_path(), 'complete_experiment': complete_history(),
        'energy_convention': 'J>0 is a declared defect cost; no mass or GeV conversion',
        'clock_convention': 'declared control/readout stages; no relativistic propagation clock claimed',
        'physical_particle_assignments': [], 'measured_particle_inputs': []}


if __name__ == '__main__':
    print(json.dumps(certificate(), ensure_ascii=False, indent=2, sort_keys=True))
