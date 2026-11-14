"""A finite central-link register with protection relative to a fixed boundary.

This classical link-memory control is separate from the static-background
coherent walk. It tests a declared move grammar; it supplies no particle mass.
"""
from collections import deque
import argparse
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path

from model.common_source import slow_geometry

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('protection_results.json')


def rank_mod2(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


class BoundaryFlux:
    """Thirty Z2 link bits and a declared region of the actual carrier faces."""

    def __init__(self, region_vertex=5):
        self.edges = tuple(tuple(sorted(e)) for e in slow_geometry.source().seams())
        index = {e: i for i, e in enumerate(self.edges)}
        self.faces = tuple(f for f in combinations(range(12), 3)
                           if all(tuple(sorted(e)) in index for e in combinations(f, 2)))
        self.face_masks = tuple(sum(1 << index[e] for e in combinations(f, 2)) for f in self.faces)
        if len(self.edges) != 30 or len(self.faces) != 20:
            raise ValueError('the link control requires the pinned icosahedral carrier')
        self.incident_faces = tuple(tuple(i for i, mask in enumerate(self.face_masks) if mask >> e & 1)
                                    for e in range(30))
        if any(len(fs) != 2 for fs in self.incident_faces):
            raise ValueError('each closed-screen link must border two faces')
        if type(region_vertex) is not int or region_vertex not in range(12):
            raise ValueError('region vertex outside the carrier')
        self.region = tuple(i for i, f in enumerate(self.faces) if region_vertex in f)
        self.region_mask = sum(1 << f for f in self.region)
        boundary = 0
        for f in self.region:
            boundary ^= self.face_masks[f]
        self.boundary = tuple(e for e in range(30) if boundary >> e & 1)
        self.boundary_mask = boundary
        self.interior = tuple(e for e, fs in enumerate(self.incident_faces) if all(f in self.region for f in fs))
        self.vertex_stars = tuple(sum(1 << e for e, ends in enumerate(self.edges) if v in ends)
                                  for v in range(12))

    @staticmethod
    def validate(state):
        if type(state) is not int or not 0 <= state < (1 << 30):
            raise ValueError('a link state is exactly thirty bits')

    def syndrome(self, state):
        self.validate(state)
        return sum(((state & mask).bit_count() & 1) << f for f, mask in enumerate(self.face_masks))

    def region_charge(self, state):
        """0 means even flux; 1 means odd flux, equivalently boundary sign minus."""
        return (self.syndrome(state) & self.region_mask).bit_count() & 1

    def flip(self, state, edge):
        """A declared interior link change; boundary and exterior bits are fixed."""
        self.validate(state)
        if type(edge) is not int or edge not in self.interior:
            raise ValueError('the protected grammar permits only interior edges')
        return state ^ (1 << edge)

    def transport(self, state, edge):
        """Move one adjacent defect; retain a no-operation record otherwise.

        This local rule is involutory: a move leaves exactly one adjacent defect.
        It is a specified finite repair event, not energy minimization.
        """
        self.validate(state)
        if type(edge) is not int or edge not in self.interior:
            raise ValueError('the transport grammar fixes boundary and exterior links')
        before = self.syndrome(state)
        pair = self.incident_faces[edge]
        count = sum(before >> f & 1 for f in pair)
        after_state = state ^ (1 << edge) if count == 1 else state
        return after_state, {'edge': list(self.edges[edge]), 'moved': count == 1,
                             'before_face_syndrome': before,
                             'after_face_syndrome': self.syndrome(after_state)}


def certificate():
    law = BoundaryFlux()
    # These identities on linear generators prove them for every 30-bit word.
    assert rank_mod2(law.vertex_stars) == 11
    assert rank_mod2(law.face_masks) == 19
    assert all(law.syndrome(star) == 0 for star in law.vertex_stars)
    assert all(law.syndrome(1 << e).bit_count() == 2 for e in range(30))
    for e in range(30):
        basis = 1 << e
        assert law.syndrome(basis).bit_count() % 2 == 0
        assert law.region_charge(basis) == ((basis & law.boundary_mask).bit_count() & 1)
    assert all(not ((1 << e) & law.boundary_mask) for e in law.interior)

    central_edge = law.edges.index((0, 1))
    initial = 1 << central_edge
    assert central_edge in law.boundary
    assert len(law.region) == len(law.interior) == len(law.boundary) == 5
    assert law.syndrome(initial).bit_count() == 2 and law.region_charge(initial) == 1
    outside = law.syndrome(initial) & ~law.region_mask
    all_states = []
    for bits in range(1 << len(law.interior)):
        state = initial
        for k, edge in enumerate(law.interior):
            if bits >> k & 1:
                state = law.flip(state, edge)
        all_states.append(state)
        assert state & law.boundary_mask == initial & law.boundary_mask
        assert law.region_charge(state) == 1
        assert law.syndrome(state) & ~law.region_mask == outside
    assert len(set(law.syndrome(state) & law.region_mask for state in all_states)) == 16

    reached, queue = {initial: []}, deque([initial])
    move_checks = hold_checks = 0
    while queue:
        state = queue.popleft()
        for edge in law.interior:
            output, record = law.transport(state, edge)
            assert law.transport(output, edge)[0] == state
            assert law.region_charge(output) == 1
            assert (law.syndrome(output) & law.region_mask).bit_count() == 1
            assert law.syndrome(output) & ~law.region_mask == outside
            move_checks += int(record['moved'])
            hold_checks += int(not record['moved'])
            if output not in reached:
                reached[output] = reached[state] + [list(law.edges[edge])]
                queue.append(output)
    assert len(reached) == 10
    destinations = {next(f for f in law.region if law.syndrome(state) >> f & 1) for state in reached}
    assert destinations == set(law.region)

    # A fixed finite scheduler executes the actual event grammar including holds.
    state, tick, visited, schedule = initial, 0, set(), []
    while (state, tick % len(law.interior)) not in visited:
        visited.add((state, tick % len(law.interior)))
        edge = law.interior[tick % len(law.interior)]
        state, record = law.transport(state, edge)
        schedule.append({'tick': tick, **record})
        tick += 1
    assert len(schedule) <= len(reached) * len(law.interior)
    assert law.region_charge(state) == 1
    assert state == initial and tick % len(law.interior) == 0
    # Permitting the selected boundary flip erases the initial pair exactly.
    assert initial ^ (1 << central_edge) == 0
    signed_trial = sum((1 - (-1 if e == central_edge else 1))**2 for e in range(30)) / Fraction(12)
    assert signed_trial == Fraction(1, 3)

    inputs = ('model/coherent_seam_repair/protection.py', 'model/common_source/slow_geometry.py',
              'evidence/native_geometric_source_20260925/spec.json',
              'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
              'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py')
    return {
        'schema': 'oph.coherent_seam_repair.boundary_flux.v1',
        'input_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs},
        'source_domain': 'declared finite classical central-link register and fixed-boundary move grammar',
        'link_memory_bits': 30, 'face_count': 20, 'vertex_count': 12,
        'region_faces': [list(law.faces[f]) for f in law.region],
        'boundary_edges': [list(law.edges[e]) for e in law.boundary],
        'interior_edges': [list(law.edges[e]) for e in law.interior],
        'initial_negative_edge': [0, 1],
        'initial_flux_faces': [list(law.faces[f]) for f in range(20) if law.syndrome(initial) >> f & 1],
        'mod2_vertex_boundary_rank': 11, 'mod2_face_boundary_rank': 19,
        'filled_sphere_first_homology_dimension': 30 - 11 - 19,
        'whole_sphere_face_flux_parity': 0,
        'region_boundary_identity_checked_on_all_link_generators': True,
        'reachable_interior_flip_states': len(all_states),
        'distinct_odd_region_syndromes': 16,
        'single_defect_transport_link_states': len(reached),
        'single_defect_positions': len(destinations),
        'involutory_transport_move_checks': move_checks,
        'retained_no_operation_checks': hold_checks,
        'round_robin_execution': schedule,
        'protection': {'fixed_boundary_odd_flux_cannot_reach_vacuum': True,
                       'arbitrary_link_updates_protect_initial_pair': False,
                       'vertex_recharting_preserves_all_face_fluxes': True,
                       'coupled_static_wave_register_is_this_dynamic_link_program': False},
        'finite_gap_bound': {'single_negative_link': 'lambda_min(L_signed) <= 4/N',
                             'actual_carrier_trial_quotient': str(signed_trial),
                             'clock_condition': 'fixed local spectral normalization during volume growth',
                             'uniform_infinite_volume_mass_lower_bound_supplied': False},
        'scope': {'physical_particle_identification': False, 'rest_mass': False,
                  'dynamical_gauge_backreaction_with_coherent_walk': False,
                  'source_selected_boundary_or_move_grammar': False},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    text = json.dumps(certificate(), sort_keys=True, indent=2) + '\n'
    if args.write:
        OUTPUT.write_text(text)
    elif not args.check.exists() or args.check.read_text() != text:
        raise ValueError('stale boundary-flux certificate')
    print('PASS finite boundary-charge protection, complete transport and mass boundary')


if __name__ == '__main__':
    main()
