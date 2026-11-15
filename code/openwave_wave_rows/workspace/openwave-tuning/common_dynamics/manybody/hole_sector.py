"""Negative-current holes in the same declared fermion-mode dynamics.

The filled reference is a supplied state, with no ground-state or observed
antiparticle interpretation. No second species or Hamiltonian is introduced.
"""
from pathlib import Path
import hashlib
import json

import numpy as np

import manybody as mb

HERE = Path(__file__).resolve().parent


def run():
    paths = [Path(__file__), Path(mb.__file__), Path(mb.common.__file__), Path(mb.field.__file__)]
    pins = {str(p.relative_to(mb.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    g = mb.field.geometry()
    B, P = g['B'], g['P']
    source_permutation = mb.common.expected_swap_permutation()
    active = [16 * mb.common.encode_position(x) + event for x in range(96) for event in range(12)]
    active_rank = {mode: i for i, mode in enumerate(active)}
    permutation = [active_rank[int(source_permutation[mode])] for mode in active]
    assert sorted(permutation) == list(range(1152))
    swaps = sum(permutation[x] != x for x in range(1152)) // 2
    assert swaps == 288
    transport_determinant = (-1) ** swaps
    reflection_determinant = (-1) ** (11 * 96)
    assert transport_determinant == reflection_determinant == 1
    # Complement-basis phases c_j|filled>=(-1)^j|filled exceptj>.
    # Check the permutation's induced exterior-complement sign by counting
    # inversions once and subtracting those involving the omitted mode.
    all_inversions = sum(permutation[i] > permutation[j]
                         for i in range(1152) for j in range(i + 1, 1152))
    assert (-1) ** all_inversions == transport_determinant
    G = np.full(96, 12, dtype=np.int64)
    G[0] -= 1
    cases = 0
    neighbours, edges = mb.common.geometry_tables()
    for x in range(96):
        for event in range(12):
            mode = 16 * mb.common.encode_position(x) + event
            rank = active_rank[mode]
            after_rank = permutation[rank]
            omitted_inversions = (sum(permutation[k] > after_rank for k in range(rank))
                                  + sum(after_rank > permutation[k] for k in range(rank + 1, 1152)))
            physical_complement_sign = (-1) ** (all_inversions - omitted_inversions)
            logical_hole_sign = (-1) ** rank * physical_complement_sign * (-1) ** after_rank
            assert logical_hole_sign == transport_determinant
            after_mode = int(source_permutation[mode])
            y = mb.common.decode_position(after_mode // 16)
            n, n_after = np.full(96, 12, dtype=np.int64), np.full(96, 12, dtype=np.int64)
            n[x] -= 1
            n_after[y] -= 1
            current = np.zeros(288, dtype=np.int64)
            if event % 2:
                edge = int(edges[x, event // 2])
                a, b = g['edges'][edge]
                current[edge] = -(int(x == a) - int(y == a))
            assert np.array_equal(n + B @ (-P[:, x]), G)
            assert np.array_equal(n_after + B @ current, n)
            cases += 1
    result = {'schema': 'oph.local12.same_source_holes.v1', 'candidate_id': mb.MODE_ID,
              'checks': {'active_W_determinant_one': True,
                         'every_hole_transport_complement_phase': True,
                         'opposite_current_from_same_occupations': True,
                         'attached_Gauss_reference': True,
                         'source_pins_unchanged': True},
              'hole_basis_atoms': cases, 'filled_physical_modes': 1152,
              'source_formula': 'Gamma(W)c_j|filled>=sum_i conjugate(W_ij)c_i|filled>, because detW=1',
              'same_real_W': 'Ideal LOCAL12 R and S are real, so the one-hole matrix is the same real W on the physical sector.',
              'native_current': 'j_hole = -(delta_x(a)-delta_y(a)); obtained from actual occupied-mode difference, not an added negative species',
              'reference': {'G': '12 at every port except11 atroot0', 'matter_total': 1151,
                            'particular_field_for_hole_at_x': '-P_x'},
              'scope': 'Supplied completely filled physical-mode reference and one-hole sector of the same proposed fermion circuit. No energy-ground, stability, spin or observed-antiparticle claim; filled-state and affine-field preparation/reset are not replayed.',
              'source_sha256': pins, 'official_rows_changed': []}
    assert pins == {str(p.relative_to(mb.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE / 'HOLE_SECTOR.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'hole_basis_atoms', 'reference')}))


if __name__ == '__main__':
    run()
