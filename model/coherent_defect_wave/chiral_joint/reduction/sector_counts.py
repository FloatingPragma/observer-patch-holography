"""Count diagonal Gauss sectors of the uniform two-cell lepton control.

No Hamiltonian, weak-singlet projector, ground state or mass is computed.
Weights come from the same native, scalar and quantum-link current matrices
as the frozen chiral construction.
"""
import argparse
from collections import Counter, defaultdict
import json

import sympy as s

from ..model import active_terms
from ...chiral_higgs.model import higgs_memory
from ...chiral_source.model import link_generator, link_labels
from ...physical_attachment.chiral import representation


def _integer_weight(value, multiplier):
    result = s.simplify(multiplier * value)
    if result.is_Integer is not True:
        raise ValueError('current does not have the stated exact integer weights')
    return int(result)


def sector_counts():
    """Return every (L,2Sz) block in each admitted diagonal Gauss sector."""
    indices, _, _, active_links = active_terms()
    currents = representation()[2]
    _, scalar_currents, _, _ = higgs_memory()
    # The actual weak Cartan is the unique nonzero diagonal current on the
    # selected leptons, excluding the final native hypercharge generator.
    weak_cartans = [
        k for k, t in enumerate(currents[:-1])
        if t == s.diag(*t.diagonal())
        and t.extract(indices, indices) != s.zeros(len(indices))
    ]
    if len(weak_cartans) != 1:
        raise ValueError('the retained native lepton module has changed')
    weak_cartan = weak_cartans[0]
    y = [_integer_weight(currents[-1][i, i], 6) for i in indices] * 2
    z = [_integer_weight(currents[weak_cartan][i, i], 2) for i in indices] * 2
    scalar_y = [_integer_weight(q, 6) for q in scalar_currents[-1].diagonal()]
    scalar_z = [_integer_weight(q, 2) for q in scalar_currents[weak_cartan].diagonal()]

    # L=N_L-N_ec; both chiral spin components have these internal charges.
    # The one z-directed edge and spin-singlet Yukawa preserve 2Sz.
    lepton_number = (-1, 1, 1) * 2
    twice_spin = (1, 1, 1, -1, -1, -1)
    local = defaultdict(Counter)
    for mask in range(64):
        occupied = [j for j in range(6) if mask >> j & 1]
        fermion_y = sum(y[j] for j in occupied)
        fermion_z = sum(z[j] for j in occupied)
        ell = sum(lepton_number[j] for j in occupied)
        spin = sum(twice_spin[j] for j in occupied)
        for scalar in range(5):
            local[fermion_y + scalar_y[scalar], fermion_z + scalar_z[scalar]][ell, spin] += 1

    # Pull endpoint weights from the actual sixty-state quantum link before
    # selecting the blank, e^c and four weak-doublet active link states.
    all_links = link_labels()
    active_positions = [0] + [
        all_links.index((indices[a], indices[b])) for a, b in active_links[1:]
    ]
    link_y = [link_generator(currents[-1], endpoint) for endpoint in (0, 1)]
    link_z = [link_generator(currents[weak_cartan], endpoint) for endpoint in (0, 1)]
    sectors = []
    for boundary_charge in (0, -1):
        blocks = Counter()
        for position in active_positions:
            ly, ry = [_integer_weight(q[position, position], 6) for q in link_y]
            lz, rz = [_integer_weight(q[position, position], 2) for q in link_z]
            left = local[-ly, -lz]
            right = local[-6 * boundary_charge - ry, -rz]
            for (ell_left, spin_left), multiplicity_left in left.items():
                for (ell_right, spin_right), multiplicity_right in right.items():
                    blocks[ell_left + ell_right, spin_left + spin_right] += multiplicity_left * multiplicity_right
        largest = max(blocks.values())
        sectors.append({
            'fixed_right_boundary_hypercharge': boundary_charge,
            'active_source_dimension_before_Gauss_constraints': 64 ** 2 * len(active_positions) * 5 ** 2,
            'diagonal_Gauss_support_dimension': sum(blocks.values()),
            'conserved_block_count': len(blocks),
            'largest_diagonal_block_dimension': largest,
            'largest_block_labels': [
                {'L': ell, 'two_Sz': spin}
                for (ell, spin), count in sorted(blocks.items()) if count == largest
            ],
            'blocks': [
                {'L': ell, 'two_Sz': spin, 'dimension': count}
                for (ell, spin), count in sorted(blocks.items())
            ],
        })
    return {
        'control': 'two cells, one z-directed quantum link, both quantum lepton Yukawa terms',
        'parameters': {'left_y_e': '1', 'right_y_e': '1', 'y_u': '0', 'y_d': '0'},
        'native_internal_indices': list(indices),
        'fermion_modes': 12,
        'scalar_dimension_per_cell': 5,
        'active_quantum_link_dimension': len(active_positions),
        'fixed_boundary_charge_is_not_a_dynamical_species': True,
        'diagonal_constraints': ['Y_left=0', 'T3_left=0', 'Y_right+q_boundary=0', 'T3_right=0'],
        'block_labels': 'L=N_L-N_ec and twice the source spin projection along the edge',
        'sectors': sectors,
        'full_weak_Gauss_projection_computed': False,
        'Hamiltonian_or_ground_computed': False,
        'physical_particle_or_mass_identified': False,
        'next_projection': 'joint kernel of actual weak raising maps into adjacent weight sectors; do not compress their image to zero weight',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', action='store_true', help='print counts without the complete block table')
    args = parser.parse_args()
    result = sector_counts()
    if args.summary:
        for sector in result['sectors']:
            print('boundary charge {fixed_right_boundary_hypercharge}: full={active_source_dimension_before_Gauss_constraints}, '
                  'diagonal Gauss={diagonal_Gauss_support_dimension}, blocks={conserved_block_count}, '
                  'largest={largest_diagonal_block_dimension}'.format(**sector))
        print('Full weak-Gauss projection, Hamiltonian and ground remain uncomputed.')
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
