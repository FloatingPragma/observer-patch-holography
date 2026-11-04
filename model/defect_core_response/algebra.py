"""Exact current checks on the original source links and full repair domain."""
import numpy as np
import sympy as s

from ..native_mass_channels.algebra import SIGMA
from ..defect_field_bridge import full_domain


def link_currents():
    nx = s.diag(1, 1, 0, 0)
    for sigma in SIGMA:
        link = s.I*sigma/2
        h = s.BlockMatrix([[s.zeros(2), link], [link.H, s.zeros(2)]]).as_explicit()
        current = s.I*s.BlockMatrix([[s.zeros(2), link], [-link.H, s.zeros(2)]]).as_explicit()
        if s.I*(h*nx-nx*h)+current != s.zeros(4):
            raise ArithmeticError('source link number-current orientation changed')
        if current.H != current or s.I*(h*(s.eye(4)-nx)-(s.eye(4)-nx)*h)-current != s.zeros(4):
            raise ArithmeticError('source current is not a Hermitian conserved bond current')
    # A full onsite CAR basis, with its centered quartic interaction.
    up, down = s.zeros(4), s.zeros(4)
    for mask in range(4):
        if mask & 1:
            up[mask ^ 1, mask] = 1
        if mask & 2:
            down[mask ^ 2, mask] = (-1)**(mask & 1)
    number = up.H*up+down.H*down
    contact = (up.H*up-s.eye(4)/2)*(down.H*down-s.eye(4)/2)
    if contact*number != number*contact:
        raise ArithmeticError('onsite contact changes the local number')
    for creator in (up.H, down.H):
        if number*creator-creator*number != creator:
            raise ArithmeticError('source creator has an incorrect relative charge')
    return {
        'original_three_Weyl_link_directions_checked': 3,
        'full_onsite_CAR_states': 4,
        'relative_addition_removal_charges': [1, -1],
        'number_current': 'i*(c_x^dagger*h_xy*c_y-c_y^dagger*h_yx*c_x)',
        'energy_current': 'i*[e_x,e_y], retaining every overlapping local-energy pair',
        'repair_energy_injection_may_be_nonzero': True,
        'Lorentz_stress_tensor_identified': False,
    }


def repair_currents(sites, law):
    dim = 1 << (2*sites)
    seam = (0, sites-1)
    rows = full_domain.finite_branches(sites=sites, flavors=2, law=law, seam=seam)
    occupations = [np.array([((mask >> (2*x)) & 3).bit_count() for mask in range(dim)])
                   for x in range(sites)]
    total = sum(occupations)
    completeness = np.zeros(dim, dtype=np.int64)
    active = 0
    for (_, raw), row in rows.items():
        current = 1 if raw in (1, 5) else -1 if raw in (2, 6) else 0
        for origin, output, sign in row:
            if (occupations[seam[0]][output]-occupations[seam[0]][origin] != -current
                    or occupations[seam[1]][output]-occupations[seam[1]][origin] != current
                    or total[output] != total[origin] or sign not in (-1, 1)):
                raise ArithmeticError('a complete native branch violates its recorded current')
            completeness[origin] += 1
            active += int(bool(current))
    if not np.all(completeness == 2) or not active:
        raise ArithmeticError('native coin branches or transported records were lost')
    # Actual asymmetry: a hole/single seam is active, a double/single holds.
    hole = full_domain.branch_action(1 << 2, law, 1, sites=2, flavors=2, seam=(0, 1))
    double = full_domain.branch_action(3 | (1 << 2), law, 0, sites=2, flavors=2, seam=(0, 1))
    if hole[3] == 0 or double[3] != 0 or double[0] != full_domain.COLLISION_RAW:
        raise ArithmeticError('native particle-hole countercontrol changed')
    return {'sites': sites, 'full_Fock_dimension': dim, 'law': law,
            'coin_raw_slots': len(rows), 'every_branch_number_and_current_checked': True,
            'collision_and_padded_slots_retained': True,
            'hole_moves_and_double_collision_holds': True,
            'Hamiltonian_particle_hole_symmetry_imposed_on_repair': False}


def certificate():
    return {'Hamiltonian_currents': link_currents(),
            'complete_native_record_currents': [
                repair_currents(sites, law) for sites in (2, 3, 4) for law in ('C', 'E')]}
