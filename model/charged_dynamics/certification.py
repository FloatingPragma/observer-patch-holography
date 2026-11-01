"""Exact finite certificates for the complete joined-chain parameter interval.

Floating eigensolves propose dyadic brackets only. Integer symmetric inertia
rechecks every bracket; a rational Lipschitz bound covers the gaps between
grid points. No floating spectral scan is used as an isolation proof.
"""
from __future__ import annotations

from fractions import Fraction as F
import hashlib
import json
import math

import numpy as np

from . import chain


GRID = 64
DENOMINATOR = 65536
LINK_NORM = F(1, 8)
CHARGED_GAP = F(1, 256)
GROUND_GAP = F(1, 128)


def structural_certificate():
    """Exact full-carrier spin/particle-hole identities used by the reduction.

    Every N sector is assembled. SU(2) then puts each odd-N multiplet in
    Sz=1/2, and every non-singlet even-N multiplet in Sz=1 as well as Sz=0.
    The chemical-potential-centered contact term and bipartite hopping are
    particle-hole invariant. These facts justify the small inertia blocks.
    """
    base = np.zeros((256, 256), dtype=np.int64)
    link = np.zeros_like(base)
    covered = []
    for up in range(5):
        for down in range(5):
            addresses, first, second = chain.hubbard_integer_parts(up, down)
            if any(type(value) is not int or abs(value) > 32
                   for matrix in (first, second) for row in matrix for value in row):
                raise ArithmeticError('bounded literal Hubbard integer entries required')
            covered.extend(addresses)
            base[np.ix_(addresses, addresses)] = first
            link[np.ix_(addresses, addresses)] = second
    if sorted(covered) != list(range(256)):
        raise ArithmeticError('complete Fock sector allocation required')
    raising = np.zeros_like(base)
    twice_z = np.zeros(256, dtype=np.int64)
    conjugation = np.zeros_like(base)
    spin_links = [np.zeros_like(base), np.zeros_like(base)]
    for address in range(256):
        twice_z[address] = sum(((address >> (2*j)) & 1)-((address >> (2*j+1)) & 1) for j in range(4))
        for site in range(4):
            if (address >> (2*site)) & 3 == 2:
                raising[address ^ (3 << (2*site)), address] = 1
        complement = 255 ^ address
        phase = (-1)**sum((address & ((1 << j)-1)).bit_count() for j in range(8))
        phase *= math.prod(-(-1)**(j//2) for j in range(8) if complement >> j & 1)
        conjugation[complement, address] = phase
        for spin in (0, 1):
            a, b = 2+spin, 4+spin
            if ((address >> a) & 1) != ((address >> b) & 1):
                spin_links[spin][address ^ (1 << a) ^ (1 << b), address] = (-1)**sum(
                    address >> j & 1 for j in range(a+1, b))
    identity = np.eye(256, dtype=np.int64)
    if not np.array_equal(raising@raising.T-raising.T@raising, np.diag(twice_z)):
        raise ArithmeticError('complete CAR spin representation failed')
    if not np.array_equal(conjugation.T@conjugation, identity):
        raise ArithmeticError('complete particle-hole unitary failed')
    number = np.diag([j.bit_count() for j in range(256)])
    if not np.array_equal(conjugation@number@conjugation.T, 8*identity-number):
        raise ArithmeticError('particle-hole number complement failed')
    for matrix in (base, link):
        if (not np.array_equal(matrix@raising, raising@matrix)
                or not np.array_equal(conjugation@matrix@conjugation.T, matrix)):
            raise ArithmeticError('full spin or particle-hole generator identity failed')
    first, second = spin_links
    if (not np.array_equal(link, first+second) or not np.array_equal(first@second, second@first)
            or any(not np.array_equal(hop@hop@hop, hop) for hop in spin_links)):
        raise ArithmeticError('exact middle-link norm decomposition failed')
    _, one_base, one_link = chain.hubbard_integer_parts(1, 0)
    hopping = np.array(one_base, dtype=np.int64)+np.array(one_link, dtype=np.int64)+8*np.eye(4, dtype=np.int64)
    if (max(np.abs(hopping).sum(axis=1)) > 2
            or base[85, 85] != -32 or link[85, 85] != 0):
        raise ArithmeticError('global number-sector lower bound or variational trial failed')
    # The monomial gauge uses only 0,+/-1,+/-i. All scaled matrix entries
    # below are small Gaussian integers, exactly represented in these arrays.
    # Equality at both affine endpoints binds the whole link interval.
    for axis in range(3):
        transform = chain.gauge(axis)
        local = chain._local(axis).generator
        separate = 16*(np.kron(np.eye(16), local)+np.kron(local, np.eye(16)))
        joined = 16*(chain.FullFock(8).hamiltonian(chain.geometry(axis)).toarray()-number/2)
        if (not np.array_equal(transform.conj().T@transform, identity)
                or not np.array_equal(transform@base@transform.conj().T, separate)
                or not np.array_equal(transform@(base+link)@transform.conj().T, joined)):
            raise ArithmeticError('exact gauge does not bind the existing complete Weyl/contact action')
    return {'full_Fock_columns': 256, 'all_number_component_blocks': 25,
        'exact_monomial_source_gauge_axes': 3,
        'integer_SU2_commutators_checked': True,
        'particle_hole_full_generator_and_number_identity_checked': True,
        'middle_link_commuting_spin_hops_each_norm_at_most': '1/16',
        'middle_link_norm_upper': '1/8',
        'one_particle_hopping_norm_upper_all_links': '1/8',
        'nonnegative_contact_term_used': True,
        'N_at_most_3_generator_lower_bound': '-15/8',
        'number_four_occupation_trial_variational_upper_bound': '-2',
        'product_dimer_ground_variational_upper_bound': '-2',
        'other_number_sectors_follow_by_particle_hole': True}


def _integer(value):
    return type(value) is int


def inertia(matrix):
    """Exact symmetric fraction-free congruence, with 1x1 and 2x2 pivots.

    Returns negative, zero and positive eigenvalue counts. A zero diagonal
    does not imply a kernel: an off-diagonal pivot contributes one of each
    sign. Every Bareiss division is checked exactly.
    """
    if (type(matrix) not in (tuple, list) or not matrix or len(matrix) > 70
            or any(type(row) not in (tuple, list) or len(row) != len(matrix)
                   or any(not _integer(value) for value in row) for row in matrix)
            or any(matrix[i][j] != matrix[j][i] for i in range(len(matrix)) for j in range(i))):
        raise ValueError('bounded symmetric Python-integer matrix required')
    value = [list(row) for row in matrix]
    size, position, previous = len(value), 0, 1
    negative = positive = 0

    def divide(numerator, denominator):
        quotient, remainder = divmod(numerator, denominator)
        if remainder:
            raise ArithmeticError('fraction-free inertia lost exact divisibility')
        return quotient

    def swap(a, b):
        value[a], value[b] = value[b], value[a]
        for row in value:
            row[a], row[b] = row[b], row[a]

    while position < size:
        selected = next((j for j in range(position, size) if value[j][j]), None)
        if selected is not None:
            swap(position, selected)
            pivot = value[position][position]
            negative += pivot*previous < 0
            positive += pivot*previous > 0
            for i in range(position+1, size):
                for j in range(i, size):
                    entry = divide(pivot*value[i][j]-value[i][position]*value[position][j], previous)
                    value[i][j] = value[j][i] = entry
            previous = pivot
            position += 1
            continue
        pair = next(((i, j) for i in range(position, size) for j in range(i+1, size)
                     if value[i][j]), None)
        if pair is None:
            return int(negative), size-position, int(positive)
        swap(position, pair[0])
        swap(position+1, pair[1])
        off = value[position][position+1]
        for i in range(position+2, size):
            for j in range(i, size):
                entry = divide(-off*off*value[i][j]+off*(value[i][position]*value[position+1][j]
                    +value[i][position+1]*value[position][j]), previous*previous)
                value[i][j] = value[j][i] = entry
        previous = divide(-off*off, previous)
        negative += 1; positive += 1; position += 2
    return int(negative), 0, int(positive)


def shifted_integer_matrix(up, down, link_strength, threshold):
    link_strength = chain._link(link_strength)
    if type(threshold) not in (int, F):
        raise ValueError('exact spectral threshold required')
    numerator, denominator = (threshold, 1) if type(threshold) is int else (threshold.numerator, threshold.denominator)
    if type(numerator) is not int or type(denominator) is not int or denominator <= 0:
        raise ValueError('canonical rational threshold components required')
    threshold = F(numerator, denominator)
    _, base, link = chain.hubbard_integer_parts(up, down)
    scale = math.lcm(link_strength.denominator, threshold.denominator)
    coefficient = link_strength.numerator*(scale//link_strength.denominator)
    shift = 16*threshold.numerator*(scale//threshold.denominator)
    return tuple(tuple(scale*a+coefficient*b-(shift if i == j else 0)
        for j, (a, b) in enumerate(zip(row, link[i], strict=True))) for i, row in enumerate(base))


def count_below(up, down, link_strength, threshold):
    return inertia(shifted_integer_matrix(up, down, link_strength, threshold))


def _floor_eigenvalue(up, down, coupling, index, hint):
    """Find the exact dyadic floor; the floating hint has no proof authority."""
    proposed = math.floor(float(hint)*DENOMINATOR)
    while count_below(up, down, coupling, F(proposed, DENOMINATOR))[0] > index:
        proposed -= 1
    while count_below(up, down, coupling, F(proposed+1, DENOMINATOR))[0] <= index:
        proposed += 1
    return proposed


def build_grid():
    rows = []
    for point in range(GRID+1):
        coupling = F(point, GRID)
        bounds = []
        for up, down, first, second in ((2, 1, 1, 2), (2, 2, 0, 1)):
            _, base, link = chain.hubbard_integer_parts(up, down)
            values = np.linalg.eigvalsh((np.array(base)+float(coupling)*np.array(link))/16)
            lower = _floor_eigenvalue(up, down, coupling, first, values[first])
            upper = _floor_eigenvalue(up, down, coupling, second, values[second])
            bounds.extend((lower+1, upper))
        rows.append([point, *bounds])
    return rows


def replay_grid(rows):
    """Verify an untrusted compact witness, including every interval cell."""
    if (type(rows) is not list or len(rows) != GRID+1
            or any(type(row) is not list or len(row) != 5
                   or any(type(x) is not int for x in row) for row in rows)
            or [row[0] for row in rows] != list(range(GRID+1))):
        raise ValueError('complete ordered exact interval grid required')
    structural_certificate()
    penalty = LINK_NORM/GRID  # two eigenvalues, nearest-grid distance 1/(2 GRID)
    charged, ground = [], []
    for point, charged_high, charged_next, ground_high, ground_next in rows:
        coupling = F(point, GRID)
        ch, cn, gh, gn = (F(value, DENOMINATOR) for value in
                           (charged_high, charged_next, ground_high, ground_next))
        if not (-3 < ch < cn < 0 and -3 < gh < gn < 0):
            raise ArithmeticError('ordered charged/reference spectral intervals required')
        if (count_below(2, 1, coupling, ch)[0] != 2
                or count_below(2, 1, coupling, cn)[0] != 2
                or count_below(3, 0, coupling, cn)[0] != 0
                or count_below(2, 2, coupling, gh)[0] != 1
                or count_below(2, 2, coupling, gn)[0] != 1
                or count_below(3, 1, coupling, gn)[0] != 0):
            raise ArithmeticError('exact full-spin spectral inertia failed')
        charged.append(cn-ch-penalty)
        ground.append(gn-gh-penalty)
    if min(charged) < CHARGED_GAP or min(ground) < GROUND_GAP:
        raise ArithmeticError('grid intervals do not certify the complete parameter path')
    return {'parameter_interval': ['0', '1'], 'grid_intervals': GRID,
        'dyadic_energy_denominator': DENOMINATOR,
        'middle_link_operator_norm_upper': str(LINK_NORM),
        'nearest_grid_gap_penalty': str(penalty),
        'charged_rank_in_each_of_N3_N5': 4,
        'charged_complement_domain': 'the same full N=3 or N=5 number sector',
        'charged_complement_gap_uniform_lower': str(CHARGED_GAP),
        'half_filled_ground_rank': 1,
        'half_filled_excitation_gap_uniform_lower': str(GROUND_GAP),
        'global_ground_excitation_gap_uniform_lower': str(GROUND_GAP),
        'charged_Sz_half_block_dimension': 24,
        'charged_Sz_three_halves_block_dimension': 4,
        'half_filled_Sz_zero_block_dimension': 36,
        'half_filled_Sz_one_block_dimension': 16,
        'other_number_sectors_lower_bound': '-15/8',
        'product_reference_variational_upper_bound': '-2',
        'other_number_sector_separation_lower': '1/8',
        'all_spin_multiplets_and_number_sectors_retained': True,
        'floating_eigenvalues_are_only_proposals': True}


def short_time_transfer():
    """Exact positive lower bound at the uniform coupling and duration four.

    H0 preserves block parity and V toggles it. Only odd Dyson terms can
    reach the opposite bare code. Their first term is i*j*T times identity
    on internal spin; all further terms are bounded by sinh(v*T)-v*T.
    """
    # (412311/100000)^2 > 17, so 1/sqrt(17) > 100000/412311.
    radical_upper = F(412311, 100000)
    if radical_upper*radical_upper <= 17:
        raise ArithmeticError('strict radical enclosure failed')
    j_lower = (1+1/radical_upper)/64
    x = F(1, 2)
    term, remainder = x**3/math.factorial(3), F()
    for order in range(3, 19, 2):
        remainder += term
        term *= x*x/((order+1)*(order+2))
    # The remaining odd-series ratios decrease monotonically.
    remainder += term/(1-x*x/(20*21))
    amplitude = 4*j_lower-remainder
    if amplitude <= F(1, 20):
        raise ArithmeticError('positive uniform-link transfer bound failed')
    return {'duration': '4', 'link_strength': '1', 'link_norm_upper': '1/8',
        'projected_hopping_lower': str(j_lower), 'odd_Dyson_remainder_upper': str(remainder),
        'opposite_bare_code_amplitude_lower': str(amplitude),
        'opposite_bare_code_probability_lower': str(amplitude*amplitude),
        'arbitrary_internal_spin_and_untouched_reference': True,
        'band_projection_or_postselection_used': False}


def certificate():
    structure = structural_certificate()
    rows = build_grid()
    proof = replay_grid(rows)
    digest = hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode('ascii')).hexdigest()
    return {'structure': structure, 'path': proof, 'exact_inertia_grid_sha256': digest,
            'exact_inertia_grid_rows': len(rows), 'short_time_transfer': short_time_transfer()}


report = certificate
