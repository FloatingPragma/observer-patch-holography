"""Exact endpoint trial moments and certified residual upper bounds.

The two-dimensional trial space is not assumed invariant. Its ground upper
bound is for the full physical source; its Ritz levels are not mass poles.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

from model.particle_common_contract import contract
from model.particle_primary_spectral_scale import scale
from model.particle_local_wilson_density import density
from model.relativistic_sector_audit.sector import sqrt_interval

ROOT = Path(__file__).resolve().parents[2]
SOURCE_INPUTS = (
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_spectral_scale/scale.py',
    'model/particle_primary_spectral_scale/PROOF.md',
    'model/particle_primary_survival/PROOF.md',
    'model/particle_primary_locality/PROOF.md',
    'model/particle_primary_channel_labels/PROOF.md',
    'model/relativistic_sector_audit/sector.py',
    'model/observer_particles/results/representation_identification_operators.json',
)


def multiply(left, right):
    out = Counter()
    for a, x in left.items():
        for b, y in right.items():
            out[tuple(i+j for i, j in zip(a, b))] += x*y
    return {k: v for k, v in out.items() if v}


@lru_cache(maxsize=1)
def character_control():
    """Independent Weyl constant terms from the actual fifteen exterior weights."""
    raw = json.loads((ROOT/SOURCE_INPUTS[-1]).read_text())
    basis = [raw['basis_subsets'][i] for i in raw['even_indices']]
    fundamental = ((1, 0, 0, -2), (0, 1, 0, -2), (-1, -1, 0, -2),
                   (0, 0, 1, 3), (0, 0, -1, 3))
    weights = [tuple(sum(fundamental[i][k] for i in subset) for k in range(4)) for subset in basis]
    if len(weights) != 15:
        raise ArithmeticError('the actual native carrier changed')
    char = Counter(weights)
    char.update(tuple(-x for x in weight) for weight in weights)
    zero = (0, 0, 0, 0)
    weyl = {zero: 1}
    roots = [(2, 1, 0, 0), (1, 2, 0, 0), (1, -1, 0, 0), (0, 0, 2, 0)]
    for root in roots:
        for sign in (-1, 1):
            weyl = multiply(weyl, {zero: 1, tuple(sign*x for x in root): -1})
    power = {zero: 1}; integrals = []
    for degree in range(1, 4):
        power = multiply(power, char)
        value = sum(weight*weyl.get(tuple(-x for x in exponent), 0)
                    for exponent, weight in power.items())
        integrals.append(F(value, 12))
    if integrals != [F(0), F(10), F(54)]:
        raise ArithmeticError('native character moments differ from the exact Weyl control')
    # A second calculation uses only fundamental/antifundamental and doublet
    # tensor singlets, rather than the weight polynomial.
    reps = [(-4, -1, 0), (-3, 0, 1), (1, 1, 1), (2, -1, 0), (6, 0, 0)]
    reps += [(-q, -color, weak) for q, color, weak in reps]
    triples = []
    block_order = (-4, -3, 1, 2, 6)
    locate = {abs(q): j for j, q in enumerate(block_order)}
    singlet_matrix = [[0]*5 for _ in range(5)]
    for rows in product(reps, repeat=3):
        if (sum(r[0] for r in rows) == 0 and sum(r[2] for r in rows) in (0, 2)
                and sorted(r[1] for r in rows) in ([-1, -1, -1], [-1, 0, 1], [0, 0, 0], [1, 1, 1])):
            triples.append(tuple(sorted(r[0] for r in rows)))
            singlet_matrix[locate[abs(rows[0][0])]][locate[abs(rows[2][0])]] += 1
    if len(triples) != 54:
        raise ArithmeticError('tensor-singlet and Weyl character controls differ')
    chars = []
    for q in block_order:
        values = [w for w in weights if w[-1] == q]
        entry = Counter(values)
        entry.update(tuple(-x for x in w) for w in values)
        chars.append(entry)
    for i, left in enumerate(chars):
        for j, right in enumerate(chars):
            polynomial = multiply(multiply(left, char), right)
            value = sum(weight*weyl.get(tuple(-x for x in exponent), 0)
                        for exponent, weight in polynomial.items())
            if F(value, 12) != singlet_matrix[i][j]:
                raise ArithmeticError('projected loop matrix differs between exact character methods')
    return {'native_weights': weights, 'Weyl_group_order': 12,
            'real_character_sum_moments_1_2_3': list(map(str, integrals)),
            'normalized_real_plaquette_second_moment': '1/90',
            'normalized_real_plaquette_third_moment': '1/500',
            'ordered_tensor_singlets': len(triples),
            'block_order_Q6': list(block_order),
            'cubic_loop_singlet_matrix': singlet_matrix,
            'singlet_charge_multisets': [{'charges': list(q), 'ordered_multiplicity': count}
                                      for q, count in sorted(Counter(triples).items())]}


def local_fourth_bound(total_term_norm, maximum_overlap_norm):
    """Product-vector centered even terms; disjoint supports are independent."""
    total = F(total_term_norm); overlap = F(maximum_overlap_norm)
    if total < 0 or overlap < 0:
        raise ValueError('nonnegative exact norm bounds required')
    return 3*total*total*overlap*overlap+16*total*overlap**3


def orthogonal_packet_fourth(total_term_norm, maximum_overlap_norm, variance):
    """Sharper disconnected part when distinct local packet vectors are orthogonal."""
    total, overlap, second = map(F, (total_term_norm, maximum_overlap_norm, variance))
    if min(total, overlap, second) < 0:
        raise ValueError('nonnegative exact packet bounds required')
    return {'disconnected_pairings_upper': 3*second*second,
            'connected_quadruples_upper': 16*total*overlap**3,
            'fourth_moment_upper': 3*second*second+16*total*overlap**3}


def ritz_trial(variance, first_packet_moment, second_packet_norm, fourth_perturbation_upper,
               cubic_perturbation_moment, bits=80):
    """A rational near-Ritz trial and a full-operator residual bound.

Inputs are exact moments on Omega and v=V Omega. The asserted fourth
moment bound is an input to this general arithmetic helper.
"""
    z, m, h2, v4, v3 = map(F, (variance, first_packet_moment, second_packet_norm,
                               fourth_perturbation_upper, cubic_perturbation_moment))
    if z <= 0 or m <= 0:
        raise ValueError('positive nonzero trial variance and diagonal moment required')
    b = m/z
    rad_lo, rad_hi = sqrt_interval(b*b+4*z, bits=bits)
    # Rationalize the small root to avoid subtractive loss.
    delta_lo, delta_hi = 2*z/(rad_hi+b), 2*z/(rad_lo+b)
    t = delta_lo/z
    norm = 1+t*t*z
    energy = (-2*t*z+t*t*m)/norm
    hfirst = m-v3
    hvar = h2-hfirst*hfirst/z
    vvar = v4-z*z-v3*v3/z
    if min(hvar, vvar) < 0:
        raise ValueError('packet moments do not admit the declared residual bound')
    cross = sqrt_interval(hvar*vvar, bits=bits)[1]
    outside = hvar+vvar+2*cross
    residual2 = ((-t*z-energy)**2+z*(1-t*b+energy*t)**2+t*t*outside)/norm
    radius = sqrt_interval(residual2, bits=bits)[1]
    if not energy < 0 or not -delta_hi <= energy:
        raise ArithmeticError('rational trial failed the Ritz variational enclosure')
    return {'packet_variance': z, 'packet_diagonal': b,
            'lower_Ritz_shift_interval': (-delta_hi, -delta_lo),
            'upper_Ritz_shift_interval': (b+delta_lo, b+delta_hi),
            'rational_trial_t': t, 'trial_norm_squared_before_normalization': norm,
            'trial_energy_shift': energy,
            'trial_minus_Ritz_energy_upper': energy+delta_hi,
            'outside_trial_space_squared_norm_upper': outside,
            'trial_residual_squared_upper': residual2,
            'trial_residual_upper': radius,
            'trial_spectral_interval_relative_to_anchor': (energy-radius, energy+radius),
            'trial_subspace_claimed_invariant': False}


def solve(matrix, rhs):
    """Exact rational positive-matrix solve; no floating eigenvalue enters a bound."""
    size = len(rhs)
    rows = [[F(x) for x in row]+[F(rhs[i])] for i, row in enumerate(matrix)]
    for i in range(size):
        pivot = rows[i][i]
        if pivot <= 0:
            raise ArithmeticError('the packet resolvent lost a positive pivot')
        rows[i] = [x/pivot for x in rows[i]]
        for j in range(i+1, size):
            scale = rows[j][i]
            rows[j] = [a-scale*b for a, b in zip(rows[j], rows[i])]
    answer = [F(0)]*size
    for i in reversed(range(size)):
        answer[i] = rows[i][-1]-sum(rows[i][j]*answer[j] for j in range(i+1, size))
    return answer


def packet_ritz(bonds, plaquettes, packets, bits=60):
    """Exact 11-dimensional compression at coupling one, with complement error."""
    order = character_control()['block_order_Q6']
    lookup = {row['integer_Q6_charge']: row for row in packets}
    rows = [lookup[q] for q in order]
    dims = [row['dimension'] for row in rows]
    delta_h = [F(row['bond_H0_excitation']) for row in rows]
    delta_l = [F(row['loop_H0_excitation']) for row in rows]
    nh = [F(bonds*d) for d in dims]; nl = F(plaquettes, 22050)
    counts = character_control()['cubic_loop_singlet_matrix']
    correction = [[-F(x, 420) for x in row] for row in counts]
    loop = [[correction[i][j]+(delta_l[i] if i == j else 0) for j in range(5)] for i in range(5)]
    if any(loop[i][i] <= sum(abs(loop[i][j]) for j in range(5) if j != i) for i in range(5)):
        raise ArithmeticError('the exact loop packet matrix is not positive diagonally dominant')
    def coefficients(d):
        h = [1/(x+d) for x in delta_h]
        matrix = [[loop[i][j]+(d if i == j else 0) for j in range(5)] for i in range(5)]
        return h, solve(matrix, [1]*5)
    def response(d):
        h, l = coefficients(d)
        return sum(x*t for x, t in zip(nh, h))+nl*sum(l)
    lo, hi = F(0), sqrt_interval(sum(nh)+5*nl)[1]
    if response(hi) > hi:
        raise ArithmeticError('the packet Schur root escaped its variance bound')
    target = F(1, 2**bits)
    while hi-lo > target:
        middle = (hi+lo)/2
        if response(middle) >= middle:
            lo = middle
        else:
            hi = middle
    th, tl = coefficients(lo)
    if any(t <= 0 for t in th+tl):
        raise ArithmeticError('packet trial coefficients are not positive')
    q = sum(x*t for x, t in zip(nh, th))+nl*sum(tl)
    square = sum(x*t*t for x, t in zip(nh, th))+nl*sum(t*t for t in tl)
    norm = 1+square
    energy = (-q-lo*square)/norm
    if not -hi <= energy < 0:
        raise ArithmeticError('rational packet trial violates the exact Ritz bracket')
    # Parent V and weighted V_t commute. Tensor cutoff gives a contraction on
    # V_parent V_t Omega, so parent fourth moments bound the finite complement.
    bond_norm = 2*sum(F(d)*t for d, t in zip(dims, th))
    loop_norm = sum(F(d)*t for d, t in zip(dims, tl))/105
    total = 30*bonds+F(plaquettes, 7)
    weighted_total = bonds*bond_norm+plaquettes*loop_norm
    fourth_parts = orthogonal_packet_fourth(total, 3*F(1272, 7), sum(nh)+5*nl)
    weighted_parts = orthogonal_packet_fourth(weighted_total, 3*(6*bond_norm+12*loop_norm), square)
    fourth = fourth_parts['fourth_moment_upper']
    weighted_fourth = weighted_parts['fourth_moment_upper']
    local_mixed = sqrt_interval(fourth*weighted_fourth)[1]
    global_mixed = total*total*square
    mixed = min(local_mixed, global_mixed)
    projected_loop = [sum(correction[i][j]*tl[j] for j in range(5)) for i in range(5)]
    outside = mixed-q*q-nl*sum(x*x for x in projected_loop)
    if outside < 0:
        raise ArithmeticError('the full-space fourth-moment bound is inconsistent')
    residual2 = ((q+energy)**2+(lo+energy)**2*square+outside)/norm
    return {'packet_dimension': 11, 'block_order_Q6': order,
            'hopping_packet_norms_squared': nh, 'real_loop_packet_norm_squared': nl,
            'hopping_diagonal': delta_h, 'loop_packet_matrix': loop,
            'lower_Ritz_shift_interval': (-hi, -lo),
            'rational_hopping_coefficients': th, 'rational_loop_coefficients': tl,
            'trial_norm_squared_before_normalization': norm,
            'trial_energy_shift': energy,
            'trial_minus_Ritz_energy_upper': energy+hi,
            'parent_V_fourth_moment_upper': fourth,
            'weighted_parent_V_fourth_moment_upper': weighted_fourth,
            'parent_fourth_moment_components': fourth_parts,
            'weighted_fourth_moment_components': weighted_parts,
            'local_mixed_packet_norm_squared_upper': local_mixed,
            'global_mixed_packet_norm_squared_upper': global_mixed,
            'selected_complement_bound': 'local covariance and connected remainder' if local_mixed < global_mixed else 'global operator norm',
            'mixed_packet_norm_squared_upper': mixed,
            'outside_trial_space_squared_norm_upper': outside,
            'trial_residual_squared_upper': residual2,
            'trial_residual_upper': sqrt_interval(residual2)[1],
            'trial_subspace_claimed_invariant': False}


def selected_bounds(member):
    contract.verify_member(member)
    character_control()
    source = member['source']; a = F(source['spacing']); volume = source['side']**3
    data = scale.selected_expansion(member)
    bonds = data['actual_interior_bonds']; plaquettes = data['actual_interior_plaquettes']
    packets = scale.vacuum_channels()['all_native_blocks']
    z = F(data['dimensionless_primary_H_variance_on_anchor'])
    hfirst = sum(bonds*F(row['bond_vacuum_transition_norm_squared'])*F(row['bond_H0_excitation'])
                 +plaquettes*F(row['both_loop_orientations_norm_squared'])*F(row['loop_H0_excitation'])
                 for row in packets)
    h2 = sum(bonds*F(row['bond_vacuum_transition_norm_squared'])*F(row['bond_H0_excitation'])**2
             +plaquettes*F(row['both_loop_orientations_norm_squared'])*F(row['loop_H0_excitation'])**2
             for row in packets)
    v3 = -F(plaquettes, 171500)
    if hfirst != 720*bonds+F(11, 350)*plaquettes:
        raise ArithmeticError('the actual packet first moment changed')
    total = F(data['dimensionless_global_perturbation_norm_upper'])
    overlap = 3*F(1272, 7)
    fourth = min(orthogonal_packet_fourth(total, overlap, z)['fourth_moment_upper'], total*total*z)
    two = ritz_trial(z, hfirst+v3, h2, fourth, v3)
    trial = packet_ritz(bonds, plaquettes, packets)
    if trial['trial_energy_shift'] >= two['trial_energy_shift']:
        raise ArithmeticError('resolving the exact packets did not improve the two-vector control')
    anchor = -90*volume
    upper = anchor+trial['trial_energy_shift']
    lower = F(member['preparation_plan']['source']['spectral_lower'])
    radius = trial['trial_residual_upper']/a
    center = upper/a
    return contract.serial({'schema': 'oph.primary_variational_endpoint.v1',
        'member_id': member['member_id'], **member['ids'],
        'coupling_evaluated': 1, 'dimensionless_H': 'a times the unchanged primary H',
        'trial_anchor': 'filled negative-beta physical product vector Omega0',
        'trial_packet': 'five global hopping packets and five real loop packets, with Omega0',
        'trial_is_selected_polynomial_preparation': False,
        'ground_energy_upper': center,
        'dimensionless_ground_energy_upper': upper,
        'inherited_full_source_spectral_lower': lower,
        'ground_energy_interval': (lower, center),
        'product_trial_ground_energy_upper': F(anchor)/a,
        'strict_ground_upper_improvement': -trial['trial_energy_shift']/a,
        'actual_bonds': bonds, 'actual_plaquettes': plaquettes,
        'variance': z, 'H0_packet_first_moment': hfirst,
        'H0_packet_second_moment': h2, 'cubic_perturbation_moment': v3,
        'perturbation_total_term_norm_upper': total,
        'perturbation_overlap_norm_upper': overlap,
        'perturbation_fourth_moment_upper': fourth,
        'dimensionless_trial': trial,
        'dimensionless_two_vector_control': two,
        'trial_rayleigh_quotient': center,
        'trial_residual_upper': radius,
        'some_full_physical_eigenvalue_interval': (center-radius, center+radius),
        'residual_interval_identified_as_ground_or_excitation': False,
        'trial_spectral_weight_outside_radius_rule': 'at most residual_squared/radius_squared',
        'original_primary_probe_residue_computed': False,
        'first_excitation_gap_or_rest_mass': None,
        'trial_is_topological_defect_or_single_charged_particle': False,
        'source_reference_state_and_clock_replaced': False,
        'full_spatial_matrix_executed': False})


def bind_member(member):
    result = selected_bounds(member)
    return {'result': result, 'result_sha256': contract.digest(result)}


def verify_bound(member, value):
    if value != bind_member(member):
        raise ValueError('the variational bound differs from the actual source and trial')


def certificate():
    return {'schema': 'oph.primary_variational_family.v1', 'character_control': character_control(),
            'members': [bind_member(contract.member(n)) for n in (0, 1, 2)],
            'physical_mass_table_GeV': []}
