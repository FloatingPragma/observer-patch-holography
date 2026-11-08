"""Same-source integer spin charge and charge-resolved local creation."""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations

import numpy as np

from model.native_ordered_ground import source, certificate
from model.centered_particle_spectrum.certificate import decimal_bounds
from model.defect_field_bridge import full_domain
from model.protected_spin_limit.experiment import full_hamiltonian
from model.strong_coupling_refinement.algebra import geometry

ROOT = Path(__file__).resolve().parents[2]


def spin_flips(mask, sites):
    """Columns of twice the full-Fock total internal Sx, including singlets."""
    return tuple(mask ^ (3 << (2*x)) for x in range(sites)
                 if ((mask >> (2*x)) & 3) in (1, 2))


def full_source_charge():
    columns = 0
    for sites in (2, 3, 4):
        for seam in combinations(range(sites), 2):
            for law in ('C', 'E'):
                for coin in (0, 1):
                    for mask in range(1 << (2*sites)):
                        left, right = {}, {}
                        for changed in spin_flips(mask, sites):
                            raw, out, sign, current = full_domain.branch_action(
                                changed, law, coin, sites=sites, flavors=2, seam=seam)
                            key = (raw, out, current)
                            left[key] = left.get(key, 0)+sign
                        raw, out, sign, current = full_domain.branch_action(
                            mask, law, coin, sites=sites, flavors=2, seam=seam)
                        for changed in spin_flips(out, sites):
                            key = (raw, changed, current)
                            right[key] = right.get(key, 0)+sign
                        if {k:v for k,v in left.items() if v} != {k:v for k,v in right.items() if v}:
                            raise ArithmeticError('a complete source branch changes total spin charge')
                        columns += 1
    _, edges = geometry((2,2,1))
    h = full_hamiltonian(4, edges, 2, 8)
    if not np.array_equal(h, h.astype(np.int64)):
        raise ArithmeticError('integer full-source Hamiltonian required')
    h = h.astype(np.int64)
    q2 = np.zeros((256,256), dtype=np.int64)
    for mask in range(256):
        for out in spin_flips(mask, 4):
            q2[out, mask] += 1
    if not np.array_equal(h @ q2, q2 @ h):
        raise ArithmeticError('signed full source does not conserve total spin charge')
    # Pull the full-Fock generator through the actual spin inclusion and
    # Marshall diagonal. This determines the staggered sign, rather than
    # choosing a convenient spin generator after seeing the spectrum.
    odd = set(source.cube()[1]); induced = np.zeros((256,256), dtype=np.int64)
    for word in range(256):
        mask = sum(1 << (2*x+((word >> x) & 1)) for x in range(8))
        for out in spin_flips(mask, 8):
            other = sum(((out >> (2*x+1)) & 1) << x for x in range(8))
            marshall = (-1)**sum(((word >> x) & 1)+((other >> x) & 1) for x in odd)
            induced[other, word] += marshall
    if not np.array_equal(induced, operators()['twice_charge']):
        raise ArithmeticError('charge does not intertwine through the actual spin inclusion')
    return {'complete_all_seam_full_Fock_columns':columns,
            'full_signed_source_charge_commutator_exactly_zero':True,
            'every_C_E_fine_record_commutes_with_full_source_charge':True,
            'source_charge_pulled_through_spin_and_Marshall_maps':True,
            'full_source_charge':'sum_x(c_up_dagger c_down+c_down_dagger c_up)/2',
            'coherent_controls_are_all_required_to_preserve_charge':False}


def operators():
    h, scale = source.matrix()
    odd = set(source.cube()[1])
    signs = tuple(-1 if x in odd else 1 for x in range(8))
    indices = np.arange(256)
    twice_charge = np.zeros((256, 256), dtype=np.int64)
    for x, sign in enumerate(signs):
        twice_charge[indices ^ (1 << x), indices] += sign
    if not np.array_equal(twice_charge @ h, h @ twice_charge):
        raise ArithmeticError('actual source does not conserve the stated U1 spin charge')
    permutation = np.array([sum(((word >> x) & 1) << (x ^ 4) for x in range(8))
                            for word in range(256)])
    if (not np.array_equal(h[np.ix_(permutation, permutation)], h)
            or not np.array_equal(twice_charge[np.ix_(permutation, permutation)], -twice_charge)):
        raise ArithmeticError('source symmetry does not reverse the conserved charge')
    z = np.diag(np.array([1-2*(word & 1) for word in range(256)], dtype=np.int64))
    xz = z[indices ^ 1, :]
    creation = {1:z+signs[0]*xz, -1:z-signs[0]*xz}
    for q, a in creation.items():
        if not np.array_equal(twice_charge @ a-a @ twice_charge, 2*q*a):
            raise ArithmeticError('local source creation is not charge homogeneous')
    if not np.array_equal(creation[1]+creation[-1], 2*z):
        raise ArithmeticError('complete local creation decomposition lost a channel')
    # W/16 is the exact tensor-X basis. All256 outcomes are retained.
    walsh = np.array([[(-1)**((word & label).bit_count()) for label in range(256)]
                      for word in range(256)], dtype=np.int64)
    charges = tuple(sum(signs[x]*(1-2*((label >> x) & 1)) for x in range(8))//2
                    for label in range(256))
    if (not np.array_equal(walsh.T @ walsh, 256*np.eye(256, dtype=np.int64))
            or not np.array_equal(twice_charge @ walsh, walsh*np.array(charges)[None, :]*2)):
        raise ArithmeticError('complete X readout does not resolve the charge spectral effects')
    parity = np.eye(256, dtype=np.int64)[indices ^ 255]
    if not np.array_equal(parity @ walsh, walsh*np.array([1 if q % 2 == 0 else -1 for q in charges])[None, :]):
        raise ArithmeticError('global parity is not the charge-parity subgroup')
    return {'H':h, 'scale':scale, 'twice_charge':twice_charge,
            'twice_creation':creation, 'local_Z':z, 'walsh':walsh,
            'charges':charges, 'signs':signs, 'charge_reversing_permutation':permutation}


def ground_input():
    with np.load(ROOT/'model/native_ordered_ground/witness.npz', allow_pickle=False) as archive:
        witness = {key:archive[key] for key in archive.files}
    checked = certificate.verify(witness)
    return witness['vectors'][:, 0].astype(object), F(checked['normalized_ground_vector_error']), checked


def report():
    op = operators()
    vector, error, checked = ground_input()
    norm = sum(int(x)**2 for x in vector)
    ref_coefficients = op['walsh'].T.astype(object) @ vector
    created_coefficients = op['walsh'].T.astype(object) @ (op['local_Z'].astype(object) @ vector)
    def weights(coefficients):
        return {q:F(sum(int(coefficients[i])**2 for i, charge in enumerate(op['charges']) if charge == q),
                    256*norm) for q in range(-4, 5)}
    reference, created = weights(ref_coefficients), weights(created_coefficients)
    if sum(reference.values()) != 1 or sum(created.values()) != 1:
        raise ArithmeticError('complete charge readout lost probability')
    if sum(p for q, p in reference.items() if q != 0) > error**2:
        raise ArithmeticError('proposal disagrees with exact unique-ground zero charge')
    if sum(p for q, p in created.items() if abs(q) != 1) > error**2:
        raise ArithmeticError('local-created proposal leaks beyond the certified charge error')
    rows = []
    for q, twice_a in op['twice_creation'].items():
        image = twice_a.astype(object) @ vector
        proposal_weight = F(sum(int(x)**2 for x in image), 4*norm)
        lower, upper = max(F(0), proposal_weight-2*error), min(F(1), proposal_weight+2*error)
        if lower <= F(4, 25):
            raise ArithmeticError('a nonzero source charge channel lost its certified weight')
        rows.append({'charge':q, 'complete_creation_probability':decimal_bounds(lower, upper, 12),
                     'proposal_probability_exact':str(proposal_weight),
                     'same_H_B_energy_at_least_ground_gap':checked['finite_field_gap_interval'][0],
                     'physical_electric_charge_identified':False, 'physical_mass':None})
    return {'source_charge':'Q=sum_x epsilon_x X_x/2 in the actual Marshall frame',
            'full_source_charge_attachment':full_source_charge(),
            'source_Hamiltonian_commutator_exactly_zero':True,
            'global_X_equals_exp_i_pi_Q':True,
            'unique_selected_ground_charge':0,
            'ground_charge_proof':'nondegenerate positive ground and exact charge-reversing source automorphism',
            'charge_spectral_multiplicities':{str(q):op['charges'].count(q) for q in range(-4, 5)},
            'complete_X_record_count':256,
            'ideal_created_charge_support':[-1, 1],
            'creation_channels':sorted(rows, key=lambda r:r['charge']),
            'ground_proposal_vector_error':str(error),
            'proposal_reference_outside_Q0_probability':str(sum(p for q,p in reference.items() if q)),
            'proposal_created_outside_Qplusminus1_probability':str(sum(p for q,p in created.items() if abs(q)!=1)),
            'complete_ideal_charge_instrument_constructed':True,
            'coherent_creation_instrument':'A_q=Pi_(q epsilon_0 X_0) Z_0; two local pointer outcomes',
            'coherent_creation_complete_effect_sum_is_identity':True,
            'coherent_creation_native_parent':'model/defect_localization_readout/results.json',
            'fine_X_outcome_grouping_is_coherent_charge_preparation':False,
            'full_native_256_record_charge_tape_executed':False,
            'native_pointer_composition_bound':'at most sum of the eight complete local-pointer channel errors; all failures retained',
            'two_charge_channels_claimed_as_two_physical_particles':False}
