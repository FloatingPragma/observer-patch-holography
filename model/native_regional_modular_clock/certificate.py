"""Exact regional checks against the original complete H66 source counts."""
from fractions import Fraction as F

import numpy as np

from model.primitive_observer import transfer
from . import regional as r


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def certificate():
    source = transfer.native_counts()
    transfer.exact_audit(source)
    edges = source.seams
    require(len(edges) == 30, 'the native seam inventory changed')
    basis, rho, normalization = r.reference(12, edges, 2)
    require(normalization == 56, 'the native transfer trace changed')
    actual_masks = tuple(sum(int(x) << i for i, x in enumerate(state))
                         for state in source.states)
    where = {mask: i for i, mask in enumerate(actual_masks)}
    require(all(rho[i, j] == F(source.counts[where[x]][where[y]], 60 * 56)
                for i, x in enumerate(basis) for j, y in enumerate(basis)),
            'the regional reference is not the literal native transfer reference')
    mixture = r.coherent_reference_mixture(12, edges, 2)
    require(mixture['accepted_probability'] == F(28, 33),
            'the conditional sampler acceptance changed')
    rows = []
    for region in edges:
        direct = r.partial_trace(12, basis, rho, region)
        for number in range(3):
            local, block, _ = r.regional_block(12, edges, 2, region, number)
            require(local == direct[number][0] and np.array_equal(block, direct[number][1]),
                    'the exact local reduction disagrees with a direct partial trace')
        crossings = r.split(12, edges, region)[2]
        occupation = np.array([sum((x >> i) & 1 for i in region) for x in basis], dtype=object)
        commutator = rho * (occupation[None, :] - occupation[:, None])
        hs_squared = sum(x * x for x in commutator.flat)
        require(hs_squared == r.cross_sector_hs_squared(12, 30, 2, len(crossings))
                == F(1, 70560), 'the regional coherence obstruction changed')
        rows.append({
            'region': region, 'crossing_edges': len(crossings),
            'regional_density_blocks': {
                str(k): [[str(x) for x in line] for line in block]
                for k, (_, block) in direct.items()},
            'number_commutator_HS_squared': str(hs_squared)})
    comparisons = []
    for level in range(1, 9):
        population, seams = 240 * 4**level, 630 * 4**level
        bounds = [r.modulus_bound(population, seams, population // 2, 12, 30, 3, k, F(1, 2))
                  for k in range(13)]
        comparisons.append({
            'level': level, 'ports': population, 'seams': seams, 'local_sectors': 13,
            'common_density': '1/2', 'scaled_modular_target': '(4/3) G_cell',
            'maximum_generator_error': str(max(item['generator_error'] for item in bounds)),
            'log_gate_passed': all(item['z'] < 1 for item in bounds)})
    return {
        'native_source': {'dimension': 66, 'P_trace': str(normalization),
                          'all_transfer_entries_authenticated': 66**2,
                          'regional_seams_checked': 30,
                          'all_local_occupation_blocks_kept': True},
        'regional_seams': rows,
        'conditional_sampler': {
            'accepted_probability': str(mixture['accepted_probability']),
            'rejected_probability': str(mixture['rejected_probability']),
            'canonical_sampler_and_addressed_coherent_mixer_are_premises': True,
            'growing_native_service_compilation_supplied': False},
        'finite_federation_density_bounds': comparisons,
        'exact_overlap_obstruction': {key: str(value) for key, value in r.overlap_witness().items()}}
