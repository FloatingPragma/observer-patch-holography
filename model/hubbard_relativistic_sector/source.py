"""Exact centered-law particle-hole duality and momentum convention."""
from functools import lru_cache

import numpy as np
from scipy import sparse

from ..field_symmetry import symmetry
from ..field_stress.energy import GaussianInteger as GI
from ..mass_generation_audit import chain


@lru_cache(maxsize=1)
def certificate():
    chain.certificate()
    length, modes, dimension = 4, 8, 256
    annihilators = [symmetry._annihilator(modes, j) for j in range(modes)]
    transformed = [a if j % 2 == 0 else symmetry._scale(a.adjoint(), (-1)**(j//2))
                   for j, a in enumerate(annihilators)]
    # Its vacuum is the old filled-down state. Construct every CAR basis
    # column from the transformed creation operators with exact signs.
    vacuum = np.zeros(dimension, dtype=np.int64)
    vacuum[sum(1 << (2*x+1) for x in range(length))] = 1
    columns = []
    for mask in range(dimension):
        vector = vacuum.copy()
        for j in reversed(range(modes)):
            if mask >> j & 1: vector = transformed[j].adjoint().real@vector
        columns.append(vector)
    shiba = GI.from_exact(np.stack(columns, axis=1))
    identity = symmetry._identity(dimension)
    if not (shiba@shiba.adjoint()-identity).is_zero():
        raise ArithmeticError('particle-hole CAR map is not unitary')
    for old, new in zip(annihilators, transformed, strict=True):
        if not (shiba@old@shiba.adjoint()-new).is_zero():
            raise ArithmeticError('staggered particle-hole operator rule failed')
    _, hopping, _ = chain.one_body(length)
    kinetic = symmetry._lift(hopping, modes)
    centered_contact = GI.from_exact(sparse.diags([
        4*sum((2*((mask >> (2*x)) & 1)-1)*(2*((mask >> (2*x+1)) & 1)-1)
              for x in range(length)) for mask in range(dimension)], format='csr', dtype=np.int64))
    if not (shiba@kinetic@shiba.adjoint()-kinetic).is_zero():
        raise ArithmeticError('particle-hole map changed the hopping')
    if not (shiba@centered_contact@shiba.adjoint()+centered_contact).is_zero():
        raise ArithmeticError('centered interaction failed to change sign')
    # T c_j^dag T^*=c_(j+1)^dag, including the CAR permutation sign.
    rows, signs = [], []
    for mask in range(dimension):
        destinations = [(j+2) % modes for j in range(modes) if mask >> j & 1]
        rows.append(sum(1 << j for j in destinations))
        signs.append((-1)**sum(a > b for i, a in enumerate(destinations) for b in destinations[i+1:]))
    translation = GI.from_exact(sparse.coo_matrix((signs, (rows, range(dimension))), shape=(dimension, dimension)))
    parity_down = GI.from_exact(sparse.diags([(-1)**sum((mask >> (2*x+1)) & 1 for x in range(length))
                                             for mask in range(dimension)], format='csr', dtype=np.int64))
    if not (shiba@translation@shiba.adjoint()+parity_down@translation).is_zero():
        raise ArithmeticError('momentum carrier-offset rule failed')
    return {'original_source': 'mass_generation_audit exact unchanged Weyl-chain to Hubbard join',
            'complete_Fock_dimension': dimension, 'periodic_length': length,
            'general_domain': 'even bipartite periodic rings; original period-four frame usesL=4m',
            'particle_hole_rule': 'd_(x,up)->d_(x,up),d_(x,down)->(-1)^x d_(x,down)^dagger',
            'centered_interaction_rule': '(n_up-1/2)(n_down-1/2)->-(n_up-1/2)(n_down-1/2)',
            'Melzer_normalization': 'centered source is2t times repulsive Melzer law; positive attractive parameteru=U/(4t) after particle-hole duality',
            'translation_identity': 'P T P^dagger=-(-1)^N_down T onL=4; general evenL has a convention-dependent constant phase',
            'relative_momentum_offset': 'pi*Delta N_down modulo2pi; constant phase cancels between admissible whole states',
            'source_frame_momentum': 'ordinary Hubbard translation is source translation combined with a uniform internal SU2 rotation, conjugated by the exact onsite frame',
            'continuum_p_convention': 'momentum relative to the chosen low-energy charge branch; carrier offsets are subtracted',
            'finite_volume_selection_rule_retained': True,
            'single_local_CAR_quasiparticle_state_constructed': False}
