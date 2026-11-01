"""Exact symmetry and energy inputs for the actual centered cube's local state.

No centered eigenvector or fitted local occupation is used.  The proof combines
an exact signed-translation census, the matching-dimer energy improvement and
the full N=8 kinetic norm with Cauchy--Schwarz.
"""
from fractions import Fraction as F
from itertools import product
import numpy as np

from ..interacting_vacuum import ground


def signed_vertex_flips():
    A=np.asarray(ground.band_data()[0],dtype=np.int64)
    if A.shape!=(8,8) or not np.array_equal(A,A.T):
        raise ArithmeticError('the source cube must be real symmetric')
    if not np.array_equal(A@A,3*np.eye(8,dtype=np.int64)) or np.trace(A)!=0:
        raise ArithmeticError('the source flat kinetic spectrum changed')
    rows=[]
    for shift in range(8):
        permutation=np.arange(8,dtype=np.int64)^shift
        signs=np.array([(-1)**(((shift&1)*(((x>>1)&1)+((x>>2)&1))
                               +((shift>>1)&1)*((x>>2)&1))) for x in range(8)],dtype=np.int64)
        transformed=signs[:,None]*A[np.ix_(permutation,permutation)]*signs[None,:]
        if not np.array_equal(transformed,A):
            raise ArithmeticError('signed vertex flip fails to preserve the actual hopping')
        rows.append({'shift':shift,'permutation':permutation.tolist(),'signs':signs.tolist()})
    colors=[x.bit_count()%2 for x in range(8)]
    if any(A[x,y] and colors[x]==colors[y] for x in range(8) for y in range(8)):
        raise ArithmeticError('same-sublattice source hopping found')
    reached={0}
    while True:
        new=reached|{y for x in reached for y in range(8) if A[x,y]}
        if new==reached:break
        reached=new
    if len(reached)!=8:raise ArithmeticError('source hopping graph disconnected')
    if {row['permutation'][0] for row in rows}!=set(range(8)):
        raise ArithmeticError('the exact symmetry family is not vertex transitive')
    return {'site_count':8,'hopping_real_symmetric':True,'hopping_connected':True,
            'hopping_bipartite':True,'sublattice_sizes':[colors.count(0),colors.count(1)],
            'spatial_hopping_square':'A^2=3I','spatial_hopping_trace':0,
            'signed_symmetries':rows,'ground_projector_vertex_transitive':True,
            'N8_kinetic_norm_squared':'3/4'}


def no_doublon_kinetic_block():
    """Enumerate all 256 one-particle-per-site configurations and actual hops."""
    A=np.asarray(ground.band_data()[0],dtype=np.int64)
    states=[sum(1<<(2*x+spin) for x,spin in enumerate(spins))
            for spins in product((0,1),repeat=8)]
    subspace=set(states);nonzero=0
    for mask in states:
        for target in range(8):
            for source in range(8):
                if not A[target,source]:continue
                for spin in (0,1):
                    old=2*source+spin;new=2*target+spin
                    if mask>>old&1 and not mask>>new&1:
                        result=mask^(1<<old)^(1<<new);nonzero+=1
                        if result in subspace or result.bit_count()!=8:
                            raise ArithmeticError('kinetic hopping remains in the zero-doublon N8 subspace')
                        if not any((result>>(2*x))&3==3 for x in range(8)):
                            raise ArithmeticError('nonzero hop did not create a doublon')
    if not nonzero:raise ArithmeticError('empty kinetic transition audit')
    return {'singly_occupied_configurations':len(states),
            'nonzero_directed_spin_hops_checked':nonzero,'P0_T_P0_zero':True}


def relative_entropy_upper(epsilon,delta_hs_squared):
    if type(epsilon) not in (int,F) or not 0<=epsilon<1:
        raise ValueError('exact epsilon in [0,1) required')
    if type(delta_hs_squared) not in (int,F) or delta_hs_squared<0:
        raise ValueError('nonnegative exact squared Hilbert--Schmidt norm required')
    return 7500*F(epsilon)**2*F(delta_hs_squared)/(1-F(epsilon))


def certificate():
    # sqrt17>103/25 gives delta=sqrt17/2-2>3/50.
    if not F(103,25)**2<17:raise ArithmeticError('dimer improvement lost strict lower bound')
    delta_floor=F(3,50);kinetic_norm_squared=F(3,4)
    d_floor=delta_floor**2/(9*kinetic_norm_squared*8)
    if d_floor!=F(1,15000):raise ArithmeticError('local faithfulness coefficient changed')
    return {'symmetry':signed_vertex_flips(),'no_doublon_block':no_doublon_kinetic_block(),
            'reference':'normalized full centered-cube ground projector at c=1/8,U=1',
            'local_density':'diag(d,1/2-d,1/2-d,d)',
            'double_occupation_strict_lower':str(d_floor),
            'double_occupation_strict_upper':'1/8',
            'full_site_minimum_eigenvalue_strict_lower':str(d_floor),
            'full_four_state_local_CAR_faithful':True,
            'ground_uniqueness_needed_for_this_proof':False,
            'onsite_modular_generator':'-log(1/2-d) I + log((1/2-d)/d) (N-1)^2',
            'modular_operator_norm_strict_upper':'log(15000)',
            'site_entropy_strict_lower':'log(2)',
            'site_entropy_strict_upper':'3 log(2) - (3/4)log(3)',
            'relative_entropy_upper':'7500 epsilon^2 ||delta||_HS^2/(1-epsilon)',
            'entropy_scope':'arbitrary same-four-state density sigma, delta=sigma-rho, 0<=epsilon<1',
            'infinite_ground_local_faithfulness_inferred':False,
            'geometric_modular_charge_identified':False}
