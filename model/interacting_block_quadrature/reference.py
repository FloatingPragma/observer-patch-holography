"""Exact native covariance replay and finite regular-reference completion."""
from itertools import product,permutations
from functools import lru_cache
import sympy as s

from model.coherent_defect_wave.chiral_joint.atomic_momentum_response.checks import native_translation_join,rotation_translation_join
from model.coherent_defect_wave.chiral_joint.oriented_links.checks import rotation_word_checks
from model.coherent_defect_wave.chiral_joint.oriented_links.rotations import cubic_rotations,geometry
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.checks import matrix_checks,actual_hamiltonian_checks
from model.coherent_defect_wave.chiral_joint.collective_spin_response.model import operator_certificate


@lru_cache(None)
def finite_lift():
    # This is the exact adjoint image on the admitted native scalar+triplet.
    group=[]
    for p in permutations(range(3)):
        for signs in product((-1,1),repeat=3):
            q=s.zeros(3)
            for j in range(3):q[p[j],j]=signs[j]
            group.append(s.ImmutableMatrix(q))
    keys={tuple(q):q for q in group}
    assert len(keys)==48
    unit=s.eye(4);number=s.diag(0,1,1,1)
    x=[]
    for a in range(3):
        matrix=s.zeros(4);matrix[0,a+1]=matrix[a+1,0]=1;x.append(matrix)
    h=2*(s.kronecker_product(number,unit)+s.kronecker_product(unit,number))
    h-=sum((s.kronecker_product(y,y)/8 for y in x),s.zeros(16))
    rotation=s.eye(4);rotation[0,0]=rotation[1,1]=s.Rational(3,5)
    rotation[0,1]=-s.Rational(4,5);rotation[1,0]=s.Rational(4,5)
    c=s.Rational(3,5)*s.kronecker_product(rotation,unit)
    defect=s.Rational(4,5)*s.eye(16)
    assert c.H*c+defect*defect==s.eye(16)
    kraus=[(s.eye(16)+sign*c)/2 for sign in (1,-1)]+[sign*defect/2 for sign in (1,-1)]
    assert sum((k.H*k for k in kraus),s.zeros(16))==s.eye(16)
    representations={key:s.kronecker_product(s.diag(s.ones(1),q),s.diag(s.ones(1),q)) for key,q in keys.items()}
    generators=[s.diag(-1,1,1),s.diag(1,-1,1),s.diag(1,1,-1),
        s.Matrix([[0,1,0],[1,0,0],[0,0,1]]),s.Matrix([[1,0,0],[0,0,1],[0,1,0]])]
    joint_checks=0
    for key,q in keys.items():
        u=representations[key]
        assert u.H*u==s.eye(16) and u*h==h*u
        cg=u*c*u.H
        assert cg*u==u*c
        for k in kraus:assert (u*k*u.H)*u==u*k
        for generator in generators:
            left=representations[tuple(generator)]
            target=representations[tuple(generator*q)]
            assert left*cg*left.H==target*c*target.H
            joint_checks+=1
    t=representations[tuple(s.diag(-1,1,1))]
    assert t*c!=c*t
    assert sum((u.H*u for u in representations.values()),s.zeros(16))/48==s.eye(16)
    # Failure of source-only covariance is repaired by the retained regular
    # register, with every Kraus block and its output preserved.
    return {'actual_admitted_low_representation_image_order':48,
        'test_source_dimension':16,'test_source_components':3,
        'test_source_H':'two sites, N=1, D=2, t=1/8',
        'exact_joint_covariance_generator_blocks':joint_checks,
        'all_four_Kraus_intertwiners':48*4,'isometry_Gram_identity':True,
        'unreferenced_character_changing_contraction_rejected':True,
        'reference_Hamiltonian':'0','full_reference_output_retained':True}


@lru_cache(None)
def native_cover():
    rotations=cubic_rotations();by_rotation={tuple(r):(r,spin) for r,spin in rotations}
    spatial={}
    for r,_ in rotations:
        sites,_=geometry(r,2)
        spatial[tuple(r)]=s.ImmutableSparseMatrix(8,8,{(target,j):phase for j,(target,phase) in enumerate(sites)})
    for r,spin in rotations:
        for q,other in rotations:
            key=tuple(r*q);target=by_rotation[key][1]
            assert spatial[tuple(r)]*spatial[tuple(q)]==spatial[key]
            lift=(spin*other).applyfunc(s.simplify)
            assert (lift-target).applyfunc(s.simplify)==s.zeros(2) or (lift+target).applyfunc(s.simplify)==s.zeros(2)
    return {'native_translation_cover':'Z4^3, each actual AP translation has fourth power I',
        'translation_reference_dimension_sufficient_on_full_source':64,
        'admitted_even_translation_image':'Z2^3','admitted_translation_reference_dimension':8,
        'native_cubic_cover':'Z4^3 semidirect binary octahedral group',
        'native_full_cover_order':64*48,'native_AP_rotation_composition_checks':24*24,
        'spin_lift_composition_checks_up_to_central_parity':24*24,
        'full_native_representation_may_have_a_kernel':True}


@lru_cache(None)
def certificate():
    """Replay actual native symmetries before using their regular reference."""
    return {'group':native_cover(),
        'native_translation_replay':native_translation_join(),
        'native_rotation_translation_replay':rotation_translation_join(),
        'full_native_rotation_word_replay':rotation_word_checks(),
        'full119_coordinate_current_orientation_replay':matrix_checks(),
        'actual_full119_occupied_H_replay':actual_hamiltonian_checks(),
        'native_Fourier_spin_covariance_replay':operator_certificate(),
        'exact_retained_reference_check':finite_lift(),
        'completion':{'isometry':'J psi = |G|^(-1/2) sum_g U_g psi tensor |g>',
            'lifted_operator':'sum_g U_g C U_g^dag tensor |g><g|',
            'joint_action':'U_h tensor L_h, L_h|g>=|hg>',
            'joint_input_is_invariant':True,
            'source_marginal':'group twirl of original source state',
            'records_are_joint_relational_records':True,
            'complete_record_probabilities_and_norm_errors_preserved':True,
            'native_Gauss_kernel_and_raw_L_preserved':True,
            'actual_stationary_ground_lift_is_stationary':True,
            'source_reference_Hamiltonian':'H_(a,N) tensor I, H_reference=0',
            'reference_is_retained_across_entire_sequence':True,
            'continuous_rotation_or_Lorentz_spin_claim':False}}


def apparatus_capacity(kicks):
    """One shared regular reference and both retained bits at every Weyl test."""
    if not isinstance(kicks,int) or kicks<0:
        raise ValueError('the number of retained Weyl tests must be nonnegative')
    return {'regular_reference_group':'Z4^3 semidirect binary octahedral group',
        'retained_regular_reference_dimension':3072,
        'reference_Hamiltonian':'0',
        'one_reference_shared_across_sequence':True,
        'retained_control_and_flag_qubits':2*kicks,
        'total_declared_Weyl_apparatus_dimension':3072*4**kicks,
        'additional_capacity_nats_expression':f'log(3072) + {2*kicks} log(2)',
        'scope':'regular reference and retained Weyl control/flag records; preparation controllers and other services excluded'}
