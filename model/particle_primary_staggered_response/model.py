"""A source Clifford cancellation sharpens one actual coherent energy band."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product, permutations

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_branch_resolution import model as branch
from model.particle_primary_odd_residual import model as odd

SOURCE_INPUTS=(
    'model/particle_primary_staggered_response/__init__.py',
    'model/particle_primary_staggered_response/__main__.py',
    'model/particle_primary_staggered_response/model.py',
    'model/particle_primary_staggered_response/test_model.py',
    'model/particle_primary_staggered_response/README.md',
    'model/particle_primary_staggered_response/PROOF.md',
    'model/particle_primary_branch_resolution/PROOF.md',
    'model/particle_primary_odd_residual/PROOF.md',
    'model/particle_primary_pair_channels/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'Lean/QFT/GroundRecordTail.lean',
)


@lru_cache(maxsize=1)
def clifford_control():
    pairs=odd.response.beta_pairs()
    chart=s.Matrix.hstack(pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1])
    gamma=odd.response.density.gamma_data();beta=gamma[3]
    raising=s.zeros(4);raising[:2,2:]=s.eye(2)
    flip=s.I*(raising-raising.H);rows=[]
    for axis in range(3):
        hopping=branch.apparatus.clean(chart.H*((-beta+beta*gamma[axis])/2)*chart)
        # On opposite-parity endpoints the one-body commutator uses an
        # anticommutator, with an irrelevant overall endpoint sign.
        residual=branch.apparatus.clean(hopping*raising+raising*hopping)
        if residual[:2,2:]!=s.zeros(2) or residual[2:,:2]!=s.zeros(2):
            raise ArithmeticError('the checkerboard pair-creation block did not cancel')
        if branch.apparatus.clean(hopping*flip+flip*hopping)!=s.zeros(4):
            raise ArithmeticError('the exact staggered Hermitian hopping symmetry failed')
        if branch.apparatus.clean(hopping[:2,2:].H*hopping[:2,2:])!=s.eye(2)/4:
            raise ArithmeticError('the actual off-diagonal Wilson block norm changed')
        rows.append({'axis':axis,'pair_creation_block_zero':True,
                     'Hermitian_flip_hopping_anticommutator_zero':True})
    rotations=[]
    for perm in permutations(range(3)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        for signs in product((-1,1),repeat=3):
            if sign*s.prod(signs)!=1:continue
            ratios=set()
            for x in product((-1,1),repeat=3):
                y=tuple(signs[i]*x[perm[i]] for i in range(3))
                ratios.add(int(s.prod(y)/s.prod(x)))
            if ratios!={sign}:raise ArithmeticError('the actual cube character changed')
            rotations.append(sign)
    return {'actual_Wilson_axis_controls':rows,'proper_cube_rotation_count':len(rotations),
        'proper_rotation_character':'sign of coordinate permutation (A2)',
        'mask_inversion_parity':-1,'local_scalar_inversion_parity':-1,
        'combined_operator_type':'A2g','CAR_parity':'even',
        'native_integer_charge':6,'native_dimension':1,
        'quasiparticle_modes':32,'maximum_degree':3,
        'residual_one_body_norm_upper':F(3,2),
        'residual_squared_form_upper':'36*D6',
        'coherent_normalization_squared':512,'probe_norm_squared':F(9,64),
        'finite_spatial_mask':'(-1)^(x1+x2+x3)',
        'mask_is_conserved_physical_momentum':False}


def selected_bound(member):
    branch._member_zero(member);control=clifford_control()
    upstream=branch.bind_member(member);branch.verify_bound(member,upstream)
    D=F(upstream['result']['q6_selected_defects']['q6_defect_count_upper'])
    if D>=F(46,25):raise ArithmeticError('selected native-six defect estimate failed')
    chi=(16-D)/512;moment=F(9,128)*D
    weight=chi-moment/9
    if weight<=F(17,1280):raise ArithmeticError('checkerboard original band lost its floor')
    norm2=F(9,64);eps=F(member['state']['target_ground_excess'])
    final=branch.readout.excitation_window((3,9),weight,eps,1,norm2)
    plan=branch.spectroscopy.clock_plan(member,'1/64','1/32768','1/32768')
    loss=2*F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['two_complete_clock_error_upper'])
    if loss>F(1,8192):raise ArithmeticError('complete clock loss exceeded allowance')
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'source_control':control,'occupation_producer_sha256':contract.digest(upstream),
        'q6_defect_expectation_upper':D,'original_success_weight_lower':chi,
        'source_frequency_center':6,'selected_centered_second_moment_upper':moment,
        'selected_centered_second_moment_simple_upper':F(207,1600),
        'original_transfer_band':[3,9],'original_transfer_weight_lower':weight,
        'original_transfer_weight_simple_lower':F(17,1280),
        'final_excitation':final,'final_excitation_weight_simple_lower':F(23,5120),
        'complete_two_clock_plan':plan,'complete_clock_event_loss_upper':loss,
        'reported_transfer_interval':[F(95,32),F(289,32)],
        'reported_transfer_original_weight_lower':weight-loss,
        'reported_transfer_simple_lower':F(539,40960),
        'complete_fine_outcomes':2*plan['clock_plan']['clock_size']**2,
        'Hermitian_quadrature_normalization':16,
        'Hermitian_quadrature_signed_energy_mean_lower':F(3,128)*(16-D),
        'Hermitian_quadrature_signed_energy_mean_simple_lower':F(531,1600),
        'all_failures_and_reference_retained':True,'selected_state_changed':False,
        'large_primary_histogram_executed':False,'final_state_irrep_fixed_by_operator':False,
        'dispersion_or_rest_mass_identified':False,'physical_mass_table_GeV':[]})


def bind_member(member):
    result=selected_bound(member)
    return {'schema':'oph.primary_staggered_response.v1','result':result,
            'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('staggered response or source binding changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*odd.SOURCE_INPUTS));pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_staggered_response_certificate.v1',
        'selected_member_bound':bind_member(contract.member(0)),'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('staggered inputs changed during replay')
    return result
