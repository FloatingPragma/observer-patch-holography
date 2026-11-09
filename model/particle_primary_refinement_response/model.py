"""Exact local Schur certificates, tiled trials and original response bounds.

The selected Gauss polynomial state and source clock are unchanged. A
capacity-dependent populated band does not imply fixed-energy survival.
"""
from functools import lru_cache
from fractions import Fraction as F

from model.particle_common_contract import contract, provenance
from model.particle_primary_occupation_bound import model as occupation
from model.particle_primary_variational_readout import variational
from model.particle_primary_spectral_scale import scale
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.particle_primary_pair_channels.model import passive_window
from model.particle_primary_mass_readout import readout
from model.particle_spatial_ground_state import ground
from model.spin_packet_readout.bounds import sqrt_upper


SOURCE_INPUTS=(
    'model/particle_primary_refinement_response/__init__.py',
    'model/particle_primary_refinement_response/__main__.py',
    'model/particle_primary_refinement_response/model.py',
    'model/particle_primary_refinement_response/test_model.py',
    'model/particle_primary_refinement_response/PROOF.md',
    'model/particle_primary_occupation_bound/PROOF.md',
    'model/particle_primary_variational_readout/PROOF.md',
    'model/particle_primary_spectral_scale/PROOF.md',
    'model/particle_primary_response_moments/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/particle_primary_pair_channels/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_primary_field_algebra/PROOF.md',
    'model/particle_primary_spatial_algebra/PROOF.md',
    'Lean/QFT/GroundRecordTail.lean',
    'Lean/QFT/ColimitSelectedStateGNS.lean',
    'Lean/QFT/ColimitGNSRepresentation.lean',
)
NATIVE=((3,F(4,5),-4),(2,F(9,20),-3),(6,F(21,20),1),
        (3,F(7,10),2),(1,F(3,10),6))


@lru_cache(maxsize=1)
def native_control():
    actual=tuple((row['dimension'],F(row['trace_metric_Casimir']),row['integer_Q6_charge'])
                 for row in scale.vacuum_channels()['all_native_blocks'])
    if actual!=NATIVE:raise ArithmeticError('native source representations or electric metric changed')
    return actual


def _bilinear(mask,i,j):
    if not mask&(1<<j):return None
    sign=(-1)**((mask&((1<<j)-1)).bit_count());mask^=1<<j
    if mask&(1<<i):return None
    sign*=(-1)**((mask&((1<<i)-1)).bit_count())
    return mask|(1<<i),sign


@lru_cache(maxsize=1)
def q6_integer_blocks():
    """All 256 states, partitioned by the two conserved four-mode numbers."""
    groups={}
    for mask in range(256):
        groups.setdefault(((mask&15).bit_count(),(mask>>4).bit_count()),[]).append(mask)
    answer=[]
    for sector,states in sorted(groups.items()):
        position={s:i for i,s in enumerate(states)};n=len(states)
        D=[[0]*n for _ in states];Q=[[0]*n for _ in states]
        for column,mask in enumerate(states):
            D[column][column]=4
            Q[column][column]=sum(bool(mask&(1<<j)) for j in (0,1,4,5))
            for site in (0,4):
                for j in (0,1):
                    for i,k in ((site+j,site+j+2),(site+j+2,site+j)):
                        value=_bilinear(mask,i,k)
                        if value:D[position[value[0]]][column]+=value[1]
            for a in (0,1):
                for b in (0,1):
                    right=_bilinear(mask,4+b,4+a)
                    left=_bilinear(right[0],a,b) if right else None
                    if left:Q[position[left[0]]][column]-=2*right[1]*left[1]
        if any(D[i][j]!=D[j][i] or Q[i][j]!=Q[j][i]
               for i in range(n) for j in range(n)):
            raise ArithmeticError('local CAR blocks are not symmetric')
        answer.append((sector,tuple(states),tuple(map(tuple,D)),tuple(map(tuple,Q))))
    return tuple(answer)


def q6_positive_shift(onsite,inverse_margin,shift):
    """Exact rational positive pivots for mD-wQ+cI in every number sector."""
    m=ground._q(onsite,'onsite');w=ground._q(inverse_margin,'inverse margin')
    c=ground._q(shift,'shift')
    if min(m,w)<=0:raise ValueError('positive local coefficients required')
    for _,states,D,Q in q6_integer_blocks():
        n=len(states)
        matrix=[[m*D[i][j]-w*Q[i][j]+(c if i==j else 0)
                 for j in range(n)] for i in range(n)]
        for k in range(n):
            pivot=matrix[k][k]
            if pivot<=0:return False
            for i in range(k+1,n):
                for j in range(i,n):
                    matrix[j][i]-=matrix[i][k]*matrix[j][k]/pivot
                    matrix[i][j]=matrix[j][i]
    return True


@lru_cache(maxsize=16)
def exact_q6_loss(onsite,inverse_margin,bits=16):
    """Dyadic full-space lower shift, proved by LDL rather than eigen rounding."""
    m=ground._q(onsite,'onsite');w=ground._q(inverse_margin,'inverse margin')
    bits=ground._i(bits,'dyadic precision',1)
    coarse=occupation.occupation_loss(1,m,w)['occupation_loss_upper']
    scale=1<<bits;lo=0;hi=-(-(coarse*scale).numerator//(coarse*scale).denominator)+1
    if not q6_positive_shift(m,w,F(hi,scale)):
        raise ArithmeticError('initial exact local upper shift failed')
    while hi-lo>1:
        middle=(lo+hi)//2
        if q6_positive_shift(m,w,F(middle,scale)):hi=middle
        else:lo=middle
    upper=F(hi,scale)
    if not q6_positive_shift(m,w,upper):raise ArithmeticError('local final pivots failed')
    return {'occupation_loss_upper':upper,'coarse_occupation_loss_upper':coarse,
        'dyadic_bracket':[F(lo,scale),upper],
        'full_Fock_dimension':256,'number_sector_count':25,
        'largest_exact_block_dimension':36,'all_positive_pivots':True,
        'all_occupations_retained':True,'eigensolver_used_as_certificate':False}


@lru_cache(maxsize=4)
def local_allocation(kind='all'):
    native_control()
    if kind=='all':retained=(F(3,20),)*5
    elif kind=='q6':retained=(F(),)*4+(F(16,25),)
    else:raise ValueError('allocation must be all or q6')
    onsite=tuple((1-r)/2 for r in retained) # degree is bounded by six.
    delta=4*sum((d*(sqrt_upper(m*m+F(1,4))-m)
                 for (d,_,_),m in zip(NATIVE,onsite)),F())
    margin=F(63,4)-delta
    if margin<=0:raise ArithmeticError('full nonconstant link block is not positive')
    rows=[]
    for (d,casimir,charge),m,r in zip(NATIVE,onsite,retained):
        w=1/(F(105,2)*casimir-delta)
        if d==1:loss=exact_q6_loss(m,w)
        else:loss=occupation.occupation_loss(d,m,w)
        rows.append({'native_dimension':d,'integer_Q6_charge':charge,
            'retained_onsite_fraction':r,'onsite_per_bond_endpoint':m,
            'inverse_electric_margin':w,**loss})
    return {'allocation':kind,'maximum_degree_used':6,
        'retained_fractions':retained,'pointwise_full_Fock_shift_upper':delta,
        'nonconstant_link_margin_lower':margin,'native_blocks':rows,
        'per_bond_occupation_loss_upper':sum((r['occupation_loss_upper'] for r in rows),F())}


@lru_cache(maxsize=1)
def tile_seed():
    zero=contract.member(0);receipt=variational.bind_member(zero)
    variational.verify_bound(zero,receipt)
    upper=F(receipt['result']['ground_energy_upper'])
    return {'member_id':zero['member_id'],'producer_result_sha256':receipt['result_sha256'],
        'dimensionless_tile_energy_upper':upper,'tile_cells':8,
        'dimensionless_gain_per_cell':(upper+720)/8,
        'tile_word_degree':1,'exact_full_physical_trial':True}


def tiled_upper(member):
    contract.verify_member(member)
    side=member['source']['side'];a=F(member['source']['spacing']);V=side**3
    if side%2:raise ValueError('the source box must admit two-cell cubic tiles')
    seed=tile_seed();tiles=V//8
    return {'seed':seed,'tile_count':tiles,'dimensionless_ground_upper':
        tiles*seed['dimensionless_tile_energy_upper'],
        'physical_ground_upper':tiles*seed['dimensionless_tile_energy_upper']/a,
        'intertile_links_constant':True,'intertile_hopping_expectation':'0',
        'intertile_plaquette_expectation':'0','trial_replaces_selected_state':False}


def _defects(member,kind):
    source=member['source'];a=F(source['spacing']);V=source['side']**3
    local=local_allocation(kind);tile=tiled_upper(member)
    numerator=tile['dimensionless_ground_upper']+90*V
    numerator+=source['interior_hopping_bonds']*local['per_bond_occupation_loss_upper']
    numerator+=F(source['interior_plaquettes'],7)+a*F(member['state']['target_ground_excess'])
    coefficient=F(9,20) if kind=='all' else F(48,25)
    if numerator<0:raise ArithmeticError('selected energy bounds are inconsistent')
    upper=min(F(60*V if kind=='all' else 4*V),numerator/coefficient)
    return {'allocation':local,'defect_count_upper':upper,
        'defect_density_upper':upper/V,'retained_defect_coefficient':coefficient}


def uniform_density_bound():
    local=local_allocation('all');seed=tile_seed()
    # B/V<=3, P/V<=3 and a*epsilon/V<=1/128 at every primary level.
    value=(seed['dimensionless_gain_per_cell']+3*local['per_bond_occupation_loss_upper']
           +F(3,7)+F(1,128))/F(9,20)
    if not value<F(97,10):raise ArithmeticError('uniform defect target was not proved')
    return value


def selected_response(member):
    contract.verify_member(member)
    source=member['source'];V=source['side']**3;a=F(source['spacing'])
    eps=F(member['state']['target_ground_excess']);beta=F(member['state']['beta'])
    all_defects=_defects(member,'all');q6_defects=_defects(member,'q6')
    uniform=uniform_density_bound()
    chi=max(F(),(12-all_defects['defect_density_upper'])/42)
    if chi<F(1,20):raise ArithmeticError('collective commutator floor failed')
    lower=1/(100*beta);upper=F(1446)/a
    window=passive_window(F(1,20),beta,6/a,(F(144)/a)**2,lower,upper)
    weight=window['original_transfer_weight_lower']
    if weight!=F(3,100):raise ArithmeticError('original passive window arithmetic failed')
    final=readout.excitation_window((lower,upper),weight,eps,100*eps)
    q6_weight=max(F(),1-q6_defects['defect_density_upper']/2)
    q6_window=(passive_window(q6_weight,beta,6/a,(F(24)/a)**2,lower,F(246)/a)
               if q6_weight>F(1,50) else None)
    clock=spectroscopy.clock_plan(member,lower/4,'1/400','1/400')
    clock_loss=F(1,200)+F(clock['two_complete_clock_error_upper'])
    return contract.serial({'schema':'oph.primary_refinement_response_result.v1',
        'member_id':member['member_id'],**member['ids'],
        'tiled_ground_upper':tiled_upper(member),
        'all_block_defects':all_defects,'q6_only_defects':q6_defects,
        'all_block_defect_density_upper':all_defects['defect_density_upper'],
        'uniform_defect_density_upper':uniform,'uniform_simple_defect_density_upper':F(97,10),
        'q6_averaged_matched_pair_weight_lower':q6_weight,
        'q6_scalar_original_transfer_window':q6_window,
        'native_collective_block_dimension':6,'native_collective_Q6_charge':1,
        'collective_operator_gauge_charge':0,'collective_normalization_squared':42,
        'collective_commutator_lower':chi,'uniform_collective_commutator_lower':F(1,20),
        'collective_probe':'C_x=sum_(alpha=1,2;i=1..6) c^dagger_(+,alpha,i) c_(-,alpha,i)/sqrt(42)',
        'complete_instrument_success':'K_x=C_x/sqrt(V)',
        'complete_instrument_failure':'sqrt(I-sum_x C_x^dagger C_x/V)',
        'complete_probe_outcomes':V+1,'all_failures_and_reference_retained':True,
        'selected_input_state_unchanged':True,'probe_preserves_input_density_claimed':False,
        'source_or_clock_changed':False,'exact_Gibbs_replacement':False,
        'source_log_weight_slope_upper':beta,'centered_moment_center':6/a,
        'centered_second_moment_upper':(F(144)/a)**2,
        'collective_original_transfer_band':[lower,upper],
        'collective_original_transfer_weight_lower':weight,
        'passivity_window':window,
        'final_excitation_band':[lower,upper+100*eps],
        'final_excitation_weight_lower':F(1,50),
        'transfer_to_excitation':final,
        'two_clock_plan':clock,'complete_fine_outcomes':(V+1)*clock['clock_plan']['clock_size']**2,
        'finite_record_transfer_band':[lower/2,upper+lower/2],
        'finite_record_transfer_weight_lower':weight-clock_loss,
        'finite_record_simple_weight_lower':F(1,50),
        'clock_histogram_executed':False,'full_density_matrix_executed':False,
        'joint_reference_unchanged_by_initial_energy_record':False,
        'uniform_original_weight_proved':True,
        'fixed_finite_energy_window_uniform_weight_proved':False,
        'isolated_particle_branch_proved':False,'physical_mass_table_GeV':[]})


def bind_member(member):
    value=selected_response(member)
    return {'schema':'oph.primary_refinement_response_bound.v1','result':value,
            'result_sha256':contract.digest(value)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('refinement response or its primary input changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*variational.SOURCE_INPUTS,*occupation.SOURCE_INPUTS))
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_refinement_response_certificate.v1',
        'members':[bind_member(contract.member(n)) for n in (0,1,2)],
        'uniform_defect_density_upper':str(uniform_density_bound()),
        'source_sha256':pins,'source_state_clock_reference_changed':False,
        'fixed_energy_continuum_response_proved':False,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('refinement source changed during replay')
    return result
