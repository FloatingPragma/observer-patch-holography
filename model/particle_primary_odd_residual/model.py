"""Exact Haar/CAR second actions on the physical scalar-odd trial."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_branch_resolution import model as branch
from model.particle_primary_residual_refinement import model as moments
from model.particle_primary_response_moments import model as response
from model.particle_primary_variational_readout import variational
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
    'model/particle_primary_odd_residual/__init__.py',
    'model/particle_primary_odd_residual/__main__.py',
    'model/particle_primary_odd_residual/model.py',
    'model/particle_primary_odd_residual/test_model.py',
    'model/particle_primary_odd_residual/PROOF.md',
    'model/particle_primary_branch_resolution/PROOF.md',
    'model/particle_primary_residual_refinement/PROOF.md',
    'model/particle_primary_variational_readout/PROOF.md',
    'model/particle_primary_pair_channels/PROOF.md',
    'model/particle_primary_response_moments/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'model/particle_common_contract/PROOF.md',
    'Lean/QFT/GroundRecordTail.lean',
)


def _rotate(value,phase):
    a,b=value
    return ((a,b),(-b,a),(-a,-b),(b,-a))[phase%4]


def _add(target,key,value):
    a,b=target.get(key,(0,0));c,d=value
    target[key]=(a+c,b+d)


def _dot(left,right):
    a,b=left;c,d=right
    return (a*c+b*d,a*d-b*c)


def _car(mask,i,j):
    if not mask&(1<<j):return None
    sign=(mask&((1<<j)-1)).bit_count()%2
    mask^=1<<j
    if mask&(1<<i):return None
    sign^=(mask&((1<<i)-1)).bit_count()%2
    return mask|(1<<i),2*sign


def _act(vector,terms):
    out={}
    for (mask,character),amplitude in vector.items():
        for i,j,shift,phase in terms:
            moved=_car(mask,i,j)
            if moved:
                _add(out,(moved[0],character+shift),_rotate(amplitude,phase+moved[1]))
    return {key:value for key,value in out.items() if value!=(0,0)}


def _haar_inner(left,right):
    answer=(0,0)
    for key,value in left.items():
        if key in right:
            a,b=answer;c,d=_dot(value,right[key]);answer=(a+c,b+d)
    return answer


def _pointwise_inner(left,right):
    """Laurent coefficient table, with all matter occupations retained."""
    left_groups=defaultdict(list);right_groups=defaultdict(list)
    for (mask,e),value in left.items():left_groups[mask].append((e,value))
    for (mask,e),value in right.items():right_groups[mask].append((e,value))
    out={}
    for mask,rows in left_groups.items():
        for e,value in rows:
            for f,other in right_groups.get(mask,()):_add(out,f-e,_dot(value,other))
    return {key:value for key,value in out.items() if value!=(0,0)}


@lru_cache(maxsize=1)
def q6_vectors():
    """Integer Gaussian amplitudes for all32 CAR modes and12 link characters.

    phi is multiplied by4 and every hopping action by2. Balanced base5
    encodes the exponents after at most two hopping actions injectively.
    """
    sites=list(product(range(2),repeat=3));where={x:i for i,x in enumerate(sites)}
    pairs=response.beta_pairs()
    chart=s.Matrix.hstack(pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1])
    gamma=response.density.gamma_data();beta=gamma[3]
    phases={s.S.One:0,s.I:1,-s.S.One:2,-s.I:3};spin=[]
    for axis in range(3):
        h=branch.apparatus.clean(chart.H*(-beta+beta*gamma[axis])*chart)
        if any(x not in phases and x!=0 for x in h):
            raise ArithmeticError('the actual hopping entries lost Gaussian integer form')
        spin.append(h)
    terms=[];edges=[]
    for x in sites:
        for axis in range(3):
            if x[axis]:continue
            y=list(x);y[axis]=1;y=tuple(y);xid=where[x];yid=where[y]
            shift=5**len(edges);edges.append((xid,yid,axis))
            for i,j in product(range(4),repeat=2):
                if spin[axis][i,j]:
                    phase=phases[spin[axis][i,j]]
                    terms.extend(((4*xid+i,4*yid+j,shift,phase),
                                  (4*yid+j,4*xid+i,-shift,(-phase)%4)))
    omega=sum(1<<(4*x+j) for x in range(8) for j in (2,3))
    phi={}
    for x,j in product(range(8),range(2)):
        moved=_car(omega,4*x+j,4*x+2+j)
        _add(phi,(moved[0],0),_rotate((1,0),moved[1]))
    first=_act(phi,terms);parts={2:{},4:{}}
    positive=sum(1<<(4*x+j) for x in range(8) for j in (0,1))
    negative=omega
    for key,value in first.items():
        mask=key[0];defects=(mask&positive).bit_count()+16-(mask&negative).bit_count()
        if defects not in parts:raise ArithmeticError('odd first action left its exact packet sectors')
        parts[defects][key]=value
    second={k:_act(vector,terms) for k,vector in parts.items()}
    return phi,parts,second,tuple(terms),tuple(edges)


@lru_cache(maxsize=1)
def q6_second_action_control():
    phi,parts,second,terms,edges=q6_vectors()
    first=[];matrix=[]
    if _pointwise_inner(phi,phi)!={0:(16,0)}:
        raise ArithmeticError('scalar-odd source trial normalization changed')
    for k,wanted in ((2,192),(4,672)):
        polynomial=_pointwise_inner(parts[k],parts[k])
        if polynomial!={0:(wanted,0)} or _pointwise_inner(phi,parts[k]):
            raise ArithmeticError('the odd first-action norm is not its pointwise constant')
        first.append(F(wanted,64))
        row=[]
        for j in (2,4):
            a,b=_haar_inner(second[k],second[j])
            if b:raise ArithmeticError('odd second-action Gram has an imaginary coefficient')
            row.append(F(a,256))
        matrix.append(row)
    return {'full_q6_CAR_modes':32,'interior_bonds':len(edges),
        'independent_Haar_link_characters':12,'exact_integer_Gaussian_arithmetic':True,
        'first_action_defect_sectors':[2,4],
        'pointwise_first_action_norm_squared':first,
        'parent_second_action_Gram':matrix,
        'first_action_nonzero_occupation_character_terms':[len(parts[k]) for k in (2,4)],
        'second_action_nonzero_occupation_character_terms':[len(second[k]) for k in (2,4)],
        'full_spatial_Hilbert_matrix_allocated':False}


def odd_parent_moment(hopping_weights,loop_weights):
    """Actual member-zero parent norm, including every omitted pair sector."""
    th=tuple(hopping_weights);tl=tuple(loop_weights)
    if len(th)!=6 or len(tl)!=5 or any(type(x) not in (int,str,F) for x in th+tl):
        raise ValueError('six hopping and five loop exact rational weights required')
    th=tuple(map(F,th));tl=tuple(map(F,tl))
    control=q6_second_action_control();native=branch.refinement.native_control()
    dimensions=[d for d,_,q in native if q!=6];B=12
    geometry=moments.geometry(2);G=geometry['hopping_connected_cumulant_coefficient']
    zother=B*sum(dimensions)
    sother=B*sum(d*t*t for d,t in zip(dimensions,th[:4]))
    qother=B*sum(d*t for d,t in zip(dimensions,th[:4]))
    other_moment=zother*sother+2*qother*qother-G*sum(d*t*t for d,t in zip(dimensions,th[:4]))
    n=control['pointwise_first_action_norm_squared'];gram=control['parent_second_action_Gram']
    z6=sum(n);s6=sum(v*t*t for v,t in zip(n,th[4:]));q6=sum(v*t for v,t in zip(n,th[4:]))
    q6moment=sum(th[4+i]*gram[i][j]*th[4+j] for i,j in product(range(2),repeat=2))
    hopping=other_moment+q6moment+zother*s6+sother*z6+4*qother*q6
    loop=moments.mixed_parent_moment(2,[0]*5,tl)
    z=zother+z6;square=sother+s6;q=qother+q6
    mixed=hopping+z*loop['weighted_magnetic_variance']+square*loop['magnetic_variance']
    mixed+=4*q*loop['magnetic_covariance']+loop['exact_magnetic_mixed_fourth']
    return {'exact_other_native_hopping_fourth':other_moment,
        'exact_q6_second_action_norm_squared':q6moment,
        'exact_hopping_second_action_norm_squared':hopping,
        'exact_magnetic_mixed_fourth':loop['exact_magnetic_mixed_fourth'],
        'pointwise_hopping_variance':z,'pointwise_weighted_hopping_variance':square,
        'pointwise_hopping_covariance':q,
        'full_weighted_variance':square+loop['weighted_magnetic_variance'],
        'full_covariance':q+loop['magnetic_covariance'],
        'exact_parent_second_action_norm_squared':mixed}


def refined_bound(member,producer=None):
    branch._member_zero(member)
    producer=branch.bind_member(member) if producer is None else producer
    branch.verify_bound(member,producer)
    inherited=producer['result']['odd_endpoint_trial']
    if [row['Q6'] for row in inherited['hopping_packets']]!=[-4,-3,1,2,6,6]:
        raise ArithmeticError('the inherited twelve-packet channel order changed')
    th=list(map(F,inherited['rational_hopping_coefficients']))
    tl=list(map(F,inherited['rational_loop_coefficients']))
    moment=odd_parent_moment(th,tl)
    square=moment['full_weighted_variance'];q=moment['full_covariance']
    if 1+square!=F(inherited['trial_norm_squared']):
        raise ArithmeticError('the fixed scalar-odd trial normalization changed')
    correction=variational.character_control()['cubic_loop_singlet_matrix']
    loop_projection=[-sum(F(correction[i][j],420)*tl[j] for j in range(5)) for i in range(5)]
    g=F(inherited['loop_packet_norm_squared'])
    outside=moment['exact_parent_second_action_norm_squared']-q*q-g*sum(v*v for v in loop_projection)
    e=F(inherited['trial_energy_shift_relative_to_anchor_pair'])
    dlo=-F(inherited['relative_Ritz_root_bracket'][1])
    residual2=((q+e)**2+(dlo+e)**2*square+outside)/(1+square)
    if min(outside,residual2)<=0:raise ArithmeticError('the full omitted moment is not positive')
    radius=sqrt_upper(residual2);center=F(inherited['odd_sector_absolute_ground_upper'])
    ground_upper=F(inherited['ground_Ritz_upper_for_comparison'])
    ground_lower=F(inherited['full_physical_ground_lower'])
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'inherited_trial_receipt_sha256':contract.digest(producer),
        'exact_parent_moments':moment,'parent_omitted_norm_squared':outside,
        'source_trial_energy':center,'source_residual_squared_upper':residual2,
        'source_residual_upper':radius,
        'inherited_source_residual_upper':F(inherited['dressed_trial_full_complement_residual_upper']),
        'some_physical_A1u_eigenvalue_interval':[center-radius,center+radius],
        'full_physical_ground_upper':ground_upper,'full_physical_ground_lower':ground_lower,
        'eigenvalue_interval_lower_minus_ground_upper':center-radius-ground_upper,
        'some_A1u_excitation_interval':[max(F(),center-radius-ground_upper),center+radius-ground_lower],
        'source_state_reference_clock_unchanged':True,
        'trial_replaces_selected_polynomial_state':False,
        'full_spatial_Hilbert_matrix_allocated':False,
        'final_product_compression_used_contractively':True,
        'source_residual_squared_is_exact_after_cutoff':False,
        'isolated_eigenvalue_or_particle_mass_proved':False,
        'selected_probe_weight_in_this_eigenvalue_interval_proved':False})


@lru_cache(maxsize=1)
def coherent_commutator_control():
    """Actual beta-chart blocks of the global scalar-pair commutator."""
    pairs=response.beta_pairs();chart=s.Matrix.hstack(pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1])
    gamma=response.density.gamma_data();beta=gamma[3]
    raising=s.zeros(4);raising[:2,2:]=s.eye(2);rows=[]
    for axis in range(3):
        h=branch.apparatus.clean(chart.H*((-beta+beta*gamma[axis])/2)*chart)
        actual=branch.apparatus.clean(h*raising-raising*h)
        if h[:2,:2]!=-s.eye(2)/2 or h[2:,2:]!=s.eye(2)/2:
            raise ArithmeticError('actual diagonal hopping blocks changed')
        Bstar=h[2:,:2]
        expected=(-Bstar).row_join(-s.eye(2)).col_join(s.zeros(2).row_join(Bstar))
        if actual!=expected or branch.apparatus.clean(Bstar.H*Bstar)!=s.eye(2)/4:
            raise ArithmeticError('scalar-pair creation/number commutator split changed')
        rows.append({'axis':axis,'pair_creation_block':str(actual[:2,2:]),
                     'number_block_operator_norm_squared':F(1,4)})
    return {'actual_axis_controls':rows,'q6_quasiparticle_modes':32,
        'coherent_pair_count':16,'interior_bonds':12,'maximum_interior_degree':3,
        'pair_creation_coefficient_Hilbert_Schmidt_norm_squared':48,
        'number_preserving_one_body_norm_upper':F(3,2),
        'number_preserving_one_body_trace':0,
        'pair_creation_squared_form_upper':'48*(1+D6/2)',
        'number_preserving_squared_form_upper':'36*D6',
        'Young_parameter':F(6,7),'probe_norm_squared':F(9,64),
        'probe_normalization_squared':512}


def selected_coherent_response(member,producer=None):
    """Original scalar-operator response, distinct from the trial's measure."""
    branch._member_zero(member)
    producer=branch.bind_member(member) if producer is None else producer
    branch.verify_bound(member,producer)
    control=coherent_commutator_control()
    D=F(producer['result']['q6_selected_defects']['q6_defect_count_upper'])
    if not D<F(46,25):raise ArithmeticError('selected coherent-response defect bound failed')
    commutator=(16-D)/512;moment=(624+858*D)/3584
    lower=F(3,4);upper=F(45,4)
    weight=commutator-moment/F(21,4)**2
    simple_moment=F(13767,22400);simple_weight=F(7057,1317120)
    if not moment<simple_moment or not weight>simple_weight>F(5,1024):
        raise ArithmeticError('selected coherent positive-band estimate failed')
    norm2=F(control['probe_norm_squared']);eps=F(member['state']['target_ground_excess'])
    transport=branch.readout.excitation_window((lower,upper),weight,eps,2,norm2)
    loss=F(transport['initial_excitation_tail_weight_upper'])
    final=F(transport['original_final_excitation_weight_lower'])
    if not final>F(10151,10536960)>F(1,2048):
        raise ArithmeticError('coherent final-excitation transport failed')
    plan=branch.spectroscopy.clock_plan(member,'1/64','1/32768','1/32768')
    clockloss=2*F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['two_complete_clock_error_upper'])
    if clockloss>F(1,8192):raise ArithmeticError('coherent two-clock event budget failed')
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'inherited_trial_receipt_sha256':contract.digest(producer),
        'probe':'K00=sum_(x,alpha) B_(x,alpha)/(16*sqrt(2))',
        'probe_CAR_parity':'even','operator_cubic_irrep':'A1','operator_inversion_parity':-1,
        'probe_norm_squared':norm2,'original_success_mass_upper':norm2,
        'q6_defect_expectation_upper':D,'original_success_mass_lower':commutator,
        'source_frequency_center':6,'selected_centered_second_moment_upper':moment,
        'selected_centered_second_moment_simple_upper':simple_moment,
        'original_transfer_interval':[lower,upper],'original_transfer_weight_lower':weight,
        'original_transfer_weight_rational_lower':simple_weight,
        'original_transfer_weight_simple_lower':F(5,1024),
        'initial_excitation_cut':2,'cooling_original_weight_loss_upper':loss,
        'same_source_excitation_transport':transport,
        'original_final_excitation_interval':[lower,upper+2],
        'original_final_excitation_weight_lower':final,
        'original_final_excitation_weight_rational_lower':F(10151,10536960),
        'original_final_excitation_weight_simple_lower':F(1,2048),
        'complete_probe_outcomes':2,'complete_two_clock_plan':plan,
        'finite_reported_transfer_interval':[lower-F(1,32),upper+F(1,32)],
        'finite_reported_original_weight_lower':weight-clockloss,
        'finite_reported_original_weight_simple_lower':F(39,8192),
        'all_failure_records_and_reference_retained':True,
        'probabilities_conditioned_on_success':False,
        'large_primary_histogram_executed':False,
        'selected_state_replaced_by_odd_trial':False,
        'final_state_cubic_irrep_is_fixed_by_operator_irrep':False,
        'response_attached_to_the_trial_eigenvalue':False,
        'isolated_line_or_rest_mass_proved':False})


def bind_member(member):
    branch._member_zero(member);producer=branch.bind_member(member)
    branch.verify_bound(member,producer)
    result={'member_id':member['member_id'],**member['ids'],
        'odd_trial_residual':refined_bound(member,producer),
        'selected_coherent_response':selected_coherent_response(member,producer),
        'physical_mass_table_GeV':[]}
    return {'schema':'oph.primary_odd_residual.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('source odd residual or coherent response changed')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*branch.SOURCE_INPUTS,*moments.SOURCE_INPUTS))
    if any(not (contract.ROOT/path).is_file() for path in paths):raise ValueError('an odd residual source input is missing')
    pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_odd_residual_certificate.v1',
        'q6_exact_second_action':q6_second_action_control(),
        'coherent_commutator':coherent_commutator_control(),
        'selected_member_bound':bind_member(contract.member(0)),
        'source_sha256':pins,'physical_mass_table_GeV':[]}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('odd residual input changed during replay')
    return contract.serial(result)
