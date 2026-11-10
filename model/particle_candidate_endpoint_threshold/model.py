"""Original-state separated endpoint bounds with complete finite records."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import isqrt

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_candidate_charged import model as dressed
from model.particle_candidate_charged_propagation import model as filtered
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
 'model/particle_candidate_endpoint_threshold/__init__.py',
 'model/particle_candidate_endpoint_threshold/__main__.py',
 'model/particle_candidate_endpoint_threshold/model.py',
 'model/particle_candidate_endpoint_threshold/test_model.py',
 'model/particle_candidate_endpoint_threshold/PROOF.md',
 'model/particle_candidate_endpoint_threshold/README.md',
 *dressed.SOURCE_INPUTS,*filtered.SOURCE_INPUTS,*spectroscopy.SOURCE_INPUTS,
)


def edge(x,y):
    return tuple(sorted((tuple(x),tuple(y))))


def paths(length):
    if type(length) is not int or length not in (1,2,3):
        raise ValueError('only shortest member-zero cube lengths one through three are admitted')
    out=[]
    for x in product(range(2),repeat=3):
        for axes in permutations(range(3),length):
            nodes=[x]
            for axis in axes:
                y=list(nodes[-1]);y[axis]^=1;nodes.append(tuple(y))
            out.append(tuple(nodes))
    return tuple(out)


def plaquettes():
    out=[]
    for axes in combinations(range(3),2):
        fixed=next(i for i in range(3) if i not in axes)
        for bit in (0,1):
            x=[0,0,0];x[fixed]=bit;nodes=[tuple(x)]
            for axis in (axes[0],axes[1],axes[0],axes[1]):
                x=x.copy();x[axis]^=1;nodes.append(tuple(x))
            out.append(frozenset(edge(a,b) for a,b in zip(nodes,nodes[1:])))
    return tuple(out)


def geometry(length):
    rows=paths(length);N=len(rows)
    starts=Counter(p[0] for p in rows);ends=Counter(p[-1] for p in rows)
    links=Counter(edge(a,b) for p in rows for a,b in zip(p,p[1:]))
    touched=tuple(sum(bool(set(edge(a,b) for a,b in zip(p,p[1:]))&q) for q in plaquettes()) for p in rows)
    if set(starts.values())!={N//8} or set(ends.values())!={N//8}:
        raise ArithmeticError('endpoint uniformity failed')
    if len(links)!=12 or set(links.values())!={N*length//12} or set(touched)!={length+1}:
        raise ArithmeticError('link or plaquette incidence differs')
    return {'path_length':length,'oriented_paths':rows,'path_count':N,
        'occurrences_per_endpoint':N//8,'occurrences_per_link':N*length//12,
        'touched_plaquettes':length+1,'coordinate_separation_squared':length,
        'per_link_character_degree':1,'source_per_link_word_cutoff':1,
        'longer_paths_or_larger_capacity_bound':False}


@lru_cache(maxsize=3)
def gauss_chain_control(length):
    paths(length)
    c=tuple(dressed.regional.creator(j,2) for j in range(2));n=tuple(x*x.H for x in c)
    I=s.eye(3);E=s.diag(-1,0,1);P=s.diag(0,1,0);L=s.zeros(3);L[2,1]=1
    dim=4*3**length;identity=s.eye(dim)
    flux=[]
    for i in range(length):
        factors=[I]*length+[s.eye(4)];factors[i]=E
        flux.append(s.kronecker_product(*factors))
    Nx=s.kronecker_product(*([I]*length+[n[0]]))
    Ny=s.kronecker_product(*([I]*length+[n[1]]))
    B=s.kronecker_product(*([L]*length+[c[0]*c[1].H]))
    generators=[6*(Nx-flux[0])]+[6*(flux[i-1]-flux[i]) for i in range(1,length)]+[6*(Ny+flux[-1]-identity)]
    H0=3*(Nx-Ny)+s.Rational(63,4)*sum((x*x for x in flux),s.zeros(dim))
    center=s.Rational(6)+s.Rational(63,4)*length
    if any(G*B!=B*G for G in generators) or H0*B-B*H0!=center*B:
        raise ArithmeticError('Gauss or character energy identity failed')
    expected=s.kronecker_product(*([P]*length+[n[1]*(s.eye(4)-n[0])]))
    if B.H*B!=expected:raise ArithmeticError('complete link projector missing')
    anchor=s.kronecker_product(*([s.eye(3)[:,1]]*length+[s.eye(4)[:,2]]))
    target=B*anchor
    if (target.H*target)[0]!=1 or any(G*v!=s.zeros(dim,1) for G in generators for v in (anchor,target)):
        raise ArithmeticError('physical control support failed')
    return {'control_dimension':dim,'Gauss_generators':length+1,
        'full_control_occupations_checked':True,'endpoint_Q6_increments':[6,-6],
        'interior_divergence_zero':True,'onsite_electric_commutator_center':F(center),
        'original_source_histogram_executed':False,'center_is_endpoint_mass':False}


@lru_cache(maxsize=1)
def endpoint_columns_control():
    pairs=dressed.regional.pairs.response.beta_pairs()
    positive=s.Matrix.hstack(*(p for p,m in pairs));negative=s.Matrix.hstack(*(m for p,m in pairs))
    gamma=dressed.regional.pairs.response.density.gamma_data();beta=gamma[3]
    hops=tuple((-beta+beta*g)/2 for g in gamma[:3]);cases=0
    for length in (1,2,3):
        for path in paths(length):
            selected=set(edge(a,b) for a,b in zip(path,path[1:]))
            for endpoint,other,modes in ((path[0],path[-1],positive),(path[-1],path[0],negative)):
                gram=s.zeros(2);neighbors=[]
                for axis in range(3):
                    y=list(endpoint);y[axis]^=1;y=tuple(y)
                    if edge(endpoint,y) in selected:continue
                    T=hops[axis] if endpoint[axis]==0 else hops[axis].H
                    image=T*modes;gram+=image.H*image;neighbors.append(y)
                if dressed.clean(gram)!=s.eye(2) or other in neighbors or len(neighbors)!=2:
                    raise ArithmeticError('off-path endpoint image is not an independent isometry')
                cases+=1
    # Each endpoint term is unitarily equivalent to two independent raising
    # pairs, whose coefficient singular values are both 1/sqrt(2).
    C=tuple(dressed.regional.creator(j,4) for j in range(4))
    raise_two=(C[0]*C[2].H+C[1]*C[3].H)/s.sqrt(2)
    maximum=max(dressed.clean(raise_two.H*raise_two).eigenvals())
    if maximum!=1:raise ArithmeticError('two-pair full-Fock bound differs')
    return {'literal_endpoint_column_cases':cases,'each_endpoint_Gram':'I_2',
        'each_endpoint_full_Fock_norm_upper':1,'summed_off_path_norm_upper':2,
        'opposite_endpoint_in_image':False,'scope':'shortest cube paths'}


@lru_cache(maxsize=1)
def controlled_increment_meter_control():
    """Charge-change labels inside temporal controls preserve interference."""
    meter=dressed.regional.meter
    K=s.Matrix([[0,0],[1,0]]);failure=s.diag(0,1)
    result=meter.coherent_meter((0,6),(K,failure),(6,0),5)
    W=result['coupling'];shift=meter.cyclic_shift(5,6)
    X=s.Matrix([[0,1],[1,0]]);I=s.eye(5)
    complete=s.zeros(10);filtered_row=s.zeros(2)
    for U in (s.eye(2),X):
        V=s.kronecker_product(U,I)
        complete+=V.H*W*s.kronecker_product(K,I)*W.H*V/2
        filtered_row+=U.H*K*U/2
    if complete!=s.kronecker_product(filtered_row,shift):
        raise ArithmeticError('temporal meter leaked the time node')
    outside=W*s.kronecker_product(filtered_row,I)*W.H
    if outside==complete:raise ArithmeticError('control failed to distinguish meter placement')
    return {'finite_control_source_dimension':2,'meter_radix':5,
        'noncommuting_temporal_controls_checked':2,
        'meter_inside_each_instantaneous_row':True,
        'temporal_interference_preserved':True,
        'moving_meter_outside_temporal_conjugation_is_equivalent':False,
        'source_meter_cells':8,'all_source_charge_meter_records':5**8,
        'absolute_charge_measured':False}


def rational_sqrt_lower(value,scale=10**6):
    value=F(value)
    if value<0:raise ValueError('positive radicand required')
    return F(isqrt(value.numerator*scale*scale//value.denominator),scale)


@lru_cache(maxsize=1)
def original_inputs():
    member=contract.member(0)
    bare=dressed.regional.charged.bind_member(member)
    odd=dressed.selected_odd.bind_member(member)
    dressed.weak_hole_residual_control();dressed.compensated_control()
    dressed.regional.pairs.local_controls();endpoint_columns_control()
    w=F(bare['result']['original_success_weight_lower'])
    epsilon=F(odd['result']['weak_plus_e_defect_count_upper'])/48
    if w<=F(993,1250) or epsilon>=F(129,400)**2:
        raise ArithmeticError('authenticated original-state margins differ')
    return member,bare,odd,w,epsilon


def separation_row(length,member,w,epsilon):
    geo=geometry(length);control=gauss_chain_control(length);N=geo['path_count']
    v=(80,140,199)[length-1];radius=(20,30,45)[length-1]
    residual=sqrt_upper(F(40*length))+2+F(2*(length+1),7)
    if residual**2>=v:raise ArithmeticError('full-H moment bound failed')
    success=1-length*(1-w);simple=1-length*(1-F(993,1250))
    bare=simple-F(v,radius**2)
    amplitude=rational_sqrt_lower(F(8,9)*bare)
    if amplitude<=F(129,400):raise ArithmeticError('dressed amplitude comparison is not positive')
    dressed_lower=(amplitude-F(129,400))**2
    center=F(6)+F(63*length,4);band=(center-radius,center+radius)
    retention=filtered.triangle_band_retention(F(1,128),band[1])
    compiler=filtered.filter_compiler(member)
    compiler['local_complete_instrument_outcomes']=4*N+1
    compiler['fine_outcomes']=(4*N+1)*compiler['inner_time_cells']
    clocks=spectroscopy.clock_plan(member,'1/64','1/32768','1/32768')
    loss=2*F(clocks['clock_plan']['complete_clock_tail_upper'])+F(clocks['two_complete_clock_error_upper'])
    recorded=dressed_lower*retention-compiler['distinguished_original_event_probability_error_upper']-loss
    margin=(F(4,25),F(2,25),F(7,250))[length-1]
    if recorded<=margin or band[0]<=0:raise ArithmeticError('recorded separated positive window failed')
    J=int(clocks['clock_plan']['clock_size'])
    return {'geometry':geo,'exact_chain_control':control,
        'original_bare_success_weight_lower':success,'original_bare_simple_strict_lower':simple,
        'path_hopping_norm_upper':sqrt_upper(F(40*length)),
        'off_path_endpoint_norm_upper':2,'magnetic_norm_upper':F(2*(length+1),7),
        'full_residual_norm_upper':residual,'centered_second_moment_upper':v,
        'onsite_electric_commutator_center':center,'center_is_static_potential_or_mass':False,
        'ideal_transfer_band':band,'ideal_bare_weight_strict_lower':bare,
        'original_dressing_residual_squared_upper':epsilon,
        'scaled_bare_amplitude_lower':amplitude,'dressed_ideal_weight_strict_lower':dressed_lower,
        'triangle_band_retention_lower':retention,'filter_compiler':compiler,
        'complete_source_clock_plan':clocks,'complete_two_clock_event_loss_upper':loss,
        'reported_transfer_interval':[band[0]-F(1,32),band[1]+F(1,32)],
        'original_filtered_recorded_weight_strict_lower':recorded,
        'original_filtered_recorded_simple_strict_lower':margin,
        'uniform_three_setting_joint_simple_strict_lower':margin/3,
        'complete_probe_outcomes':4*N+1,
        'complete_fine_records_with_two_clocks':compiler['fine_outcomes']*J*J,
        'complete_fine_records_with_two_clocks_and_charge_meters':compiler['fine_outcomes']*J*J*5**8,
        'failure_and_reference_records_retained':True,
        'increment_meter_inside_controlled_instantaneous_row':True,
        'meter_records_final_absolute_charge':False,
        'weight_is_conditioned_on_success':False,
        'finite_pair_interval_is_single_endpoint_mass':False}


def bind_member(member):
    contract.verify_member(member)
    actual,bare,odd,w,epsilon=original_inputs()
    if member!=actual:raise ValueError('separated positive weights are proved only on primary member zero')
    rows=[separation_row(length,member,w,epsilon) for length in (1,2,3)]
    value=contract.serial({'member_id':member['member_id'],**member['ids'],'level':0,
        'bare_source_result_sha256':bare['result_sha256'],
        'selected_odd_source_result_sha256':odd['result_sha256'],
        'same_source_state_reference_clock':True,'endpoint_column_control':endpoint_columns_control(),
        'controlled_increment_meter':controlled_increment_meter_control(),
        'separations':rows,'uniform_setting_register_size':3,
        'complete_unfiltered_probe_outcomes_with_setting':sum(row['complete_probe_outcomes'] for row in rows),
        'closed_source_total_charge':0,'regional_factors_odd_complete_probe_even':True,
        'pair_only_mass_ambiguity':'m -> m+d; V(R) -> V(R)-2d preserves E_pair(R)=2m+V(R)',
        'charged_asymptotic_threshold_proved':False,'joint_zero_momentum_pole_proved':False,
        'physical_particle_map':[],'physical_mass_table_GeV':[],
        'larger_capacity_or_large_separation_extrapolation':False})
    return {'schema':'oph.particle_candidate_endpoint_threshold_bound.v1',
        'result':value,'result_sha256':contract.digest(value)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('endpoint response result or primary source changed')


def certificate():
    paths_=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths_)
    result={'schema':'oph.particle_candidate_endpoint_threshold.v1',
        'selected_member_bound':bind_member(contract.member(0)),'source_sha256':pins}
    if provenance.definition_hashes(paths_)!=pins:raise ValueError('endpoint response source changed during replay')
    return result
