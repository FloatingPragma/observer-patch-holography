"""Gauss-complete charged spinor pairs on the original primary source."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_charged_string import model as charged
from model.particle_primary_charged_readout import model as meter
from model.particle_primary_pair_channels import model as pairs
from model.particle_primary_charge_conjugation import model as conjugation
from model.particle_primary_geometric_cut import model as cut
from model.particle_primary_joint_spectroscopy import spectroscopy
from model.spin_packet_readout.bounds import sqrt_upper

SOURCE_INPUTS=(
    'model/particle_primary_regional_fermion/__init__.py',
    'model/particle_primary_regional_fermion/__main__.py',
    'model/particle_primary_regional_fermion/model.py',
    'model/particle_primary_regional_fermion/test_model.py',
    'model/particle_primary_regional_fermion/PROOF.md',
    'model/particle_primary_regional_fermion/README.md',
    'model/particle_physical_boundary_state/PROOF.md',
    'Lean/QFT/CPRestrictionNet.lean',
    'Lean/QFT/GroundRecordTail.lean',
    'Lean/Screen/DiscreteGauss.lean',
    *charged.SOURCE_INPUTS,*meter.SOURCE_INPUTS,*pairs.SOURCE_INPUTS,
    *conjugation.SOURCE_INPUTS,*cut.SOURCE_INPUTS,
)


def clean(A):
    return A.applyfunc(s.simplify)


def creator(mode,count):
    if type(count) is not int or not 1<=count<=8 or not 0<=mode<count:
        raise ValueError('an admitted small CAR block is required')
    return s.SparseMatrix(1<<count,1<<count,{
        (mask|(1<<mode),mask):(-1)**((mask&((1<<mode)-1)).bit_count())
        for mask in range(1<<count) if not mask>>mode&1})


@lru_cache(maxsize=1)
def local_control():
    """All endpoint occupations and an exact character-shift link control."""
    fields,defect=pairs.local_controls();I=s.eye(16)
    C=[creator(j,4) for j in range(4)];N=[A*A.H for A in C]
    P0=s.diag(0,1,0);L=s.zeros(3);L[2,1]=1;E=s.diag(-1,0,1)
    Nx=N[0]+N[1];Ny=N[2]+N[3]
    Gx=6*(s.kronecker_product(s.eye(3),Nx)-s.kronecker_product(E,I))
    Gy=6*(s.kronecker_product(s.eye(3),Ny-2*I)+s.kronecker_product(E,I))
    H0=3*s.kronecker_product(s.eye(3),Nx-Ny)+s.Rational(63,4)*s.kronecker_product(E**2,I)
    gauge=s.diag(*(int(Gx[j,j]==0 and Gy[j,j]==0) for j in range(48)))
    anchor=s.eye(48)[:,16+12]
    success=[]
    for A in fields:
        B=s.kronecker_product(L,A);effect=clean(B.H*B)
        if clean(Gx*B-B*Gx)!=s.zeros(48) or clean(Gy*B-B*Gy)!=s.zeros(48):
            raise ArithmeticError('endpoint matter and Wilson flux fail Gauss matching')
        if clean(H0*B-B*H0)!=s.Rational(87,4)*B:
            raise ArithmeticError('the literal onsite/electric center differs')
        if effect!=s.kronecker_product(P0,clean(A.H*A)):
            raise ArithmeticError('charged pair success has an omitted link factor')
        v=B*anchor
        if clean(v.H*v)[0]!=1 or gauge*v!=v:
            raise ArithmeticError('a charged spin-pair anchor is absent from physical support')
        reduced=s.Matrix(4,4,lambda a,b:sum(v[32+a+4*c]*s.conjugate(v[32+b+4*c]) for c in range(4)))
        if reduced!=s.diag(0,s.Rational(1,2),s.Rational(1,2),0):
            raise ArithmeticError('a normalized spin coefficient is not a two-spinor pair')
        success.append(B)
    raw=[C[a]*C[2+b].H for a,b in product(range(2),repeat=2)]
    if clean(sum((A.H*A for A in raw),s.zeros(16)))!=clean(sum((A.H*A for A in fields),s.zeros(16))):
        raise ArithmeticError('the full spin coefficient transform is not complete')
    effect=clean(sum((B.H*B/4 for B in success),s.zeros(48)))
    values=(s.eye(48)-effect).eigenvals()
    if any(v<0 for v in values):raise ArithmeticError('the complete charged spin menu is not positive')
    bare=s.kronecker_product(s.eye(3),C[0])
    if gauge*bare*gauge!=s.zeros(48) or bare*anchor==s.zeros(48,1):
        raise ArithmeticError('a bare charged creator was admitted without compensation')
    return {'fields':fields,'defect':defect,'success':tuple(success),
        'gauge_projector':gauge,'anchor':anchor,'Gx':Gx,'Gy':Gy,'H0':H0,
        'bare_inside_creator':bare,'physical_control_dimension':int(s.trace(gauge)),
        'endpoint_spin_reduced_anchor':s.diag(0,s.Rational(1,2),s.Rational(1,2),0)}


@lru_cache(maxsize=1)
def spin_control():
    beta_pairs=pairs.response.beta_pairs()
    positive=s.Matrix.hstack(*(p for p,m in beta_pairs))
    negative=s.Matrix.hstack(*(m for p,m in beta_pairs))
    sigma,_=pairs.pair_basis();rows=[]
    for q in pairs.channels.group.elements():
        S=pairs.channels.spin_lift(q);g=clean(positive.H*S*positive);h=clean(negative.H*S*negative)
        if clean(g.H*g)!=s.eye(2) or g!=h:
            raise ArithmeticError('positive and negative native endpoint spin charts differ')
        R=s.Matrix(pairs.channels.group.rotation(q));rep=s.diag(1,R)
        for a in range(4):
            if clean(g*sigma[a]*h.H)!=clean(sum((rep[b,a]*sigma[b] for b in range(4)),s.zeros(2))):
                raise ArithmeticError('endpoint spinor tensor is not scalar plus vector')
        character=pairs.channels.group.characters(q)['G1']
        if s.simplify(s.trace(g)-pairs.channels._sym(character))!=0:
            raise ArithmeticError('the native endpoint character is not G1')
        rows.append({'endpoint':g,'pair':rep})
    return {'rows':rows,'binary_group_order':48,'endpoint_dimension':2,
        'endpoint_irrep':'G1','compensator_irrep':'conjugate G1',
        'pair_internal_irreps':{'A1':1,'T1':3},
        'endpoint_2pi_sign':-1,'complete_pair_2pi_sign':1,
        'endpoint_native_Q6':6,'compensating_native_Q6':-6,
        'endpoint_weak_and_color_singlet':True,
        'Lorentz_spin_or_physical_hypercharge_dictionary_proved':False}


@lru_cache(maxsize=1)
def hopping_control():
    """Literal source Wilson columns and exact all-occupation raising bound."""
    beta_pairs=pairs.response.beta_pairs()
    P=s.Matrix.hstack(*(p for p,m in beta_pairs));M=s.Matrix.hstack(*(m for p,m in beta_pairs))
    gamma=pairs.response.density.gamma_data();beta=gamma[3]
    hops=[(-beta+beta*g)/2 for g in gamma[:3]]
    cases=0
    for selected in range(3):
        other=[i for i in range(3) if i!=selected]
        for signs in product((1,-1),repeat=2):
            for modes in (P,M):
                images=[(hops[i] if sign>0 else hops[i].H)*modes for i,sign in zip(other,signs)]
                if clean(sum((A.H*A for A in images),s.zeros(2)))!=s.eye(2):
                    raise ArithmeticError('actual remaining Wilson image modes are not orthonormal')
                cases+=1
    maxima=[]
    for active in range(5):
        if active==0:maxima.append(F());continue
        J=s.zeros(1<<active)
        for mask in range(1<<active):
            for j in range(active):
                if not mask>>j&1:J[mask|(1<<j),mask]=1
        eigen=(J.H*J/2).eigenvals();maximum=max(eigen)
        if maximum>3:raise ArithmeticError('the full occupation raising block exceeds three')
        maxima.append(F(maximum))
    return {'literal_Wilson_column_cases':cases,
        'remaining_mode_Gram':'I_2 at each endpoint; disjoint endpoint image supports',
        'active_raising_pair_counts':list(range(5)),
        'full_Fock_squared_norm_by_active_count':maxima,
        'other_bond_squared_norm_upper':3,'selected_bond_squared_norm_upper':40,
        'magnetic_norm_upper':F(4,7)}


@lru_cache(maxsize=1)
def spin_increment_control():
    """Actual regional doubled-spin increments and coherent pointer alphabets."""
    C=[creator(j,4) for j in range(4)];N=[A*A.H for A in C]
    Zp=N[0]-N[1];Zm=N[2]-N[3]
    raw=[C[a]*C[2+b].H for a,b in product(range(2),repeat=2)]
    shifts=[]
    for a,b in product(range(2),repeat=2):
        A=raw[2*a+b];sp=(-1)**a;sm=-(-1)**b
        if Zp*A-A*Zp!=sp*A or Zm*A-A*Zm!=sm*A:
            raise ArithmeticError('the charged raw entry lacks its actual spin increment')
        shifts.append([sp,sm])
    effect=sum((A.H*A for A in raw),s.zeros(16))
    if effect!=(2*s.eye(16)-N[0]-N[1])*(N[2]+N[3]):
        raise ArithmeticError('the complete raw-spin effect is not independently spin invariant')
    if effect*Zp!=Zp*effect or effect*Zm!=Zm*effect:
        raise ArithmeticError('complete failure does not admit zero spin increment')
    Z=s.diag(1,-1);twirl=sum((r['endpoint']*Z*r['endpoint'].H for r in spin_control()['rows']),s.zeros(2))
    if clean(twirl)!=s.zeros(2):raise ArithmeticError('actual cube invariance does not cancel summed source spin')
    return {'raw_fields':tuple(raw),'positive_doubled_spin':Zp,'negative_doubled_spin':Zm,
        'raw_doubled_spin_increments':shifts,'failure_doubled_spin_increment':[0,0],
        'spin_meter_radix':3,'spin_meter_registers':16,
        'spin_meters_commute_with_native_charge_meters':True,
        'unknown_initial_spin_cancelled_coherently':True,
        'spin_increment_is_final_spin_after_energy_clock':False,
        'cube_average_source_summed_spin_vanishes':True}


def directed_menu():
    cells=tuple(product(range(2),repeat=3));where={x:i for i,x in enumerate(cells)};rows=[]
    for x in cells:
        for axis in range(3):
            if x[axis]:continue
            y=tuple(v+int(i==axis) for i,v in enumerate(x))
            for target,origin in ((x,y),(y,x)):
                for a in range(4):
                    delta=[0]*8;delta[where[target]]=6;delta[where[origin]]=-6
                    rows.append({'target':list(target),'compensator':list(origin),
                        'coefficient':a,'internal_probe_irrep':'A1' if a==0 else 'T1',
                        'charge_increment':delta})
    if len(rows)!=96:raise ArithmeticError('the complete regional endpoint menu changed')
    return rows


def bind_member(member):
    contract.verify_member(member)
    if member['source']['level']!=0:raise ValueError('the spin-resolved charged band requires member zero')
    parent=charged.bind_member(member);charged.verify_bound(member,parent);r=parent['result']
    local_control();spin=spin_control();hopping=hopping_control();rows=directed_menu();increments=spin_increment_control()
    w=F(r['original_success_weight_lower'])
    residual=sqrt_upper(F(40))+sqrt_upper(F(3))+F(4,7)
    if residual**2>=75:raise ArithmeticError('the source charged spin residual exceeds its rational cap')
    band=(F(35,4),F(139,4));ideal=w-F(75,169)
    clock=spectroscopy.clock_plan(member,'1/16','1/128','1/128')
    loss=2*F(clock['clock_plan']['complete_clock_tail_upper'])+F(clock['two_complete_clock_error_upper'])
    recorded=ideal-loss
    if recorded<=F(5,16):raise ArithmeticError('the original charged spin energy margin failed')
    result=contract.serial({'member_id':member['member_id'],**member['ids'],
        'charged_occupation_result_sha256':parent['result_sha256'],
        'regional_cut':'one target cell; crossing link and compensating endpoint retained outside',
        'source_carrier_state_reference_clock_changed':False,
        'global_probe_CAR_parity':0,'inside_and_outside_matter_factor_CAR_parity':1,
        'inside_factor_preserves_full_Gauss_carrier_alone':False,
        'regional_endpoint_irrep':spin['endpoint_irrep'],
        'regional_endpoint_dimension':2,'endpoint_charge_increment':6,
        'compensating_charge_increment':-6,'total_charge_increment':0,
        'internal_pair_probe_irreps':spin['pair_internal_irreps'],
        'full_source_total_spin_state_irrep_inferred':False,
        'per_chosen_coefficient_original_success_lower':w,
        'per_chosen_coefficient_centered_second_moment_upper':75,
        'reference_transfer_center':F(87,4),
        'hopping_control':hopping,'full_residual_norm_upper':residual,
        'ideal_transfer_interval':band,'per_chosen_coefficient_ideal_weight_lower':ideal,
        'complete_clock_plan':clock,'complete_two_clock_event_loss_upper':loss,
        'reported_transfer_interval':[band[0]-F(1,8),band[1]+F(1,8)],
        'per_chosen_coefficient_recorded_weight_lower':recorded,
        'per_chosen_coefficient_recorded_simple_lower':F(5,16),
        'complete_menu':rows,'complete_success_outcomes':96,'complete_probe_outcomes':97,
        'complete_menu_each_coefficient_recorded_weight_lower':recorded/4,
        'complete_menu_scalar_recorded_simple_lower':F(5,64),
        'complete_menu_vector_recorded_simple_lower':F(15,64),
        'complete_joint_probe_energy_outcomes':97*int(clock['clock_plan']['clock_size'])**2,
        'coherent_charge_meter_radix':5,'coherent_charge_meter_registers':8,
        'coherent_charge_meter_dimension':5**8,
        'charge_meter_initial_charge_readout':False,
        'joint_charge_energy_bounds_require_success_record':True,
        'source_C_conjugated_apparatus_has_same_original_complete_record_law':True,
        'source_C_changes_endpoint_charge_sign':True,
        'raw_spin_apparatus':{
            'complete_success_outcomes':96,'complete_probe_outcomes':97,
            'inside_spin_projection_increments':[F(1,2),F(-1,2)],
            'outside_spin_projection_increments':[F(-1,2),F(1,2)],
            'raw_doubled_spin_increments':increments['raw_doubled_spin_increments'],
            'chosen_axis_per_endpoint_spin_meter_radix':3,'spin_meter_registers':16,
            'spin_meter_dimension':3**16,'joint_charge_spin_meter_dimension':5**8*3**16,
            'complete_joint_probe_clock_charge_spin_outcomes':97*int(clock['clock_plan']['clock_size'])**2*5**8*3**16,
            'source_cube_invariance_consumed_for_spin_deficit_bound':True,
            'per_fixed_spin_entry_original_success_lower':w,
            'per_fixed_spin_entry_centered_second_moment_upper':70,
            'reported_transfer_interval':r['reported_transfer_interval'],
            'per_fixed_spin_entry_recorded_weight_lower':F(r['reported_transfer_original_weight_lower']),
            'per_fixed_spin_entry_recorded_simple_lower':F(1,4),
            'each_of_four_joint_spin_rows_recorded_weight_lower':F(r['reported_transfer_original_weight_lower'])/4,
            'each_of_four_joint_spin_rows_recorded_simple_lower':F(1,16),
            'charge_and_spin_counters_preserve_original_complete_reference_history':True,
            'initial_absolute_spin_or_charge_measured':False,
            'simultaneous_noncommuting_spin_axes_claimed':False,
            'spin_increment_is_final_spin_after_energy_clock':False,
            'raw_spin_and_Pauli_records_are_same_joint_apparatus':False},
        'terminal_regional_spin_analyzer_retains_prior_classical_energy_records':True,
        'terminal_regional_spin_analyzer_preserves_source_reference_state':False,
        'terminal_occupation_one_in_the_energy_band_positive_weight_proved':False,
        'endpoint_statistics_are_regional_CAR_factor_statistics':True,
        'spinor_label_is_physical_Lorentz_spin':False,
        'pair_energy_is_single_endpoint_rest_energy':False,
        'pair_energy_divided_by_two_is_justified':False,
        'isolated_charged_particle_or_physical_electron_identified':False,
        'physical_mass_table_GeV':[]})
    return {'schema':'oph.primary_regional_fermion_bound.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value):
    if value!=bind_member(member):raise ValueError('regional fermion source or result differs')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_regional_fermion_certificate.v1',
        'members':[bind_member(contract.member(0))],
        'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('regional fermion inputs changed during replay')
    return result
