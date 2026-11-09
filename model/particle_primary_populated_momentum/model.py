"""Actual Wilson phase dependence with original, complete energy probabilities."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import numpy as np
import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_pair_channels import model as pairs
from model.particle_primary_staggered_response import model as staggered

branch=staggered.branch
CELLS=tuple(product(range(2),repeat=3))
SOURCE_INPUTS=tuple(dict.fromkeys((
    'model/particle_primary_populated_momentum/__init__.py',
    'model/particle_primary_populated_momentum/__main__.py',
    'model/particle_primary_populated_momentum/model.py',
    'model/particle_primary_populated_momentum/test_model.py',
    'model/particle_primary_populated_momentum/README.md',
    'model/particle_primary_populated_momentum/PROOF.md',
    'model/particle_primary_momentum/PROOF.md',
    'model/particle_primary_momentum_transfer/PROOF.md',
    'model/particle_primary_origin_momentum/PROOF.md',
    'Lean/QFT/PeriodicMomentum.lean',
    *pairs.SOURCE_INPUTS,*staggered.SOURCE_INPUTS,*staggered.odd.SOURCE_INPUTS,
)))


@lru_cache(maxsize=1)
def source_phase_control():
    """Check each actual Clifford axis with two symbolic endpoint phases."""
    beta_pairs=staggered.odd.response.beta_pairs()
    chart=s.Matrix.hstack(beta_pairs[0][0],beta_pairs[1][0],beta_pairs[0][1],beta_pairs[1][1])
    gamma=staggered.odd.response.density.gamma_data();beta=gamma[3]
    wx,wy=s.symbols('w_x w_y');raising=s.zeros(4);raising[:2,2:]=s.eye(2)
    rows=[]
    for axis in range(3):
        h=branch.apparatus.clean(chart.H*((-beta+beta*gamma[axis])/2)*chart)
        residual=(h*wy*raising-wx*raising*h).applyfunc(s.expand)
        if residual[:2,2:]!=-(wx+wy)*s.eye(2)/2 or residual[2:,:2]!=s.zeros(2):
            raise ArithmeticError('phase-dependent creation block differs from the actual Wilson block')
        if residual[:2,:2]!=-wx*h[2:,:2] or residual[2:,2:]!=wy*h[2:,:2]:
            raise ArithmeticError('phase-dependent number block differs')
        if branch.apparatus.clean(h[2:,:2].H*h[2:,:2])!=s.eye(2)/4:
            raise ArithmeticError('off-diagonal Wilson block norm differs')
        rows.append({'axis':axis,'creation_block':'-(w_x+w_y)*I_2/2',
                     'annihilation_block_zero':True,'number_block_norm_per_edge':F(1,2)})
    return {'actual_axis_identities':rows,'oriented_edges_per_axis':8,
        'spin_rows_per_edge':2,'creation_coefficient_HS_squared':'16*sum_i cos(q_i/2)^2',
        'number_remainder_squared_form_upper':'36*D6',
        'coherent_normalization_squared':512,'probe_norm_squared':F(9,64),
        'pair_commutator':'(16-D6)/512','native_integer_charge':6,
        'native_block_dimension':1,'CAR_parity':'even','global_charge_transfer':0}


@lru_cache(maxsize=4)
def analyzer(offset=F(1,16)):
    """Full32-row spatial/internal analyzer; offset is a rational multiple of pi."""
    offset=F(offset)
    base,labels=pairs.coefficient_transform('walsh')
    phases=[s.exp(-s.I*s.pi*s.Rational(offset.numerator,offset.denominator)*sum(x))
            for x in CELLS for _ in range(4)]
    value=base*s.diag(*phases)
    if value*value.H!=s.eye(32):raise ArithmeticError('the complete phase-scanned analyzer is not unitary')
    return value,labels


def centered_moment(defects,creation_cosine_sum,young):
    D=F(defects);C=F(creation_cosine_sum);t=F(young)
    if not 0<=D<=32 or not 0<=C<=3 or t<=0:raise ValueError('invalid source moment bound')
    return ((1+t)*16*C*(1+D/2)+(1+1/t)*36*D)/512


def phase_scan_rows(defects):
    """Certified pi/16 scan; rational pi<22/7 bounds every shifted component."""
    D=F(defects);chi=(16-D)/512
    small_sine_squared=F(121,12544)  # (pi/32)^2 < (11/112)^2.
    plans=((F(6,7),F(21,4)),(F(1),F(5)),(F(3,2),F(9,2)),(F(9),F(3)))
    rows=[]
    for k in CELLS:
        r=sum(k);t,radius=plans[r];C=3-r+r*small_sine_squared
        moment=centered_moment(D,C,t);weight=chi-moment/radius**2
        if weight<=0:raise ArithmeticError('a complete scanned spatial bin lost positive weight')
        rows.append({'fourier_bin':list(k),'scalar_component':0,
            'wavevector_over_pi':[F(1,16)+j for j in k],
            'principal_wavevector_over_pi':[F(1,16)-j for j in k],
            'shifted_axis_count':r,'cosine_squared_sum_upper':C,'young_parameter':t,
            'centered_second_moment_upper':moment,'original_success_weight_lower':chi,
            'original_transfer_interval':[6-radius,6+radius],
            'original_transfer_weight_lower':weight})
    return rows


def analyzer_error(phase_diagonal_error):
    value=F(phase_diagonal_error)
    if value<0:raise ValueError('phase implementation error must be nonnegative')
    return min(F(1),value)


@lru_cache(maxsize=1)
def physical_bond_control():
    """A distinct physical six-state CAR/link control, not the primary density."""
    n=[pairs.car_bilinear(4,i,i) for i in range(4)]
    I=s.eye(16);flux=s.diag(-1,0,1);shift=s.zeros(3);shift[1,0]=shift[2,1]=1
    one_body=s.Matrix([[-1,s.I],[s.I,1]])/2
    hop=sum((one_body[a,b]*s.kronecker_product(shift,pairs.car_bilinear(4,a,2+b))
             for a in range(2) for b in range(2)),s.zeros(48))
    onsite=s.kronecker_product(s.eye(3),3*(n[0]-n[1]+n[2]-n[3]))
    full=onsite+s.Rational(63,4)*s.kronecker_product(flux**2,I)+hop+hop.H
    qx=s.kronecker_product(s.eye(3),n[0]+n[1]-I)
    qy=s.kronecker_product(s.eye(3),n[2]+n[3]-I)
    electric=s.kronecker_product(flux,I)
    indices=[j for j in range(48) if qx[j,j]==electric[j,j] and qy[j,j]==-electric[j,j]]
    V=s.eye(48)[:,indices]
    if len(indices)!=6 or (s.eye(48)-V*V.H)*full*V!=s.zeros(48,6):
        raise ArithmeticError('the physical bond sector is not invariant and complete')
    H=V.H*full*V
    B=tuple(V.H*s.kronecker_product(s.eye(3),pairs.car_bilinear(4,2*j,2*j+1))*V for j in range(2))
    complements=tuple(s.eye(6)-b.H*b for b in B)
    if any(f*f!=f for f in complements):raise ArithmeticError('bond complementary effects are not projections')
    root=40*s.eye(6)-H;rho=root*root/s.trace(root*root)
    # The finite inversion combines endpoint exchange, beta signs and link reversal.
    invert=s.zeros(48);perm=(2,3,0,1)
    for e in range(3):
        for mask in range(16):
            occupied=[j for j in range(4) if mask&(1<<j)]
            mapped=[perm[j] for j in occupied]
            inversions=sum(mapped[i]>mapped[j] for i in range(len(mapped)) for j in range(i+1,len(mapped)))
            sign=(-1)**(inversions+sum(j%2 for j in occupied))
            target=sum(1<<j for j in mapped);invert[(2-e)*16+target,e*16+mask]=sign
    parity=V.H*invert*V
    if parity.H*H*parity!=H or parity.H*B[0]*parity!=-B[1]:
        raise ArithmeticError('physical bond inversion covariance failed')
    return {'H':H,'rho':rho,'root':root/s.sqrt(s.trace(root*root)),
        'probes':B,'complements':complements,'parity':parity,'physical_embedding':V,
        'full_H':full,'Gauss_x':qx-electric,'Gauss_y':qy+electric}


def bond_instrument(q):
    c=physical_bond_control();phase=np.exp(-1j*float(q))
    operators=[]
    for family in (c['probes'],c['complements']):
        A,B=(np.array(v,dtype=complex) for v in family)
        operators.extend(((A+phase*B)/2,(A-phase*B)/2))
    return tuple(operators)


def control_receipt():
    c=physical_bond_control();H=np.array(c['H'],dtype=complex);rho=np.array(c['rho'],dtype=complex)
    q=np.pi/16;K=bond_instrument(q);reverse=bond_instrument(-q)
    eig,U=np.linalg.eigh(H);p=np.diag(U.conj().T@rho@U).real
    law=lambda rows:[p[None,:]*abs(U.conj().T@b@U)**2 for b in rows]
    full=law(K);opposite=law(reverse)
    if abs(sum(np.sum(b) for b in full)-1)>1e-12:raise ArithmeticError('bond complete energy law failed')
    covariance=max(np.max(abs(a-b)) for a,b in zip(full,opposite))
    effect=K[0].conj().T@K[0]
    nonconserved=float(np.linalg.norm(H@effect-effect@H))
    nonprojective=float(np.linalg.norm(effect@effect-effect))
    mixture=sum(np.trace(rho@np.array(b,dtype=complex).conj().T@np.array(b,dtype=complex)).real
                for b in c['probes'])/4
    coherent=float(np.trace(rho@effect).real)
    if covariance>1e-12 or nonconserved<=.8 or nonprojective<=.3 or abs(coherent-mixture)<=3e-5:
        raise ArithmeticError('the independent physical bond control failed')
    return {'domain':'separate exact physical bond CAR/link control; not primary source or state',
        'physical_dimension':6,'complete_instrument_outcomes':4,'complete_energy_histories':144,
        'numerical_probability_tolerance':'1/1000000000000',
        'probability_normalization_verified':True,'inversion_energy_covariance_verified':True,
        'display_decimal_places':12,
        'success_effect_commutator_H_norm':round(nonconserved,12),
        'success_effect_idempotence_defect':round(nonprojective,12),
        'coherent_first_bin_probability':round(coherent,12),
        'incoherent_first_bin_probability':round(float(mixture),12),
        'coherent_incoherent_probability_difference':round(float(abs(coherent-mixture)),12)}


def selected_bound(member,phase_diagonal_error='0'):
    branch._member_zero(member);source=source_phase_control();analyzer()
    upstream=branch.bind_member(member);branch.verify_bound(member,upstream)
    D=F(upstream['result']['q6_selected_defects']['q6_defect_count_upper'])
    if D>=F(46,25):raise ArithmeticError('the selected-state defect bound changed')
    plan=branch.spectroscopy.clock_plan(member,'1/64','1/32768','1/32768')
    clockloss=2*F(plan['clock_plan']['complete_clock_tail_upper'])+F(plan['two_complete_clock_error_upper'])
    if clockloss>F(1,8192):raise ArithmeticError('complete clock loss exceeded allowance')
    phase_error=analyzer_error(phase_diagonal_error);rows=phase_scan_rows(D)
    for row in rows:
        lo,hi=row['original_transfer_interval']
        row['reported_transfer_interval']=[lo-F(1,32),hi+F(1,32)]
        row['reported_joint_original_weight_lower']=max(F(0),row['original_transfer_weight_lower']-clockloss-phase_error)
        row['final_excitation']=branch.readout.excitation_window((lo,hi),row['original_transfer_weight_lower'],
            F(member['state']['target_ground_excess']),2,F(9,64))
    # The successful scalar rows are disjoint events in the same64-outcome instrument.
    aggregate=sum(row['original_transfer_weight_lower'] for row in rows)
    union_final_tail=F(member['state']['target_ground_excess'])/8
    if aggregate-union_final_tail<=F(15,256) or aggregate-clockloss<=F(1,16):
        raise ArithmeticError('the unconditioned joint union lost its source-derived floor')
    uniform=F(177,6400)-F(13767,22400)/F(21,4)**2
    neighborhood=F(177,6400)-F(477,2800)/9
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'source_control':source,'occupation_producer_sha256':contract.digest(upstream),
        'q6_selected_defect_upper':D,'configured_phase_offset_over_pi':[F(1,16)]*3,
        'wavevector_alias_period_over_pi':2,'wavevector_units':'radians per source-cell spacing; a=1',
        'phase_gate_order':32,'complete_probe_outcomes':64,'scalar_spatial_rows':rows,
        'all_wavevector_uniform_band':[F(3,4),F(45,4)],
        'all_wavevector_uniform_original_weight_lower':uniform,
        'checkerboard_neighborhood_coordinate_radius':F(1,4),
        'checkerboard_neighborhood_band':[3,9],
        'checkerboard_neighborhood_original_weight_lower':neighborhood,
        'selected_scan_union_original_weight_lower':aggregate,
        'selected_scan_union_reported_original_weight_lower':max(F(0),aggregate-clockloss-phase_error),
        'scalar_subinstrument_effect_norm_upper':F(1,4),
        'union_final_excitation_initial_cut':2,
        'union_final_excitation_tail_upper':union_final_tail,
        'union_final_excitation_original_weight_lower':aggregate-union_final_tail,
        'union_final_excitation_simple_lower':F(15,256),
        'complete_two_clock_plan':plan,'complete_clock_event_loss_upper':clockloss,
        'complete_phase_implementation_half_trace_error_upper':phase_error,
        'complete_fine_energy_spatial_outcomes':64*plan['clock_plan']['clock_size']**2,
        'selected_preparation_changed':False,'source_clock_or_energy_zero_changed':False,
        'all_failures_and_quantum_reference_retained':True,
        'finite_cube_energy_covariance':'mu(q)=mu(R q)=mu(-q), with the corresponding transformed apparatus',
        'phase_implementation_is_exact_declared_algebraic_service':phase_error==0,
        'spatial_analyzer_is_conserved_momentum_PVM':False,
        'populated_relative_GNS_momentum_component_proved':False,
        'wavevector_zero_establishes_a_rest_frame':False,'refinement_survival_proved':False,
        'large_primary_matrix_or_histogram_executed':False,'physical_mass_table_GeV':[]})


def bind_member(member,phase_diagonal_error='0'):
    result=selected_bound(member,phase_diagonal_error)
    return {'schema':'oph.primary_populated_momentum_bound.v1','result':result,'result_sha256':contract.digest(result)}


def verify_bound(member,value,phase_diagonal_error='0'):
    if value!=bind_member(member,phase_diagonal_error):raise ValueError('populated spatial response or source binding changed')


def certificate():
    paths=provenance.definition_inputs(SOURCE_INPUTS);pins=provenance.definition_hashes(paths)
    result={'schema':'oph.primary_populated_momentum_certificate.v1',
        'selected_member_bound':bind_member(contract.member(0)),
        'independent_physical_control':control_receipt(),'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('populated spatial inputs changed during replay')
    return result
