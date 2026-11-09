"""Source-bound charge-transfer meters without initial charge dephasing."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

import sympy as s

from model.particle_common_contract import contract, provenance
from model.particle_primary_charged_string import model as charged
from model.particle_primary_channel_labels import channels

SOURCE_INPUTS=(
    'model/particle_primary_charged_readout/__init__.py',
    'model/particle_primary_charged_readout/__main__.py',
    'model/particle_primary_charged_readout/model.py',
    'model/particle_primary_charged_readout/test_model.py',
    'model/particle_primary_charged_readout/PROOF.md',
    'model/particle_primary_charged_readout/README.md',
    'model/particle_common_contract/PROOF.md',
    'model/particle_primary_charged_string/PROOF.md',
    'model/particle_physical_boundary_state/PROOF.md',
    'model/particle_primary_joint_spectroscopy/PROOF.md',
    'model/particle_primary_mass_readout/PROOF.md',
    'Lean/QFT/GroundRecordTail.lean',
)


def clean(matrix):
    return matrix.applyfunc(s.simplify)


def cyclic_shift(modulus, amount):
    if type(modulus) is not int or modulus<1 or type(amount) is not int:
        raise ValueError('positive integer modulus and integer shift required')
    return s.SparseMatrix(modulus,modulus,
        {((j+amount)%modulus,j):1 for j in range(modulus)})


def transfer_record(value):
    if type(value) is not int or not 0<=value<5:
        raise ValueError('a radix-five meter record is required')
    return {0:0,1:6,4:-6}.get(value)


def coherent_meter(charges, operators, shifts, modulus):
    """Execute all finite meter outcomes on a small charge-eigenbasis control."""
    charges=tuple(charges);operators=tuple(s.Matrix(k) for k in operators);shifts=tuple(shifts)
    if not charges or len(charges)>8 or any(type(q) is not int for q in charges):
        raise ValueError('one through eight integer charge basis values required')
    if type(modulus) is not int or not 1<=modulus<=41:
        raise ValueError('a small positive integer control modulus is required')
    if len(operators)!=len(shifts) or not operators or any(type(q) is not int for q in shifts):
        raise ValueError('one integer charge shift per outcome required')
    d=len(charges);Q=s.diag(*charges)
    if any(k.shape!=(d,d) for k in operators):raise ValueError('one common source carrier required')
    if clean(sum((k.H*k for k in operators),s.zeros(d)))!=s.eye(d):
        raise ValueError('the original instrument is not complete')
    if any(clean(Q*k-k*Q-shift*k)!=s.zeros(d) for k,shift in zip(operators,shifts)):
        raise ValueError('a probe violates the actual charge-shift identity')
    W=s.diag(*(cyclic_shift(modulus,q) for q in charges))
    if W.H*W!=s.eye(d*modulus):raise ArithmeticError('the charge coupling is not unitary')
    rows=[]
    for b,(k,shift) in enumerate(zip(operators,shifts)):
        full=clean(W*s.kronecker_product(k,s.eye(modulus))*W.H)
        expected=s.kronecker_product(k,cyclic_shift(modulus,shift))
        if full!=expected:raise ArithmeticError('the coherent charge counter did not factor')
        for r in range(modulus):
            out=s.Matrix(d,d,lambda i,j:full[i*modulus+r,j*modulus])
            rows.append({'probe':b,'meter':r,'K':out})
    if clean(sum((r['K'].H*r['K'] for r in rows),s.zeros(d)))!=s.eye(d):
        raise ArithmeticError('the meter lost complete original probability')
    return {'coupling':W,'rows':rows,'source_dimension':d,'meter_dimension':modulus}


@lru_cache(maxsize=1)
def native_charge_menu():
    data=channels.native_number_control();blocks=data['blocks']
    lower=sum(row['Q6']*row['CAR_modes_per_cell'] for row in blocks if row['Q6']<0)
    upper=sum(row['Q6']*row['CAR_modes_per_cell'] for row in blocks if row['Q6']>0)
    if (lower,upper)!=(-72,72):raise ArithmeticError('the primary local charge range changed')
    charge=charged.endpoint_control()['native_Q6'];alphabet=(-charge,0,charge)
    radix=next(k for k in range(1,2*abs(charge)+2) if len({q%k for q in alphabet})==3)
    decoding={q%radix:q for q in alphabet}
    return {'local_charge_range':[lower,upper],
        'local_charge_difference_range':[lower-upper,upper-lower],
        'complete_absolute_charge_labels':upper-lower+1,
        'difference_meter_radix':radix,
        'admitted_local_transfer_alphabet':list(alphabet),
        'meter_decoding':[decoding.get(k) for k in range(radix)],
        'smallest_ordinary_cyclic_native_number_shift_radix':radix,
        'matter_modes_per_cell':sum(row['CAR_modes_per_cell'] for row in blocks),
        'laboratory_electric_charge_identified':False}


def source_shift_table():
    cells=tuple(product(range(2),repeat=3));where={x:i for i,x in enumerate(cells)};rows=[]
    for x in cells:
        for axis in range(3):
            if x[axis]:continue
            y=list(x);y[axis]=1;y=tuple(y)
            for target,origin in ((x,y),(y,x)):
                for alpha in range(2):
                    delta=[0]*8;delta[where[target]]=6;delta[where[origin]]=-6
                    rows.append({'target':list(target),'origin':list(origin),
                        'matched_beta_component':alpha,'charge_shift':delta})
    if len(rows)!=48 or any(sum(r['charge_shift']) for r in rows):
        raise ArithmeticError('the primary directed-bond charge menu changed')
    rows.append({'failure':True,'charge_shift':[0]*8})
    return rows


def encoder_error(cells, per_cell_error):
    if type(cells) is not int or cells<1:raise ValueError('a positive cell count is required')
    error=F(per_cell_error)
    if error<0:raise ValueError('nonnegative operator error required')
    return min(F(1),2*cells*error)


@lru_cache(maxsize=1)
def quantum_control():
    """A physical neutral three-state control with a noncommuting source clock."""
    H=s.Matrix([[0,2,2],[2,2,0],[2,0,2]])
    Q=s.diag(0,6,-6);P=[s.diag(*(int(i==j) for i in range(3))) for j in range(3)]
    plus=s.zeros(3);plus[1,0]=1/s.sqrt(2)
    minus=s.zeros(3);minus[2,0]=1/s.sqrt(2)
    operators=(plus,minus,s.diag(0,1,1));shifts=(6,-6,0)
    meter=coherent_meter((0,6,-6),operators,shifts,5)
    root=5*s.eye(3)-H;rho=root**2/s.trace(root**2)
    if H*rho!=rho*H or s.trace(rho)!=1 or H*Q==Q*H:
        raise ArithmeticError('the control must be stationary with a noncommuting charge')
    energies=(-2,2,4);energy=[]
    for value in energies:
        e=s.eye(3)
        for other in energies:
            if other!=value:e=e*(H-other*s.eye(3))/(value-other)
        energy.append(e)
    if sum(energy,s.zeros(3))!=s.eye(3):raise ArithmeticError('the energy PVM is incomplete')
    original=[];with_terminal=[];premature=[];wrong_final_charge=0
    for i,ei in enumerate(energy):
        for b,k in enumerate(operators):
            for f,ef in enumerate(energy):
                L=ef*k*ei;prob=s.simplify(s.trace(L*rho*L.H));original.append(prob)
                terminal=[s.simplify(s.trace(p*L*rho*L.H*p)) for p in P]
                with_terminal.append(sum(terminal))
                premature.append(s.simplify(sum(s.trace(ef*k*p*ei*rho*ei*p*k.H*ef) for p in P)))
                # In this control the raising successes end at charge +/-1
                # before the final energy instrument. Its later charge record
                # need not retain that value.
                if b<2:wrong_final_charge+=sum(terminal[j] for j in range(3) if (0,6,-6)[j]!=shifts[b])
    if original!=with_terminal or sum(original)!=1:
        raise ArithmeticError('terminal charge refinement changed prior energy records')
    distance=s.simplify(sum(abs(a-b) for a,b in zip(original,premature))/2)
    if distance<=0 or wrong_final_charge<=0:
        raise ArithmeticError('the noncommuting-history controls lost their distinction')
    return {'H':H,'Q':Q,'rho':rho,'root':root/s.sqrt(s.trace(root**2)),
        'operators':operators,'energy_projectors':tuple(energy),'meter':meter,
        'original_energy_probabilities':tuple(original),
        'premature_charge_readout_record_distance':distance,
        'post_energy_charge_not_equal_probe_endpoint_probability':s.simplify(wrong_final_charge)}


def selected_bound(member, per_cell_encoder_error='0'):
    contract.verify_member(member)
    if member['source']['level']!=0:raise ValueError('the charged response requires member zero')
    parent=charged.bind_member(member);charged.verify_bound(member,parent)
    r=parent['result'];menu=native_charge_menu();rows=source_shift_table();cells=8
    if any(menu['meter_decoding'][q%menu['difference_meter_radix']]!=q
           for row in rows for q in row['charge_shift']):
        raise ArithmeticError('the actual intervention transfer cannot be decoded')
    loss=encoder_error(cells,per_cell_encoder_error)
    weight=max(F(0),F(r['reported_transfer_original_weight_lower'])-loss)
    return contract.serial({'member_id':member['member_id'],**member['ids'],
        'charged_energy_producer_result_sha256':parent['result_sha256'],
        'charge_menu':menu,'source_probe_charge_records':rows,
        'difference_meter_registers':cells,'difference_meter_dimension':menu['difference_meter_radix']**cells,
        'meter_records_distinguish_proved_probe_transfer_alphabet':True,
        'arbitrary_local_charge_differences_resolved':False,
        'one_coupling_controlled_CAR_number_shifts':cells*menu['matter_modes_per_cell'],
        'two_couplings_controlled_CAR_number_shifts':2*cells*menu['matter_modes_per_cell'],
        'meter_read_before_cancellation':False,'initial_absolute_charge_record_retained':False,
        'complete_original_source_reference_history_preserved_exactly_for_exact_gates':True,
        'full_history_preservation_applies_before_optional_terminal_charge':True,
        'charge_record_kind':'local source charge transferred by the probe intervention',
        'terminal_absolute_charge_readout_available':True,
        'terminal_charge_refinement_preserves_prior_energy_record_probabilities':True,
        'terminal_charge_readout_changes_source_reference_state':True,
        'terminal_charge_inferred_equal_to_transfer_label':False,
        'absolute_charge_instrument_outcomes':menu['complete_absolute_charge_labels']**cells,
        'complete_source_energy_probe_record_count':r['all_fine_energy_probe_records'],
        'complete_energy_probe_charge_meter_record_count':r['all_fine_energy_probe_records']*menu['difference_meter_radix']**cells,
        'all_zero_weight_meter_outcomes_retained':True,
        'per_cell_encoder_operator_error':F(per_cell_encoder_error),
        'complete_encoder_half_trace_error_upper':loss,
        'reported_transfer_interval':r['reported_transfer_interval'],
        'reported_joint_opposite_charge_transfer_original_weight_lower':weight,
        'positive_event_requires_probe_success':True,
        'positive_event_includes_failure_energy_records':False,
        'exact_gate_joint_weight_simple_lower':F(1,4),
        'source_clock_energy_zero_changed':False,'source_preparation_changed':False,
        'reference_and_all_failures_retained':True,
        'large_source_or_meter_histogram_executed':False,
        'physical_single_particle_charge_spin_or_rest_mass_identified':False,
        'physical_mass_table_GeV':[]})


def bind_member(member,**kwargs):
    value=selected_bound(member,**kwargs)
    return {'schema':'oph.primary_charged_readout_bound.v1','result':value,'result_sha256':contract.digest(value)}


def verify_bound(member,value,**kwargs):
    if value!=bind_member(member,**kwargs):raise ValueError('charged readout or source binding differs')


def certificate():
    paths=provenance.definition_inputs((*SOURCE_INPUTS,*charged.SOURCE_INPUTS));pins=provenance.definition_hashes(paths)
    q=quantum_control()
    result={'schema':'oph.primary_charged_readout_certificate.v1',
        'selected_member_bound':bind_member(contract.member(0)),
        'quantum_control':{'source_dimension':3,'meter_dimension':5,'energy_probe_histories':27,
            'full_meter_instrument_outcomes':len(q['meter']['rows']),
            'original_energy_probabilities':[str(x) for x in q['original_energy_probabilities']],
            'premature_charge_readout_record_distance':str(q['premature_charge_readout_record_distance']),
            'post_energy_charge_not_equal_probe_endpoint_probability':str(q['post_energy_charge_not_equal_probe_endpoint_probability']),
            'control_is_selected_primary_source':False},
        'source_sha256':pins}
    if provenance.definition_hashes(paths)!=pins:raise ValueError('charged-readout inputs changed during replay')
    return result
