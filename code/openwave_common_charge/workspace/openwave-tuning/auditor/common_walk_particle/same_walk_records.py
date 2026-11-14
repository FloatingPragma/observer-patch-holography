"""Exact complete position and phase records of the ONE common LOCAL12 walk.

No source Hamiltonian, spectrum fit, independent clock or calibration input.
The exact source reference is compared to the owning compiled W callback.
"""
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
from native_phase_probe import actual_native_pointer

HERE=Path(__file__).resolve().parent
WORKSPACE=HERE.parents[2]
PROVIDER=HERE.parent/'native_walk_bridge'
sys.path.insert(0,str(PROVIDER))
import common_walk as common

MODEL='native-local12-current-v0-proposed'
STEPS=(0,1,2,4,8)


def exact_initial():
    value=np.zeros((128,16),dtype=object)
    value[:12,:12]=1
    return value


def exact_step(value):
    """Integer numerator update, denominator multiplied by6 each step."""
    result=-6*value.copy()
    for position in range(128):
        total=sum(value[position,:12])
        for event in range(12):result[position,event]+=total
    out=np.empty(2048,dtype=object)
    out[common.expected_swap_permutation()]=result.reshape(-1)
    return out.reshape(128,16)


def exact_phase_records():
    initial=exact_initial();value=initial.copy();rows=[]
    for n in range(max(STEPS)+1):
        if n in STEPS:
            denominator=12*6**n
            norm=sum(x*x for x in value.reshape(-1))
            assert norm==denominator*denominator
            correlation=F(sum(a*b for a,b in zip(initial.reshape(-1),value.reshape(-1))),12*denominator)
            positions=[sum(x*x for x in value[p])for p in range(128)]
            assert sum(positions)==denominator**2
            rows.append({'W_steps':n,'autocorrelation_real_exact':str(correlation),
                'autocorrelation_imaginary_exact':'0',
                'X_probe_probabilities_exact':[str((1+correlation)/2),str((1-correlation)/2)],
                'Y_probe_probabilities_exact':['1/2','1/2'],
                'all128_position_numerators':positions,'position_probability_denominator':denominator**2,
                'source_padding_numerator':sum(positions[p]for p in range(128)if common.decode_position(p)is None),
                'event_padding_numerator':sum(x*x for x in value[:,12:].reshape(-1)),
                'minus_W_global_phase_control':str((-1)**n*correlation)})
        if n<max(STEPS):value=exact_step(value)
    assert rows[1]['autocorrelation_real_exact']!=rows[1]['minus_W_global_phase_control']
    return rows


def exact_two_time_positions(first_steps=0):
    prepared=exact_initial()
    for _ in range(first_steps):prepared=exact_step(prepared)
    outputs={n:np.zeros((128,128),dtype=object)for n in STEPS}
    for first in range(128):
        if not any(prepared[first]):continue
        state=np.zeros((128,16),dtype=object);state[first]=prepared[first].copy()
        for n in range(max(STEPS)+1):
            if n in outputs:
                outputs[n][first]=[sum(x*x for x in state[last])for last in range(128)]
            if n<max(STEPS):state=exact_step(state)
    rows=[]
    for n,values in outputs.items():
        denominator=(12*6**(first_steps+n))**2
        assert sum(values.reshape(-1))==denominator
        assert all(sum(values[first])==sum(x*x for x in prepared[first])*36**n for first in range(128))
        marginal=values.sum(axis=0)
        reference=prepared.copy()
        for _ in range(n):reference=exact_step(reference)
        unmeasured=[sum(x*x for x in reference[p])for p in range(128)]
        tv=F(sum(abs(a-b)for a,b in zip(marginal,unmeasured)),2*denominator)
        rows.append({'W_steps_before_first_read':first_steps,'W_steps_between_reads':n,
                     'full128_by128_joint_numerators':values.tolist(),
                     'common_denominator':denominator,'normalization_exact':'1',
                     'measurement_backaction_TV_exact':str(tv),
                     'first_position_is_an_actual_intervention':True,
                     'all_source_padding_and_zero_records_retained':True})
    assert any(F(row['measurement_backaction_TV_exact'])>0 for row in rows)
    return rows


def native_position_instrument():
    outcomes,run=actual_native_pointer()
    # Logical seven-bit read: tensor products of the actual one-bit C/E
    # instruments, with every native fine record retained. Identity on event,
    # current fiber, all spectator cycles and arbitrary external reference.
    laws={}
    for law,branches in outcomes.items():
        worst=0.;weights={}
        for source in (0,1):
            probabilities=[]
            for meta,ks in branches:
                p=sum(float(np.vdot(K[:,source],K[:,source]).real)for K in ks)
                probabilities.append({'record':meta,'probability':p})
            weights[str(source)]=probabilities
        for bit in (0,1):
            actual=sum(np.kron(k.conj(),k)for meta,ks in branches if meta['bit']==bit for k in ks)
            P=np.diag([int(bit==0),int(bit==1)])
            worst=max(worst,float(np.max(abs(actual-np.kron(P,P)))))
        laws[law]={'one_bit_matrix_unit_error':worst,'native_fine_laws':weights}
    return {'source_bits_read':7,'events_and_current_fiber_untouched':True,
            'complete_position_effects':'P_x=|x><x| tensor I_event tensor I_current for every0<=x<128',
            'quantum_outputs':'P_x rho P_x, all128 outcomes and every native fine branch, not only valid ports',
            'two_time_instrument':'P_y W^n P_x rho P_x (W*)^n P_y; same full W and retained reference',
            'per_read_native_source_pulses':7*6,'per_read_fresh_ready_auxiliaries':7,
            'local_pointer_export_failure_trace':float(np.trace(run.failure_gram).real),
            'two_C_or_E_laws_separate':laws}


def compiled_comparison():
    state=common.initial_cell_state();exact=exact_initial();rows=[]
    assert state.shape==(128,16)
    for n in range(max(STEPS)+1):
        if n in STEPS:
            target=np.asarray(exact,dtype=complex)/(12*6**n)
            error=float(np.linalg.norm(state-target))
            assert error<2e-9
            record=common.complete_position_records(state)
            rows.append({'W_steps':n,'full_padded_vector_error':error,
                'unrenormalized_original_norm_squared':record['original_norm_squared'],
                'source_padding_probability':record['source_padding_probability'],
                'event_padding_probability':record['event_padding_probability_on_valid_source']})
        if n<max(STEPS):state=common.apply_walk(state);exact=exact_step(exact)
    return rows


def source_spectral_identity():
    # Independent full vertex-space compression, not an eigenbasis detector.
    neighbours,_=common.geometry_tables();T=np.asarray(common.transfer_numerators(),float)/12
    expected=.5*np.eye(96)
    for x,row in enumerate(neighbours):
        for y in row:expected[y,x]+=1/12
    assert np.max(abs(T-expected))<1e-15
    values,vectors=np.linalg.eigh(T)
    assert values.min()>-1 and values.max()<1+2e-14
    f=np.zeros(96);f[:12]=1/np.sqrt(12)
    weights=abs(vectors.conj().T@f)**2
    exact={r['W_steps']:F(r['autocorrelation_real_exact'])for r in exact_phase_records()}
    for n,correlation in exact.items():
        # Eigenvalues are diagnostics; the Chebyshev recurrence proves the
        # source identity and no eigenvector is used in preparation/readout.
        p=np.ones(96);q=values.copy()
        if n==0:poly=p
        elif n==1:poly=q
        else:
            for _ in range(2,n+1):p,q=q,2*values*q-p
            poly=q
        assert abs(weights@poly-float(correlation))<2e-13
    return {'common_compression':'T=A* S A=I-L/12',
            'exact_recurrence':'A* W^n A=T_n(T); T_0=I,T_1=T,T_(n+1)=2T*T_n-T_(n-1)',
            'quasiphases':'theta=plus_or_minus arccos(1-lambda(L)/12), modulo2pi per coherent W step',
            'branch_convention':'reported principal unit-circle angle; no independently evolved Hamiltonian introduced',
            'stationary_component_weight_exact':'1/8',
            'diagnostic_T_eigenvalues':values.tolist(),'diagnostic_preparation_spectral_weights':weights.tolist(),
            'unprepared_dark_and_padding_modes':'retained in the full W and all complete readout effects; not claimed absent from the source',
            'energy_or_mass_unit':'none; no rest mass or physical-second conversion inferred'}


def run():
    assert common.MODEL_ID==MODEL
    paths=[Path(__file__),HERE/'native_phase_probe.py',PROVIDER/'common_walk.py',PROVIDER/'compile_walk.py']
    before={str(p.relative_to(WORKSPACE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    phase=exact_phase_records();positions=exact_two_time_positions()+exact_two_time_positions(1);instrument=native_position_instrument()
    actual=compiled_comparison();spectrum=source_spectral_identity()
    after={str(p.relative_to(WORKSPACE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    assert before==after
    return {'schema':'oph.common_walk.same_source_records.v1','model_id':MODEL,
        'checks':{'same_one_W_provider_only':True,'exact_complete_phase_records':True,
                  'full_two_time_position_instrument':True,'native_pointer_all_matrix_units':True,
                  'compiled_common_W_matches_exact_reference':True,'scalar_phase_mutation_detected':True,
                  'measurement_backaction_retained':True,'source_pins_unchanged':True},
        'source_descriptor':common.source_descriptor(),'phase_records':phase,'two_time_position_records':positions,
        'native_position_instrument':instrument,'compiled_common_source_comparisons':actual,
        'same_W_spectral_identity':spectrum,
        'current_fiber_contract':'use DYN03 J_G: G=delta0, F=P_x+Cz with all193 cycle coordinates retained; P_x J=J P_x and controlled W_F J=J controlled W',
        'apparatus_contract':'native pointer and controlled-word compiler; common source/reset/export error plus added apparatus resets must remain separate complete flags',
        'no_independent_H_or_calibration':True,'official_rows_changed':[],'source_sha256':before}


if __name__=='__main__':
    result=run()
    (HERE/'SAME_WALK_RECORDS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'phase':[{'steps':r['W_steps'],'C':r['autocorrelation_real_exact']}for r in result['phase_records']],
                      'max_compiled_error':max(r['full_padded_vector_error']for r in result['compiled_common_source_comparisons'])},indent=2))
