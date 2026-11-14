"""Native preparation, complete occupation readout and source binding for N=2."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
import two_charge as tc

ROOT=tc.INPUT
BASE=ROOT/'openwave-tuning/common_dynamics'
sys.path[:0]=[str(BASE/'particles'),str(BASE/'vector_field')]
import massive_coupled as massive
import composition


def words(y=1,dark=True):
    mb=tc.af.mb;compiler=mb.compiler
    address=mb.common.encode_position(y);offset=16*address
    flip=mb.architecture.primitive_routines()['flip']
    # Port0 uniform and a disjoint local orbital are orthogonal. Acting on
    # vacuum in the stated order only changes a global preparation phase.
    out=composition.preparation_word(1)
    out+=compiler.remap(flip,[offset])
    if dark:
        out+=mb.givens_word(offset,offset+2,-math.pi/4)
    else:
        for k in range(1,12):out+=mb.givens_word(offset,offset+k,math.asin(1/math.sqrt(13-k)))
    field=next(iter(tc.initial(y=y,dark=dark)))[2]
    for e,value in enumerate(field):
        for bit in range(4):
            if value>>bit&1:out+=compiler.remap(flip,[2048+4*e+bit])
    return out


def errors(y=1,dark=True,mass=0,steps=3):
    mb=tc.af.mb
    root=json.loads((ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    base=json.loads((BASE/'vector_field/COMPOSITION.json').read_text())
    massdata=json.loads((BASE/'particles/MASSIVE_COUPLED.json').read_text())
    massrow=next(r for r in massdata['experiments'] if r['mass_parameter']==str(mass))
    bound=massive.mass_coin.source_error.exact_template_bound
    flip=F(root['primitive_operator_errors_exact']['flip'])
    uniform=F(base['source_preparation_error_exact'])-flip
    if dark:
        rotation=mb.givens_word(0,2,-math.pi/4)
        template=composition.givens_template(0,2,composition.cert._pi(-1,4))
        second=flip+bound(rotation,template)
    else:second=flip+uniform
    field=next(iter(tc.initial(y=y,dark=dark)))[2]
    field_flips=sum(int(v).bit_count() for v in field)
    preparation=F(base['source_preparation_error_exact'])+second+field_flips*flip
    pointer=massive.particle_records.apparatus.pointer_bound()
    p0=F(base['native_blank_reset_success_lower_exact'])
    pr=F(base['native_ready_reset_success_lower_exact'])
    reset=3203*(1-p0)+3200*(1-pr)
    reset_controller=massive.particle_records.apparatus.configuration
    reset_circuit=massive.particle_records.apparatus.circuit
    blank_budget=reset_controller.budget(4,32,full=True,mode='reset',target=reset_circuit.decode(0))
    ready_budget=reset_controller.budget(4,32,full=True,mode='reset',target=reset_circuit.decode(massive.particle_records.apparatus.pointer.READY_CODE))
    evolution=steps*F(massrow['native_complete_error']['whole_step_operator_error_exact'])
    total=preparation+reset+evolution+3200*pointer
    field_data=json.loads((BASE/'vector_field/NATIVE_FIELD.json').read_text())
    pulses_per_step=(base['source_pulses_per_matter_current_layer']+
                    field_data['oneface_field_step_native_pulses']+
                    massrow['native_complete_error']['additional_native_pulses_per_step'])
    return {'preparation_operator_error_exact':str(preparation),
            'whole_step_operator_error_exact':massrow['native_complete_error']['whole_step_operator_error_exact'],
            'blank_resets':3203,'ready_pointer_resets':3200,
            'reset_timeout_union_upper_exact':str(reset),
            'reset_native_attempts_max_each_blank':blank_budget['native_attempts_max'],
            'reset_native_attempts_max_each_ready':ready_budget['native_attempts_max'],
            'all_reset_native_attempts_max':3203*blank_budget['native_attempts_max']+3200*ready_budget['native_attempts_max'],
            'reset_controller':{'seed_scans':4,'retries':32,'full78_input_domain':True,'native_attempts_are_not_seconds':True},
            'complete_native_TV_upper_exact':str(total),'complete_native_TV_upper':math.nextafter(float(total),math.inf),
            'source_logical_binary_wires':3203,'readout_native_pointers':3200,
            'retained_matter_occupation_bits':2048,'retained_field_bits':1152,
            'source_clean_work_bits':[3200,3201,3202],
            'preparation_native_pulses':len(words(y,dark)),
            'evolution_native_pulses_per_step':pulses_per_step,
            'complete_readout_native_pulses':19200,
            'total_native_pulses_excluding_resets':len(words(y,dark))+steps*pulses_per_step+19200,
            'field_basis_preparation_bit_flips':field_flips,
            'reset_semantics':'Worst-conditional lower success bounds and union bound; no factorization of correlated initial reset probabilities.',
            'domain':'3203 blank source carriers plus3200 fresh READY pointers. Coherent source services, all source/export complements, timeout and fine pointer outcomes retained. No success normalization.',
            'clock':'One application is one source macrostep at formal t=n*h*tau0. This consumer uses no clock/probe apparatus and awards no physical-clock claim.'}


def source_checks():
    """Replay actual finite native pulses, with vacuum and padding retained."""
    mb=tc.af.mb; compiler=mb.compiler
    emitted=json.loads((BASE/'manybody/MANYBODY_PROGRAMS.json').read_text())['mode']['onsite_reflection_word']
    assert emitted==mb.fock_reflection_word()
    kick=json.loads((BASE/'particles/MASSIVE_KICK_PROGRAM.json').read_text())['onsite16_mode_word']
    assert kick==massive.fock_kick_word(F(1,8),F(1))
    pairs=list(combinations(range(16),2));rng=np.random.default_rng(1402)
    amplitude=rng.normal(size=120)+1j*rng.normal(size=120);amplitude/=np.linalg.norm(amplitude)
    antisym=np.zeros((16,16),complex)
    source=np.zeros((1<<16,4),complex);source[0,0]=1;source[1<<15,1]=1;source[-1,3]=1
    for (i,j),value in zip(pairs,amplitude):
        antisym[i,j]=value;antisym[j,i]=-value;source[(1<<i)|(1<<j),2]=value
    actual=compiler.execute(kick+emitted,16,source)
    k=tc.coin(1); transformed=k@antisym@k.T
    target=np.zeros_like(source);target[0,0]=1
    target[1<<15,1]=-np.exp(1j/8)
    target[-1,3]=-np.exp(14j/8)
    for i,j in pairs:target[(1<<i)|(1<<j),2]=transformed[i,j]
    residual=float(np.linalg.norm(actual-target))
    assert residual<2e-10,residual
    # Entire exterior-square matrix checks the rank-one identity.
    u=np.zeros(16);u[:12]=1/np.sqrt(12);p=np.outer(u,u)
    r=2*p-np.eye(16)
    identity_error=0.
    for i,j in pairs:
        a=np.zeros((16,16));a[i,j]=1;a[j,i]=-1
        identity_error=max(identity_error,float(np.max(abs(r@a@r.T-(a-2*(p@a+a@p))))))
    assert identity_error<2e-14
    # Check a dark preparation using actual source pulses on its local modes.
    v=np.zeros(8,complex);v[1]=1
    got=compiler.execute(mb.givens_word(0,2,-math.pi/4),3,v)
    want=np.zeros(8,complex);want[1]=1/np.sqrt(2);want[4]=-1/np.sqrt(2)
    dark_error=float(np.linalg.norm(got-want));assert dark_error<2e-11
    # All2048 modes: compare fast current/CAR routing with the frozen map.
    transport_error=0
    examples=[(i,j) for i,j in ((0,1),(1,17),(2046,2047),(1900,2))]
    examples += [tuple(sorted(map(int,rng.choice(2048,2,replace=False)))) for _ in range(64)]
    for left,right in examples:
        left,right=sorted((left,right));a,i=divmod(left,16);b,j=divmod(right,16)
        x,y,k,l,sign,current=tc.pair_move(a,b,i,j)
        out,f,s,c=tc.af.mode_field.mode_transport((left,right),np.zeros(288,dtype=np.int64),16)
        expected=(16*x+k,16*y+l)
        assert expected==out and sign==s
        for e,delta in current:transport_error=max(transport_error,abs(int(c[e])-delta))
        assert bytes(f.astype(np.uint8))==tc.shift_field(bytes(288),current)
    return {'frozen_emitted_reflection_word_equal':True,'frozen_emitted_mass_word_equal':True,
            'full16_mode_actual_native_pulse_replay_residual':residual,
            'native_replay_columns':['vacuum','padded_single','random_complex_two_particle_all120_pairs','full16_mode_occupation'],
            'all120_exterior_square_rank_one_identity_residual':identity_error,
            'dark_preparation_native_residual':dark_error,'frozen_transport_examples':len(examples),
            'transport_current_residual':transport_error,'native_pulses_replayed':len(kick)+len(emitted)}


def run():
    checks=source_checks(); rows=[];programs=[]
    for y in (1,2,5):
        for mass in (0,1):
            rows.append({'second_vertex':y,'mass':mass,'dark':True,**errors(y,True,mass,3)})
    rows.append({'second_vertex':1,'mass':0,'dark':False,**errors(1,False,0,2)})
    for y,dark in ((1,True),(2,True),(5,True),(1,False)):
        programs.append({'second_vertex':y,'dark':dark,'source_preparation_word':words(y,dark)})
    data={'schema':'oph.two_charge.native_contract.v1','candidate_id':tc.MODEL,
          'checks':checks,'cases':rows,'readout':'Read every2048matter occupancy and1152field bits with the frozen native_bit_pointer_word; decode pair positions, events and last current classically. All sectors, padding and failures retained.',
          'one_particle_parity_decoder_used':False,'native_complete_instrument_bound':True,
          'ideal_numerical_law_scope':'Native TV bounds concern the exact ideal circuit, including arbitrary spectator references. Binary64 decimal tables have residual checks and no interval-enclosed probability certificate.',
          'source_reset_domain':'Bound concerns terminal particle/current/field marginal and original untouched quantum reference. Fine reset tapes retained by instrument but not asserted independent of source/reference.',
          'programs':programs}
    (tc.HERE/'NATIVE_CONTRACT.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'checks':checks,'TV_bounds':[r['complete_native_TV_upper'] for r in rows]}),flush=True)


if __name__=='__main__':run()
