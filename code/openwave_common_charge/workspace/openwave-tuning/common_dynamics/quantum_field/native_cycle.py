"""Native preparation, clock and field readout for the exact cycle sector."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

import quantum_state as qs

HERE=Path(__file__).resolve().parent
ROOT=qs.ROOT
nf=qs.nf
compiler=nf.compiler
cert=nf.cert
iv=nf.iv
sys.path.insert(0,str(HERE.parent/'clock'))
import phase_clock
from model.reference_recovery import recovery
from model.native_clock import configuration
from model.observer.source_discrimination import circuit
from model.native_interface import pointer


def gray_preparation_bound(word, values):
    """Outward proof for the complete native Gray amplitude preparation."""
    size=len(values);bits=(size-1).bit_length()
    if size != 1<<bits or bits>6 or not any(values):
        raise ValueError('complete nonzero integer target through six bits required')
    gray=[i^(i>>1) for i in range(size)]
    offset=0;error=F();rows=[]
    for row in range(1,size):
        a,b=gray[row-1],gray[row];difference=a^b
        target=difference.bit_length()-1
        orientation=(-1)**((a>>target)&1)
        variables=[]
        for mask in range(size):
            if mask&difference:continue
            label=''.join('Y' if j==target else 'Z' if mask>>j&1 else 'I' for j in range(bits))
            offset,fixed,angles=cert.consume(word,offset,cert.pauli(label))
            assert len(angles)==1
            variables.append((mask,angles[0]));error+=fixed
        theta=orientation*(1<<(bits-1))*variables[0][1]
        assert all(v==orientation*theta*(-1)**((a&mask).bit_count())/(1<<(bits-1)) for mask,v in variables)
        left=values[a]**2;right=sum(values[gray[j]]**2 for j in range(row,size))
        assert values[a]>=0
        error+=cert.two_vector_error(theta,left,right,False)
        rows.append({'plane':[a,b],'left_weight':left,'right_weight':right})
    assert offset==len(word)
    return error+cert.controls_bound(word,{}),rows


def state_preparation_word():
    state,values,phases,data=qs.gaussian_cycle(16,F(1,3),8,cap=7)
    amplitude=np.asarray(values,dtype=float)/math.sqrt(sum(v*v for v in values))
    real_word=compiler.compile_real_state(amplitude)
    real_error,planes=gray_preparation_bound(real_word,values)
    phase_word=compiler.compile_diagonal(phases)
    phase_error=nf.phase_bound(phase_word,[iv.Interval.point(F(float(p))) for p in phases])
    return real_word+phase_word,real_error+phase_error,state,values,phases,data,planes


def encoder_blocks():
    """Map a pivot coordinate into every physical edge of the native face."""
    face=qs.af.data()[2][:,0]
    pivot=int(np.flatnonzero(face)[0]);sign=int(face[pivot])
    blocks=[]
    for edge in np.flatnonzero(face):
        if edge!=pivot:
            blocks.append({'source_edge':pivot,'target_edge':int(edge),
                           'add_sign':int(face[edge]*sign)})
    for z in range(16):
        expected=(z*face*sign)%16
        actual=np.zeros(288,dtype=int);actual[pivot]=z
        for block in blocks:actual[block['target_edge']]=(block['add_sign']*z)%16
        assert np.array_equal(actual,expected)
        assert not np.any((qs.af.data()[1]@actual)%16)
    return pivot,blocks


def encoder_word(inverse=False):
    pivot,blocks=encoder_blocks();out=[]
    for block in blocks:
        gates=nf.add_program() if block['add_sign']==1 else list(reversed(nf.add_program()))
        mapping={j:2048+4*pivot+j for j in range(4)}
        mapping.update({j+4:2048+4*block['target_edge']+j for j in range(4)})
        mapping.update({8+j:nf.SCRATCH[j] for j in range(3)})
        out.extend(nf.arithmetic_word(gates,mapping))
    return compiler.inverse(out) if inverse else out


def compressed_word(setting='X'):
    prep,prep_error,state,values,phases,data,planes=state_preparation_word()
    magnetic,mag_error=nf.plaquette_word(data['beta'])
    electric,elec_error=nf.electric_word(data['alpha']*qs.CYCLE_NORM2)
    H=compiler.remap(nf.H(0),[4])
    word=prep+H
    for _ in range(8):
        word+=magnetic+electric+[compiler.pulse('Z',-float(data['delta'])/16,[4])]
    word+=compiler.inverse(nf.qft_word())
    if setting=='Y':word+=[compiler.pulse('Z',math.pi/4,[4])]
    elif setting!='X':raise ValueError('X or Y native clock readout required')
    word+=H
    qft_error=nf.root_certificate.exact_template_bound(nf.qft_word(),nf.qft_template())
    error=prep_error+8*(mag_error+elec_error)+qft_error+phase_clock.exact_word_bound(8,setting,data['delta'])
    return word,error,{'preparation_error':prep_error,'magnetic_error':mag_error,
                       'electric_compressed_error':elec_error,'qft_error':qft_error,
                       'initial_state':state,'data':data,'values':values,'phases':phases,'planes':planes}


def full_word(setting='X'):
    """Pure-gauge specialization of the actual full source, no dense array."""
    prep,*_=state_preparation_word();pivot,blocks=encoder_blocks()
    yield from compiler.remap(prep,[2048+4*pivot+j for j in range(4)])
    yield from encoder_word()
    yield from compiler.remap(nf.H(0),[3203])
    par=qs.parameters(F(1,3),8)
    couplings=[F()]*208;couplings[0]=par['beta']
    for _ in range(8):
        yield from nf.field_word(par['alpha'],couplings)
        yield compiler.pulse('Z',-float(par['delta'])/16,[3203])
    yield from encoder_word(True)
    yield from compiler.remap(compiler.inverse(nf.qft_word()),[2048+4*pivot+j for j in range(4)])
    if setting=='Y':yield compiler.pulse('Z',math.pi/4,[3203])
    elif setting!='X':raise ValueError('declared clock quadrature required')
    yield from compiler.remap(nf.H(0),[3203])


def joint_effect(effect_rows, labels):
    """Complete native fine/coarse record probability uses this tensor effect."""
    effect=np.ones((1,1),complex)
    for label in reversed(labels):effect=np.kron(effect,effect_rows[label])
    return effect


def run():
    paths=[Path(__file__),Path(qs.__file__),Path(nf.__file__),Path(compiler.__file__),
           Path(cert.__file__),Path(configuration.__file__),Path(phase_clock.__file__),
           HERE.parent/'vector_field/RESET_DOMAIN.json',HERE.parent/'vector_field/MANIFEST.json',
           ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    instruments,_=phase_clock.actual_native_pointer()
    rows=[];programs={};largest=0.;resources=[]
    root=json.loads((ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    primitive={k:F(v) for k,v in root['primitive_operator_errors_exact'].items()}
    counts=Counter(k for k,_ in nf.add_program())
    adder=sum((number*primitive[k] for k,number in counts.items()),F())
    pivot,blocks=encoder_blocks()
    p0=min(map(F,configuration.success_probabilities(4,32,full=True,mode='reset',target=circuit.decode(0))))
    pr=min(map(F,configuration.success_probabilities(4,32,full=True,mode='reset',target=circuit.decode(pointer.READY_CODE))))
    reset=3204*(1-p0)+5*(1-pr)
    for setting in ('X','Y'):
        word,error,parts=compressed_word(setting)
        actual=recovery.execute_complete(word,5,np.eye(32,dtype=complex)[:,[0]])
        v=actual.success[:,0]
        target=parts['initial_state'].copy()
        for _ in range(8):target=qs.source_step(target,parts['data']['alpha'],parts['data']['beta'])
        target=np.fft.fft(target,norm='ortho')
        clock=compiler.execute(phase_clock.local_words(8,setting,F(1,8)),1,np.array([1,0],complex))
        expected=np.kron(clock,target)
        discrepancy=float(np.linalg.norm(v-expected));largest=max(largest,discrepancy)
        assert discrepancy<float(error)+2e-12
        laws={}
        for law,branches in instruments.items():
            effects=[np.zeros((2,2),complex),np.zeros((2,2),complex)];fine=[]
            for metadata,kraus in branches:
                effect=sum((K.conj().T@K for K in kraus),np.zeros((2,2),complex))
                effects[metadata['bit']]+=effect
                fine.append({'record':metadata,'effect_real':effect.real.tolist(),'effect_imag':effect.imag.tolist()})
            probabilities=[]
            for outcome in range(32):
                effect=joint_effect(effects,[(outcome>>j)&1 for j in range(5)])
                probabilities.append(float(np.vdot(v,effect@v).real))
            assert abs(sum(probabilities)-float(np.vdot(v,v).real))<2e-11
            laws[law]={'joint_angle_bin_clock_probabilities':probabilities,
                       'fine_one_bit_effects':fine,
                       'complete_five_bit_fine_probability_rule':'v* (tensor in descending bit order of the five retained fine effects) v; first native export/reset failures are additional outcomes',
                       'success_normalization':False}
        electric_full,ee=nf.electric_word(parts['data']['alpha'])
        full_per_step=288*ee+parts['magnetic_error']+2*len(blocks)*adder
        clock_and_qft=parts['qft_error']+phase_clock.exact_word_bound(8,setting,F(1,8))
        # exact_word_bound includes the clock pointer once; the four field
        # coordinate pointers are counted here. Encoder and decoder each
        # have one adder per nonpivot face edge.
        bound=reset+parts['preparation_error']+2*len(blocks)*adder+8*full_per_step+clock_and_qft+4*phase_clock.pointer_bound()
        rows.append({'setting':setting,'native_word_pulses':len(word),'native_state_error':discrepancy,
                     'native_export_failure_probability':float(actual.failure_gram[0,0].real),
                     'compressed_native_error_exact':str(error),'native_pointer_laws':laws,
                     'full_physical_program_record_TV_bound_exact':str(bound),
                     'full_physical_program_record_TV_bound_upper':math.nextafter(float(bound),math.inf)})
        programs[setting]=word
        count=0;max_wire=-1
        for pulse in full_word(setting):
            count+=1;max_wire=max(max_wire,*pulse['targets'])
        resources.append({'setting':setting,'full_word_pulses':count,'max_wire':max_wire})
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE/'NATIVE_CYCLE_PROGRAM.json').write_text(json.dumps({'compressed_five_bit_words':programs,
        'physical_pivot':pivot,'encoder_blocks':blocks,'full_source_generator':'native_cycle.full_word',
        'full_source_wires':3204,'pure_gauge_reference':'all2048fermionmode occupations zero; all field/work bits reset before Gaussian seed',
        'clock_wire':3203,'field_pointer_bits':[2048+4*pivot+j for j in range(4)],
        'blank_resets':3204,'ready_pointer_resets':5,
        'reset_success_lower_exact':{'blank':str(p0),'ready':str(pr)},
        'complete_failure_rule':'Retain the first reset timeout, every export complement and every native pointer fine record. No success conditioning or independence assumption.'},indent=2)+'\n')
    result={'schema':'oph.local12.native_quantum_cycle.v1','checks':{
       'integer_amplitude_preparation_outward_certificate':True,
       'actual_native_cosine_and_electric_words':True,
       'same_clock_on_fresh_wire':True,'both_native_pointer_laws':True,
       'all32_joint_angle_clock_outcomes_retained':True,'exact_full_field_cycle_encoder':True,
       'factorized_full_source_error_and_correlated_resets':True,
       'full_physical_word_streamed':True,'source_bytes_unchanged':True},
       'experiments':rows,'full_resources':resources,'maximum_numeric_state_error':largest,
       'full_initialization_reset_union_bound_exact':str(reset),
       'input_domain':'arbitrary coherent/correlated total-two native carriers before reset; exact conclusive rank-one Kraus argument in vector_field/RESET_DOMAIN; arbitrary reference retained through initialized word',
       'scope':'Exact pure-gauge one-face reducing code, with complete physical field encoder and native finite instrument. The numerical full-source execution is factorized through its proven16-dimensional cycle code, not a2^3204 matrix. No photon, matter vacuum, or Maxwell continuum claim.',
       'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'NATIVE_CYCLE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'largest_native_residual':largest,'full_resources':resources,
         'full_record_bounds':[r['full_physical_program_record_TV_bound_upper'] for r in rows]}),flush=True)


if __name__=='__main__':run()
