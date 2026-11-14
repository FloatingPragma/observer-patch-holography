"""Native32-mode sea and193-cycle Gauss preparation, complete failures.

This program initializes a state. It does not import the one-particle
Dirac compiler as a full Fock evolution or interacting-vacuum certificate.
"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

import reference as ref

HERE=Path(__file__).resolve().parent
ROOT=ref.ROOT
sys.path.insert(0,str(HERE.parent/'quantum_field'))
import native_cycle
from model.native_clock import configuration
from model.observer.source_discrimination import circuit
from model.native_interface import pointer

nf=ref.nf;compiler=nf.compiler;cert=nf.cert
compiled=ref.dp.spin.compiled
FIELD_BASE=32
SCRATCH=(1184,1185,1186)
TOTAL_WIRES=1187


def validated_remap(word,targets):
    """Exact address lift after the source word and distinct map are checked."""
    return [{'kind':p['kind'],'angle':p['angle'],
             'targets':[int(targets[q]) for q in p['targets']]} for p in word]


def increment_program():
    """Fieldbits1..4 += occupationbit0; scratch5..7 is restored."""
    gates=[]
    for target in reversed(range(4)):
        controls={0:1};controls.update({1+j:1 for j in range(target)})
        if len(controls)==1:gates.append(('cx',(0,1+target)))
        else:gates+=compiled.multi_control(controls,1+target,(5,6,7))
    return gates


def prepare_blocks(kind='vacuum',particle=0,hole=1,base=16):
    G=ref.reference_charge(kind);geometry=ref.dp.field.current_walk.geometry()
    pattern=ref.sector_pattern(kind,particle,hole)
    seed=[]
    for occupancy,(a,b,sign) in zip(pattern,ref.pairs()):
        if occupancy:seed.append(a)
        if occupancy==2:seed.append(b)
    tree_offsets=(-geometry['P']@G)%16
    offsets=[FIELD_BASE+4*e+bit for e,value in enumerate(tree_offsets)
             for bit in range(4) if (int(value)>>bit)&1]
    shear=[]
    for e in geometry['tree']:
        for j in np.flatnonzero(geometry['C'][e]):
            shear.append({'source_edge':int(geometry['chords'][int(j)]),'target_edge':int(e),
                          'sign':int(geometry['C'][e,j])})
    controlled=[]
    for mode in range(32):
        for e in np.flatnonzero(geometry['P'][:,12*(mode//4)]):
            controlled.append({'occupation':mode,'target_edge':int(e),
                               'sign':int(geometry['P'][e,12*(mode//4)])})
    return {'kind':kind,'particle':particle,'hole':hole,'rational_gaussian_base':base,
            'matter_seed_occupied_modes':seed,'negative_orbital_pairs':ref.pairs(),
            'chord_registers':geometry['chords'],'tree_offset_flip_wires':offsets,
            'cycle_shears':shear,'occupation_field_additions':controlled,
            'data_wires':1184,'clean_scratch':list(SCRATCH),'blank_resets':1187,
            'terminal_data_pointer_reads':1184,'pointer_auxiliary_qubits':[1187,1188],
            'clock_or_spectral_probe_included':False,'reference_G':G.tolist()}


def word(program,gaussian_word):
    primitive=compiled.controls()['primitive_words']
    mb=ref.dp.field.mb
    for mode in program['matter_seed_occupied_modes']:
        yield from validated_remap(primitive['flip'],[mode])
    for a,b,sign in ref.pairs():
        yield from mb.givens_word(a,b,sign*math.pi/4)
    for edge in program['chord_registers']:
        yield from validated_remap(gaussian_word,[FIELD_BASE+4*edge+j for j in range(4)])
    for wire in program['tree_offset_flip_wires']:
        yield from validated_remap(primitive['flip'],[wire])
    for block in program['cycle_shears']:
        gates=nf.add_program() if block['sign']==1 else list(reversed(nf.add_program()))
        mapping={j:FIELD_BASE+4*block['source_edge']+j for j in range(4)}
        mapping.update({4+j:FIELD_BASE+4*block['target_edge']+j for j in range(4)})
        mapping.update({8+j:SCRATCH[j] for j in range(3)})
        for kind,targets in gates:yield from validated_remap(primitive[kind],[mapping[q] for q in targets])
    for block in program['occupation_field_additions']:
        gates=increment_program() if block['sign']==1 else list(reversed(increment_program()))
        mapping={0:block['occupation']}
        mapping.update({1+j:FIELD_BASE+4*block['target_edge']+j for j in range(4)})
        mapping.update({5+j:SCRATCH[j] for j in range(3)})
        for kind,targets in gates:yield from validated_remap(primitive[kind],[mapping[q] for q in targets])


def classical_encoder(mask,z,program):
    F0=np.zeros(288,dtype=np.int64)
    F0[program['chord_registers']]=z
    for wire in program['tree_offset_flip_wires']:
        index=wire-FIELD_BASE;F0[index//4]^=1<<(index%4)
    for block in program['cycle_shears']:
        F0[block['target_edge']]=(F0[block['target_edge']]+block['sign']*F0[block['source_edge']])%16
    for block in program['occupation_field_additions']:
        F0[block['target_edge']]=(F0[block['target_edge']]+block['sign']*((mask>>block['occupation'])&1))%16
    return F0


def run():
    paths=[Path(__file__),Path(ref.__file__),Path(native_cycle.__file__),Path(compiler.__file__),
           Path(configuration.__file__),Path(cert.__file__),HERE.parent/'vector_field/RESET_DOMAIN.json',
           ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    values=ref.rational_gaussian_amplitudes(16)
    target=np.array(values,dtype=float)/math.sqrt(sum(v*v for v in values))
    gaussian_word=compiler.compile_real_state(target)
    gaussian_error,planes=native_cycle.gray_preparation_bound(gaussian_word,values)
    recovery=ref.dp.spin.recovery
    executed=recovery.execute_complete(gaussian_word,4,np.eye(16,dtype=complex))
    gaussian_residual=float(np.linalg.norm(executed.success[:,0]-target))
    assert gaussian_residual<float(gaussian_error)+2e-13
    program=prepare_blocks()
    geometry=ref.dp.field.current_walk.geometry()
    encoded_C=np.zeros((288,193),dtype=np.int64)
    for j,e in enumerate(program['chord_registers']):encoded_C[e,j]=1
    chord_index={e:j for j,e in enumerate(program['chord_registers'])}
    for block in program['cycle_shears']:
        encoded_C[block['target_edge'],chord_index[block['source_edge']]]+=block['sign']
    assert np.array_equal(encoded_C,geometry['C'])
    encoded_P=np.zeros((288,32),dtype=np.int64)
    for block in program['occupation_field_additions']:
        encoded_P[block['target_edge'],block['occupation']]+=block['sign']
    assert np.array_equal(encoded_P,geometry['P'][:,np.repeat(np.arange(8)*12,4)])
    root=json.loads((ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    primitive_error={k:F(v) for k,v in root['primitive_operator_errors_exact'].items()}
    adder_counts=Counter(k for k,t in nf.add_program())
    increment_counts=Counter(k for k,t in increment_program())
    counts=Counter()
    counts['flip']=len(program['matter_seed_occupied_modes'])+len(program['tree_offset_flip_wires'])
    for k,n in adder_counts.items():counts[k]+=n*len(program['cycle_shears'])
    for k,n in increment_counts.items():counts[k]+=n*len(program['occupation_field_additions'])
    error=193*gaussian_error+sum((n*primitive_error[k] for k,n in counts.items()),F())
    # Exact pi templates retain every Jordan-Wigner spectator and scalar.
    slater_error=F();slater_pulses=0
    mb=ref.dp.field.mb
    for a,b,sign in ref.pairs():
        first=['I']*(b+1);second=first.copy()
        first[a],first[b]='X','Y';second[a],second[b]='Y','X'
        for k in range(a+1,b):first[k]=second[k]='Z'
        template=cert.pauli(''.join(first),cert._pi(-sign,8))+cert.pauli(''.join(second),cert._pi(sign,8))
        native=mb.givens_word(a,b,sign*math.pi/4)
        slater_error+=nf.root_certificate.exact_template_bound(native,template)
        slater_pulses+=len(native)
    local=compiler.execute(mb.givens_word(0,1,math.pi/4),2)
    assert np.linalg.norm(local[:,1]-np.array([0,1,1,0])/np.sqrt(2))<2e-13
    primitives=compiled.controls()['primitive_words']
    for kind,targets in (('flip',[91]),('cx',[17,48]),('ccx',[17,48,1184])):
        assert validated_remap(primitives[kind],targets)==compiler.remap(primitives[kind],targets)
    assert validated_remap(gaussian_word,[32,33,34,35])==compiler.remap(gaussian_word,[32,33,34,35])
    error+=slater_error
    # Exhaustively test this new controlled increment's clean/dirty extension.
    actual=compiled.execute_reversible(increment_program(),np.arange(256,dtype=np.int64))
    expected=np.array([(x&1)+((((x>>1)&15)+(x&1))%16)*2 for x in range(32)])
    assert np.array_equal(actual[:32],expected)
    assert len(np.unique(actual))==256 and np.array_equal(actual>>5,np.arange(256)>>5)
    assert np.array_equal(compiled.execute_reversible(list(reversed(increment_program())),actual),np.arange(256))
    # All field basis columns and selected full Fock records bind the whole
    # reversible encoder to the exact Gauss formula, without dense2^1187.
    rng=np.random.default_rng(141413)
    laws,_=ref.complete_occupation_law()
    encoder_tests=0
    for index in (0,1,17,511,4095,65535):
        mask=laws[index][0]
        for _ in range(3):
            z=rng.integers(0,16,size=193)
            assert np.array_equal(classical_encoder(mask,z,program),ref.dressed_field(mask,z,ref.reference_charge()))
            encoder_tests+=1
    # The generator emits the full native source word; the count check
    # detects omitted/remapped blocks without storing a multi-million tape.
    count=maximum=0;digest=hashlib.sha256()
    for pulse in word(program,gaussian_word):
        count+=1;maximum=max(maximum,max(pulse['targets']))
        digest.update(json.dumps(pulse,sort_keys=True,separators=(',',':')).encode())
        digest.update(b'\n')
    primitives=compiled.controls()['primitive_words']
    expected_count=slater_pulses+193*len(gaussian_word)+sum(n*len(primitives[k]) for k,n in counts.items())
    assert count==expected_count and maximum<=1186
    p0=min(map(F,configuration.success_probabilities(4,32,full=True,mode='reset',target=circuit.decode(0))))
    pr=min(map(F,configuration.success_probabilities(4,32,full=True,mode='reset',target=circuit.decode(pointer.READY_CODE))))
    reset=1187*(1-p0)+1184*(1-pr)
    pointer_error=nf.root_certificate.exact_template_bound(
        pointer._local_word(),cert.cnot(2,0)+cert.cnot(0,1))
    total=min(F(1),error+reset+1184*pointer_error)
    program.update(gaussian_word=gaussian_word,gaussian_amplitude_integers=values,
                   increment_local_gates=increment_program(),adder_local_gates=nf.add_program(),
                   primitive_words=primitives,complete_word_sha256=digest.hexdigest())
    (HERE/'PREPARATION_PROGRAM.json').write_text(json.dumps(program,indent=2)+'\n')
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result={'schema':'oph.dirac.gauss_reference_native_preparation.v1',
      'checks':{'native_Gaussian_all16columns_and_export':True,'outward_integer_amplitude_certificate':True,
               'all16native_Slater_orbitals_with_CAR_strings':True,'controlled_increment_clean_dirty_columns':True,
               'complete_shared193cycle_affine_encoder':True,'full_native_preparation_word_streamed':True,
               'correlated_reset_and_complete_data_readout_bound':True,'source_bytes_unchanged':True},
      'allocation':{'Fock_modes':32,'field_qubits':1152,'clean_scratch':3,'source_total':1187,
                    'pointer_auxiliaries':[1187,1188],'clock_probe_included':False},
      'native_preparation_pulses':count,'maximum_source_wire':maximum,
      'full_word_sha256':digest.hexdigest(),'gaussian_pulses_per_cycle':len(gaussian_word),
      'gaussian_native_residual':gaussian_residual,'gaussian_source_error_exact':str(gaussian_error),
      'full_preparation_source_error_exact':str(error),'full_preparation_source_error_upper':float(error),
      'native_blank_reset_success_lower_exact':str(p0),'native_ready_reset_success_lower_exact':str(pr),
      'complete_reset_timeout_upper_exact':str(reset),'pointer_per_bit_error_exact':str(pointer_error),
      'complete_preparation_data_record_TV_bound_exact':str(total),
      'complete_preparation_data_record_TV_bound_upper':float(total),'encoder_basis_tests':encoder_tests,
      'complete_data_law':'For every32occupation outcome and288Z16field outcome retain|Slater(n)|² product_j|g(z_j)|² whenF=P(n−G)+Cz; all otherbasis outcomes have zero idealweight. Native reset/export/pointer failure and fine records are additional outcomes. No success normalization.',
      'scope':'Actual native preparation and full terminal occupation/field instrument on a distinct32modeFock layout. FullFock Dirac evolution and an interacting ground/vacuum are not certified by this preparation.',
      'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'NATIVE_PREPARATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'pulses':count,'source_error':float(error),'complete_TV':float(total)}))


if __name__=='__main__':run()
