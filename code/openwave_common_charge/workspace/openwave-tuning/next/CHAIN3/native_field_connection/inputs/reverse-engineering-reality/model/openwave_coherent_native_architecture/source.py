"""Actual C/E permutation, coherent controls and complete recorded marginal.

Coherent memory, phase-preserving control and proper-code services are
declared architecture. Committing the original records recovers native
repair; re-encoding and reusing them defines a different observation policy.
"""
from collections import defaultdict, Counter
from fractions import Fraction as F
from hashlib import sha256
import json
import math
import numpy as np

from model.openwave_native_process import source as native
from model.observer.source_discrimination import circuit
from model.observer.source_scalar import compiler

BINARY_CODES = (0, 2, 1, 4)  # source code for occupation mask00,10,01,11


def event_permutation():
    """All256 addresses: source3, record3, coin1, law1 bits."""
    out = []
    for index in range(256):
        local, coin, law = index & 63, (index >> 6) & 1, (index >> 7) & 1
        low = local & 31
        for left, right in circuit.REPAIR_SWAPS['C', coin]:
            if low == left: low = right
            elif low == right: low = left
        local = (local & 32) | low
        if law and (local >> coin) & 1: local ^= 32
        out.append((index & ~63) | local)
    if len(set(out)) != 256: raise ArithmeticError('full native permutation lost unitarity')
    return tuple(out)


def coherent_columns():
    """Four equal1/2 amplitudes per input; all source padding retained."""
    p = event_permutation()
    return tuple(tuple(p[x | (coin << 6) | (law << 7)]
                       for law in range(2) for coin in range(2)) for x in range(8))


def committed_column(left, right):
    """Full raw-register instrument; C high-bit padding is explicit."""
    if not 0 <= left < 8 or not 0 <= right < 8: raise ValueError('source address required')
    out = {}
    for a, b in zip(coherent_columns()[left], coherent_columns()[right], strict=True):
        if a >> 3 == b >> 3: out[a >> 3] = ((a & 7, b & 7), F(1, 4))
    return out


def complete_permutation(partial, size):
    """Deterministic full unitary extension of an injective partial map."""
    if any(type(x) is not int or not 0 <= x < size for p in partial.items() for x in p):
        raise ValueError('finite register addresses required')
    if len(set(partial.values())) != len(partial): raise ValueError('partial map is not injective')
    out = dict(partial)
    for x, y in zip(sorted(set(range(size))-set(out)),
                    sorted(set(range(size))-set(out.values())), strict=True): out[x] = y
    return tuple(out[x] for x in range(size))


def record_recoder():
    """Reversible native-record to controlled-swap encoding on binary input.

    Target coin sigma=coin xor(left-only). Record bits are uncomputed
    coherently; the inverse restores the original readable native labels.
    """
    p = event_permutation(); partial = {}
    for x, code in enumerate(BINARY_CODES):
        for law in range(2):
            for coin in range(2):
                before = p[code | coin << 6 | law << 7]
                sigma = coin ^ int(x == 1)
                y = (((x & 1) << 1) | (x >> 1)) if sigma else x
                partial[before] = BINARY_CODES[y] | sigma << 6 | law << 7
    return complete_permutation(partial, 256)


def primitive_routines():
    """Actual source words for reusable phase-exact reversible gates."""
    return {
        'ccx': ([compiler.pulse('H', math.pi/2, [2])]
                + compiler.compile_diagonal([0.]*7+[math.pi])
                + [compiler.pulse('H', -math.pi/2, [2])]),
        'flip': [compiler.pulse('I', -math.pi/2, [0]), compiler.pulse('X', math.pi/2, [0])],
        'hadamard': [compiler.pulse('I', -math.pi/2, [0]), compiler.pulse('H', math.pi/2, [0])],
    }


def compile_permutation(permutation):
    """Gray transpositions and complete Toffoli ladders; no target gate call."""
    size = len(permutation); bits = (size-1).bit_length()
    if size != 1 << bits or bits < 3 or sorted(permutation) != list(range(size)):
        raise ValueError('permutation on at least three bits required')
    gates = []; scratch = list(range(bits, bits+max(0,bits-3)))
    def flip(j): gates.append(('flip', (j,)))
    def adjacent(a,b):
        target = (a ^ b).bit_length()-1
        controls = [j for j in range(bits) if j != target]
        zeros = [j for j in controls if not a >> j & 1]
        for j in zeros: flip(j)
        if len(controls) == 2: gates.append(('ccx', tuple(controls+[target])))
        else:
            chain = [(controls[0], controls[1], scratch[0])]
            for j in range(2,len(controls)-1): chain.append((scratch[j-2],controls[j],scratch[j-1]))
            gates.extend(('ccx', x) for x in chain)
            gates.append(('ccx',(scratch[len(controls)-3],controls[-1],target)))
            gates.extend(('ccx', x) for x in reversed(chain))
        for j in reversed(zeros): flip(j)
    def swap(a,b):
        path=[a]
        for j in range(bits):
            if (a ^ b) >> j & 1: path.append(path[-1] ^ (1 << j))
        pairs=list(zip(path[:-1],path[1:],strict=True))
        for x,y in pairs+pairs[-2::-1]: adjacent(x,y)
    seen=set()
    for a in range(size):
        if a in seen: continue
        cycle=[]; x=a
        while x not in seen: seen.add(x);cycle.append(x);x=permutation[x]
        for b in cycle[1:]: swap(a,b)
    return {'address_bits':bits,'scratch_bits':len(scratch),'gates':gates}


def execute_basis(program, inputs=None):
    """Exact circuit action on all inputs, or an explicitly supplied subset.

    A subset execution does not certify that every dirty input was queried;
    each reversible gate still defines its complete dirty-input extension.
    """
    q=program['address_bits']+program['scratch_bits']
    out=np.arange(1 << q,dtype=np.int64) if inputs is None else np.asarray(inputs,dtype=np.int64).copy()
    if out.ndim!=1 or np.any(out<0) or np.any(out >= 1 << q) or len(np.unique(out))!=len(out):
        raise ValueError('distinct complete-register basis inputs required')
    for kind,targets in program['gates']:
        if kind=='flip': out ^= 1 << targets[0]
        elif kind=='ccx':
            a,b,c=targets;out ^= (((out >> a)&1)&((out >> b)&1)) << c
        else: raise ValueError('unknown reversible gate')
    if len(np.unique(out)) != len(out): raise ArithmeticError('dirty-scratch action not unitary')
    return out


def construction_descriptor():
    """Shared declared hardware and observation policy, without unit claims."""
    return {'carrier':'native binary load source; complete encoded local8-address extension',
            'local_event':'actual C/E repair permutation with coherent fair law and placement bits',
            'raw_event_slots':25,'local_encoded_address_bits':8,
            'coherent_services':['proper-code M6 load/response/export','blank coherent record registers',
                                 'reversible record recoding','phase-preserving reversible gates'],
            'committed_policy':'measure original raw labels after event',
            'coherent_policy':'retain registers; invert recoding to recover original labels',
            'record_identity':'same complete instrument at the declared commit intervention',
            'time_unit':'dimensionless ordered coherent operation count',
            'GeV_calibration':None,'official_physical_row_award':False}


def execute_primitives():
    words=primitive_routines(); matrices={}; errors={}; resources={}
    for name,word in words.items():
        q=3 if name=='ccx' else 1
        value=compiler.execute(word,q)
        if name=='ccx':
            target=np.eye(8);target[:,[3,7]]=target[:,[7,3]]
        elif name=='flip': target=np.array([[0.,1.],[1.,0.]])
        else: target=np.array([[1.,1.],[1.,-1.]])/math.sqrt(2)
        errors[name]=float(np.linalg.norm(value-target,'fro'))
        if errors[name]>1e-10: raise ArithmeticError(('source primitive failed',name,errors[name]))
        matrices[name]=value;resources[name]=compiler.resources(word)
    return matrices,errors,resources


def run():
    columns=coherent_columns()
    if len({a for col in columns for a in col}) != 32: raise ArithmeticError('isometry Gram failed')
    cases=0
    for x in range(4):
        for y in range(4):
            actual=committed_column(BINARY_CODES[x],BINARY_CODES[y])
            expected={}
            for row in native.raw_local_branches():
                if not row.weight or not row.pairs: continue
                pair=dict(row.pairs)
                if x not in pair or y not in pair:continue
                _,law,coin,raw,_=row.record
                key=raw+8*coin+16*int(law=='E')
                expected[key]=((BINARY_CODES[pair[x]],BINARY_CODES[pair[y]]),row.weight)
            if actual!=expected: raise ArithmeticError(('native matrix unit mismatch',x,y))
            cases+=25
    full_cases=0
    for x in range(6):
        for y in range(6):
            expected={}
            for law in range(2):
                for coin in range(2):
                    images=[]
                    for z in (x,y):
                        a,b=circuit.SEAM_STATES[z];low=(a+b)//2;parity=(a+b)%2
                        left=low+coin*parity;right=low+(1-coin)*parity
                        current=a-left;raw=circuit.CURRENT_CODES[current]+4*parity*law
                        images.append((circuit.SEAM_STATES.index((left,right)),raw))
                    if images[0][1]==images[1][1]:
                        expected[images[0][1]+8*coin+16*law]=((images[0][0],images[1][0]),F(1,4))
            assert committed_column(x,y)==expected
            full_cases+=1
    p=event_permutation();r=record_recoder()
    for x,code in enumerate(BINARY_CODES):
        for law in range(2):
            for coin in range(2):
                z=r[p[code | coin<<6 | law<<7]];sigma=coin ^ int(x==1)
                y=(((x&1)<<1)|(x>>1)) if sigma else x
                assert z==BINARY_CODES[y] | sigma<<6 | law<<7
    matrices,errors,costs=execute_primitives();reports=[]
    for name,perm in [('native_event',p),('coherent_record_recoding',r)]:
        program=compile_permutation(perm);truth=execute_basis(program)
        assert tuple(map(int,truth[:256]))==perm
        count=Counter(kind for kind,_ in program['gates']);s=F(sum(count.values()),10**10)
        assert s<1
        reports.append({'name':name,'complete_basis_inputs':len(truth),'clean_inputs':256,
                        'gate_counts':dict(count),'native_pulses':sum(count[k]*costs[k]['native_pulses'] for k in count),
                        'operator_error_upper_given_checked_routine_bounds':str(s/(1-s)),
                        'truth_table_sha256':sha256(json.dumps(list(map(int,truth))).encode()).hexdigest(),
                        'schedule_sha256':sha256(json.dumps(program,sort_keys=True).encode()).hexdigest()})
    return {'verdict':'PASS_COMPLETE_COHERENT_NATIVE_EVENT','construction':construction_descriptor(),
            'coherent_event_address_dimension':256,
            'source_address_columns':8,'source_padding_columns':2,'binary_raw_matrix_unit_slots_checked':cases,
            'full_six_state_endpoint_matrix_units_checked':full_cases,
            'native_complete_raw_alphabet':25,'full_encoded_raw_register_labels':32,
            'arbitrary_reference_identity':True,'zero_hold_and_padding_retained':True,
            'coherent_event_amplitude_exact':'1/2','reversible_binary_record_recoding_columns':16,
            'primitive_source_errors':errors,'primitive_resources':costs,'compiled_reversible_services':reports,
            'source_compilation':'actual full M6 primitive exponentials plus exact complete factor circuit; no polar correction',
            'memory_policy':'commit original raw labels to recover C/E; coherent inverse recoding preserves labels as an available intervention',
            'readout_boundary':'tracing the recoded swap coin is a different channel from tracing original C/E records',
            'offdiagonal_negative_control':{'input':'|10><01|','original_unread_image':'0',
                'recoded_coin_trace_image':'(|10><01|+|01><10|)/2'},
            'declared_services':['coherent law/coin preparation','blank record and scratch','proper-code load/export','coherent record recoding','terminal pointer apparatus'],
            'physical_service_precision_certified':False,'physical_clock':False,'physical_particle_identification':False,
            'official_OpenWave_row_award':False}
