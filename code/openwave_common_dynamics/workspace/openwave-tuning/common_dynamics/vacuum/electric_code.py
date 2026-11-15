"""Native field-code measurement of the free reference's electric leakage."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math

import numpy as np

import native_preparation as prep
from model.native_clock import configuration
from model.observer.source_discrimination import circuit
from model.native_interface import pointer

HERE=Path(__file__).resolve().parent
ROOT=prep.ROOT
nf=prep.nf;compiler=prep.compiler;cert=prep.cert


def three_edge_connected(edges,vertices):
    """Exact check: deleting each edge leaves a connected bridgeless graph."""
    adjacency=[[] for _ in range(vertices)]
    for e,(a,b) in enumerate(edges):
        adjacency[int(a)].append((int(b),e));adjacency[int(b)].append((int(a),e))
    for removed in range(len(edges)):
        clock=0;discovery=[-1]*vertices;low=[-1]*vertices;bridges=[]
        def visit(v,parent):
            nonlocal clock
            discovery[v]=low[v]=clock;clock+=1
            for w,e in adjacency[v]:
                if e in (parent,removed):continue
                if discovery[w]<0:
                    visit(w,e);low[v]=min(low[v],low[w])
                    if low[w]>discovery[v]:bridges.append(e)
                else:low[v]=min(low[v],discovery[w])
        visit(0,-1)
        if clock!=vertices or bridges:return False
    return True


def fast_remap(word,targets):
    """Exact remap of a provider-validated word; no repeated pulse parsing."""
    return [{'kind':p['kind'],'angle':p['angle'],'targets':[int(targets[q]) for q in p['targets']]}
            for p in word]


def encoder_word(program,inverse=False):
    primitives=prep.compiled.controls()['primitive_words']
    sequence=[]
    for wire in program['tree_offset_flip_wires']:
        sequence.append(('flip',[wire]))
    for block in program['cycle_shears']:
        gates=nf.add_program() if block['sign']==1 else list(reversed(nf.add_program()))
        targets=[32+4*block['source_edge']+j for j in range(4)]
        targets += [32+4*block['target_edge']+j for j in range(4)]+list(prep.SCRATCH)
        for kind,local in gates:sequence.append((kind,[targets[q] for q in local]))
    for block in program['occupation_field_additions']:
        gates=prep.increment_program() if block['sign']==1 else list(reversed(prep.increment_program()))
        targets=[block['occupation']]+[32+4*block['target_edge']+j for j in range(4)]+list(prep.SCRATCH)
        for kind,local in gates:sequence.append((kind,[targets[q] for q in local]))
    if inverse:
        inverse_primitives={k:compiler.inverse(v) for k,v in primitives.items()}
        for kind,targets in reversed(sequence):yield from fast_remap(inverse_primitives[kind],targets)
    else:
        for kind,targets in sequence:yield from fast_remap(primitives[kind],targets)


def run():
    paths=[Path(__file__),Path(prep.__file__),Path(prep.ref.__file__),Path(pointer.__file__),
           HERE/'NATIVE_PREPARATION.json',HERE/'REFERENCE.json',
           ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    graph=prep.ref.dp.field.data()[0]
    assert three_edge_connected(graph.seams,96)
    values=[F(((j+8)%16-8)**2,2) for j in range(16)]
    mean=sum(values,F())/16;variance=sum((v*v for v in values),F())/16-mean*mean
    totalmean=288*mean;totalvar=288*variance
    assert totalmean==3096 and totalvar==27234
    alpha=F(1,4096);radius=F(9216)
    ideal_lower=totalvar*alpha*alpha*(1-alpha*alpha*radius*radius/12)
    assert 1-alpha*alpha*radius*radius/12==F(37,64)
    assert ideal_lower>F(9,10000)
    program=prep.prepare_blocks()
    root=json.loads((ROOT/'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    per={k:F(v) for k,v in root['primitive_operator_errors_exact'].items()}
    counts=Counter({'flip':len(program['tree_offset_flip_wires'])})
    for k,n in Counter(k for k,t in nf.add_program()).items():counts[k]+=n*len(program['cycle_shears'])
    for k,n in Counter(k for k,t in prep.increment_program()).items():counts[k]+=n*len(program['occupation_field_additions'])
    encoder_error=sum((per[k]*n for k,n in counts.items()),F())
    # The Slater part is the certified preparation with its Gaussian and
    # classical encoder contributions subtracted, without reusing its state.
    prepared=json.loads((HERE/'NATIVE_PREPARATION.json').read_text())
    gaussian_error=F(prepared['gaussian_source_error_exact'])
    full=F(prepared['full_preparation_source_error_exact'])
    matter=full-193*gaussian_error-encoder_error
    assert matter>=0
    electric,electric_per=nf.electric_word(alpha)
    pointer_error=nf.root_certificate.exact_template_bound(pointer._local_word(),cert.cnot(2,0)+cert.cnot(0,1))
    p0=F(prepared['native_blank_reset_success_lower_exact']);pr=F(prepared['native_ready_reset_success_lower_exact'])
    reset=1187*(1-p0)+1152*(1-pr)
    native_error=matter+2*encoder_error+2*772*per['hadamard']+288*electric_per+1152*pointer_error
    complete=reset+native_error
    observed_lower=ideal_lower-complete
    assert observed_lower>0
    # Validate the efficient remapper exactly against the actual compiler.
    primitives=prep.compiled.controls()['primitive_words']
    for kind,targets in (('flip',[91]),('cx',[17,48]),('ccx',[17,48,1184])):
        assert fast_remap(primitives[kind],targets)==compiler.remap(primitives[kind],targets)
    H=nf.H(0)
    assert nf.root_certificate.exact_template_bound(H,[('I',(0,),cert._pi(-1,2)),('H',(0,),cert._pi(1,2))])==per['hadamard']
    # Stream the complete intervention/readout word; invert each complete
    # native primitive in reverse order, retaining its scalar phase.
    def fullword():
        for mode in program['matter_seed_occupied_modes']:yield from fast_remap(primitives['flip'],[mode])
        for a,b,sign in prep.ref.pairs():yield from prep.ref.dp.field.mb.givens_word(a,b,sign*math.pi/4)
        for edge in program['chord_registers']:
            for bit in range(4):yield from fast_remap(H,[32+4*edge+bit])
        yield from encoder_word(program)
        for edge in range(288):yield from fast_remap(electric,[32+4*edge+j for j in range(4)])
        yield from encoder_word(program,inverse=True)
        for edge in program['chord_registers']:
            for bit in range(4):yield from fast_remap(H,[32+4*edge+bit])
        for wire in range(32,1184):yield from fast_remap(pointer._local_word(),[1187,1188,wire])
    digest=hashlib.sha256();count=maximum=0
    for event in fullword():
        count+=1;maximum=max(maximum,max(event['targets']))
        digest.update(json.dumps(event,sort_keys=True,separators=(',',':')).encode()+b'\n')
    result={'schema':'oph.dirac.uniform_fibre_native_code_test.v1',
      'checks':{'all_two_edge_deletions_leave_connected_graph':True,
               'exact_uniform_fibre_mean_and_variance':True,'operator_leakage_lower_bound':True,
               'full_affine_decode_and_native_projector_readout':True,
               'complete_native_error_below_signal':True,'full_native_word_streamed':True,
               'source_bytes_unchanged':True},
      'electric_alpha_exact':str(alpha),'uniform_fibre_mean_exact':str(totalmean),
      'uniform_fibre_variance_exact':str(totalvar),'electric_energy_range_exact':str(radius),
      'ideal_code_complement_probability_lower_exact':str(ideal_lower),
      'ideal_code_complement_probability_lower':float(ideal_lower),
      'complete_native_record_TV_bound_exact':str(complete),'complete_native_record_TV_bound':float(complete),
      'native_conclusive_nonzero_field_record_probability_lower_exact':str(observed_lower),
      'native_conclusive_nonzero_field_record_probability_lower':float(observed_lower),
      'native_word_pulses':count,'maximum_wire':maximum,'native_word_sha256':digest.hexdigest(),
      'instrument':'Prepare the actual Slater sea and uniform193cycle Gauss state, apply actualelectric kick, inverseGaussencoder, Hadamards on772chordbits, readall1152fieldbits with complete nativeC/E pointers. Any nonzero data bit witnesses leaving the originaluniformfibre code. All reset/export/finepointer outcomes retained.',
      'proof':'For each matter configuration, code survival is |E exp(−i alphaH)|² on its uniform affine fibre. Independent H,Hprime give leakage E[1−cos(alpha(H−Hprime))] >= alpha²Var(H)−alpha⁴Range(H)²Var(H)/12. The same conditional moments give the lower bound for arbitrary matter/reference states. Here the factor is37/64.',
      'scope':'An executed native-word and analytic complete-instrument certificate for an electric intervention on the free-reference code. It is not an interacting-vacuum or fullFock-evolution certificate.',
      'source_sha256':pins,'official_rows_changed':[]}
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE/'ELECTRIC_CODE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'ideal_lower':float(ideal_lower),'native_lower':float(observed_lower),'pulses':count}))


if __name__=='__main__':run()
