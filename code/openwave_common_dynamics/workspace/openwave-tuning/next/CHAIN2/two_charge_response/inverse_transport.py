"""Literal native terminal inverse of the frozen full shared-current S."""
from collections import Counter
import hashlib
import json
import two_charge as t


def native_inverse_transport_word():
    mb=t.af.mb;compiler=mb.compiler
    program=mb.mode_program();block=t.af.mode_field.field_block(16)
    words=mb.compiled.controls()['primitive_words']
    for a,b,edge in reversed(program['transport_pairs']):
        # The forward native_field.complete_word emits current then swap.
        yield from compiler.inverse(mb.fock_swap_word(a,b))
        mapping={q:2048+4*edge+q for q in range(4)}
        mapping.update({4:a,5:b,6:3200,7:3201,8:3202})
        for kind,targets in reversed(block['program']['gates']):
            yield from compiler.inverse(compiler.remap(words[kind],[mapping[q]for q in targets]))


def run():
    digest=hashlib.sha256();count=0;largest=-1;kinds=Counter()
    for event in native_inverse_transport_word():
        line=(json.dumps(event,sort_keys=True,separators=(',',':'))+'\n').encode()
        digest.update(line);count+=1;largest=max(largest,*event['targets']);kinds[event['kind']]+=1
    base=json.loads((t.INPUT/'openwave-tuning/common_dynamics/vector_field/COMPOSITION.json').read_text())
    expected=base['source_pulses_per_matter_current_layer']-128*len(t.af.mb.fock_reflection_word())
    assert count==expected and largest==3202
    result={'schema':'oph.two_charge.native_inverse_transport.v1','candidate_id':t.MODEL,
            'generator':'inverse_transport.native_inverse_transport_word',
            'literal_JSONL_encoding':'json.dumps(event,sort_keys=True,separators=(comma,colon))+newline in chronological order',
            'native_pulse_count':count,'native_pulse_kind_counts':dict(kinds),
            'largest_source_wire':largest,'source_wires':3203,
            'emitted_word_sha256':digest.hexdigest(),
            'construction':'Reverse all288frozen transport pairs. For each pair emit inverse full-CAR swap, then reverse and invert every native primitive of its shared Z16current block. This is the literal inverse of native_field.complete_word transport segment, with all phase pulses retained.',
            'error_binding':'At most the frozen COMPOSITION complete_matter_current_layer_error_exact, which conservatively includes the omitted onsite reflection error.',
            'quantum_reference':'Full-Fock source domain and arbitrary spectator references; clean source work and all complements as frozen composition.',
            'clock':'No clock word or clock work wires used.'}
    (t.HERE/'INVERSE_TRANSPORT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':run()
