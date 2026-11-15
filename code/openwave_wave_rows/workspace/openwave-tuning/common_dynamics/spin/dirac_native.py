"""Factorized native word for the finite encoded Dirac route macro.

Each edge increment is checked on its full clean and dirty register domain.
This compiles the macro's complete encoded channel, not physical locality of
every gate prefix. Full routing through empty physical transit modes is a
separate locality implementation obligation.
"""
from collections import Counter
from pathlib import Path
import hashlib
import json
import numpy as np

import dirac_proposal as dp

HERE=Path(__file__).resolve().parent
SCRATCH=tuple(range(9,14))


def edge_program(axis,edge,charge=1):
    gates=[]
    for cell in range(8):
        for rail,direction in ((0,1),(1,-1)):
            current=charge*int(dp.route(cell,axis,direction)[2][edge])
            assert current in (-1,0,1)
            if not current:continue
            values={q:(cell>>q)&1 for q in range(3)};values[4]=rail
            for target in reversed(range(4)):
                condition=dict(values)
                condition.update({5+q:int(current>0) for q in range(target)})
                gates+=dp.spin.compiled.multi_control(condition,5+target,SCRATCH)
    return gates


def program(charge=1):
    primitives=dp.spin.compiled.controls()['primitive_words']
    axes=[];counts=Counter();clean_checks=0;dirty_checks=0
    for axis in range(3):
        blocks=[]
        for edge in range(288):
            gates=edge_program(axis,edge,charge)
            if not gates:continue
            actual=dp.spin.compiled.execute_reversible(gates,np.arange(1<<14,dtype=np.int64))
            expected=np.arange(1<<9,dtype=np.int64)
            for index in range(1<<9):
                cell=index&7;rail=(index>>4)&1;f=(index>>5)&15
                current=charge*int(dp.route(cell,axis,1 if rail==0 else -1)[2][edge])
                expected[index]=(index&31)+(((f+current)%16)<<5)
            assert np.array_equal(actual[:1<<9],expected)
            assert len(np.unique(actual))==1<<14
            # The ladder uncomputes scratch even when it starts dirty.
            assert np.array_equal(actual>>9,np.arange(1<<14)>>9)
            assert np.array_equal(actual&31,np.arange(1<<14)&31)
            clean_checks+=1<<9;dirty_checks+=1<<14
            counts.update(kind for kind,targets in gates)
            blocks.append(dict(edge=edge,local_gates=gates,
                               local_wires={'cell':[0,1,2],'spin':[3,4],'field':[5,6,7,8],'scratch':list(SCRATCH)}))
        axes.append(dict(axis=axis,basis_word=dp.diagonalizer_word(axis),field_blocks=blocks,
                         cell_flip_bit=2-axis,restore_basis_word=dp.spin.compiler.inverse(dp.diagonalizer_word(axis))))
    pulses=sum(number*len(primitives[kind]) for kind,number in counts.items())
    pulses+=sum(2*len(row['basis_word'])+len(primitives['flip']) for row in axes)
    return dict(charge=charge,axes=axes,primitive_words=primitives,
                data_wires={'cell':[0,1,2],'spin':[3,4],'field_edge_e':'5+4*e through8+4*e'},
                fresh_scratch=list(range(1157,1162)),
                clean_truth_columns=clean_checks,dirty_truth_columns=dirty_checks,
                classical_gate_counts=dict(counts),native_pulses_excluding_mass=pulses)


def native_word(h,m,charge=1):
    """Yield the complete native matter macro without materializing its tape."""
    yield from dp.spin.compiler.remap(dp.mass_word(h,m),[3,4])
    primitives=dp.spin.compiled.controls()['primitive_words']
    for axis in range(3):
        basis=dp.spin.compiler.remap(dp.diagonalizer_word(axis),[3,4])
        yield from basis
        for edge in range(288):
            mapping={q:q for q in range(5)}
            mapping.update({5+q:5+4*edge+q for q in range(4)})
            mapping.update({9+q:1157+q for q in range(5)})
            for kind,targets in edge_program(axis,edge,charge):
                yield from dp.spin.compiler.remap(primitives[kind],[mapping[q] for q in targets])
        yield from dp.spin.compiler.remap(primitives['flip'],[2-axis])
        yield from dp.spin.compiler.inverse(basis)


def run():
    sourcepins={str(path.relative_to(dp.spin.ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (Path(__file__),Path(dp.__file__),Path(dp.spin.compiled.__file__))}
    out=program(1)
    count=sum(1 for _ in native_word(.1,.4,1))
    assert count==out['native_pulses_excluding_mass']+len(dp.mass_word(.1,.4))
    receipt=dict(schema='oph.local12.dirac_native_macro.v1',model_id=dp.MODEL_ID,
                 checks={'all_edge_clean_columns':True,'all_edge_dirty_columns_and_scratch_return':True,
                         'native_word_complete_factorization':True,'source_current_route_match':True},
                 clean_columns=out['clean_truth_columns'],dirty_columns=out['dirty_truth_columns'],
                 matter_native_pulses=count,field_edges=288,encoded_data_qubits=1157,fresh_scratch_qubits=5,
                 scope='Exact reversible truth tables plus source-native primitive words. Full1157-qubit exponential execution, coherent preparation/reset, aggregate physical error certificate and locally Gauss-preserving gate-prefix routing are not supplied by this factorization.',
                 source_sha256=sourcepins,official_rows_changed=[])
    assert sourcepins=={str(path.relative_to(dp.spin.ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in (Path(__file__),Path(dp.__file__),Path(dp.spin.compiled.__file__))}
    (HERE/'DIRAC_NATIVE_PROGRAM.json').write_text(json.dumps(out,indent=2)+'\n')
    (HERE/'DIRAC_NATIVE.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))


if __name__=='__main__':run()
