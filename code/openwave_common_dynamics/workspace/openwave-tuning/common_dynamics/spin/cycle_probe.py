"""A native spin readout of an actual Gauss-preserving route-cycle character."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import numpy as np

import dirac_proposal as dp

HERE=Path(__file__).resolve().parent


def cycle_word():
    gates=[]
    # Cell bit2; eigen-rail bit1; Z16 cycle bits3..6; clean scratch7..8.
    # The same signed current is cell-rail. Separate controlled +/-1
    # additions cancel when both controls are one and need one less scratch.
    for control,increment in ((2,1),(1,-1)):
        for target in reversed(range(4)):
            values={control:1}
            values.update({3+q:int(increment>0) for q in range(target)})
            gates+=([('cx',(control,3))] if target==0 else
                    dp.spin.compiled.multi_control(values,3+target,(7,8)))
    gates += [('flip',(2,))]
    actual=dp.spin.compiled.execute_reversible(gates,np.arange(512,dtype=np.int64))
    expected=[]
    for index in range(128):
        cell=(index>>2)&1;rail=(index>>1)&1;t=(index>>3)&15
        increment=int(cell==1 and rail==0)-int(cell==0 and rail==1)
        expected.append((index&3)+((1-cell)<<2)+(((t+increment)%16)<<3))
    assert np.array_equal(actual[:128],expected)
    assert len(np.unique(actual))==512 and np.array_equal(actual>>7,np.arange(512)>>7)
    primitive=dp.spin.compiled.controls()['primitive_words']
    word=dp.diagonalizer_word(0)
    for kind,targets in gates:word+=dp.spin.compiler.remap(primitive[kind],list(targets))
    word+=dp.spin.compiler.inverse(dp.diagonalizer_word(0))
    return gates,word


def embedded_state(theta):
    other,path,jplus=dp.route(0,0,1)
    other2,path2,jminus=dp.route(0,0,-1)
    assert other==other2
    D=jplus-jminus;B=dp.field.data()[1]
    assert not np.any(B@D) and max(abs(D))==1
    assert np.array_equal(dp.route(other,0,1)[2],-jminus)
    assert np.array_equal(dp.route(other,0,-1)[2],-jplus)
    sparse={}
    for t in range(16):
        F=bytes((t*D%16).astype(np.uint8))
        sparse[0,F]=np.array([np.exp(1j*theta*t)/4,0,0,0],complex)
    return sparse,other,jplus,D


def run():
    paths=[Path(__file__),Path(dp.__file__),Path(dp.spin.__file__),Path(dp.spin.compiled.__file__)]
    pins={str(p.relative_to(dp.spin.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    gates,word=cycle_word()
    initial=np.zeros((512,2),complex)
    rows=[]
    for column,theta in enumerate((0.,np.pi)):
        state,other,jplus,D=embedded_state(theta)
        out=dp.axis_step(state,0)
        dp.records(out)
        logical=np.zeros(128,complex)
        for t in range(16):
            F=bytes(((jplus+t*D)%16).astype(np.uint8))
            logical[8*t+4:8*t+8]=out[other,F]
            initial[8*t,column]=np.exp(1j*theta*t)/4
        assert abs(np.vdot(logical,logical)-1)<1e-13
        expected_spin=0 if column==0 else 3
        assert abs(sum(abs(logical[8*t+4+expected_spin])**2 for t in range(16))-1)<1e-13
        rows.append(dict(cycle_character=theta,expected_spin_bits=expected_spin,logical=logical))
    executed=dp.spin.recovery.execute_complete(word,9,initial)
    instruments,pointer=dp.spin.actual_native_pointer()
    for column,row in enumerate(rows):
        target=np.zeros(512,complex);target[:128]=row.pop('logical')
        error=float(np.linalg.norm(executed.success[:,column]-target))
        assert error<2e-11
        values=executed.success[:,column].reshape(256,2)
        rho=values.T@values.conj()
        laws={}
        for law,outcomes in instruments.items():
            bits=[0.,0.];fine=[]
            for meta,kraus in outcomes:
                p=sum(float(np.trace(K@rho@K.conj().T).real) for K in kraus)
                bits[meta['bit']]+=p;fine.append(dict(record=meta,probability=p))
            assert abs(bits[column]-1)<2e-10
            laws[law]=dict(bits=bits,fine_records=fine)
        row.update(native_state_error=error,native_pointer_laws=laws)
    state,other,jplus,D=embedded_state(0.)
    receipt=dict(schema='oph.local12.dirac_native_cycle_probe.v1',model_id=dp.MODEL_ID,
                 checks={'actual_native_closed_route_cycle':True,'exact_Gauss_onecycle_embedding':True,
                         'complete_128_column_clean_route_truth':True,'complete_512_column_dirty_route_truth':True,
                         'native_route_spin_words_executed':True,'native_field_character_changes_spin_readout':True},
                 cycle_edge_current=[(int(e),int(D[e])) for e in np.flatnonzero(D)],
                 position_roots=[0,12*other],cycle_modulus=16,native_qubits=9,native_pulses=len(word),
                 scratch_bits=[7,8],native_export_failure_gram=executed.failure_gram.real.tolist(),
                 pointer_export_failure_trace=float(np.trace(pointer.failure_gram).real),experiments=rows,
                 scope='The invariant two-root/one-cycle restriction is exact for one axis of this declared Dirac matter source. Native dynamics and C/E pointer execute on supplied cycle-character source preparations. A source-native preparation/reset word and whole-macro magnetic-moment measurement are not inferred.',
                 source_sha256=pins,official_rows_changed=[])
    assert pins=={str(p.relative_to(dp.spin.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE/'CYCLE_PROBE_PROGRAM.json').write_text(json.dumps(dict(gates=gates,native_word=word),indent=2)+'\n')
    (HERE/'CYCLE_PROBE.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('checks','native_pulses','native_export_failure_gram')}))


if __name__=='__main__':run()
