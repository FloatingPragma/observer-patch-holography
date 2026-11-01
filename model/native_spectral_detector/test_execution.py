import numpy as np
import pytest
from dataclasses import replace

from model.protected_particle_readout import execute
from . import contraction,execution as e,circuits as c


def test_fused_contraction_matches_independent_full_tensor_arithmetic():
    rng=np.random.default_rng(911)
    for bits in (3,5,7):
        for targets in ((0,),(2,0),(1,2,0)):
            for columns in (1,5,32):
                state=np.ascontiguousarray(rng.normal(size=(1<<bits,columns))+1j*rng.normal(size=(1<<bits,columns)))
                d=1<<len(targets);matrix=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
                value,axes,shape=execute.arranged(state,targets,bits);flat=value.reshape(d,-1)
                density=flat@flat.conj().T
                target=(matrix@flat).reshape(shape).transpose(np.argsort(axes)).reshape(state.shape)
                actual=contraction.apply_inplace(state,matrix,targets,bits)
                assert np.allclose(actual,density,atol=4e-12,rtol=2e-14)
                assert np.allclose(state,target,atol=1e-13,rtol=2e-14)


def test_complete_column_partition_keeps_all_failure_mass():
    rng=np.random.default_rng(23);state=rng.normal(size=(8,9))+1j*rng.normal(size=(8,9));state/=np.linalg.norm(state)
    blocks=((c.charge_blocks()[0][0],(0,1,2)),)
    whole=e.complete(blocks,3,state)
    parts=[e.complete(blocks,3,state[:,start:start+3]) for start in (0,3,6)]
    assert np.allclose(whole.success,np.column_stack([part.success for part in parts]),atol=1e-14)
    assert np.allclose(whole.failure_weights,sum(part.failure_weights for part in parts),atol=1e-43)
    assert whole.normalization_error<1e-12
    # A retained first-failure vector reconstructs its stored slot weight.
    failed=whole.failure(0,0,0)
    assert abs(np.linalg.norm(failed)**2-whole.failure_weights[0,0])<1e-43


def test_reject_masked_or_incomplete_or_aliasing_layouts():
    with pytest.raises(ValueError):e.complete((),3,np.zeros((7,1)))
    with pytest.raises(ValueError):e.complete((),3,np.ma.array(np.zeros((8,1))))
    with pytest.raises(ValueError):contraction.layout(3,(1,1))
    with pytest.raises(ValueError):contraction.apply_inplace(np.zeros((8,1)),np.eye(2),(0,),3)
    contraction.layout(3,(1,))
    with pytest.raises(ValueError):contraction.layout(3,(True,))
    with pytest.raises(ValueError):contraction.layout(3,[1])


def test_nonfinite_native_matrix_and_effects_fail_closed():
    state=np.zeros((8,1),complex);state[0]=1
    bad=np.eye(2,dtype=complex);bad[0,0]=np.nan
    with pytest.raises(ValueError,match='finite local matrix'):
        contraction.apply_inplace(state,bad,(0,),3)
    instrument=c.charge_blocks()[0][0]
    effects=np.array(instrument.effects);effects[0,0,0,0]=np.nan
    with pytest.raises(ArithmeticError,match='nonfinite native failure'):
        e.complete(((replace(instrument,effects=effects),(0,1,2)),),3,state)
