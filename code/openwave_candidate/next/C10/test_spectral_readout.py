import numpy as np
from joint_source import JointNativeSource
from spectral_readout import GroupedEnergyInstrument


def test_grouped_complex_coherence_and_continuation():
    model=JointNativeSource();instrument=GroupedEnergyInstrument(model)
    group=min((x for x in instrument.groups if x['source_number']==3 and x['degeneracy']==2),key=lambda x:abs(x['energy']+.7613052172833566))
    ids=instrument.arrays[group['position_indices_archive_key']];basis=instrument.arrays[group['basis_archive_key']]
    v=basis@np.array([1,1j])/np.sqrt(2);rho=np.zeros((model.n,model.n),complex);rho[np.ix_(ids,ids)]=np.outer(v,v.conj())
    rows=instrument.measure(rho,'complex_control');assert abs(sum(r['probability'] for r in rows)-1)<1e-12
    row=next(r for r in rows if r['probability']>.99);branch=instrument.reconstruct(row)
    assert np.max(abs(branch-rho))<1e-12
    coefficient=instrument.arrays[row['coefficient_archive_key']];assert np.max(abs(coefficient.imag))>.49
    fine=np.zeros_like(rho);fine[np.ix_(ids,ids)]=basis@np.diag(np.diag(coefficient))@basis.T
    assert abs(np.linalg.norm(branch-fine,ord='nuc')/2-.5)<1e-12
    out,ledger=instrument.continue_branch(branch)
    assert abs(np.trace(out)-1)<1e-12 and abs(sum(r['joint_probability'] for r in ledger)-1)<1e-12
    fineout,_=instrument.continue_branch(fine)
    assert np.linalg.norm(out-fineout,ord='nuc')/2>.1


def test_invalid_group_tolerances_rejected_before_apparatus_creation():
    import pytest
    for value in (float("nan"),float("inf"),0,-1,True,None):
        with pytest.raises(ValueError):GroupedEnergyInstrument(None,value)
