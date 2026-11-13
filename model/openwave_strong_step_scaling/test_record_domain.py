"""Exact carrier support retains the full numerical distribution unchanged."""
import numpy as np
import pytest
from model.openwave_nuclear_distance.spatial import carrier
from .certificate import compatible_native_occupations,complete_success_records


@pytest.mark.parametrize('length',[2,3,4,5])
def test_exact_domain_covers_all_block_vectors_and_zero_tiny_outcomes(length):
    d=carrier(length)
    nu=d['counts'][:,[j for j,o in enumerate(d['orbitals']) if o[2]==0]].sum(axis=1)
    ns=d['counts'][:,[j for j,o in enumerate(d['orbitals']) if o[1]==0]].sum(axis=1)
    ids=np.flatnonzero((nu==1)&(ns==2))
    support=compatible_native_occupations(d,ids)
    assert len(support)==9*length**3
    outside=np.ones(len(d['full']),dtype=bool);outside[support]=False
    B=d['B'][:,ids].tocsr()
    assert B[outside].nnz==0  # Exact structural zero for every state and time.
    state=np.zeros(len(ids),complex);state[0]=1;state[-1]=1e-25j
    full=np.asarray(B@state).ravel();dimension=len(d['labels'])
    rows=complete_success_records(full,support,dimension)
    dense=np.zeros(len(full))
    for i,p in rows:dense[i]=p
    expected=np.array([float(abs(x)**2/dimension) for x in full])
    assert np.array_equal(dense,expected)
    assert len(rows)==len(support) and any(p==0 for i,p in rows)
    assert any(0<p<1e-40 for i,p in rows)
    # No change in the complete complementary failure distribution.
    diagonal=np.asarray(abs(d['B']).power(2).sum(axis=1)).ravel()/dimension
    assert np.array_equal(diagonal-dense,diagonal-expected)
