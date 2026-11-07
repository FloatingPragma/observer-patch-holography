"""Source exhaustion, complete response, and numerical-scope regression checks."""
from fractions import Fraction
import json
from pathlib import Path
import shutil

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from scipy import sparse

from .interacting_carrier import (HERE, ROOT, NativeCarrier, build_join,
    load_operators, primary_family_matrices, produce_operators)
from .precision_audit import matrix_rounding_error
from .run_joined_response import (addition_data, block_eigh, maxabs,
    operator_checks, schur_checks)


@pytest.fixture(scope='module')
def replayed(tmp_path_factory):
    root=tmp_path_factory.mktemp('joined-source')
    receipt=produce_operators(ROOT,root,3)
    return root,receipt


def test_complete_source_rebuild_preserves_exact_basis_and_operators(replayed):
    root,receipt=replayed
    assert receipt['complete_Gauss_dimensions']=={'0':438,'1':1596}
    assert receipt['complete_operator_leakage_columns_checked']==126108
    for name in ('joined_basis.json','joined_exact_operators.json'):
        assert (root/name).read_bytes()==(HERE/name).read_bytes()
    with np.load(root/'joined_operators.npz',allow_pickle=False) as new, \
         np.load(HERE/'joined_operators.npz',allow_pickle=False) as saved:
        assert set(new.files)==set(saved.files)
        for name in new.files:
            np.testing.assert_allclose(new[name],saved[name],rtol=0,atol=1e-14)


def test_prepared_family_expressions_follow_exact_native_records():
    from ..dependencies.source_reference import exact_source_packet
    source=exact_source_packet(ROOT)
    expressions=json.loads((HERE.parent/'dependencies'/'exact_family_operators.json').read_text())
    for name in ('Y','D','J'):
        expected=sp.Matrix(expressions[name])
        assert all(sp.simplify(x)==0 for x in expected-source[name])


def _representation_dimension(families,q):
    # Independent tensor-product counting: charge, color type, weak dimension.
    irreps=((0,0,1),(-4,-1,1),(-3,0,2),(1,1,2),(2,-1,1),(6,0,1))
    from math import comb
    total=0
    for nh in (0,1):
        ne=q+nh
        for e0 in range(ne+1):
            e1=ne-e0
            for lsite in ((-1,) if nh==0 else (0,1)):
                charge=6*e0-3*(lsite==0)-3*nh-6*q
                multiplicity=0
                for a in irreps:
                    for b in irreps:
                        if a[0]-b[0]!=charge or a[1]!=b[1]:continue
                        if lsite==1:
                            multiplicity+=a[2]!=b[2]
                        elif a[2]==b[2]:
                            multiplicity+=a[2] if nh else 1
                total+=multiplicity*comb(2*families,e0)*comb(2*families,e1)*(2*families if nh else 1)
    return total


@pytest.mark.parametrize('families',range(1,7))
def test_representation_count_includes_full_link_complement(families):
    assert _representation_dimension(families,0)==6+48*families**2
    assert _representation_dimension(families,1)==16*families+64*families**3-20*families**2


def test_every_retained_basis_vector_has_correct_charge_norm_and_orthogonality():
    packet=json.loads((HERE/'joined_basis.json').read_text())
    for q in (0,1):
        sector=packet[str(q)]
        assert len(sector['labels'])==len(sector['basis'])==sector['dimension']
        reverse={}
        for col,(basis,label) in enumerate(zip(sector['basis'],sector['labels'],strict=True)):
            norm=Fraction(0)
            for h,raw,l0,l1,c in basis['entries']:
                mask=int(raw);coefficient=Fraction(c)
                assert 0<=l0<60 and 0<=l1<60
                assert mask>=0 and mask.bit_length()<=180
                assert mask.bit_count()-2*(h>=0)==q
                norm+=coefficient**2
                reverse.setdefault((h,mask,l0,l1),[]).append((col,coefficient))
            assert norm==Fraction(sector['raw_norms'][col])>0
            assert label['N_F']-2*label['N_H']==q
        cross={}
        for row in reverse.values():
            for i,a in row:
                for j,b in row:
                    if i<j:cross[i,j]=cross.get((i,j),Fraction(0))+a*b
        assert all(x==0 for x in cross.values())


def test_complete_addition_probes_keep_pauli_blocking_and_reference():
    matrices,_=load_operators(HERE)
    packet=json.loads((HERE/'joined_basis.json').read_text())['0']
    for f in range(3):
        for spin in (0,1):
            b=matrices[f'addition_{f}_{spin}']
            m=((9*3+f)*2+spin)
            occupation=[]
            for row,norm in zip(packet['basis'],packet['raw_norms'],strict=True):
                value=sum(Fraction(c)**2 for h,mask,l0,l1,c in row['entries'] if int(mask)&(1<<m))
                occupation.append(float(value/Fraction(norm)))
            np.testing.assert_allclose((b.getH()@b).toarray(),np.diag(1-np.asarray(occupation)),rtol=0,atol=2e-15)


def test_primary_interaction_and_calibration_on_complete_carrier():
    y,sigma=primary_family_matrices()
    make=lambda **kw:build_join(HERE,family0=y,family1=sigma,yukawa=y,**kw)
    model=make(phase0=.17,phase1=-.23,higgs_phase=.11)
    for q in (0,1):
        checks=operator_checks(model['sectors'][q])
        assert checks['Higgs_hopping_commutator_frobenius']>1
    base=make();shift=make(fermion_energy=.375,scalar_energy=.25,reference_energy=-.375)
    for q in (0,1):
        assert maxabs(base['sectors'][q]['H']-shift['sectors'][q]['H'])<2e-15
    check=schur_checks(base['sectors'][1],{},'test')
    assert check['coupling_gram_min_eigenvalue']>1.8
    assert check['offdiagonal_is_literal_Higgs_max']==0
    assert abs(check['resolvent_checks'][0]['self_energy_operator_norm']-1.793879153064)<2e-12
    damaged=dict(base['sectors'][1])
    damaged['H_Higgs']=0*damaged['H_Higgs']
    with pytest.raises(ArithmeticError,match='Higgs interaction'):
        schur_checks(damaged,{},'damaged')


def test_near_degenerate_preparation_retains_each_initial_energy():
    # Both initial energies are below the numerical grouping tolerance, but
    # at long times they are not an exactly degenerate preparation.
    eps=1e-15;e0=np.array([0.,eps]);e1=np.array([3*eps,8*eps])
    dummy=lambda e:{'H':sparse.diags(e,format='csr',dtype=complex),
        'observables':{'twice_Sz':sparse.csr_matrix((2,2),dtype=complex)}}
    initial=block_eigh(dummy(e0));final=block_eigh(dummy(e1))
    assert initial[3]['ground_numerical_multiplicity']==2
    model={'addition':{(0,0):sparse.eye(2,format='csr')}};arrays={}
    receipt=addition_data(model,{0:initial,1:final},'test',arrays)
    assert receipt['all_initial_final_pairs_retained']==4
    t=1/eps
    actual=np.einsum('n,nab->ab',np.exp(-1j*t*arrays['test_addition_frequencies']),arrays['test_addition_residues'])[0,0]
    expected=(np.exp(-3j)+np.exp(-7j))/2
    assert abs(actual-expected)<2e-15
    incorrect=(np.exp(-3j)+np.exp(-8j))/2
    assert abs(actual-incorrect)>.4


def test_rounding_diagnostic_includes_floating_entries_at_exact_zeros():
    with mp.workdps(50):
        exact={(0,0):mp.mpf(1),(1,1):mp.mpf(2)}
        floating=np.array([[1.,3e-8],[4e-8,2.]],complex)
        assert abs(float(matrix_rounding_error(exact,[0,1],floating))-5e-8)<1e-23


def test_precision_audit_restores_context_when_exact_input_is_missing(tmp_path):
    from .precision_audit import run
    y,sigma=primary_family_matrices()
    np.savez(tmp_path/'joined_response_arrays.npz',source_Y=y,source_sigma=sigma)
    with mp.workdps(37):
        with pytest.raises(FileNotFoundError):
            run(HERE,tmp_path/'absent_exact_references.json',response_root=tmp_path)
        assert mp.mp.dps==37


def test_eigensolver_rejects_nonhermitian_input():
    model={'H':sparse.csr_matrix([[1.,1.],[0.,2.]]),
           'observables':{'twice_Sz':sparse.csr_matrix((2,2))}}
    with pytest.raises(ValueError,match='Hermitian'):
        block_eigh(model)


def test_basis_corruption_is_rejected_before_readout(tmp_path):
    for name in ('joined_operator_receipt.json','joined_operators.npz','joined_basis.json','joined_exact_operators.json'):
        shutil.copy2(HERE/name,tmp_path/name)
    with (tmp_path/'joined_basis.json').open('a') as file:file.write(' ')
    with pytest.raises(ArithmeticError,match='checksum'):
        load_operators(tmp_path)
