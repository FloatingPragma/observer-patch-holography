"""Corruptions cannot change source scope or bypass strict bridge verification."""
from copy import deepcopy
import json
import subprocess
import sys
import pytest
import sympy as sy
from model.ccf_lt2.bridges import verify
from model.ccf_lt2.bridges.gravity import coframe
from model.ccf_lt2.bridges.reference.check import multiplicities, gauss_rank, stationarity_boundary


def test_exact_Gauss_rank_agrees_with_literal_small_carrier():
    import itertools
    charges=(-1,1)
    counts=multiplicities(charges)
    B=sy.Matrix([[1],[-1]])
    rank,_=gauss_rank(B,1,counts)
    occupations=[sum(charges[i] for i in range(2) if mask>>i&1) for mask in range(4)]
    exact=sum(e==q and -e==r for q,r,e in itertools.product(occupations,occupations,range(-1,2)))
    assert rank==exact==6
    assert gauss_rank(B,1,{0:4})[0]!=exact


def test_stationarity_requires_actual_coherent_matrix_entry():
    assert stationarity_boundary({'local_capacity_entry_exact':'9/32'})['entropy_plus_area_first_variation_exact']=='9/16384'
    with pytest.raises(AssertionError):stationarity_boundary({'local_capacity_entry_exact':'0'})


@pytest.mark.parametrize('frame', [sy.zeros(3),sy.diag(-1,1,1),sy.diag(sy.I,-sy.I,1)])
def test_invalid_coframe_is_not_an_accepted_background(frame):
    with pytest.raises(ValueError):coframe.expression({}, {'edges':[],'normals':[],'static':{}},frame)


def test_optimized_python_cannot_skip_checks():
    run=subprocess.run([sys.executable,'-O','-m','model.ccf_lt2.bridges.verify'],cwd=verify.ROOT,
                       capture_output=True,text=True)
    assert run.returncode!=0 and 'requires enabled assertions' in run.stderr


def test_saved_bridge_receipt_binds_current_source():
    assert verify.validate_saved()['all_passed'] is True


@pytest.mark.parametrize('damage',['binding','candidate','physical','passed'])
def test_saved_bridge_rejects_mutation(tmp_path,damage):
    receipt=deepcopy(verify.validate_saved())
    if damage=='binding':receipt['scientific_inputs_sha256'].pop(next(iter(receipt['scientific_inputs_sha256'])))
    if damage=='candidate':receipt['candidate']='other-law'
    if damage=='physical':receipt['physical_rows_awarded']=[26]
    if damage=='passed':receipt['all_passed']=False
    path=tmp_path/'damaged.json';path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError):verify.validate_saved(path)


@pytest.mark.parametrize('damage',['omitted_theorem','axiom','source'])
def test_Lean_receipt_rejects_missing_or_unchecked_claims(tmp_path,monkeypatch,damage):
    from model.ccf_lt2.bridges import verify_lean
    data=deepcopy(verify_lean.validate_saved())
    if damage=='source':data['source_sha256'].pop(next(iter(data['source_sha256'])))
    else:
        module=next(iter(data['modules']));theorem=next(iter(data['modules'][module]))
        if damage=='omitted_theorem':del data['modules'][module][theorem]
        else:data['modules'][module][theorem].append('sorryAx')
    path=tmp_path/'lean.json';path.write_text(json.dumps(data));monkeypatch.setattr(verify_lean,'OUTPUT',path)
    with pytest.raises(ValueError):verify_lean.validate_saved()
