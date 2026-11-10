"""Reject altered proposal custody and replay without a numerical eigensolver."""
import json

import numpy as np
import pytest
import scipy.linalg as la
import sympy as s

from ..eigenbasis_witness import (BITS,SCALE,WitnessError,block_key,load_witness,
                                 write_witnesses)
from .spectrum import exact_dyadic_eigenbasis_error


@pytest.fixture
def proposal(tmp_path):
    path=tmp_path/'witnesses.npz';key=block_key((0,0),2,0)
    k=s.ImmutableSparseMatrix(s.diag(1,2));g=s.ImmutableSparseMatrix(s.eye(2))
    d=np.array([SCALE,2*SCALE],dtype=np.int64)
    z=np.eye(2,dtype=np.int64)*SCALE
    write_witnesses(path,'test',( (key,k,g,d,z), ))
    return path,key,k,g,d,z


def test_dyadic_archive_round_trip_and_exact_full_basis(proposal):
    path,key,k,g,d,z=proposal
    actual_d,actual_z=load_witness(path,'test',key,k,g)
    assert np.array_equal(d,actual_d) and np.array_equal(z,actual_z)
    assert not actual_d.flags.writeable and not actual_z.flags.writeable
    assert exact_dyadic_eigenbasis_error(k,g,actual_d,actual_z)==(0,0)
    assert np.array_equal(np.asarray(actual_z,dtype=float)*2**-BITS*SCALE,z)


def test_wrong_domain_key_matrix_metric_and_dimension_are_rejected(proposal):
    path,key,k,g,d,z=proposal
    with pytest.raises(WitnessError,match='domain mismatch'):
        load_witness(path,'other',key,k,g)
    for wrong_key in (block_key((0,1),2,0),block_key((0,0),2,1),
                      block_key((0,0),2,0,1),block_key((0,0),2,0,component=(0,1))):
        with pytest.raises(WitnessError,match='missing eigenbasis witness block'):
            load_witness(path,'test',wrong_key,k,g)
    for wrong_k,wrong_g in ((k+s.eye(2),g),(k,2*g),(s.eye(3),s.eye(3))):
        with pytest.raises(WitnessError,match='exact K/G mismatch'):
            load_witness(path,'test',key,wrong_k,wrong_g)


def test_changed_integer_payload_is_rejected_even_when_shapes_match(proposal):
    path,key,k,g,d,z=proposal
    with np.load(path,allow_pickle=False) as archive:
        arrays={name:archive[name] for name in archive.files}
    name=next(name for name in arrays if name.endswith('_z'))
    arrays[name][0,0]+=1
    np.savez_compressed(path,**arrays)
    with pytest.raises(WitnessError,match='integer digest mismatch'):
        load_witness(path,'test',key,k,g)


def test_corrupt_manifest_missing_archive_and_missing_arrays_are_rejected(proposal,tmp_path):
    path,key,k,g,d,z=proposal
    with pytest.raises(WitnessError,match='missing eigenbasis witnesses'):
        load_witness(tmp_path/'absent.npz','test',key,k,g)
    with np.load(path,allow_pickle=False) as archive:
        arrays={name:archive[name] for name in archive.files}
    manifest=json.loads(arrays['manifest'].tobytes())
    manifest['bits']=43
    corrupt=dict(arrays,manifest=np.frombuffer(json.dumps(manifest).encode(),dtype=np.uint8))
    np.savez_compressed(path,**corrupt)
    with pytest.raises(WitnessError,match='precision mismatch'):
        load_witness(path,'test',key,k,g)
    arrays.pop(next(name for name in arrays if name.endswith('_d')))
    np.savez_compressed(path,**arrays)
    with pytest.raises(WitnessError,match='missing eigenbasis witness integer array'):
        load_witness(path,'test',key,k,g)


def test_hash_valid_bad_basis_still_fails_exact_rank_check(proposal):
    path,key,k,g,d,z=proposal
    z[:,1]=z[:,0]
    write_witnesses(path,'test',((key,k,g,d,z),))
    actual_d,actual_z=load_witness(path,'test',key,k,g)
    with pytest.raises(AssertionError,match='Gram bound does not prove full rank'):
        exact_dyadic_eigenbasis_error(k,g,actual_d,actual_z)


def test_hash_valid_wrong_energies_still_fail_exact_enclosure(proposal,monkeypatch):
    from . import spectrum
    path,key,k,g,d,z=proposal
    d[1]=4*SCALE
    write_witnesses(path,'test',((key,k,g,d,z),))
    actual_d,actual_z=load_witness(path,'test',key,k,g)
    error,gram=exact_dyadic_eigenbasis_error(k,g,actual_d,actual_z)
    assert gram==0 and error>s.Rational(1,1<<26)
    monkeypatch.setattr(spectrum,'WITNESS_PATH',path)
    monkeypatch.setattr(spectrum,'WITNESS_DOMAIN','test')
    monkeypatch.setattr(spectrum,'physical_block',lambda *args: ((0,1),g,g,k,k))
    spectrum.rational_block.cache_clear();spectrum.certified_block.cache_clear()
    try:
        with pytest.raises(AssertionError,match='exact eigenbasis error exceeds energy enclosure'):
            spectrum.certified_block(0,2,0)
    finally:
        spectrum.rational_block.cache_clear();spectrum.certified_block.cache_clear()


def test_source_spectral_replay_never_calls_eigh(monkeypatch):
    from . import spectrum as open_source
    from ..periodic_spectrum import spectrum as periodic_source
    from . import readout as open_readout
    from ..periodic_spectrum import readout as periodic_readout
    def forbidden(*args,**kwargs):
        pytest.fail('certificate replay called a numerical eigensolver')
    monkeypatch.setattr(la,'eigh',forbidden)
    for source,readout,label in ((open_source,open_readout,(0,0,0)),
                                (periodic_source,periodic_readout,(0,0,0,1))):
        source.rational_block.cache_clear();source.numerical_block.cache_clear()
        source.certified_block.cache_clear()
        values,vectors=source.numerical_block(*label)
        result=source.certified_block(*label)
        assert len(values)==vectors.shape[0]==len(result['energy_intervals'])>0
        readout.approximation.cache_clear()
        z,d,residuals,gram=readout.approximation(label)
        assert z.shape==vectors.shape and len(d)==len(residuals)==len(values)
        assert gram<s.Rational(1,1<<23) and all(r>=0 for r in residuals)


def test_receipt_failure_reports_changed_json_paths_without_writing(tmp_path):
    from ..receipt_replay import check_receipt
    path=tmp_path/'results.json'
    expected={'spectra':[[1,2]],'input_sha256':{'a/b.py':'old'}}
    actual={'spectra':[[1,3]],'input_sha256':{'a/b.py':'changed'}}
    stored=json.dumps(expected)+'\n';path.write_text(stored)
    check_receipt(path,stored,'test')
    with pytest.raises(SystemExit) as failure:
        check_receipt(path,json.dumps(actual)+'\n','test')
    report=json.loads(str(failure.value))
    assert report['difference_paths']==['/input_sha256/a~1b.py','/spectra/0/1']
    assert not report['serialization_only'] and path.read_text()==stored
    with pytest.raises(SystemExit) as formatting:
        check_receipt(path,json.dumps(expected,indent=2)+'\n','test')
    assert json.loads(str(formatting.value))['serialization_only']
