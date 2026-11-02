from fractions import Fraction as F
import importlib.util
import io
from pathlib import Path
import types
import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('ordered_readout',Path(__file__).with_name('read_ordered_spin_energies.py'))
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
ROOT=Path.cwd()

@pytest.fixture(scope='module')
def data():
    source=types.ModuleType('source')
    inputs=r.checked_inputs(ROOT)
    exec(compile(inputs['model/native_ordered_ground/source.py'],'source.py','exec'),source.__dict__)
    h,scale=source.matrix()
    with np.load(io.BytesIO(inputs['model/native_ordered_ground/witness.npz']),allow_pickle=False) as archive:
        witness={k:archive[k] for k in archive.files}
    return witness,h,scale

@pytest.fixture(scope='module')
def proof(data):return r.verify_witness(*data)


def test_independent_certificate_matches_frozen_source(proof):
    assert r.canonical(proof[0])==r.canonical(r.load_json(ROOT/'model/native_ordered_ground/results.json')['finite_spectrum'])


def test_all_256_levels_and_isolated_clusters_retained(proof):
    _,energy,radius,ground=proof
    rows,clusters=r.extract(energy,radius,ground)
    assert len(rows)==256 and len(clusters)==68
    assert sum(c['certified_state_count'] for c in clusters)==256
    assert rows[1]['exact_ratio_interval']==['1','1']
    assert all(F(x['exact_gap_interval'][0])>0 for x in rows[1:])
    assert all(F(x['exact_ratio_interval'][0])>=1 for x in rows[1:])
    assert all(x['exact_degeneracy_claimed'] is False for x in clusters)
    for index,row in enumerate(rows[1:],1):
        lo,hi=map(F,row['exact_ratio_interval'])
        center=(energy[index]-energy[0])/(energy[1]-energy[0])
        assert lo<=center<=hi and hi-lo<F(1,10**8)


def test_outward_rounding_on_both_sides_of_zero():
    assert r.decimal_interval(-F(1,3),F(2,3),3)==['-0.334','0.667']
    with pytest.raises(ValueError):r.decimal_interval(F(2),F(1))


def test_literal_corruption_fails_before_acceptance(data):
    witness,h,scale=data
    broken={k:v.copy() for k,v in witness.items()};broken['vectors'][0,0]+=1<<48
    with pytest.raises(ValueError,match='Gram'):r.verify_witness(broken,h,scale)


def test_reordered_energies_rejected(data):
    witness,h,scale=data;broken={k:v.copy() for k,v in witness.items()}
    broken['energies'][0],broken['energies'][-1]=broken['energies'][-1],broken['energies'][0]
    with pytest.raises(ValueError,match='ordered'):r.verify_witness(broken,h,scale)


def test_nonliteral_float_witness_rejected(data):
    witness,h,scale=data;broken=dict(witness);broken['vectors']=broken['vectors'].astype(float)
    with pytest.raises(ValueError,match='int64'):r.verify_witness(broken,h,scale)


def test_lost_first_gap_rejected(proof):
    _,energy,radius,ground=proof
    with pytest.raises(ValueError,match='enclosures'):r.extract(energy,F(100),ground)


def test_changed_current_source_rejected(tmp_path):
    p=tmp_path/'model/native_ordered_ground/source.py';p.parent.mkdir(parents=True);p.write_text('changed')
    with pytest.raises(ValueError,match='source pin changed'):r.checked_inputs(tmp_path)


def test_json_scope_alias_and_duplicate_keys_rejected(tmp_path):
    assert r.canonical({'relativistic_masses_identified':False})!=r.canonical({'relativistic_masses_identified':0})
    p=tmp_path/'bad.json';p.write_text('{"a":1,"a":2}')
    with pytest.raises(ValueError,match='duplicate'):r.load_json(p)


def test_shared_ground_corner_extrema_and_interior_ratios():
    lo,hi=r.shared_ground_ratio((3,4),(1,2),(-2,-1),True)
    assert (lo,hi)==(F(5,4),F(5,2))
    for e in [F(3),F(7,2),F(4)]:
        for ref in [F(1),F(3,2),F(2)]:
            for ground in [F(-2),F(-3,2),F(-1)]:
                assert lo<=(e-ground)/(ref-ground)<=hi
    assert r.shared_ground_ratio((103,104),(101,102),(98,99),True)==(lo,hi)
    with pytest.raises(ValueError,match='separated'):
        r.shared_ground_ratio((3,4),(1,2),(1,2))


def test_declared_hamiltonian_norm_is_independently_checked(data):
    witness,h,scale=data
    broken=h.copy();broken[0,0]=1000
    with pytest.raises(ValueError,match='norm bound'):
        r.verify_witness(witness,broken,scale)


def test_tool_change_during_readout_rejected(monkeypatch):
    original=Path.read_bytes;own=Path(r.__file__).resolve();calls=0
    def changed(path):
        nonlocal calls
        raw=original(path)
        if path.resolve()==own:
            calls+=1
            if calls>1:return raw+b'changed during execution'
        return raw
    monkeypatch.setattr(Path,'read_bytes',changed)
    with pytest.raises(ValueError,match='tool changed'):
        r.report(ROOT)
