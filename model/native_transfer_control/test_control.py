from copy import deepcopy
from fractions import Fraction as F
import json
import sys
import numpy as np
import pytest
from . import spectral as sp, control as c, packet


@pytest.fixture(scope='module')
def source():return c.response.data()


def test_exact_source_interval_and_ground(source):
    rows,pi=sp.source_intervals(source)
    for i,row in enumerate(rows):
        for j,lo,hi in row:
            wanted=float(source['P'][i,j])*np.sqrt(float(pi[j]/pi[i]))
            assert float(lo)-1e-15<=wanted<=float(hi)+1e-15
    P=np.array(source['P'],float);p=np.array(pi,float)
    S=P*np.sqrt(p)[None,:]/np.sqrt(p)[:,None]
    assert np.linalg.norm(S@np.sqrt(p)-np.sqrt(p))<1e-14
    assert np.linalg.norm(S-S.T)<1e-14
    assert np.linalg.eigvalsh(S)[0]>float(c.transfer.LOWER)


def test_residual_certificate_without_numerical_eigenvalue_assumption(source):
    assert sp.certify(data=source)['spectral_eigenvalue_gaps_required'] is False
    columns,values=sp.load();bad=[list(x) for x in columns];bad[0][0]+=F(1,1000)
    with pytest.raises(ArithmeticError):sp.certify(tuple(map(tuple,bad)),values,source)
    badvalues=list(values);badvalues[0]+=F(1,100)
    with pytest.raises(ArithmeticError):sp.certify(columns,tuple(badvalues),source)


def test_nonreversible_source_mutation_rejected(source):
    bad=dict(source);P=source['P'].as_mutable();P[0,0]+=F(1,1000);bad['P']=P
    with pytest.raises(ArithmeticError):sp.source_intervals(bad)
    bad=dict(source);pi=source['pi'].as_mutable();pi[0]+=F(1,1000);bad['pi']=pi
    with pytest.raises(ArithmeticError):sp.source_intervals(bad)


@pytest.mark.parametrize('value',[True,0,.1,F(2),F(1,2**65)])
def test_transfer_parameter_rejects_alias_and_unbounded_inputs(value):
    with pytest.raises(ValueError):sp.time_parameter(value)


@pytest.mark.parametrize('time',[F(-1),F(),F(1,10),F(1)])
def test_log_grid_bounds_cover_all_declared_endpoints(time):
    _,values=sp.load()
    for v in values:
        lo,hi=sp.log_interval(v);n=sp.phase_grid(v,time)
        a,b=sorted((-time*lo,-time*hi))
        assert max(abs(F(n,c.iv.AQ)-a),abs(F(n,c.iv.AQ)-b))<=F(1,c.iv.AQ)
    assert c.protocol_budget(time)['complete_halfdiamond_strict_upper']=='1/11000'


def test_proposal_schema_and_scalar_types_reject_mutations(tmp_path):
    original=json.loads(sp.PROPOSAL.read_text())
    for key,value in [('denominator',True),('schema','untyped'),('columns',original['columns'][:-1])]:
        data=deepcopy(original);data[key]=value;p=tmp_path/'bad.json';p.write_text(json.dumps(data))
        with pytest.raises(ValueError):sp.load(p)
    data=deepcopy(original);data['columns'][0][0]=True;p.write_text(json.dumps(data))
    with pytest.raises(ValueError):sp.load(p)


def test_rational_polar_comparison_against_independent_numeric_diagnostic(source):
    cols,values=sp.load();V=np.array(cols,float).T;lam=np.array(values,float)
    eig,U=np.linalg.eigh(V.T@V);W=V@(U*(1/np.sqrt(eig)))@U.T
    pi=np.array(source['pi'],float).ravel();S=np.array(source['P'],float)*np.sqrt(pi)[None,:]/np.sqrt(pi)[:,None]
    assert np.linalg.norm(W-V,2)<float(sp.GRAM_GATE)
    assert np.linalg.norm(S-(W*lam)@W.T,2)<float(sp.RESIDUAL_GATE+2*sp.GRAM_GATE)


def test_ground_replacement_preserves_arbitrary_reference_and_keeps_archive():
    rng=np.random.default_rng(2027);x=rng.normal(size=(66,2))+1j*rng.normal(size=(66,2));x/=np.linalg.norm(x)
    run=c.detector.run(7,(0,)*11);out,records=c.ideal_ground_branch(run['slots'],x)
    weight,mapping,_=c.detector.branch(run['slots']);i=next(iter(mapping))
    assert np.allclose(out.conj().T@out,float(weight)*np.outer(x[i].conj(),x[i]))
    assert records==run['records']
    zero=list(run['slots']);zero[0]=next(k for k,r in enumerate(c.detector.kernel(c.detector.tree()[0])) if not r.pairs)
    assert not c.ideal_ground_branch(tuple(zero),x)[0].any()


def test_native_readout_scores_initial_not_disturbed_terminal():
    found=False
    for i in range(66):
        for bits in (0,2047):
            run=c.detector.run(i,tuple(bits>>k&1 for k in range(11)))
            assert c.readout_score(run['slots'])==int(i==c.TARGET)
            found|=(i==c.TARGET)!=(run['final']==c.TARGET)
    assert found


def test_ground_word_and_phase_use_actual_compiler():
    w=c.ground_word(65);assert len(w)<=c.GROUND_PULSES;c.syn.certify_word(w)
    x=np.zeros((66,1),complex);x[65]=1;y,abort=c.menus.successful_word(w,x)
    _,pi=sp.source_intervals();omega=np.sqrt(np.array(pi,float))
    assert np.linalg.norm(y[:,0]-omega)<1e-8
    assert np.trace(abort).real<1e-8
    w=c.intervention_word();c.syn.certify_word(w);z,_=c.menus.successful_word(w,y)
    wanted=y.copy();wanted[65]*=1j
    assert np.linalg.norm(z-wanted)<1e-8
    stopped=c.syn.interpreter.execute_word(c.syn.interpreter.source_word(w,9),x)
    assert stopped['status']=='RESOURCE_STOP' and stopped['pulses_executed']==0


def test_complete_failure_and_padding_addresses():
    vector=tuple((F(i==0),F()) for i in range(66));word=c.menus.phase_word(vector,F(),1)
    x=np.eye(66,2,dtype=complex)/np.sqrt(2)
    event=word[0];dimension=1<<len(event['targets'])
    for complement in range(dimension,6):
        output=c.complete_word_failure((),word,0,complement,x)
        assert output.shape==(512,2) and np.isfinite(output).all()
    address=next(i for i in range(512) if i not in c.syn.tomo.codes())
    row=c.complete_padding_row((),word,address,x)
    assert row.shape==(2,)
    with pytest.raises(ValueError):c.complete_padding_row((),word,True,x)
    with pytest.raises(ValueError):c.complete_padding_row((),word,c.syn.tomo.codes()[0],x)


def test_witness_includes_unequal_trace_error():
    w=c.witness()
    assert F(w['ideal_intervention_probability_increment_lower'])>F(1,1000)
    assert F(w['stored_map_score_increment_lower'])==F(w['ideal_intervention_probability_increment_lower'])-F(4,11000)
    assert F(w['stored_map_score_increment_lower'])>F(1,2000)
    assert c.outcomes()['preparation_tree']['formal_raw_words']==25**11
    assert c.protocol_budget()['native_pulses_upper']==3877451


def test_numeric_stored_primitive_mutation_is_rejected(monkeypatch):
    original=c.syn.compiler.native_pulse
    def bad(kind,angle):
        out=original(kind,angle).copy();out[0,0]+=1e-3;return out
    monkeypatch.setattr(c.syn.compiler,'native_pulse',bad);c.syn.certify_primitive.cache_clear()
    try:
        with pytest.raises(ArithmeticError):c.syn.certify_word(c.intervention_word())
    finally:c.syn.certify_primitive.cache_clear()


def test_canonical_receipt_rejects_bool_alias_and_duplicate_keys(tmp_path,monkeypatch):
    original={'finite':1,'scope':False};monkeypatch.setattr(packet,'produce',lambda:original)
    p=tmp_path/'receipt.json';p.write_text(packet.canonical({'finite':True,'scope':0}))
    with pytest.raises(ValueError):packet.verify(p)
    p.write_text('{"finite":1,"finite":1,"scope":false}')
    with pytest.raises(ValueError):packet.verify(p)


def test_all_loaded_local_imports_are_bound():
    import subprocess
    code = """
import sys
from pathlib import Path
from model.native_transfer_control import packet
bound=set(packet.source_bindings());missing=[]
for module in tuple(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if not path:continue
    try:relative=str(Path(path).resolve().relative_to(packet.ROOT))
    except ValueError:continue
    if relative.startswith(('model/','tools/')) and relative.endswith('.py') and relative not in bound:
        missing.append(relative)
assert not missing,missing
"""
    subprocess.run([sys.executable,'-c',code],cwd=packet.ROOT,check=True)


def test_frozen_receipt_source_pins_match():
    value=packet.strict_json.load(packet.OUTPUT)
    assert value['source_bindings']==packet.source_bindings()
    assert len(value['actual_compilation']['transfer_phases'])==66
    assert value['ground_response']['complete_increment_strict_lower']=='1/2000'
