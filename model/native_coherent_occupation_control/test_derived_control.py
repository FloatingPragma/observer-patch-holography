from copy import deepcopy
from fractions import Fraction as F
import numpy as np
import pytest
from . import derived_control as d


def test_exact_four_native_generator_identity_and_full_complement():
    value=d.exact_identity()
    nonzero={(i,j):x for i,row in enumerate(value)for j,x in enumerate(row)if x}
    assert nonzero=={(1,2):1,(2,1):-1,(4,5):1,(5,4):-1}
    result=d.execute(d.declare(F(25,32)))
    assert all(x==0 for row in result['complement_rows']for x in row)
    assert result['full_complement_rows'][0][5]!=0
    assert result['full_complement_rows'][1][4]!=0
    assert result['full_complement_rows'][0][4]!=1
    assert abs(result['full_processor'][1][2]+result['full_processor'][2][1])<=2*d.GATE
    assert result['literal_exponential_error_diagnostic_upper']<=d.GATE


def test_exact_stored_return_and_complete_native_boundary_handoff():
    for level in (0,1):
        result=d.coupled_bell(level)
        assert result['number_leakage_mass']==0
        assert result['native_handoff_executed']
        assert result['native_raw_slots_retained']==25
        assert result['proper_code_complement_rows_retained']==4
        assert result['ideal_success_interval'][0]>F(4999,10000)
        assert result['stored_success_mass']>F(4999,10000)
        assert result['complete_halfdiamond_error_upper']<F(1,10**12)
        assert abs(result['stored_complete_output_mass']-1)<2*result['complete_halfdiamond_error_upper']
        assert not result['normalization_or_number_projection_inserted']


@pytest.mark.parametrize('angle',(F(4,5),F(1,2**1074),F(2),True,0.25))
def test_angle_domain_retains_exact_binary64_control_payloads(angle):
    with pytest.raises(ValueError):d.declare(angle)


def test_subnormal_half_angle_boundary_and_zero_services():
    admitted=d.declare(F(1,2**1073))
    assert F(admitted['integrated_controls'][6])==F(1,2**1074)
    zero=d.execute(d.declare(F()))
    assert zero['success']==tuple(tuple(F(i==j)for j in range(4))for i in range(4))
    assert zero['declaration']['resources']['native_response_pulses']==1


@pytest.mark.parametrize('field',('source_sha256','integrated_controls','resources','cell'))
def test_complete_declaration_cannot_omit_or_mutate_evidence(field):
    value=d.declare(F(25,32));changed=deepcopy(value)
    if field=='cell':changed[field]=False
    elif field=='integrated_controls':changed[field][4]='0'
    else:changed[field]={}
    with pytest.raises(ValueError):d.execute(changed)


def test_wrong_source_generator_identity_is_rejected(monkeypatch):
    original=d.source.specification
    spec=deepcopy(original())
    rows=list(spec[4]);row=list(rows[0]);cell=list(row[0]);part=list(cell[0]);part[0]+=1
    cell[0]=tuple(part);row[0]=tuple(cell);rows[0]=tuple(row);spec=list(spec);spec[4]=tuple(rows)
    monkeypatch.setattr(d.source,'specification',lambda:tuple(spec))
    with pytest.raises(ArithmeticError):d.exact_identity()


def test_structural_number_leak_is_refused_before_native_handoff(monkeypatch):
    actual=d.expm
    def changed(generator):
        matrix=actual(generator);matrix[0,1]=2**-100
        return matrix
    monkeypatch.setattr(d,'expm',changed)
    with pytest.raises(ArithmeticError):d.execute(d.declare(F(25,32)))


def test_existing_frozen_word_alphabet_is_unchanged():
    with pytest.raises(ValueError):
        d.native.source_word([{'kind':'occupation_exchange','targets':[0,1],'angle':.5}],12)
    assert d.compiler.KINDS==('X','Y','Z','I','H','J')


def test_literal_addressed_load_export_preserves_spectators_and_reference():
    declaration=d.declare(F(25,32),(1,5),cell=7,level=1)
    a,b=85,89;spectator=3
    x=(1<<a)|(1<<spectator);y=(1<<b)|(1<<spectator)
    rho={(x,8,x,8):F(1,2),(x,8,y,9):F(1,2),
         (y,9,x,8):F(1,2),(y,9,y,9):F(1,2)}
    carrier={'kind':'federation','level':1,'total':2}
    lifted=d.lift(declaration,rho,carrier)
    assert lifted['absolute_targets']==(a,b)
    assert len(lifted['complements'])==2
    assert all(not row['density']for row in lifted['complements'])
    assert {r for _,r,_,_ in lifted['success']}=={8,9}
    assert all((xx>>spectator)&1 and (yy>>spectator)&1
               for xx,_,yy,_ in lifted['success'])
    assert all(xx.bit_count()==yy.bit_count()==2 for xx,_,yy,_ in lifted['success'])
    assert any(r!=s for _,r,_,s in lifted['success'])
    with pytest.raises(ValueError):d.lift(declaration,rho,{**carrier,'level':0})
    with pytest.raises(ValueError):d.lift(declaration,rho,{**carrier,'total':1})


def test_declaration_charges_exact_payload_and_actual_geometry():
    declaration=d.declare(F(25,32),(1,11),cell=319,level=2)
    resources=declaration['resources']
    assert resources['angle_numerator_bits']==5
    assert resources['angle_denominator_bits']==6
    assert resources['cell_address_bits']==9
    assert resources['integrated_control_payload_bits']==768
    assert resources['canonical_declaration_bytes']==len(d.native.canonical(declaration).encode())
    for args in ({'cell':320,'level':2},{'cell':20,'level':0},{'level':3},{'cell':True}):
        with pytest.raises(ValueError):d.declare(F(25,32),**args)
    changed=deepcopy(declaration);changed['federation_level']=1
    with pytest.raises(ValueError):d.authenticate(changed)


def test_public_native_admission_and_fifo_are_really_executed(monkeypatch):
    original=d.native.execute;calls=[]
    def tracked(program,density,**kwargs):
        calls.append(program)
        return original(program,density,**kwargs)
    monkeypatch.setattr(d.native,'execute',tracked)
    result=d.coupled_bell(0)
    assert len(calls)==1 and result['public_native_status']=='COMPLETE'
    assert result['public_native_history_resources']=={
        'attempts_executed':1,'complete_expanded_blocks':25,'maximum_FIFO_occupied':1,
        'maximum_exported_records':0,'terminal_status_records_per_path':1}


def test_full_complement_rotation_cannot_be_erased(monkeypatch):
    actual=d.expm
    def erased(generator):
        matrix=actual(generator)
        matrix[4:6,4:6]=np.eye(2)
        return matrix
    monkeypatch.setattr(d,'expm',erased)
    with pytest.raises(ArithmeticError):d.execute(d.declare(F(25,32)))
