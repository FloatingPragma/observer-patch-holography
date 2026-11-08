from fractions import Fraction as F
from copy import deepcopy
import pytest
from . import diagonal_control as d


def test_exact_full_source_diagonal_directions():
    coefficients,identity=d.exact_identity('conditional_phase')
    assert coefficients==(F(1),)*12
    assert identity==((F(1),F()),)*3+((F(),F()),)*3
    coefficients,identity=d.exact_identity('rail_phase')
    assert coefficients==(F(1),)*4+(F(),)*8
    assert tuple(b-a for a,b in zip(identity[1],identity[2]))==(F(1),F(1))


@pytest.mark.parametrize('kind',d.KINDS)
def test_complete_source_phase_and_exact_structure(kind):
    result=d.execute(d.declare(kind,F(1,2),(1,3),cell=10,level=0))
    assert result['exact_stored_number_preservation']
    assert result['literal_exponential_error_diagnostic_upper']<=d.beam.GATE
    assert result['full_complement_diagonal']==((F(1),F()),)*2
    assert len(result['complement_rows'])==2
    declaration=result['declaration']
    assert declaration['resources']['canonical_declaration_bytes']==len(d.beam.native.canonical(declaration).encode())


def test_conditional_phase_entangles_two_dual_rails_in_one_cell():
    # Two logical qubits on (0,1),(2,3), each logical1 at its first port.
    masks=tuple((1<<(0 if a else 1))|(1<<(2 if b else 3))for a,b in ((0,0),(0,1),(1,0),(1,1)))
    rho={(x,0,y,0):F(1,4)for x in masks for y in masks}
    value=d.declare('conditional_phase',F(1,2),(0,2))
    result=d.lift(value,(rho,{}),{'kind':'federation','level':0,'total':2})
    rr,ii=result['success'];assert ii
    # Pure output amplitudes have a nonzero two-qubit determinant. For the
    # ideal law its modulus is |sin(theta/2)|/2; actual stored det stays nonzero.
    entries=d.execute(value)['stored_diagonal']
    def multiply(x,y):return x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]
    aa=multiply(entries[0],entries[3]);bb=multiply(entries[1],entries[2])
    determinant=tuple((a-b)/4 for a,b in zip(aa,bb,strict=True))
    assert sum(x*x for x in determinant)>F(3,25)**2
    assert all(x.bit_count()==y.bit_count()==2 for x,_,y,_ in set(rr)|set(ii))
    assert len(result['complements'])==2
    assert all(not p for row in result['complements']for p in row['parts'])
    # Public native repair accepts both exact parts separately.
    program=d.beam.native.program(federation={'level':0,'total':2},biased=False,
        commands=[{'op':'repair','selector':{'kind':'fixed','seam':[0,1]},'law':'E','repeat':1}])
    assert d.beam.native.execute(program,rr).status=='COMPLETE'
    assert d.beam.native.execute(program,ii).status=='COMPLETE'


@pytest.mark.parametrize('field',('kind','schema','integrated_controls','resources'))
def test_typed_evidence_bound_diagonal_admission(field):
    value=deepcopy(d.declare('conditional_phase',F(1,2)))
    value[field]='bad'
    with pytest.raises((ValueError,TypeError)):d.execute(value)


def test_phase_lift_keeps_reference_indices_and_imaginary_part():
    x,y=1,2
    rho={(x,7,x,7):F(1,2),(x,7,y,9):F(1,2),(y,9,x,7):F(1,2),(y,9,y,9):F(1,2)}
    value=d.declare('rail_phase',F(1,4))
    result=d.lift(value,(rho,{}),{'kind':'federation','level':0,'total':1})
    rr,ii=result['success']
    assert set(rr)==set(rho)
    assert set(ii)=={(x,7,y,9),(y,9,x,7)}
    assert ii[x,7,y,9]==-ii[y,9,x,7]
    with pytest.raises(ValueError):d.lift(value,(rho,{}),{'kind':'federation','level':1,'total':1})


def test_any_stored_offdiagonal_phase_is_rejected(monkeypatch):
    original=d.expm
    def mutation(value):
        matrix=original(value);matrix[0,1]=2**-100
        return matrix
    monkeypatch.setattr(d,'expm',mutation)
    with pytest.raises(ArithmeticError):d.execute(d.declare('rail_phase',F(1,2)))


def test_conditional_phase_uses_the_existing_literal_J_request():
    declaration=d.declare('conditional_phase',F(1,2))
    actual=d.beam.source.literal_controls('J',.5)
    assert tuple(map(F,declaration['integrated_controls']))==actual==(F(1,2),)*12
