from fractions import Fraction as F
import pytest
from . import source as s

def test_odd_and_even_charge_domains():
    for q in(-3,-1,1,3):assert s.equal_integer_halves(q)is None
    assert s.equal_integer_halves(2)==(1,1)
    with pytest.raises(ValueError):s.equal_integer_halves(F(1,2))

def test_joint_cut_records_keep_all_weights():
    w=s.joint.step(s.joint.prepare(),electric_phase=True)
    result=s.cut_records(w,(0,))
    rows=result['complete_joint_charge_flux_records']
    assert sum(F(r['probability'])for r in rows)==1
    assert all(r['charge_in_q_star']==r['outward_flux_in_q_star']for r in rows)
    assert not result['postselection_or_renormalization']

def test_undressed_current_fails_gauss_record():
    w=s.joint.prepare();w.dressed=False
    w=s.joint.step(w,electric_phase=True)
    with pytest.raises(ArithmeticError):s.cut_records(w,(0,))

def test_duplicate_cut_rejected():
    with pytest.raises(ValueError):s.cut_records(s.joint.prepare(),(0,0))
