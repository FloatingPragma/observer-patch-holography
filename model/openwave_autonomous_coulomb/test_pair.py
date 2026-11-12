import numpy as np
import sympy as s
from .source import basis,matrices,column,REFERENCE_MASK,ingredients
from .experiment import current_matrices
def test_all_invariant_rows_and_loops():
    h,v=matrices();assert h.shape==(18,18)and h==h.H
    assert len({f for _,_,_,f in basis()})==9
def test_exact_native_continuity():assert len(current_matrices())==4
def test_reference_is_a_full_operator_argument():
    key=basis()[0][0];a=column(key);b=column(key,reference_mask=0)
    assert a[key]!=b[key]
def test_position_pinning_changes_source():
    h,_=matrices();P=s.diag(*[int(y==1)for _,_,y,_ in basis()])
    assert h*P!=P*h
def test_native_coordinate_has_no_discarded_output():
    U,labels,_,_=ingredients();outside=[i for i in range(119)if i not in labels]
    assert U.extract(outside,labels)==s.zeros(116,3)
    assert U.H.extract(outside,labels)==s.zeros(116,3)
