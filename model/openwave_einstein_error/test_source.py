import sympy as s
from itertools import product
from . import source as q

def test_every_symmetric_coordinate_and_metric_nullspace():
    for i,j in product(range(4),repeat=2):
        B=s.zeros(4);B[i,j]=B[j,i]=1
        assert q.reconstruct(B)==B-q.metric_coefficient(B)*q.ETA
    assert q.charges(123*q.ETA)==s.zeros(9,1)

def test_error_constants_are_sharp_and_all_vertices_are_bounded():
    _,_,R=q.design();C=2+s.sqrt(3)
    for signs in product((-1,1),repeat=9):
        error=q.tracefree_matrix(R*s.Matrix(signs))
        assert all(abs(v)<=C for v in error)
    result=q.run();assert result['sharp_maximum_entry_error_constant']=='2 + sqrt(3)'

def test_central_connectivity_counterexample():
    values={i:(-1)**i for i in range(8)}
    assert all(values[(i+1)%8]-values[(i-1)%8]==0 for i in values)
    assert any(values[(i+1)%8]!=values[i] for i in values)
