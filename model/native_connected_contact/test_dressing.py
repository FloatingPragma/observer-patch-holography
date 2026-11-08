from fractions import Fraction as F
from itertools import combinations

import pytest
import sympy as sp

from . import dressing as d


def control():
    edges = ((0,1),(1,2),(0,2),(0,3),(1,4),(2,5),(3,4),(4,5),(3,5))
    neighbours = tuple(tuple(sorted(b if a == i else a for a, b in edges if i in (a,b))) for i in range(6))
    return d.Problem(neighbours, (0,1,2))


def dense(problem):
    pairs = tuple(combinations(range(len(problem.neighbours)), 2)); where = {p:i for i,p in enumerate(pairs)}
    result = sp.zeros(len(pairs))
    # Independent edge-transposition construction, not apply().
    edges = [(a,b) for a,row in enumerate(problem.neighbours) for b in row if a < b]
    for j,pair in enumerate(pairs):
        for a,b in edges:
            image = tuple(sorted(b if x == a else a if x == b else x for x in pair))
            result[j,j] += 1; result[where[image],j] -= 1
    return pairs,result


def psd(matrix):
    assert matrix == matrix.T
    for n in range(1,matrix.rows+1):
        for indices in combinations(range(matrix.rows),n):
            assert matrix.extract(indices,indices).det() >= 0


@pytest.mark.parametrize('shift',[F(1,2),F(3),F(24)])
def test_exact_full_resolvent_and_monotone_exterior_error(shift):
    p=control();pairs,g=dense(p);ids=[i for i,pair in enumerate(pairs) if set(pair)<=set(p.core)]
    exterior=[i for i in range(len(pairs)) if i not in ids]
    a=g.extract(ids,ids);b=g.extract(exterior,ids);dd=g.extract(exterior,exterior)
    sigma=b.T*(dd+sp.Rational(shift)*sp.eye(len(exterior))).inv()*b
    exact=(g+sp.Rational(shift)*sp.eye(len(pairs))).inv().extract(ids,ids)
    assert exact==(a+sp.Rational(shift)*sp.eye(len(ids))-sigma).inv()
    previous=sp.zeros(len(ids))
    for degree in range(4):
        result=d.approximation(p,shift,degree);approx=result['schur'].inv()
        assert result['moments'][degree]==b.T*(24*sp.eye(len(exterior))-dd)**degree*b
        psd(approx-previous);psd(exact-approx)
        psd(sp.Rational(result['core_resolvent_error_upper'])*sp.eye(len(ids))-(exact-approx))
        previous=approx


def test_actual_source_core_compression_and_no_local_primitive_reprojection():
    p=d.problem(0);data=d.moments(p,2);q=d.contact.data()['primitive']
    assert data['core_action']==d.contact.data()['two']+data['coupling_gram']
    assert data['coupling_gram'].rank()==30
    assert data['reachable_exterior_counts']==(33,204,425)
    assert q*data['core_action'] != data['core_action']*q
    assert q*data['moments'][0] != data['moments'][0]*q


def test_connected_dressing_is_distinct_from_isolated_response():
    isolated=-24*d.contact.response(-24)
    intervals=[d.response(d.problem(0),24,j)['normalized_Abel_interval'] for j in range(3)]
    assert intervals[0][0]<intervals[1][0]<intervals[2][0]
    assert intervals[0][1]>intervals[1][1]>intervals[2][1]
    assert isolated-intervals[-1][1]>F(1,120)
    assert intervals[-1][1]-intervals[-1][0]==F(1,2304)


def test_colliding_outward_images_are_rejected():
    edges=((0,1),(1,2),(0,2),(0,3),(1,3))
    rows=tuple(tuple(sorted(b if a==i else a for a,b in edges if i in (a,b))) for i in range(4))
    with pytest.raises(ValueError,match='distinct outward'):
        d.blocks(d.Problem(rows,(0,1,2)))
    edges=((0,1),(1,2),(0,2),(0,3),(0,4))
    rows=tuple(tuple(sorted(b if a==i else a for a,b in edges if i in (a,b))) for i in range(5))
    with pytest.raises(ValueError,match='one outward seam'):
        d.blocks(d.Problem(rows,(0,1,2)))


def test_mutated_outward_amplitude_fails_source_coupling(monkeypatch):
    original=d.apply
    def wrong(problem,values,exterior=False):
        result=original(problem,values,exterior)
        if not exterior:
            for pair in result:
                if not set(pair)<=set(problem.core): result[pair]*=2
        return result
    monkeypatch.setattr(d,'apply',wrong)
    with pytest.raises(ValueError,match='distinct outward'):
        d.blocks(control())


@pytest.mark.parametrize('shift',[0,-1,True,1.0,'24'])
def test_off_spectrum_shift_domain(shift):
    with pytest.raises(ValueError): d.errors(shift,2)


@pytest.mark.parametrize('degree',[True,-1,2.0,7])
def test_finite_exact_depth_domain(degree):
    with pytest.raises(ValueError): d.moments(control(),degree)


def test_topology_and_level_type_guards():
    with pytest.raises(ValueError):d.problem(True)
    with pytest.raises(ValueError):d.problem(3)
    with pytest.raises(ValueError):d.problem(0,True)
    with pytest.raises(ValueError):d.validate(d.Problem(((1,),(0,True)),(0,1)))
    with pytest.raises(ValueError):d.validate(d.Problem(((1,),()),(0,1)))
