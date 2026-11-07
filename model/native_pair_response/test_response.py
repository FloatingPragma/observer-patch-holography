from fractions import Fraction as F
from itertools import combinations,product
import sympy as sp
import pytest
from . import response as r,records


def family(V):return r.s.source.Family((1,)*V,tuple((i,i+1) for i in range(V-1)),(1,)*(V-1))


@pytest.mark.parametrize('V',(4,5,6))
def test_independent_one_walker_counts_give_killed_compression_and_complete_collision(V):
    f=family(V);pairs=tuple(combinations(range(V),2));M=len(f.seams)
    for pair in pairs:
        free={};collision=0
        # Independent construction: edge and walker, with no native coin formula.
        for edge,walker in product(f.seams,(0,1)):
            after=list(pair);a,b=edge
            if after[walker]==a:after[walker]=b
            elif after[walker]==b:after[walker]=a
            if after[0]==after[1]:collision+=1
            else:
                target=tuple(sorted(after));free[target]=free.get(target,0)+1
        actual=r.killed_counts(f,{pair:1})
        assert actual==free and sum(actual.values())+collision==2*M
        assert collision==2*int(pair in f.seams)
        native=r.s.counts_step(f,{pair:1},2)
        assert {x:native.get(x,0)-actual.get(x,0) for x in pairs if native.get(x,0)!=actual.get(x,0)}==(
            {pair:2} if pair in f.seams else {})


@pytest.mark.parametrize('T',(0,1,2,3,6))
def test_deterministic_kernel_telescoping_entrywise_control(T):
    f=family(4);pairs=tuple(combinations(range(4),2));M=len(f.seams)
    P=sp.Matrix([[F(r.s.counts_step(f,{x:1},2).get(y,0),2*M) for x in pairs] for y in pairs])
    R=sp.Matrix([[F(r.killed_counts(f,{x:1}).get(y,0),2*M) for x in pairs] for y in pairs])
    A=sp.diag(*(int(x in f.seams) for x in pairs));D=P**T-R**T
    assert D==sum((P**(T-1-k)*A*R**k/M for k in range(T)),sp.zeros(len(pairs)))
    assert all(x>=0 for x in D)
    for col,pair in enumerate(pairs):
        law=r.collision_law(f,{pair:F(1)},T)
        assert sum(law['surviving_pairs'].values())+law['collision']==1
        assert sum(D[:,col])==law['collision']


def test_collision_is_not_success_conditioned_or_simultaneous_two_walk_clock():
    f=family(4);a=((0,1),(2,3))
    law=r.collision_law(f,{(1,2):F(1)},1)
    assert law['collision']==F(1,3) and sum(law['surviving_pairs'].values())==F(2,3)
    native=r.response_matrix(f,(a,),1)[1][0][0]
    independent=r.response_matrix(f,(a,),1,killed=True)[1][0][0]
    assert native-independent==F(1,12)
    # Separate P1 clocks can move both labels in one step; the comparison cannot.
    assert (1,2) not in r.killed_counts(f,{(0,3):1})
    one0=r.s.counts_step(f,{(0,):1},1);one3=r.s.counts_step(f,{(3,):1},1)
    assert one0[(1,)]*one3[(2,)]>0


@pytest.mark.parametrize('q',(2,3,4))
@pytest.mark.parametrize('T',(0,1,2))
def test_cross_matrix_matches_independent_complete_native_readout(q,T):
    f=family(6);probes=(((0,1),(2,3)),((1,2),(4,5)))
    matrix=r.response_matrix(f,probes,T)[T];v=4*r.s.norm_factor(6,2,q)
    for i,j in product(range(2),repeat=2):
        assert records.direct_cross(f,q,T,probes[j],probes[i])==v*matrix[i][j]
        error=r.comparison_error(f,probes[j],probes[i],T)
        killed=r.response_matrix(f,probes,T,killed=True)[T]
        assert abs(matrix[i][j]-killed[i][j])<=error['normalized_cross_error_upper']


@pytest.mark.parametrize('level',(0,1,2))
def test_actual_spatial_matrices_have_positive_gram_and_bounded_local_output(level):
    f,probes=r.federation_menu(level);rows=r.response_matrix(f,probes,3)
    for matrix in rows:
        S=sp.Matrix(matrix);assert S==S.T
        for size in range(1,4):
            for indices in combinations(range(3),size):assert S.extract(indices,indices).det()>=0
    assert rows[0][0][2]==0
    local=r.local_error(f,len(f.weights)//3,probes[0],probes[2],3,1)
    assert set(x for e in probes[0]+probes[2] for x in e)<=set(local['roots'])<=set(local['window'])
    assert local['normalized_response_error_upper']<=2*local['complete_local_halfdiamond_upper']/local['normalization_variance']


@pytest.mark.parametrize('T',(0,1,2,3,5,8))
def test_separated_initial_pairs_retain_exact_finite_collision_cone(T):
    f=family(6);probe=((0,1),(4,5));value=r.comparison_error(f,probe,probe,T)
    assert value['minimum_initial_pair_distance']==3
    assert value['exact_collision_probability']<=value['collision_upper']
    if T<3:assert value['exact_collision_probability']==0


@pytest.mark.parametrize('call',(
    lambda:r.probe(family(4),((0,1),(1,2))),
    lambda:r.probe(family(4),((0,1),(0,3))),
    lambda:r.probe(family(4),((False,1),(2,3))),
    lambda:r.response_matrix(family(4),(((0,1),(2,3)),),True),
    lambda:r.response_matrix(family(4),(((0,1),(2,3)),),13),
    lambda:r.response_matrix(family(4),(((0,1),(2,3)),),1,killed=1),
    lambda:r.collision_law(family(4),{(0,1):F(1,2)},1),
    lambda:r.collision_law(family(4),{(0,1):F(2),(2,3):F(-1)},1),
    lambda:r.collision_law(family(4),{(0,1):1},65),
    lambda:r.local_error(family(4),1,((0,1),(2,3)),((0,1),(2,3)),1,1),
))
def test_exact_source_domains_fail_closed(call):
    with pytest.raises(ValueError):call()
