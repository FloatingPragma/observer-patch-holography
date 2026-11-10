import sympy as s
from . import model as m

def test_independent_Taylor_and_exact_rational_inertia():
    row=m.audit()
    assert row['all_13_compression_intervals_independently_certified']
    assert len(row['exact_rational_inertia_endpoint_counts'])==13
    for j,(lo,hi) in enumerate(row['exact_rational_inertia_endpoint_counts'],1):assert lo<j<=hi
    assert m.F(row['actual_ground_nonvacuum_probability_positive_lower'])>m.F('1e-12')

def test_two_photon_invariant_contraction_with_arbitrary_real_direction():
    f=(s.Rational(3,5),s.Rational(4,5),s.Integer(0))
    assert sum(x*x for x in f)==1
    # In the normalized two-particle basis, double occupations have f_j^2
    # and distinct-mode occupations have sqrt2 f_i f_j.
    diagonal=[x*x for x in f]
    mixed=[s.sqrt(2)*f[i]*f[j] for i in range(3) for j in range(i+1,3)]
    assert sum(x*x for x in diagonal+mixed)==1
    assert sum(diagonal)/s.sqrt(3)==1/s.sqrt(3)

def test_weighted_adjacency_quotient_controls_complete_graph():
    wr,wb,r=s.symbols('wr wb r',positive=True)
    edges=m.candidate.mobile.actual_edges();A=s.zeros(13)
    for e,(u,v) in enumerate(edges):A[u,v]=A[v,u]=wr if e<12 else wb
    vector=s.Matrix([12*wr/r]+[s.Integer(1)]*12)
    residual=A*vector-r*vector
    assert residual[0]==0
    assert all(s.expand(r*value)==12*wr**2+5*wb*r-r**2 for value in residual[1:])


def test_public_replay_and_scientific_comparison_fail_closed():
    import copy,pytest
    from .certificate import certificate,compare_research
    value=m.audit()
    assert compare_research(value)['all_scientific_fields_equal']
    wrong=copy.deepcopy(value)
    wrong['all_13_compression_intervals_independently_certified']=False
    with pytest.raises(ValueError):compare_research(wrong)
    public=certificate()
    assert public['independent_research_comparison']['fresh_public_calculation_replayed']
