import numpy as np
from .coupling import tree_loop,tree_coefficients,creutz,jackknife,loops,gauge_control,RECTANGLES,K


def test_tree_plaquette_has_exact_finite_normalization():
    for side in (3,4,6,8):
        assert abs(tree_loop(side,1,1)-(1-side**-4)/3)<2e-15


def test_tree_covariance_against_independent_position_inverse():
    side=3;coords=np.indices((side,)*4).reshape(4,-1).T;volume=len(coords)
    index={tuple(x):i for i,x in enumerate(coords)};lap=np.zeros((volume,volume))
    for i,x in enumerate(coords):
        for mu in range(4):
            y=x.copy();y[mu]=(y[mu]+1)%side;j=index[tuple(y)]
            lap[i,i]+=2;lap[i,j]-=1;lap[j,i]-=1
    green=np.linalg.inv(lap+np.ones_like(lap)/volume)-np.ones_like(lap)/volume
    for r,t in [(1,1),(2,1),(2,2)]:
        current=np.zeros((4,volume));x=np.zeros(4,dtype=int)
        for mu,length,sign in [(0,r,1),(1,t,1),(0,r,-1),(1,t,-1)]:
            for _ in range(length):
                if sign<0:x[mu]=(x[mu]-1)%side
                current[mu,index[tuple(x)]]+=sign
                if sign>0:x[mu]=(x[mu]+1)%side
        assert not x.any()
        independent=2/3*sum(j@green@j for j in current)
        assert abs(independent-tree_loop(side,r,t))<3e-14


def test_normalization_removes_free_geometry_scale_dependence():
    for side in (6,8):
        w=np.exp(-1.2*np.array([tree_loop(side,r,t) for r,t in RECTANGLES]))
        assert np.allclose(creutz(w)/tree_coefficients(side),1.2,atol=3e-14)
        assert np.ptp(creutz(w))>.1  # Geometry alone would counterfeit running without k_R.
        result=jackknife(np.tile(w,(80,1)),side,10)
        assert max(result['standard_errors'])<1e-13


def test_actual_closed_loops_are_gauge_invariant_and_match_plaquette():
    from model.colour_deconfinement.su3 import random_su3
    side=4;u=random_su3(side**4*4,np.random.default_rng(81)).reshape(side**4,4,3,3)
    up,_=K.neighbours(*(side,)*4)
    assert abs(loops(u,up)[0]-sum(K.plaquette(u,up))/2)<3e-15
    assert gauge_control(u,side)<3e-15


def test_nonpositive_means_cannot_hide_in_squared_denominators():
    import pytest
    for row in [[.6,-.2,.1,-.05,.03],[.6,.2,.1,0,.03],
                [.6,.2,.1,float('nan'),.03],[.6,.2,.1,float('inf'),.03],
                [.6,.2,.1,.03],[[.6,.2,.1,.05,.03]]]:
        with pytest.raises(ValueError):creutz(row)
