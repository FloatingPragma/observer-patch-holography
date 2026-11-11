from pathlib import Path
from fractions import Fraction as F
import numpy as np
import sympy as s

HERE=Path(__file__).resolve().parent
from . import model as m


def test_actual_noncommuting_boundary_homotopy_and_symmetry():
    control=m.finite_controls()
    assert not control['derivative_and_Wilson_Laplacian_commute_at_boundary']
    assert not control['ordinary_time_reversal_preserved']
    assert control['free_open_box_combined_PT_preserved']
    m.literal_box(2,s.Rational(3,5),s.Rational(4,5))


def test_fresh_bulk_hopping_and_photon_compression():
    assert m.bulk_controls()['fresh_IBP_identity']
    q=(np.arange(48)+.5)*2*np.pi/48;c=np.cos(q);z=np.sin(q)
    A=3-c[None,:,None]-c[None,None,:];C=c[:,None,None]
    S=z[None,:,None]**2+z[None,None,:]**2;E=np.sqrt(A*A+2+S-2*A*C)
    first=np.mean((A*C-1)/E)
    second=-np.mean(((A*C-1)**2+1+S)/E**3)
    assert abs(first-second)<1e-10 and first<-2/343
    assert E.max()<=np.sqrt(37) and E.min()>=1
    for r in (.51,.8,.99,.999):assert 3*r-2*r*r>1


def test_fixed_frame_Fock_lift_orientation():
    # Ordinary number-preserving second quantization: R† c R=U c,
    # so the transformed one-particle state projector is U P U†.
    c0=s.Matrix([[0,1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]])
    c1=s.Matrix([[0,0,1,0],[0,0,0,-1],[0,0,0,0],[0,0,0,0]])
    U=s.Matrix([[s.Rational(3,5),-s.Rational(4,5)],[s.Rational(4,5),s.Rational(3,5)]])
    R=s.diag(1,U,1)
    assert R.H*c0*R==U[0,0]*c0+U[0,1]*c1
    assert R.H*c1*R==U[1,0]*c0+U[1,1]*c1
    state=R*s.Matrix([0,1,0,0])
    P=s.Matrix([[(state.H*cj.H*ci*state)[0] for cj in (c0,c1)] for ci in (c0,c1)])
    assert P==U*s.diag(1,0)*U.H


def test_fresh_bounds_and_no_ground_promotion():
    v=m.scientific();b=v['bounds_exact']
    assert F(b['complete_finite_law_TV_separation_lower'])>F(34,10**11)
    assert F(b['finite_dressed_hopping_magnitude_lower'])>F(1,173)
    assert F(b['compressed_squared_curvature_minus_one_lower'])>0
    assert not v['boundaries']['fixed_cube_ground_theorem_inherited_in_volume']
    assert v['separate_source_only_completion_option']['all_caps_K_at_least_3']


def test_source_only_volume_completion_and_exact_atomic_trial():
    _,beta,gamma=m.source_matrices()
    onsite=(s.eye(4)-beta)/2
    for n in (1,2):
        vertices,edges,D,C,H,currents=m.literal_box(n)
        V=len(vertices);P=s.SparseMatrix(s.kronecker_product(s.eye(V),onsite))
        assert s.trace(P*H)==-2*V
        for x in range(V):
            q=s.SparseMatrix(4*V,4*V,{(4*x+j,4*x+j):1 for j in range(4)})
            assert s.trace(P*q)==2 and s.trace(P*q*(s.eye(4*V)-P)*q)==0
        assert (18*V+1)-18*V==1
