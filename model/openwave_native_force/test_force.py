"""Independent native-matrix, energy, ground and finite-window checks."""
from fractions import Fraction as F
from itertools import product
import mpmath as mp
import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import spsolve
import sympy as s

from model.native_equilibration import vacuum_repair as native
from . import force


def test_native_matrix_and_derivative_without_target_action_formula():
    r,matrix,_,charge=force.source_matrices()
    for value in (F(9,10),F(1),F(11,10)):
        actual=native.hamiltonian(2,2,value)/np.sqrt(14)
        np.testing.assert_allclose(np.asarray(matrix.subs(r,s.Rational(value)),float),actual,rtol=0,atol=2e-15)
    errors=[]
    for step in (F(1,100),F(1,200),F(1,400)):
        derivative=(native.hamiltonian(2,2,1+step)-native.hamiltonian(2,2,1-step))/(2*float(step)*np.sqrt(14))
        errors.append(np.linalg.norm(derivative-np.asarray(charge,float)))
    assert all(3.9<errors[i]/errors[i+1]<4.1 for i in (0,1))


def test_true_ground_charge_signs_with_independent_high_precision_eigenvectors():
    result=force.ground_counterexample()
    _,_,m,j=force.source_matrices()
    with mp.workdps(80):
        mm=mp.matrix([[mp.mpf(str(x.p))/int(x.q) for x in row] for row in m.tolist()])
        jj=mp.matrix([[mp.mpf(str(x.p))/int(x.q) for x in row] for row in j.tolist()])
        eig,vec=mp.eigsy(mm)
        g=vec[:,0];jground=(g.T*jj*g)[0]
        for row in result['opposite_charge_positive_energy_states']:
            candidate=mp.matrix([mp.mpf(F(x).numerator)/F(x).denominator for x in row['binary_rational_state']])
            norm=(candidate.T*candidate)[0]
            energy=(candidate.T*mm*candidate)[0]/norm-eig[0]
            charge=(candidate.T*jj*candidate)[0]/norm-jground
            for number,key in ((energy,'excess_energy_interval'),(charge,'geometry_charge_interval')):
                lo,hi=[mp.mpf(F(x).numerator)/F(x).denominator for x in row[key]]
                assert lo<number<hi
            assert energy>0
            assert charge<0 if row['sign']<0 else charge>0


def test_elastic_solution_matches_independent_sparse_variation():
    n=5;eye=sparse.eye(n)
    one=sparse.diags((-np.ones(n-1),2*np.ones(n),-np.ones(n-1)),(-1,0,1))
    lap=(sparse.kron(sparse.kron(one,eye),eye)
         +sparse.kron(sparse.kron(eye,one),eye)
         +sparse.kron(sparse.kron(eye,eye),one)).tocsr()
    j=np.zeros((n,n,n));j[2,2,2]=9/140000;j[3,2,2]=9/140000
    direct=-spsolve(lap,j.ravel()).reshape(n,n,n)
    got=force.green_solve(j)
    np.testing.assert_allclose(got,direct,rtol=0,atol=2e-20)
    perturb=np.random.default_rng(1757).normal(size=got.shape)*1e-5
    energy=lambda q:float(q.ravel()@lap@q.ravel()/2+j.ravel()@q.ravel())
    assert energy(got+perturb)>energy(got)
    assert np.linalg.norm(force.laplacian(-got)+j)>1e-4


def test_exact_minimum_enclosure_and_finite_carrier_window():
    r=force.finite_energy_certificate(5,((2,2,2),))
    lo,hi=map(F,r['minimum_energy_interval'])
    assert lo<hi<0
    assert F(r['geometry_max_absolute'])<F(r['finite_geometry_window'])
    # Fully enumerated rational two-site analogue challenges the residual bound.
    a=s.Matrix([[2,-1],[-1,2]]);j=s.Matrix([s.Rational(1,4),s.Rational(1,8)])
    exact_min=-(j.T*a.inv()*j)[0]/2
    points=[s.Matrix([s.Rational(x,16),s.Rational(y,16)]) for x,y in product(range(-8,9),repeat=2)]
    energies=[(q.T*a*q)[0]/2+(j.T*q)[0] for q in points]
    q=points[energies.index(min(energies))];b=a*q+j
    assert energies[energies.index(min(energies))]-(b.T*b)[0]/2<=exact_min<=min(energies)
    with pytest.raises(ValueError):force.finite_energy_certificate(4,((1,1,1),))
    with pytest.raises(ValueError):force.finite_energy_certificate(5,((2,2,2),),geometry_window=F(1,2))


def test_nonlinear_remainder_is_positive_over_whole_declared_window():
    r,matrix,_,j=force.source_matrices()
    window=s.Rational(1,50000)
    bound=(7-3*window)/(14*(1-window))
    for k in range(-10,11):
        q=window*k/10
        actual=(matrix[1,1]-matrix[0,0]).subs(r,1+q)-1-s.Rational(9,14)*q
        assert 0<=actual<=bound*q*q
        assert s.factor(actual-q*q*(7+3*q)/(14*(1+q)))==0
    # Reversing the response charge leaves a first-order error, not the proved quadratic error.
    q=window
    wrong=(matrix[1,1]-matrix[0,0]).subs(r,1+q)-1+s.Rational(9,14)*q
    assert wrong>1000*bound*q*q


def test_force_and_range_controls_do_not_hide_extension_choices():
    value=force.force_run()
    assert all(F(x)>0 for x in value['full_nonlinear_attractive_force_interval_at_six'])
    assert F(value['full_nonlinear_native_cross_energy_interval'][1])<0
    assert value['strength_changes_under_stiffness_change']
    assert value['positive_geometry_mass_screens_long_range_response']
    assert not value['physical_Newton_constant_or_universal_mass_coupling_derived']
