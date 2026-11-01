from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import momentum as m
from ..coarse_energy_ward import ward as w
from ..full_fock_field.field import FieldGeometry,FullFock


@pytest.mark.parametrize('sides',[(3,1,1),(1,3,1),(1,1,3),(2,2,1)])
def test_independent_complete_matrix_force_and_stress(sides):
    geometry=FieldGeometry(sides)
    data=m.construct(geometry)
    space=FullFock(geometry.modes)
    h=space.hamiltonian(geometry).toarray()
    potential=np.diag([sum(((mask>>(2*x))&1)*((mask>>(2*x+1))&1)
                           for x in range(geometry.volume)) for mask in range(space.dimension)])
    np.testing.assert_array_equal(data['H'].matrix(geometry.modes),h)
    for column in data['columns']:
        divergence=w.divergence(column['stress'],geometry.volume)
        for x,p in enumerate(column['local']):
            pm=p.matrix(geometry.modes)
            fm=column['force'][x].matrix(geometry.modes)
            np.testing.assert_allclose(1j*(potential@pm-pm@potential),fm,atol=0,rtol=0)
            np.testing.assert_allclose(1j*(h@pm-pm@h)+divergence[x].matrix(geometry.modes),fm,atol=0,rtol=0)
            assert np.linalg.norm(pm,2)<=1+1e-14
            assert np.linalg.norm(fm,2)<=1+1e-14
        for current in column['stress'].values():
            assert np.linalg.norm(current.matrix(geometry.modes),2)<=2*geometry.speed+1e-14


def test_all_three_cartesian_differences_commute_on_open_box():
    geometry=FieldGeometry((3,2,2))
    h,p,_=m.one_particle(geometry)
    for left,right in product(p,repeat=2):
        assert left*right==right*left
    sigma=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    expected=sum((-s.Rational(1,8)*s.kronecker_product(s.eye(geometry.volume),spin)*q
                  for spin,q in zip(sigma,p)),s.zeros(geometry.modes))
    assert h==expected


def test_three_dimensional_exact_forces_and_blocks():
    geometry=FieldGeometry((3,3,3))
    data=m.construct(geometry)
    assert [len(z['stress']) for z in data['columns']]==[54]*3
    assert m.force_gram(data['columns'])==tuple(tuple(F(9,4) if i==j else F() for j in range(3)) for i in range(3))
    for axis,column in enumerate(data['columns']):
        assert w.commutator(data['H0'],column['total']).is_zero()
        assert not w.commutator(data['H'],column['total']).is_zero()
        assert all(mask.bit_count()==4 for value in column['force'] for mask in value.terms)
        assert all(mask.bit_count()==2 for value in column['stress'].values() for mask in value.terms)
        for (x,y),current in column['stress'].items():
            assert len(current.support())<=6
            assert all(max(abs(a-b) for a,b in zip(geometry.sites[z],geometry.sites[x]))<=2 for z in current.support())
        labels=(0,)*geometry.volume
        assert not w.aggregate_flows(column['stress'],labels)
        assert w.aggregate(column['force'],labels)==(column['total_force'],)


@pytest.mark.parametrize('sides',[(2,1,1),(3,2,1),(2,2,2)])
def test_force_gram_counts_edges_at_arbitrary_rational_contact(sides):
    geometry=FieldGeometry(sides,speed=.25,interaction=.5)
    gram=m.force_gram(m.construct(geometry)['columns'])
    for i in range(3):
        count=(sides[i]-1)*np.prod([sides[j] for j in range(3) if j!=i])
        for j in range(3):
            assert gram[i][j]==(F(int(count),32) if i==j else 0)


def test_force_vanishes_only_on_free_or_one_particle_controls():
    free=m.construct(FieldGeometry((3,1,1),interaction=0.))
    assert all(column['total_force'].is_zero() for column in free['columns'])
    interacting=m.construct(FieldGeometry((3,1,1)))
    force=interacting['columns'][0]['total_force'].matrix(6)
    low=[0]+[1<<p for p in range(6)]
    assert np.count_nonzero(force[:,low])==0
    assert np.linalg.norm(force,2)>0


def test_exact_doublon_witness_recovers_all_offsite_coefficients():
    for sites in (2,3,6):
        witnesses=m.no_go_matrix_elements(sites)
        assert len(witnesses)==4*sites*(sites-1)
        assert len({(p,q) for p,q,_,_,_ in witnesses})==len(witnesses)
        assert {sign for *_,sign in witnesses}<={-1,1}


def test_hidden_spin_and_number_exhaust_onsite_transport_dimension():
    geometry=FieldGeometry((2,2,2))
    basis=m.conserved_quadratics(geometry)
    assert len(basis)==4
    potential=m.contact(geometry)
    for q in basis:
        assert m.w.commutator(potential,m.lift(q)).is_zero()
    # Connected bonds transport a freely chosen 2x2 root matrix. Verify the
    # full finite linear constraint rank independently on this small cube.
    h,_,_=m.one_particle(geometry)
    columns=[]
    for x in range(geometry.volume):
        for a,b in product(range(2),repeat=2):
            q=s.zeros(geometry.modes)
            q[2*x+a,2*x+b]=1
            columns.append(s.Matrix(h*q-q*h).reshape(geometry.modes**2,1))
    constraints=s.Matrix.hstack(*columns)
    assert len(columns)-constraints.rank()==4


def test_force_integral_is_required_in_finite_cell():
    from .packet import time_cell_diagnostic
    result=time_cell_diagnostic()
    assert result['forced_balance_residual']<1e-12
    assert result['omitted_force_residual']>.01
    assert result['point_force_sample_residual']>1e-4


@pytest.mark.parametrize('kwargs',[dict(radius=2),dict(duration=-1),dict(duration=.1),
                                   dict(steps=0),dict(samples=0),dict(speed=.125),
                                   dict(interaction=-1),dict(duration=100)])
def test_reject_invalid_error_domain(kwargs):
    args=dict(radius=12,duration=F(1,8),steps=10000,samples=100000)
    args.update(kwargs)
    with pytest.raises(ValueError):m.cell_error(**args)


def test_error_bound_converges_and_retains_force():
    coarse=m.cell_error(12,F(1,8),10000,samples=100000)
    fine=m.cell_error(18,F(1,8),10000000,samples=10**12)
    assert fine['forced_cell_residual_bound']<F(1,10**10)
    assert fine['forced_cell_residual_bound']<coarse['forced_cell_residual_bound']
    assert fine['integrated_force_error']>0
    assert m.cell_error(12,0,100)['forced_cell_residual_bound']==0


def test_retained_receipt_is_exactly_source_bound():
    from tools import strict_json
    from model.observer.scalar_response.experiment import compare_receipt
    from . import packet
    compare_receipt(strict_json.loads(packet.OUTPUT.read_text()),packet.run())


def test_malformed_supplied_receipt_path_is_rejected(tmp_path):
    import subprocess
    import sys
    path=tmp_path/'bad.json'
    path.write_text('{"all_pass":true}')
    result=subprocess.run([sys.executable,'-m','model.lattice_momentum_balance','--check',str(path)],
                          capture_output=True,text=True)
    assert result.returncode!=0
