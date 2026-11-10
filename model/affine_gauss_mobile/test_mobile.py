"""Literal finite-matrix controls for the coupled source estimates."""
from fractions import Fraction as F
from pathlib import Path
import importlib.util
import itertools
import math
import numpy as np
import scipy.linalg as la
import sympy as s

from . import model as source


def test_actual_parent_and_global_ground_bound():
    parent,sha=source.source()
    row=source.universal_bounds()
    assert row['unique_actual_physical_ground_unchanged']
    assert row['actual_physical_ground_gap']==s.sqrt(4-s.sqrt(5))
    assert row['degree']==[12]+[6]*12
    assert row['neutral_field_records_unchanged_exactly']
    assert not row['static_exact_charged_energy_line_law_transferred']
    assert parent['architecture']['physical_source_dimension']==44695552


def test_rational_uniform_arrival_and_field_occupation_margins():
    rows=source.universal_bounds()
    arrival=rows['neighbor_arrival'];joint=rows['neighbor_arrival_and_transverse_occupation']
    assert F(arrival['unconditional_probability_lower'])==F(2209,339738624)>0
    assert F(joint['unconditional_probability_lower'])==F(169,347892350976)>0
    assert len(arrival['complete_records'])==3 and len(joint['complete_records'])==5
    for row in rows['energy_records']:
        assert F(row['unconditional_pump_and_energy_bin_hit_lower'])>=F(1,6)


def test_four_level_Peierls_norm_is_exact_without_Weyl_relations():
    z=s.zeros(4)
    for j in range(1,4):z[j-1,j]=z[j,j-1]=s.sqrt(j)
    x=s.Symbol('x')
    assert z.charpoly(x).as_expr()==x**4-6*x*x+3
    assert z**4-6*z*z+3*s.eye(4)==s.zeros(4)
    for edge in (0,12):
        frequencies,f2,variance,pt=source.exact_peierls_data(edge)
        assert variance>0 and variance<s.Rational(1,2)
        assert (3+s.sqrt(6))*variance<s.Rational(11,4)
        assert pt/s.Integer(72)>s.Rational(1,108)


def test_first_integral_intervals_cover_literal_165_state_matrix_control():
    edge=12;y=1;x=2;time=F(1,4)
    record=source.first_dyson_interval(edge,y,x,time)
    frequencies,f2,variance,pt=source.exact_peierls_data(edge)
    basis=list(source.total_occupations(8,3));ids={n:j for j,n in enumerate(basis)}
    A=np.zeros((165,165));w=np.array([float(s.N(v,30)) for v in frequencies])
    for col,n in enumerate(basis):
        for j,value in enumerate(n):
            if value:
                target=list(n);target[j]-=1
                amplitude=math.sqrt(value*float(s.N(f2[j],30)))
                A[ids[tuple(target)],col]+=amplitude
                A[col,ids[tuple(target)]]+=amplitude
    vacuum=np.zeros(165);vacuum[0]=1
    state=la.expm(1j*A)@vacuum
    charges,_=source.parent.charge_geometry()
    delta=float(charges[x]['static_electric_energy']-charges[y]['static_electric_energy'])+np.array([np.dot(n,w) for n in basis])
    response=float(time)*np.sinc(delta*float(time)/(2*np.pi))
    first=float(source.ETA)*response*state
    norm=np.linalg.norm(first);photon=np.linalg.norm(first[1:])
    for actual,key in ((norm,'first_arrival_norm_interval'),(photon,'first_arrival_nonvacuum_norm_interval')):
        low,high=map(F,record[key])
        assert float(low)-1e-13<=actual<=float(high)+1e-13
    assert F(record['unconditional_arrival_nonvacuum_probability_lower'])>0


def test_literal_full_three_mode_CAR_hopping_ground_residual_and_complete_law():
    # Independent three-vertex cycle control, one transverse oscillator with
    # total cutoff3. This tests the same algebraic architecture at small size.
    sites=3;photon=4;eta=float(source.ETA);omega=math.sqrt(3)
    annih=[]
    for y in range(sites):
        c=np.zeros((8,8))
        for mask in range(8):
            if mask>>y&1:c[mask^(1<<y),mask]=(-1)**((mask&((1<<y)-1)).bit_count())
        annih.append(c)
    number=sum(c.T@c for c in annih)
    a=np.diag(np.sqrt(np.arange(1,photon)),1);q=(a+a.T)/math.sqrt(2*omega)
    hp=np.diag(omega*np.arange(photon))
    G=np.eye(3)/3-np.ones((3,3))/9
    charge_energy=[]
    for mask in range(8):
        n=np.array([mask>>j&1 for j in range(3)]);rho=n.copy();rho[0]-=sum(n)
        charge_energy.append(rho@G@rho/2)
    H0=np.kron(2*number+np.diag(charge_energy),np.eye(photon))+np.kron(np.eye(8),hp)
    H=H0.astype(complex)
    for u,v,sign in ((0,1,1),(1,2,1),(0,2,-1)):
        W=la.expm(1j*sign*q/math.sqrt(3));hop=np.kron(annih[v].T@annih[u],W)
        H-=eta*(hop+hop.conj().T)
    vacuum=np.zeros(32);vacuum[0]=1
    assert np.linalg.norm(H@vacuum)==0
    assert np.linalg.eigvalsh(H)[1]>1
    C=np.kron(annih[0].T/2,np.eye(photon))
    initial=2*C@vacuum
    energy=float(np.vdot(initial,H@initial).real)
    assert abs(energy-2)<1e-14
    assert abs(np.linalg.norm(H@initial-energy*initial)**2-2*eta**2)<1e-14
    U=la.expm(-.25j*H);P=np.kron(annih[1].T@annih[1],np.eye(photon))
    hit=float(np.linalg.norm(P@U@C@vacuum)**2)
    assert hit>float(F(2209,339738624))
    fail=np.diag(np.sqrt(np.diag(np.eye(32)-C.conj().T@C)))
    kraus=(fail,P@U@C,(np.eye(32)-P)@U@C)
    assert np.max(abs(sum(k.conj().T@k for k in kraus)-np.eye(32)))<1e-14
    assert abs(sum(np.linalg.norm(k@vacuum)**2 for k in kraus)-1)<1e-14


def test_literal_number_current_with_unrelated_finite_link_phases():
    # No mutual commutation of different link matrices is used by continuity.
    c0=np.array([[0,1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]],complex)
    c1=np.array([[0,0,1,0],[0,0,0,-1],[0,0,0,0],[0,0,0,0]],complex)
    W=la.expm(1j*np.array([[.2,.3],[.3,-.1]]))
    T=np.kron(c1.conj().T@c0,W)
    H=-float(source.ETA)*(T+T.conj().T)
    J=1j*float(source.ETA)*(T-T.conj().T)
    n1=np.kron(c1.conj().T@c1,np.eye(2))
    n0=np.kron(c0.conj().T@c0,np.eye(2))
    assert np.max(abs(1j*(H@n1-n1@H)-J))<1e-15
    assert np.max(abs(1j*(H@n0-n0@H)+J))<1e-15



def test_fresh_scientific_calculation_matches_preserved_independent_receipt():
    from .certificate import compare_research,scientific_fields
    value=source.calculate()
    assert compare_research(value)['all_scientific_fields_equal']
    import copy
    import pytest
    corrupted=copy.deepcopy(value)
    key='architecture' if 'architecture' in value else 'source'
    corrupted[key]['clock_mutation']='a different source clock'
    with pytest.raises(ValueError):compare_research(corrupted)
