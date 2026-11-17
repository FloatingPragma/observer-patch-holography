"""Outward-interval inertia and full-carrier Schur bounds; no eigensolver."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json

import mpmath as mp
import sympy as s

from model.affine_gauss_mobile import model as mobile
from model.affine_gauss_mobile.certificate import certificate as parent_certificate
from model.affine_gauss_mobile.certificate import validate_pins

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PARENT='model/affine_gauss_mobile/results.json'
PARENT_SHA='46942e200db305e6aa281b836ba2a85b411dc419baa1a60730ae9868e9bdc61e'


def ivf(value):
    value=F(value)
    return mp.iv.mpf(value.numerator)/value.denominator


def endpoint(value, upper=False):
    return mobile.endpoint(value,upper)


def bounds(value, places=36):
    return mobile.rounded(value,places)


@lru_cache(None)
def admitted_parent():
    raw=(ROOT/PARENT).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PARENT_SHA:
        raise ValueError('the actual mobile source changed')
    saved=json.loads(raw)
    validate_pins(saved['input_sha256'])
    if saved!=parent_certificate():
        raise ValueError('mobile source fresh replay differs')
    return saved


@lru_cache(None)
def compression():
    admitted_parent()
    mp.iv.dps=90
    iv=mp.iv;eta=ivf(F(1,24));sq6=iv.sqrt(6)
    # All spectral projector diagonals are identical within each actual edge
    # orbit. This proves the two overlap values apply to all42 source edges.
    for power in mobile.parent.geometry()[8]:
        assert all(power[j,j]==power[0 if j<12 else 12,0 if j<12 else 12] for j in range(42))
    aa=iv.sqrt(3-sq6);bb=iv.sqrt(3+sq6)
    overlaps=[];second=[];variances=[]
    for edge in (0,12):
        variance=mobile.ivalue(mobile.exact_peierls_data(edge)[2])
        sigma=iv.sqrt(variance)
        c2=(iv.cos(sigma*bb)-iv.cos(sigma*aa))/(2*sq6)
        c0=iv.cos(sigma*aa)-(3-sq6)*c2+c2
        two=iv.sqrt(2)*c2
        assert 0<endpoint(c0)<endpoint(c0,True)<1
        assert endpoint(two,True)<0
        overlaps.append(c0);second.append(two);variances.append(variance)
    # Verify the actual cone adjacency, its degree sequence and full spectrum.
    adjacency=s.zeros(13)
    for u,v in mobile.actual_edges():adjacency[u,v]=adjacency[v,u]=1
    x=s.Symbol('x')
    assert adjacency.charpoly(x).as_poly()==s.Poly((x*x-5*x-12)*(x*x-5)**3*(x+1)**5,x)
    assert [sum(adjacency[i,j] for j in range(13)) for i in range(13)]==[12]+[6]*12
    A=[[iv.mpf(0) for _ in range(13)] for _ in range(13)]
    charges,_=mobile.parent.charge_geometry()
    for j in range(13):A[j][j]=mobile.ivalue(charges[j]['total_energy_above_actual_neutral_ground'])
    for e,(u,v) in enumerate(mobile.actual_edges()):A[u][v]=A[v][u]=-eta*overlaps[e>=12]
    wr,wb=[iv.sqrt(1-k*k) for k in overlaps]
    v=(5+iv.sqrt(73))/48
    r=eta*(5*wb+iv.sqrt(25*wb*wb+48*wr*wr))/2
    q_lower=2+iv.sqrt(4-iv.sqrt(5))-v
    # <B|2_f>=1/sqrt(30), with |B>=sum_j|2_j>/sqrt(30).
    dmin=min(-endpoint(a,True) for a in second)
    dressing=F(1,24)**2*dmin*dmin/65
    assert endpoint(r,True)<1 and endpoint(q_lower)>3
    return A,overlaps,second,variances,v,r,q_lower,dressing


def inertia(A, shift):
    """Certified unpivoted LDL inertia whenever every pivot excludes zero."""
    n=len(A);L=[[mp.iv.mpf(0) for _ in range(n)] for _ in range(n)];D=[]
    for i in range(n):
        L[i][i]=mp.iv.mpf(1)
        pivot=A[i][i]-ivf(shift)-sum((L[i][j]*L[i][j]*D[j] for j in range(i)),mp.iv.mpf(0))
        if endpoint(pivot)<=0<=endpoint(pivot,True):
            raise ArithmeticError('uncertain LDL pivot; no inertia claimed')
        D.append(pivot)
        for row in range(i+1,n):
            L[row][i]=(A[row][i]-sum((L[row][j]*L[i][j]*D[j] for j in range(i)),mp.iv.mpf(0)))/pivot
    return sum(endpoint(p,True)<0 for p in D)


def eigen_enclosure(A, ordinal, bits=110):
    lo,hi=F(1),F(3)
    assert inertia(A,lo)==0 and inertia(A,hi)==len(A)
    for _ in range(bits):
        mid=(lo+hi)/2
        if inertia(A,mid)<ordinal:lo=mid
        else:hi=mid
    assert inertia(A,lo)<ordinal<=inertia(A,hi)
    return lo,hi


def schur_lower(a,c,r):
    a,c,r=map(ivf,(a,c,r))
    return endpoint((a+c-mp.iv.sqrt((c-a)**2+4*r*r))/2)


@lru_cache(None)
def spectral_bounds():
    A,overlaps,second,variances,v,r,c,dressing=compression()
    c_lo=endpoint(c);r_hi=endpoint(r,True)
    rows=[]
    for j in range(1,14):
        low,high=eigen_enclosure(A,j)
        full_low=schur_lower(low,c_lo,r_hi)
        full_high=high-dressing if j==1 else high
        assert 0<full_low<full_high<c_lo
        rows.append({'ordinal_in_one_charge_sector':j,
            'zero_photon_compression_energy_interval':[str(low),str(high)],
            'actual_H_energy_interval':[str(full_low),str(full_high)]})
    low1,high1=map(F,rows[0]['actual_H_energy_interval'])
    low2=F(rows[1]['actual_H_energy_interval'][0])
    trial_high=F(rows[0]['zero_photon_compression_energy_interval'][1])
    assert low2>high1
    fidelity=(low2-trial_high)/(low2-low1)
    photon_upper=r_hi*r_hi/((c_lo-trial_high)**2+r_hi*r_hi)
    photon_lower=dressing*dressing/(4*r_hi*r_hi)
    assert 0<fidelity<1 and 0<photon_upper<1
    assert 0<photon_lower<photon_upper
    return {'low_charged_levels':rows,
        'fourteenth_and_higher_energy_lower':str(c_lo),
        'one_charge_ground_is_unique':True,
        'one_charge_ground_gap_lower':str(low2-high1),
        'strict_ground_energy_reduction_by_photon_dressing_lower':str(dressing),
        'ground_nonvacuum_transverse_probability_upper':str(photon_upper),
        'ground_nonvacuum_transverse_probability_lower':str(photon_lower),
        'ground_preparation':{
            'trial':'positive normalized ground of the actual zero-photon compression, tensored with photon vacuum',
            'pump':'D_u/2; compensated creation in the declared normalized spatial mode u',
            'source_vacuum_is_actual_interacting_ground':True,
            'conditional_trial_fidelity_lower':str(fidelity),
            'pump_failure_probability':'3/4',
            'ground_window_hit_probability_interval':[str(fidelity/4),'1/4'],
            'energy_complement_probability_interval':['0',str((1-fidelity)/4)],
            'exact_law':['3/4','<u,0|P_ground(H_eta)|u,0>/4','(1-<u,0|P_ground(H_eta)|u,0>)/4'],
            'ground_window':rows[0]['actual_H_energy_interval'],
            'success_state':'the unique actual one-charge ground',
            'finite_spectral_preparation_and_detector_services_declared':True,
            'efficient_native_circuit_claimed':False},
        'reference_energy_retained':'1',
        'source_energy_after_subtracting_fixed_reference':[str(low1-1),str(high1-1)],
        'operator_bounds':{'hopping_norm_upper':bounds(v),'off_diagonal_block_norm_upper':bounds(r),
            'nonvacuum_block_energy_lower':bounds(c),'radial_and_surface_vacuum_overlaps':[bounds(a) for a in overlaps],
            'radial_and_surface_two_photon_amplitudes':[bounds(a) for a in second],
            'radial_and_surface_field_variances':[bounds(a) for a in variances]}}


def certificate():
    parent=admitted_parent()
    result={'schema':'oph.affine_gauss_charged_spectrum.v1',
        'parent_receipt':PARENT,'parent_receipt_sha256':PARENT_SHA,
        'same_H_state_reference_and_clock':True,
        'source':'the exact finite mobile affine-Gauss source; eta=1/24; fixed compensator at vertex1',
        'one_charge_dimension':70928,'full_physical_dimension':44695552,
        'numerical_eigensolver_used_for_certification':False,
        'method':'outward interval LDL inertia for the13-state compression; analytic full-complement Schur bounds',
        'spectral':spectral_bounds(),
        'scope':{'finite_boundary_reference_problem':True,'relativistic_rest_mass_identified':False,
                 'physical_spin_identified':False,'continuum_limit_claimed':False,
                 'physical_particle_map':[],'mass_predictions_GeV':[]}}
    pins=dict(parent['input_sha256']);pins[PARENT]=PARENT_SHA
    pins.update({p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')})
    result['input_sha256']=pins
    return result
