"""Finite free Dirac sea and a shared Gauss dressing on the actual routes.

The energy is an explicitly chosen principal one-body phase lift. No
interacting ground state or physical species identification is asserted.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'spin'),str(HERE.parent/'vector_field')]
import dirac_proposal as dp
import native_field as nf


def pairs():
    """One occupied negative-band orbital per antipodal/spin pair."""
    return [(4*c+s,4*(7-c)+s,-1 if s%2==0 else 1)
            for c in range(4) for s in range(4)]


def reference_charge(kind='vacuum'):
    G=np.zeros(96,dtype=np.int64);G[::12]=2
    if kind=='particle':G[0]+=1
    elif kind=='hole':G[0]-=1
    elif kind not in ('vacuum','neutral_pair'):raise ValueError('named reference sector required')
    return G


def sector_pattern(kind='vacuum',particle=0,hole=1):
    if kind not in ('vacuum','particle','hole','neutral_pair'):
        raise ValueError('named finite excitation required')
    if any(type(j) is not int or not 0<=j<16 for j in (particle,hole)):
        raise ValueError('actual orbital index0..15 required')
    if kind=='neutral_pair' and particle==hole:
        raise ValueError('distinct particle/hole orbital slots required')
    pattern=[1]*16
    if kind in ('particle','neutral_pair'):pattern[particle]=2
    if kind in ('hole','neutral_pair'):pattern[hole]=0
    return pattern


def occupations(mask):
    if type(mask) is not int or not 0<=mask<1<<32:raise ValueError('32 mode bits required')
    n=np.zeros(96,dtype=np.int64)
    for cell in range(8):n[12*cell]=((mask>>(4*cell))&15).bit_count()
    return n


def complete_occupation_law(kind='vacuum',particle=0,hole=1):
    """Every original outcome; exact integer numerator and common denominator."""
    pattern=sector_pattern(kind,particle,hole)
    free=[j for j,v in enumerate(pattern) if v==1]
    answer=[]
    for choice in range(1<<len(free)):
        occupied=[];coefficient=1;cursor=0
        for j,(a,b,sign) in enumerate(pairs()):
            if pattern[j]==0:continue
            if pattern[j]==2:occupied.extend((a,b));continue
            bit=(choice>>cursor)&1;cursor+=1
            occupied.append(b if bit else a)
            if bit:coefficient*=sign
        inversions=sum(occupied[i]>occupied[j] for i in range(len(occupied)) for j in range(i+1,len(occupied)))
        coefficient*=(-1)**inversions
        answer.append((sum(1<<mode for mode in occupied),coefficient))
    assert len({mask for mask,sign in answer})==len(answer)
    return answer,len(answer)


def particular(mask,G):
    n=occupations(mask)
    if int(n.sum())!=int(np.sum(G)):raise ValueError('torus Gauss total mismatch')
    P=dp.field.current_walk.geometry()['P']
    F0=P@(n-G)
    assert np.array_equal(dp.field.data()[1]@F0,G-n)
    return F0


def dressed_field(mask,z,G,modulus=16):
    C=dp.field.current_walk.geometry()['C']
    z=np.asarray(z,dtype=np.int64)
    if z.shape!=(193,):raise ValueError('all193 shared cycle coordinates required')
    return (particular(mask,G)+C@z)%modulus


def rational_gaussian_amplitudes(base=16):
    """Positive exact integer amplitudes proportional to((base-1)/base)^(z²)."""
    if type(base) is not int or not 2<=base<=64:raise ValueError('base2..64 required')
    centered=[(z+8)%16-8 for z in range(16)]
    return [base**(64-z*z)*(base-1)**(z*z) for z in centered]


def cycle_overlap(shift,amplitudes):
    Z=sum(a*a for a in amplitudes);q=len(amplitudes)
    return F(sum(amplitudes[z]*amplitudes[(z+shift)%q] for z in range(q)),Z)


def fibre_overlap(shift,amplitudes):
    result=F(1)
    for h in shift:
        if h:result*=cycle_overlap(int(h),amplitudes)
    return result


def gaussian_energy(mask,G,amplitudes):
    """Complete finite diagonal expectation via independent cycle marginals."""
    g=dp.field.current_walk.geometry();C=g['C'];P0=particular(mask,G)
    weights=np.array([float(a*a) for a in amplitudes]);weights/=weights.sum()
    distribution_cache={}
    energy=0.
    for e,row in enumerate(C):
        pos=int(np.count_nonzero(row==1));neg=int(np.count_nonzero(row==-1))
        key=(pos,neg)
        if key not in distribution_cache:
            distribution=np.zeros(16);distribution[0]=1.
            for sign,count in ((1,pos),(-1,neg)):
                kernel=weights[np.array([(sign*z)%16 for z in range(16)])]
                for _ in range(count):
                    distribution=np.array([sum(distribution[j]*kernel[(k-j)%16] for j in range(16)) for k in range(16)])
            distribution_cache[key]=distribution
        values=(np.arange(16)+int(P0[e])+8)%16-8
        energy+=float(distribution_cache[key]@(values*values))/2
    return energy


def run():
    paths=[Path(__file__),Path(dp.__file__),Path(dp.field.current_walk.__file__),
           HERE.parent/'spin/MANIFEST.json',HERE.parent/'vector_field/MANIFEST.json',
           HERE.parent/'units/MANIFEST.json']
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    mu=.04;h=.1;m=.4
    zero=np.zeros(288)
    U=np.kron(np.eye(8),np.diag(np.exp(-1j*mu*np.diag(dp.BETA).real)))
    for axis in range(3):U=dp.background_matrix(axis,zero)@U
    direct=np.zeros((32,32),complex)
    for cell in range(8):direct[4*(7-cell):4*(7-cell)+4,4*cell:4*cell+4]=np.diag(np.exp(-1j*mu*np.diag(dp.BETA).real))
    assert np.linalg.norm(U-direct)<2e-13
    minus=np.zeros((32,16));plus=minus.copy();negative=[];positive=[]
    for j,(a,b,sign) in enumerate(pairs()):
        minus[a,j]=plus[a,j]=1/np.sqrt(2)
        minus[b,j]=sign/np.sqrt(2);plus[b,j]=-sign/np.sqrt(2)
        beta=1 if a%2==0 else -1
        energy_minus=-((math.pi-mu)/h if beta==1 else mu/h)
        energy_plus=(mu/h if beta==1 else (math.pi-mu)/h)
        negative.append(energy_minus);positive.append(energy_plus)
        assert np.linalg.norm(U@minus[:,j]-np.exp(-1j*h*energy_minus)*minus[:,j])<2e-13
        assert np.linalg.norm(U@plus[:,j]-np.exp(-1j*h*energy_plus)*plus[:,j])<2e-13
    V=np.column_stack((minus,plus));assert np.linalg.norm(V.T@V-np.eye(32))<2e-13
    P=minus@minus.T
    assert np.linalg.norm(P@P-P)<2e-13 and np.max(abs(np.diag(P)-.5))<2e-13
    law,den=complete_occupation_law();lookup=dict(law)
    G=reference_charge();sample=[]
    for i,(mask,coefficient) in enumerate(law):
        selected=[j for j in range(32) if (mask>>j)&1]
        permuted=[4*(7-j//4)+j%4 for j in selected]
        inversions=sum(permuted[a]>permuted[b] for a in range(16) for b in range(a+1,16))
        target=sum(1<<j for j in permuted)
        assert lookup[target]==coefficient*(-1)**inversions
        assert sum(1 if j%2==0 else -1 for j in selected)==0
        if i in (0,1,17,511,4095,65535):sample.append(mask)
    # General Gauss code checked on distinct complete-record configurations.
    rng=np.random.default_rng(1414);geometry=dp.field.current_walk.geometry()
    for mask in sample:
        for _ in range(4):
            z=rng.integers(0,16,size=193)
            field=dressed_field(mask,z,G)
            assert not np.any((occupations(mask)+geometry['B']@field-G)%16)
    # Every actual Dirac route gives a closed correction in the same fibre.
    routes=[]
    for cell,axis,direction in product(range(8),range(3),(-1,1)):
        other,path,current=dp.route(cell,axis,direction)
        correction=geometry['P'][:,12*cell]+current-geometry['P'][:,12*other]
        shift=correction[geometry['chords']]
        assert np.array_equal(correction,geometry['C']@shift)
        routes.append(shift)
    gaussians=[]
    for base in (4,8,16,32,64):
        amplitudes=rational_gaussian_amplitudes(base)
        overlaps=[fibre_overlap(v,amplitudes) for v in routes]
        # Orthogonal rail projectors and the route permutation make the
        # one-axis isometry defect the maximum branch translation defect.
        error_by_axis=[]
        for axis in range(3):
            worst=min(overlaps[j] for j,(c,a,d) in enumerate(product(range(8),range(3),(-1,1))) if a==axis)
            error_by_axis.append(math.sqrt(float(2*(1-worst))))
        gaussians.append({'base':base,'cycle_shift1_overlap_exact':str(cycle_overlap(1,amplitudes)),
             'worst_route_overlap':float(min(overlaps)), 'one_macro_free_transport_norm_bound':min(2.,sum(error_by_axis)),
             'energy_on_six_named_occupation_records':[gaussian_energy(mask,G,amplitudes) for mask in sample]})
    # Generic electric phases do not preserve the uniform affine code.
    face=dp.field.data()[2][:,0]
    assert np.count_nonzero(face)==3 and np.max(abs(face))==1
    electric_phase_contrast=abs(1-np.exp(-1j*(1/16)*3/2))
    assert electric_phase_contrast>.09
    outcomes={}
    for kind in ('vacuum','particle','hole','neutral_pair'):
        records,d=complete_occupation_law(kind)
        expected=16+(kind=='particle')-(kind=='hole')
        assert all(mask.bit_count()==expected for mask,sign in records)
        assert d==len(records)
        sectorG=reference_charge(kind)
        for mask,sign in records[:3]:particular(mask,sectorG)
        outcomes[kind]={'outcomes':d,'probability_each_exact':str(F(1,d)),
                        'total_occupation':expected,'reference_total':int(sectorG.sum())}
    # Fixed G cannot silently admit a lone extra charged particle on a torus.
    rejected=False
    try:particular(complete_occupation_law('particle')[0][0][0],G)
    except ValueError:rejected=True
    assert rejected
    assert pins=={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result={'schema':'oph.dirac.shared_gauss_reference.v1',
      'checks':{'actual_zero_character32mode_source':True,'all_negative_and_positive_orbitals':True,
               'full65536vacuum_occupation_coefficients':True,'fermionic_spectator_signs':True,
               'shared_reference_Gauss_dressing':True,'all48route_cycle_corrections':True,
               'finite_gaussian_preparations_and_translation_bounds':True,
               'uniform_fibre_electric_kick_control':True,'complete_excitation_occupation_laws':True,
               'wrong_charged_reference_rejected':True,'source_bytes_unchanged':True},
      'candidate':'explicit CAR lift of native-local12-dirac-routed-v3-proposed',
      'energy_lift':'K=-(1/h) Arg(U) on the principal one-body branch; H_F=dGamma(K)-E_sea',
      'parameters':{'h':h,'m':m,'mu':mu,'modulus':16},
      'negative_onebody_energies':negative,'positive_onebody_energies':positive,
      'sea_energy':sum(negative),'vacuum_macro_eigenphase':'1; exp(i8pi)=1',
      'positive_excitation_gap_in_declared_lift':m,'neutral_pair_threshold_in_declared_lift':2*m,
      'low_excitation_degeneracy_each_charge':8,'zero_and_pi_valleys_retained':True,
      'occupation_laws':outcomes,'vacuum_reference':'G=2 on all8root ports,0 elsewhere',
      'gaussian_dressing':gaussians,
      'uniform_fibre_electric_energy_exact':str(F(288*(16**2+2),24)),
      'uniform_fibre_energy_formula':'288*(m_field^2+2)/24; every edge marginal is uniform since graph has no bridges',
      'electric_phase_noninvariance_contrast':electric_phase_contrast,
      'dressing_identity':'J_g|n>=|n> sum_z g(z)|P(n-G)+Cz>; distinct n remain orthogonal; BC=0 and BP=delta0−I',
      'translation_bound':'Each axis intertwining defect <= max_route sqrt(2−2 Re product_j <g_j,T_hj g_j>); sum the three axis bounds per macro and over macros. A full Fock bound must sum all transported rail currents, rather than reuse the one-particle bound.',
      'scope':['Free Gamma(U) reference and its exact Gauss dressing at zero electric coefficient; homogeneous plaquette kicks act as a retained scalar on a uniform zero-character fibre.',
               'Generic electric phases destroy uniform-fibre invariance. No interacting eigenstate, vacuum, pole or physical species is proved.',
               'The additive Fock energy lift is declared; discrete-time Gamma(U) alone does not choose this logarithm or its vacuum.',
               'Gaussian product cycles are finite-energy preparations with an explicit finite free-transfer error, not claimed ground states.',
               'The new32mode Fock allocation differs from the encodedoneparticle Dirac allocation. Native preparation is provided separately; fullFock evolution/error is not inherited from that compiler.'],
      'source_sha256':pins,'official_rows_changed':[]}
    (HERE/'REFERENCE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'occupations':outcomes,'gaussians':gaussians,
                      'uniform_energy':result['uniform_fibre_electric_energy_exact']}))


if __name__=='__main__':run()
