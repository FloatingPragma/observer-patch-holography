"""Finite native collective matrices, Gaussian comparison and dilute KG bounds."""
from functools import lru_cache
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
RER=HERE.parents[1]
sys.path.insert(0,str(RER))
import mpmath as mp
import numpy as np
import sympy as sy
from scipy import sparse
from scipy.sparse.linalg import eigsh

from .reference import apparatus_capacity, certificate as reference_certificate

PARENT='model/interacting_block_propagation/results.json'


def source():
    value=json.loads((RER/PARENT).read_text())
    pins=dict(value['input_sha256'])
    pins[PARENT]=hashlib.sha256((RER/PARENT).read_bytes()).hexdigest()
    for path,expected in pins.items():
        assert hashlib.sha256((RER/path).read_bytes()).hexdigest()==expected,path
    return value,pins


def occupations(N,colors=3):
    return [n for n in itertools.product(range(N+1),repeat=colors) if sum(n)<=N]


def exact_collective(N):
    """Compile the real replicated |g><e_alpha| units, then compress exactly."""
    basis=occupations(N)
    ids={n:i for i,n in enumerate(basis)}
    raw=list(itertools.product(range(4),repeat=N))
    W=sy.zeros(4**N,len(basis))
    multiplicity={n:math.factorial(N)//(math.factorial(N-sum(n))*math.prod(math.factorial(k) for k in n)) for n in basis}
    for i,word in enumerate(raw):
        n=tuple(word.count(a) for a in (1,2,3))
        W[i,ids[n]]=1/sy.sqrt(multiplicity[n])
    assert W.H*W==sy.eye(len(basis))
    lowers=[]
    for color in range(3):
        primitive=sy.zeros(4);primitive[0,color+1]=1
        collective=sy.zeros(4**N)
        for i in range(N):
            collective+=sy.kronecker_product(*[primitive if j==i else sy.eye(4) for j in range(N)])/sy.sqrt(N)
        expected=sy.zeros(len(basis))
        for j,n in enumerate(basis):
            if n[color]:
                target=list(n);target[color]-=1
                expected[ids[tuple(target)],j]=sy.sqrt(n[color])*sy.sqrt(1-sy.Rational(sum(n)-1,N))
        assert sy.simplify(W.H*collective*W-expected)==sy.zeros(len(basis))
        assert sy.simplify(collective*W-W*expected)==sy.zeros(4**N,len(basis))
        lowers.append(expected)
    anticommutator=sum((b.H*b+b*b.H for b in lowers),sy.zeros(len(basis)))
    expected=sy.diag(*[2*sum(n)+3-sy.Rational(2*sum(n)**2+2*sum(n),N) for n in basis])
    assert sy.simplify(anticommutator-expected)==sy.zeros(len(basis))
    return {'constituent_blocks_per_cell':N,'raw_low_dimension':4**N,
            'symmetric_dimension':len(basis),'exact_collective_columns':3*len(basis),
            'symmetric_sector_reducing':True,'HP_identity_exact':True,
            'quadrature_occupation_bound_exact':True}


def sparse_collective(N):
    basis=occupations(N);ids={n:i for i,n in enumerate(basis)}
    lowers=[]
    for color in range(3):
        row=[];col=[];data=[]
        for j,n in enumerate(basis):
            if n[color]:
                target=list(n);target[color]-=1
                row.append(ids[tuple(target)]);col.append(j)
                data.append(math.sqrt(n[color]*(1-(sum(n)-1)/N)))
        lowers.append(sparse.csr_matrix((data,(row,col)),shape=(len(basis),len(basis))))
    return basis,lowers


def gaussian_number_moments(occupation,one_particle=False):
    """Exact squeezed-vacuum cumulants; largest mode bounds one excitation."""
    ns=[mp.mpf(n) for n in occupation]
    k1=sum(ns);k2=sum(2*n*(n+1) for n in ns)
    k3=sum(4*n*(n+1)*(2*n+1) for n in ns)
    k4=sum(8*n*(n+1)*(6*n*n+6*n+1) for n in ns)
    if one_particle:
        n=max(ns)
        k1+=1+2*n;k2+=4*n*(n+1)
        k3+=8*n*(n+1)*(2*n+1)
        k4+=16*n*(n+1)*(6*n*n+6*n+1)
    m2=k1*k1+k2
    center=k1+2
    m4=center**4+6*center*center*k2+3*k2*k2+4*center*k3+k4
    return {'mean':k1,'second_moment':m2,'n_plus_two_fourth_moment':m4}


def free_parameters(a,side):
    a=mp.mpf(a);D=1/a;t=(1-a*a)/(12*a)
    frequencies=[]
    for menu in itertools.product(range(side),repeat=3):
        cosine=sum(mp.cos(2*mp.pi*n/side) for n in menu)
        frequencies.append(mp.sqrt(D*D-4*D*t*cosine))
    occ=[(D/w+w/D-2)/4 for w in frequencies for _ in range(3)]
    return D,t,frequencies,occ


def rational_moment_bounds(a,side,one_particle=False):
    a=F(a);modes=3*side**3
    u=(a+1/a-2)/4
    k1=modes*u;k2=modes*2*u*(u+1)
    k3=modes*4*u*(u+1)*(2*u+1)
    k4=modes*8*u*(u+1)*(6*u*u+6*u+1)
    if one_particle:
        k1+=1+2*u;k2+=4*u*(u+1)
        k3+=8*u*(u+1)*(2*u+1)
        k4+=16*u*(u+1)*(6*u*u+6*u+1)
    center=k1+2
    return {'mean':k1,'second_moment':k1*k1+k2,
            'n_plus_two_fourth_moment':center**4+6*center*center*k2+3*k2*k2+4*center*k3+k4}


def ceil_fraction(value):
    value=F(value)
    return -(-value.numerator//value.denominator)


def ceil_sqrt_fraction(value):
    n=math.isqrt(ceil_fraction(value))
    return n if n*n>=value else n+1


def dec(value):
    value=F(value)
    return mp.mpf(value.numerator)/value.denominator


def controlled_resources(a,side):
    """Every resource decision uses exact rational upper bounds and margins."""
    a=F(a);D=1/a;t=(1-a*a)/(12*a)
    V=side**3;z=6;edges=z*V//2
    vac=rational_moment_bounds(a,side)
    one=rational_moment_bounds(a,side,True)
    M2=vac['second_moment']
    kicks=2;segments=3
    displacement_sum=F(1,4)
    R=displacement_sum/a
    displacement_factor=1+8*R+6*R**2+4*R**3+R**4
    M4=displacement_factor**2*one['n_plus_two_fourth_moment']
    tail_count=2*segments+4*kicks
    energy_tolerance=F(1,10**8)
    dynamics_tolerance=F(1,1000)
    C=D+12*t*edges
    K=max(2,ceil_sqrt_fraction(64*C*C*M4/dynamics_tolerance**2),ceil_sqrt_fraction(2*M2))
    def cutoff(k):
        q=M2/(k+1)**2
        return ((D+2*z*t)*M2/(k+1)+6*z*t*V*M2/(k+1)**2)/(1-q)
    while cutoff(K)>energy_tolerance/2:K*=2
    hp_target=min(energy_tolerance/2,dynamics_tolerance/4)
    N=max(K+2,ceil_fraction(24*z*V*t*(K+1)*(K+2)/hp_target))
    hp=12*z*V*t*(K+1)*(K+2)/N
    cut=cutoff(K)
    assert cut+hp<=energy_tolerance
    initial_reference_error=F(1,5000)
    assert 2*(cut+hp)<=initial_reference_error**2
    assert hp<=dynamics_tolerance/4
    assert 4*C*C*M4/(K+3)**2 <= (dynamics_tolerance/4)**2
    assert tail_count**2*M4/(N+3)**4 <= (dynamics_tolerance/4)**2
    dynamics=dec(hp)+2*dec(C)*mp.sqrt(dec(M4))/(K+3)+tail_count*mp.sqrt(dec(M4))/(N+3)**2
    density=3*z*t/(N*(D-2*z*t))
    return {'spacing':str(a),'side':side,'macro_sites':V,
            'physical_box_side':str(a*side),
            'time_window_Delta_times_T':'1','momentum_window':'|p_i|<=1',
            'smallest_nonzero_physical_momentum_decimal':str(2*mp.pi/dec(a*side)),
            'D_over_Delta':str(D),'t_over_Delta':str(t),
            'minimum_free_frequency_over_Delta':'1',
            'free_vacuum_moment_rational_upper_bounds':{k:str(v) for k,v in vac.items()},
            'single_normal_mode_moment_rational_upper_bounds':{k:str(v) for k,v in one.items()},
            'complete_sequence':{'Weyl_kicks':kicks,'evolution_segments':segments,
                'total_time_Delta_times_T':'1','canonical_displacement_sum_upper':str(displacement_sum),
                'all_four_records_per_Weyl_test_retained':True,
                'displacement_amplification_R':str(R),
                'fourth_moment_upper_exact':str(M4),'projection_tail_multiplicity':tail_count,
                'uncontrolled_Halmos_output_error_factor':'sqrt(2)',
                'normalized_controlled_Hadamard_output_error_factor':'1'},
            'total_number_cutoff_K':K,'actual_constituent_blocks_per_macrocell_N':N,
            'physical_source_cells':8*N*V,
            'kinematic_source_capacity_nats_expression':str(N*V)+' * (240 log(2) + 8 log(5) + 24 log(119))',
            'kinematic_source_capacity_nats_decimal':str(mp.mpf(N*V)*(240*mp.log(2)+8*mp.log(5)+24*mp.log(119))),
            'retained_Weyl_apparatus':apparatus_capacity(kicks),
            'source_and_declared_Weyl_apparatus_capacity_nats_expression':str(N*V)+' * (240 log(2) + 8 log(5) + 24 log(119)) + log(3072) + '+str(2*kicks)+' log(2)',
            'source_and_declared_Weyl_apparatus_capacity_nats_decimal':str(mp.mpf(N*V)*(240*mp.log(2)+8*mp.log(5)+24*mp.log(119))+mp.log(3072)+2*kicks*mp.log(2)),
            'capacity_scope':'log dimension of full logical block carriers; upper bound on Gauss-physical source capacity; controllers and apparatus excluded',
            'states_were_not_enumerated_at_this_capacity':True,
            'ground_energy_error_over_Delta_upper_exact':str(cut+hp),
            'actual_symmetric_global_ground_infidelity_upper_decimal':str(dec(cut+hp)),
            'cutoff_energy_error_over_Delta_upper_exact':str(cut),
            'hardcore_Hamiltonian_error_over_Delta_upper_exact':str(hp),
            'one_excitation_evolution_vector_error_upper_exact':str(3*dynamics_tolerance/4),
            'complete_two_Weyl_three_segment_vector_error_upper_exact':str(3*dynamics_tolerance/4),
            'complete_two_Weyl_three_segment_vector_error_sharper_decimal':str(dynamics),
            'initial_actual_ground_vector_error_upper_exact':str(initial_reference_error),
            'complete_records_including_actual_ground_error_upper_exact':str(3*dynamics_tolerance/4+initial_reference_error),
            'ground_occupation_per_constituent_upper_exact':str(density),
            'targets_are_record_accuracy_not_measured_physics':True}


def gaussian_pair_coefficients(D,t,basis):
    """Actual two-site, three-color Gaussian Fock coefficients, no cutoff fit."""
    wp=math.sqrt(D*(D-2*t));wm=math.sqrt(D*(D+2*t))
    zp=(D-wp)/(D+wp);zm=(D-wm)/(D+wm)
    Z=np.zeros((6,6))
    for a in range(3):
        Z[a,a]=Z[a+3,a+3]=(zp+zm)/2
        Z[a,a+3]=Z[a+3,a]=(zp-zm)/2
    vacuum=float(np.linalg.det(np.eye(6)-Z@Z)**.25)
    @lru_cache(None)
    def coefficient(n):
        if not any(n):return vacuum
        i=next(i for i in range(6) if n[i])
        reduced=list(n);reduced[i]-=1
        value=0.
        for j in range(6):
            if reduced[j]:
                final=reduced.copy();final[j]-=1
                value+=Z[i,j]*math.sqrt(reduced[j])*coefficient(tuple(final))
        return value/math.sqrt(n[i])
    return np.array([coefficient(left+right) for left in basis for right in basis]),wp,wm


def finite_hardcore_check(N):
    # Two-site calculation checks the full three-polarization algebra.
    # Its graph degree is one. Cubic volume results use degree six above.
    D=2.;t=.125
    basis,lowers=sparse_collective(N);dim=len(basis)
    number=sparse.diags([sum(n) for n in basis],format='csr',dtype=float)
    ident=sparse.eye(dim,format='csr')
    H=D*(sparse.kron(number,ident)+sparse.kron(ident,number))
    for b in lowers:
        x=b+b.T;H-=t*sparse.kron(x,x)
    energy,vector=eigsh(H.tocsr(),k=1,which='SA',tol=1e-13,v0=np.ones(dim*dim))
    ground=vector[:,0]
    if ground[0]<0:ground=-ground
    assert min(ground)>-1e-10
    free,wp,wm=gaussian_pair_coefficients(D,t,basis)
    Efree=1.5*(wp+wm-2*D)
    fidelity=float(abs(free@ground)**2)
    error=float(energy[0]-Efree)
    assert error>=-1e-12
    assert 1-fidelity<=error/wp+1e-10
    return {'N':N,'exact_symmetric_H_dimension':dim*dim,'degree':1,
            'all_three_native_polarizations_retained':True,
            'ground_energy':float(energy[0]),'free_ground_energy':Efree,
            'ground_energy_error':error,'free_ground_overlap_squared':fidelity,
            'ground_infidelity_upper_from_actual_energy':error/wp,
            'free_Gaussian_weight_inside_local_capacity':float(free@free)}


def dispersion_controls():
    rows=[]
    for a in (1.,.5,.25,.125,.0625):
        maximum=0.;ratio=0.;covmax=0.
        for p in itertools.product(np.linspace(-1,1,9),repeat=3):
            p=np.asarray(p)
            w=math.sqrt(1+(1-a*a)/(3*a*a)*sum(1-np.cos(a*p)))
            w0=math.sqrt(1+sum(p*p)/6)
            bound=a*a*(sum(p*p)/12+sum(p**4)/144)
            error=abs(w-w0)
            assert error<=bound+1e-14
            cov=abs(np.exp(-1j*w)/(2*w)-np.exp(-1j*w0)/(2*w0))
            assert cov<=error+1e-14
            maximum=max(maximum,error);covmax=max(covmax,cov)
            if bound>0:ratio=max(ratio,error/bound)
        rows.append({'spacing':a,'box_p_max':1,'time_Delta_times_T':1,
                     'frequency_error_max_over_Delta':maximum,
                     'frequency_error_uniform_upper_over_Delta':a*a*(3/12+3/144),
                     'error_to_bound_max_ratio':ratio,'covariance_error_max_times_Delta':covmax})
    return rows


def calculate():
    mp.mp.dps=60
    parent,pins=source()
    exact=[exact_collective(N) for N in (1,2,3)]
    finite=[finite_hardcore_check(N) for N in (1,2,3,4,6)]
    resources=[controlled_resources(a,side) for a,side in ((F(1,2),16),(F(1,4),32),(F(1,8),64))]
    result={'schema':'oph.codex.collective_quadrature_refinement.v1',
            'source_gap_interval':parent['source_bounds']['block_gap_interval'],
            'construction':{'same_source_native_matrix_unit':'Pg O_alpha Pt/(60 sqrt(r))',
                            'constituent_source':'actual eight-cell H0+2^-34 T',
                            'full_physical_complement_retained':True,
                            'gauge_invariant_even_raw_L_preserving_couplers':True,
                            'clock':'exp(-it H_(a,N)) with h/a at each constituent',
                            'declared_hardware':'N replica blocks per macrocell plus quadrature exchange',
                            'measured_inputs':[],'limiting_rest_frequency_in_source_units':'actual Delta',
                            'limiting_speed_in_declared_length_clock_units':'Delta/sqrt(6)'},
            'retained_native_symmetry_reference':reference_certificate(),
            'exact_replicated_native_algebra':exact,
            'finite_full_triplet_ground_checks':finite,
            'finite_momentum_covariance_controls':dispersion_controls(),
            'explicit_source_only_capacity_choices':resources,
            'theorems':{'actual_stationary_symmetric_global_ground_exists':True,
                        'higher_block_sectors_cannot_minimize_global_energy':True,
                        'ground_energy_and_overlap_error_explicit':True,
                        'finite_time_full_hardcore_record_error_explicit':True,
                        'massive_free_KG_limit_in_increasing_capacity_subsector':True,
                        'three_components_are_internal_scalar_labels':True,
                        'finite_N_hardcore_interactions_vanish_on_controlled_domain':True,
                        'fixed_N_relativistic_limit_proved':False,
                        'nontrivial_interacting_relativistic_limit_proved':False,
                        'Lorentz_vector_or_Proca_identification_proved':False,
                        'full_finite_capacity_unique_ground_or_uniform_gap_proved':False},
            'input_sha256':pins}
    return result
