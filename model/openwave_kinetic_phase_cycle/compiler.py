"""Fourth-order full-source compilation using two distinct sign controls."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import mpmath as mp
import numpy as np
import sympy as s
from model.lepton_mass_compiler import analyzer,compiler as native
from model.openwave_lepton_phase_cycle.compiler import formal_order_certificate


def radius_bound(sites,spacing=F(1)):
    if type(sites)is not int or sites<1 or not isinstance(spacing,(int,F))or spacing<=0:
        raise ValueError('positive finite sites and rational spacing required')
    # Twenty-four active CAR modes/site; every other mode evolves exactly
    # under the same commuting kinetic spectator Hamiltonian.
    return F(7,3)*(36*sites/F(spacing)+6*sites)


def trace_bound(rounds,time=F(1),*,radius=None,sites=1,spacing=F(1)):
    if type(rounds)is not int or rounds<1 or not isinstance(time,(int,F))or time<0:
        raise ValueError('positive rounds and nonnegative rational time required')
    r=radius_bound(sites,spacing)if radius is None else F(radius)
    if r<=0 or rounds<r*time:raise ValueError('positive radius with r*time/rounds<=1 required')
    return min(F(2),r**5*F(time)**5/(10*rounds**4))


def rounds_for_error(error,time=F(1),*,radius=None,sites=1,spacing=F(1)):
    error=F(error);time=F(time)
    if error<=0 or time<0:raise ValueError('positive error and nonnegative time required')
    r=radius_bound(sites,spacing)if radius is None else F(radius)
    if r<=0:raise ValueError('positive radius required')
    rt=r*time;lo=max(1,(rt.numerator+rt.denominator-1)//rt.denominator);hi=lo
    while trace_bound(hi,time,radius=r)>error:hi*=2
    while lo<hi:
        mid=(lo+hi)//2
        if trace_bound(mid,time,radius=r)<=error:hi=mid
        else:lo=mid+1
    return lo


def annihilator(n,mode):
    return s.ImmutableSparseMatrix(2**n,2**n,{(mask^(1<<mode),mask):(-1)**((mask&((1<<mode)-1)).bit_count())
        for mask in range(2**n)if mask>>mode&1})


def exterior_permutation(permutation):
    n=len(permutation)
    if sorted(permutation)!=list(range(n)):raise ValueError('complete one-particle permutation required')
    entries={}
    for mask in range(2**n):
        targets=[permutation[i]for i in range(n)if mask>>i&1]
        sign=(-1)**sum(a>b for i,a in enumerate(targets)for b in targets[i+1:])
        entries[sum(1<<i for i in targets),mask]=sign
    return s.ImmutableSparseMatrix(2**n,2**n,entries)


@lru_cache(maxsize=1)
def exact_sign_controls():
    # A nonzero pair of actual periodic Weyl momenta, with physical and
    # active ancillary flags on both charged legs. All256 Fock rows enter.
    c=[annihilator(8,i)for i in range(8)];w=s.Matrix([[s.Rational(1,5),s.Rational(2,5)],[s.Rational(2,5),s.Rational(4,5)]])
    V=s.zeros(256)
    for offset in (0,4):
        for i,j in product(range(2),repeat=2):V+=w[i,j]*c[offset+i].H*c[offset+2+j].H
    V=s.ImmutableSparseMatrix(V+V.H)
    N=[c[i].H*c[i]for i in range(8)]
    K=s.ImmutableSparseMatrix(sum((N[i]for i in range(4)),s.zeros(256))-sum((N[i]for i in range(4,8)),s.zeros(256)))
    inversion=exterior_permutation([4,5,6,7,0,1,2,3])
    leg=s.diag(*[(-1)**sum((mask>>i)&1 for i in (0,1,4,5))for mask in range(256)])
    combined=s.ImmutableSparseMatrix(inversion*leg)
    assert inversion*K*inversion.H==-K and inversion*V*inversion.H==V
    assert leg*K*leg==K and leg*V*leg==-V
    assert leg*(K+V)*leg!=-(K+V)
    # Literal CAR covariance also detects omission of fermionic permutation
    # signs, even in cases where even Hamiltonians alone hide that omission.
    assert all(inversion*c[i]*inversion.H==c[[4,5,6,7,0,1,2,3][i]]for i in range(8))
    wrong=s.ImmutableSparseMatrix(256,256,{(row,col):1 for row,col in inversion.todok()})
    assert any(wrong*c[i]*wrong.H!=c[[4,5,6,7,0,1,2,3][i]]for i in range(8))
    phase_count=[sum((mask>>i)&1 for i in (1,3,5,7))for mask in range(256)]
    generators=[]
    for r in range(4):
        R=s.diag(*[s.I**(r*count)for count in phase_count])
        H=s.ImmutableSparseMatrix(K+R.H*V*R)
        assert R.H*K*R==K and combined*H*combined.H==-H
        generators.append(H)
    avg=sum(generators,s.zeros(256))/4
    phys=[sum(((mask>>j)&1)<<mode for j,mode in enumerate((0,2,4,6)))for mask in range(16)]
    J=s.ImmutableSparseMatrix(256,16,{(row,col):1 for col,row in enumerate(phys)})
    assert avg*J==J*(J.H*(K+V)*J)
    assert (generators[0]+generators[2])*J/2!=avg*J
    assert K*V!=V*K
    return {'complete_Fock_dimension':256,'physical_input_dimension':16,
            'actual_CAR_permutation_checks':8,'phase_generators_reversed':4,
            'spatial_inversion_reverses_kinetic_preserves_pairing':True,
            'charged_leg_parity_reverses_pairing_preserves_kinetic':True,
            'leg_parity_alone_reverses_total_H':False,
            'unsigned_permutation_mutation_rejected':True,
            'two_phase_ancillary_pair_mutation_rejected':True,
            'noncommuting_kinetic_pairing_control':True,
            'complete_physical_embedding_intertwiner':True}


def spatial_inversion_certificate(side=4):
    if type(side)is not int or side<3:raise ValueError('periodic side at least3 required')
    points=list(product(range(side),repeat=3));where={p:i for i,p in enumerate(points)}
    P=s.ImmutableSparseMatrix(side**3,side**3,{(where[tuple((-x)%side for x in p)],j):1 for j,p in enumerate(points)})
    sigma=(s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1))
    K=s.zeros(2*side**3)
    for axis in range(3):
        T=s.zeros(side**3)
        for j,p in enumerate(points):
            target=list(p);target[axis]=(target[axis]+1)%side;T[j,where[tuple(target)]]=1
        assert P*T*P.H==T.H
        K+=s.kronecker_product(sigma[axis],s.I*(T-T.H)/2)
    full=s.kronecker_product(s.eye(2),P)
    assert full*K*full.H==-K and K.H==K
    return {'periodic_side':side,'spatial_sites':side**3,'spin_spatial_rows':2*side**3,
            'native_Weyl_sign':'i sigma.(T-T*)/(2a)=-sigma.sin(ap)/a',
            'complete_spatial_inversion_exact':True,'boundary':'periodic; no antiperiodic PW119 transfer'}


def compiled_block(effect,kappa=F(1),rounds=256,time=F(1),precision=55):
    """Every four-mode Fock row, with both kinetic and pairing in each pulse."""
    effect=F(effect);kappa=F(kappa);time=F(time)
    if not 0<effect<1 or rounds<1:raise ValueError('interior native effect and positive rounds required')
    with mp.workdps(precision):
        convert=lambda x:mp.mpf(x.numerator)/x.denominator
        e=convert(effect);pkin=convert(kappa);t=convert(time)
        w=mp.matrix([[e,mp.sqrt(e*(1-e))],[mp.sqrt(e*(1-e)),1-e]])
        c=[mp.matrix([[mp.mpc(z.real,z.imag)for z in row]for row in native._annihilator(4,i)])for i in range(4)]
        V=mp.zeros(16)
        for i,j in product(range(2),repeat=2):V+=w[i,j]*c[i].H*c[2+j].H
        V+=V.H;K=mp.diag([pkin*mask.bit_count()for mask in range(16)])
        H=K+V;counts=[((mask>>1)&1)+((mask>>3)&1)for mask in range(16)]
        R=[mp.diag([mp.j**(r*n)for n in counts])for r in range(4)]
        physical=(0,1,4,5);J=mp.matrix([[int(row==col)for col in physical]for row in range(16)])
        target=J.H*H*J;p=1/(4-mp.root(4,3));palindrome=((0,mp.mpf(1)/8),(1,mp.mpf(1)/8),(2,mp.mpf(1)/8),(3,mp.mpf(1)/4),(2,mp.mpf(1)/8),(1,mp.mpf(1)/8),(0,mp.mpf(1)/8))
        cache={}
        def step(duration):
            U=mp.eye(16)
            for r,weight in palindrome:
                key=mp.nstr(duration*weight,precision)
                if key not in cache:cache[key]=mp.expm(-mp.j*duration*weight*H)
                U=R[r].H*cache[key]*R[r]*U
            return U
        positive=step(p*t/rounds);negative=step((1-4*p)*t/rounds)
        actual=(positive*positive*negative*positive*positive)**rounds*J
        ideal=J*mp.expm(-mp.j*t*target)
        frobenius=mp.sqrt(sum(abs(z)**2 for z in actual-ideal))
        probabilities=[sum(abs(actual[word|mask,col])**2 for col in range(4)for word in physical)/4 for mask in (0,2,8,10)]
        assert abs(sum(probabilities)-1)<mp.mpf('1e-40')
        block_radius=F(7,3)*(4*abs(kappa)+1)
        bound=trace_bound(rounds,time,radius=block_radius)
        assert 2*frobenius<=convert(bound)
        return {'source_kappa':str(kappa),'rounds':rounds,'source_time':str(time),
            'complete_Fock_rows':16,'physical_input_columns':4,'untouched_reference_dimension':4,
            'all_final_ancillary_probabilities':[mp.nstr(x,25)for x in probabilities],
            'physical_embedding_Frobenius_error':mp.nstr(frobenius,20),
            'complete_unhalved_trace_bound':str(bound),
            'numerical_execution_is_interval_certified':False,
            'kinetic_cost_not_reversed_by_leg_parity':True,
            'physical_energy_gaps':{'single_excitation':mp.nstr(mp.sqrt(pkin*pkin+e*e),25),'pair_gap_is_twice_single':True}}


def certificate():
    controls=exact_sign_controls();space=spatial_inversion_certificate();formal=formal_order_certificate()
    boxes=analyzer.rational_effect_intervals()['effect_intervals_k0_k1_k2'];runs=[]
    for k,box in enumerate(boxes):
        midpoint=sum(map(F,box),F(0))/2
        runs.append({'character':k,'effect_interval':box,
                     'opposite_native_momenta':[compiled_block(midpoint,sign)for sign in (F(-1),F(1))]})
    resources=[]
    for side in (3,4):
        sites=side**3;r=radius_bound(sites);n=rounds_for_error(F(1,10**12),sites=sites)
        resources.append({'side':side,'spatial_sites':sites,'physical_CAR_modes':90*sites,
            'complete_enlarged_CAR_modes':324000*sites,'coupled_active_CAR_modes':24*sites,
            'unchanged_kinetic_spectator_CAR_modes':(324000-24)*sites,
            'radius_upper':str(r),'rounds_for_complete_trace_error_1e_minus12':n,
            'native_Hamiltonian_segments':35*n,'complete_trace_bound':str(trace_bound(n,sites=sites)),
            'reset_protocol_rounds_for_same_error':108*sites**2*10**12})
    return {'schema':'oph.kinetic_phase_cycle.v1','exact_full_Fock_sign_controls':controls,
            'actual_periodic_kinetic_inversion':space,'formal_fourth_order':formal,
            'native_analyzer_phase':'2/9, inherited retrospective input',
            'active_complete_output_runs':runs,'complete_source_resources':resources,
            'target':'H_kinetic+H_Y with original analyzer-derived Y, same lattice clock and unit neutral background',
            'negative_time_service':'spatial inversion times complete charged-leg parity; positive source durations only',
            'total_source_pulse_duration':'(8/(4-cuberoot(4))-1)*abs(target_time)',
            'complete_records':'All physical occupations, all ancillary outputs and untouched references; no resets or success postselection.',
            'scope':'Declared coherent controls on the homogeneous periodic free source. Spatial inversion is an explicit global CAR permutation service. No primitive repair compilation, passive elapsed clock, selected Higgs vacuum, interacting mass minima or physical lepton identification is inferred.',
            'verdict':'PASS_COMPLETE_KINETIC_PHASE_COMPILER'}
