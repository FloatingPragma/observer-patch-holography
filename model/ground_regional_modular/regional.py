"""Full-output CAR traces of the actual five-state ground mixture.

Exact matrices are pairs (A,B)/d representing (A+sqrt(3) B)/d.
All regional occupation states, including the kernel, are retained.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb, factorial, lcm

import numpy as np

from ..interacting_vacuum import ground
from ..spin_resolved_spectrum import spin


def sites_checked(sites):
    if (type(sites) is not tuple or not sites or
            any(type(x) is not int or not 0 <= x < 8 for x in sites) or
            tuple(sorted(set(sites))) != sites):
        raise ValueError('a nonempty sorted tuple of distinct cube sites is required')
    return sites


def split_mask(mask, sites):
    """Move retained creation modes to the left, keeping both orders."""
    sites_checked(sites)
    kept=tuple(2*x+s for x in sites for s in (0,1))
    outside=tuple(x for x in range(16) if x not in kept)
    inside_mask=sum(((mask>>x)&1)<<j for j,x in enumerate(kept))
    exterior=sum(((mask>>x)&1)<<j for j,x in enumerate(outside))
    crossings=sum(1 for x in kept if mask&(1<<x)
                  for y in outside if y<x and mask&(1<<y))
    return inside_mask,exterior,(-1)**crossings


def density(sites):
    sites_checked(sites)
    if len(sites)>3:
        raise ValueError('the dense exact replay retains at most three sites')
    return _density(sites)


@lru_cache(maxsize=None)
def _density(sites):
    common=lcm(*spin.NORMS)
    denominator=5*1296*common
    dim=1<<(2*len(sites))
    a=np.zeros((dim,dim),dtype=object);b=np.zeros_like(a)
    for norm,state in zip(spin.NORMS,spin.ground_family(),strict=True):
        groups={}
        for mask,(x,y) in state.items():
            inside,exterior,sign=split_mask(mask,sites)
            groups.setdefault(exterior,[]).append((inside,sign*x,sign*y))
        for group in groups.values():
            for i,x,y in group:
                for j,u,v in group:
                    a[i,j]+=(common//norm)*(x*u+3*y*v)
                    b[i,j]+=(common//norm)*(x*v+y*u)
    # A common gcd reduction keeps the retained exact arrays small.
    from math import gcd
    divisor=denominator
    for value in (*a.flat,*b.flat):divisor=gcd(divisor,int(value))
    a//=divisor;b//=divisor;denominator//=divisor
    if np.trace(a)!=denominator or np.trace(b)!=0:
        raise ArithmeticError('signed regional trace lost total probability')
    if not np.array_equal(a,a.T) or not np.array_equal(b,b.T):
        raise ArithmeticError('regional trace is not self-adjoint')
    for value in (a,b):value.setflags(write=False)
    return a,b,denominator


def floating(value):
    a,b,d=value
    return (np.asarray(a,float)+np.sqrt(3)*np.asarray(b,float))/d


def creator(modes,mode):
    if type(modes) is not int or type(mode) is not int or not 0<=mode<modes:
        raise ValueError('valid local CAR mode required')
    a=np.zeros((1<<modes,1<<modes),dtype=np.int64)
    for mask in range(1<<modes):
        item=ground._create(mask,mode)
        if item is not None:a[item[0],mask]=item[1]
    return a


def pair_support():
    """Twice P: vacuum, all one-particle states, and spin-triplet N=2."""
    p=np.zeros((16,16),dtype=np.int64)
    for mask in (0,1,2,4,8,5,10):p[mask,mask]=2
    # S^-|up,up> = |up,down>+|down,up> in canonical site order.
    p[6,6]=p[9,9]=p[6,9]=p[9,6]=1
    return p


def pair_hopping(sign):
    if sign not in (-1,0,1):raise ValueError('the source bond is signed or absent')
    cs=[creator(4,j) for j in range(4)]
    return sign*sum((cs[s]@cs[2+s].T+cs[2+s]@cs[s].T for s in (0,1)),
                    np.zeros((16,16),dtype=np.int64))


def pair_certificate(sites):
    sites_checked(sites)
    if len(sites)!=2:raise ValueError('two sites required')
    a,b,d=density(sites);p2=pair_support()
    if not np.array_equal(p2@p2,2*p2):raise ArithmeticError('triplet support is not a projection')
    if not np.array_equal(p2@a,2*a) or not np.array_equal(p2@b,2*b):
        raise ArithmeticError('ground density leaks outside its claimed support')
    sign=int(ground.band_data()[0][sites[0],sites[1]])
    t=pair_hopping(sign)
    expected_a=np.zeros((16,16),dtype=object)
    expected_b=np.zeros_like(expected_a)
    # Exact density from Bernoulli natural-orbital weights, after SU(2) twirl.
    if sign:
        expected_a[0,0]=F(1,6)
        for mask in (1,2,4,8):expected_a[mask,mask]=F(1,6)
        for i in range(16):
            for j in range(16):
                if i.bit_count()==1 and j.bit_count()==1:expected_b[i,j]=F(-int(t[i,j]),12)
        for i in range(16):
            for j in range(16):
                if i.bit_count()==2 and j.bit_count()==2:expected_a[i,j]=F(int(p2[i,j]),36)
    else:
        expected_a[0,0]=F(1,4)
        for mask in (1,2,4,8):expected_a[mask,mask]=F(1,8)
        for i in range(16):
            for j in range(16):
                if i.bit_count()==2 and j.bit_count()==2:expected_a[i,j]=F(int(p2[i,j]),24)
    if any(F(int(x),d)!=y for x,y in zip(a.flat,expected_a.flat,strict=True)):
        raise ArithmeticError('actual regional rational density differs from twirl formula')
    if any(F(int(x),d)!=y for x,y in zip(b.flat,expected_b.flat,strict=True)):
        raise ArithmeticError('actual regional irrational density differs from twirl formula')
    if not np.array_equal(p2@t,t@p2):raise ArithmeticError('regional kinetic action fails to preserve support')
    if sign:
        n1=np.diag([int(x.bit_count()==1) for x in range(16)])
        if not np.array_equal(t@t,n1):
            # H also acts on the discarded N=2 singlet and N=3 blocks;
            # only the actual support is used in the modular identity.
            if not np.array_equal(p2@(t@t),2*n1):
                raise ArithmeticError('supported kinetic eigenvalues changed')
    return {'sites':list(sites),'bond_sign':sign,'full_dimension':16,'support_rank':8,
            'kernel_rank':8,'trace':'1','exact_signed_source_trace':True,
            'number_probabilities':['1/6','2/3','1/6'] if sign else ['1/4','1/2','1/4'],
            'nonzero_eigenvalues':[
                {'value':'1/6','multiplicity':1},
                {'value':'(2-sqrt3)/12','multiplicity':2},
                {'value':'(2+sqrt3)/12','multiplicity':2},
                {'value':'1/18','multiplicity':3}] if sign else [
                {'value':'1/4','multiplicity':1},{'value':'1/8','multiplicity':4},
                {'value':'1/12','multiplicity':3}],
            'supported_modular_generator':'log(6) I + 16 log(2+sqrt3) H_A + log(N_A+1)' if sign else
                'log(4) I + log(N_A+1)',
            'geometric_cap_identification':False}


def spin_twirl_weights(particles):
    if type(particles) is not int or particles<0:raise ValueError('nonnegative particle number required')
    return tuple(F(comb(particles,k)*factorial(k)*factorial(particles-k),factorial(particles+1))
                 for k in range(particles+1))


def generic_spectrum(occupations):
    """Numerical oracle for the proved regional natural-orbital formula."""
    values=np.asarray(occupations,float)
    if values.ndim!=1 or not np.isfinite(values).all() or np.any(values<0) or np.any(values>1):
        raise ValueError('natural occupations must lie in [0,1]')
    eig=[]
    for mask in range(1<<len(values)):
        k=mask.bit_count()
        weight=float(np.prod([x if mask&(1<<j) else 1-x for j,x in enumerate(values)]))
        eig.extend([weight/(k+1)]*(k+1))
    eig.extend([0.0]*((1<<(2*len(values)))-len(eig)))
    return np.sort(eig)


def certificate():
    pairs=[pair_certificate(sites) for sites in combinations(range(8),2)]
    singles=[]
    for site in range(8):
        a,b,d=density((site,))
        if not np.array_equal(4*a,d*np.diag([2,1,1,0])) or np.any(b):
            raise ArithmeticError('single-site complete ground density changed')
        singles.append({'site':site,'eigenvalues':['1/2','1/4','1/4','0']})
    if sum(bool(p['bond_sign']) for p in pairs)!=12:raise ArithmeticError('cube bond census changed')
    return {'source':'P_N4_ground/5','single_sites':singles,'all_site_pairs':pairs,
            'spin_twirl_weights':{str(k):[str(x) for x in spin_twirl_weights(k)] for k in range(5)},
            'general_eigenvalue':'w_I/(|I|+1), multiplicity |I|+1; retain the entire zero eigenspace',
            'entropy':'sum_j h_2(lambda_j) + sum_I w_I log(|I|+1)',
            'faithful_on_full_regional_CAR':False,'support_first_law':True,
            'full_CAR_first_law_inferred':False,'geometric_cap_or_Einstein_identified':False,
            'physical_particles_or_masses_identified':False}
