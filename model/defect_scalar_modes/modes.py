"""Complete scalar stiffness, frequency and vacuum spectral response.

The texture-dependent scalar stiffness is a declared local architecture.
The positive matrix square root defines the frequency in the occupation
cap. Its blank reference is distinct from a joint interacting ground.
"""
from functools import lru_cache
from itertools import product

import numpy as np
from scipy import linalg as la, optimize, sparse as sp
from scipy.sparse.linalg import spsolve

from model.compact_defect_pair import topology
from model.defect_mode_overlap import overlap


def graph_laplacian(points):
    points=tuple(map(tuple,points));address={p:i for i,p in enumerate(points)}
    matrix=np.zeros((len(points),len(points)))
    for p,i in address.items():
        for a in range(3):
            q=tuple(x+int(a==j) for j,x in enumerate(p))
            if q in address:
                k=address[q]
                matrix[i,i]+=1;matrix[k,k]+=1
                matrix[i,k]-=1;matrix[k,i]-=1
    return matrix


def stiffness(points, *, mass=64, half_separation=1, coupled=True):
    """Canonical q-K-q stiffness, not the scalar occupation frequency."""
    if type(mass) is not int or mass<8:
        raise ValueError('integer texture amplitude at least eight required')
    if type(coupled) is not bool:
        raise ValueError('the defect-amplitude coupling control must be explicit')
    points=tuple(map(tuple,points))
    amplitude_squared=np.array([int(any(topology.direction(p,half_separation))) for p in points],float)
    laplacian=graph_laplacian(points)
    return np.eye(len(points))+laplacian+(mass*mass*np.diag(amplitude_squared) if coupled else 0),laplacian,amplitude_squared


def _bound_profiles(k, core):
    """Resolve both core parities by full-exterior Schur inversion.

    The graph and scalar potential have an exact core-exchange reflection.
    This avoids numerical mixing of two close but nondegenerate roots.
    """
    n=len(k);outer=np.setdiff1d(np.arange(n),core)
    d=k[np.ix_(outer,outer)];b=k[np.ix_(outer,core)]
    de,du=la.eigh(d)
    result=[]
    for sign in (1,-1):
        u=np.array([1.,float(sign)])/np.sqrt(2)
        rotated=du.T@(b@u)
        def equation(energy):
            return 7.-float(np.sum(rotated**2/(de-energy)))-energy
        energy=optimize.brentq(equation,1.,7.,xtol=2e-14)
        f=np.zeros(n);f[core]=u
        f[outer]=-spsolve(sp.csr_matrix(d-energy*np.eye(len(d))),b@u)
        f/=la.norm(f)
        residual=la.norm(k@f-energy*f)
        if residual>2e-10:
            raise ArithmeticError('complete exterior scalar profile failed')
        result.append((energy,f,residual))
    if not result[0][0]<result[1][0] or np.min(result[0][1])<=0:
        raise ArithmeticError('the selected Perron mode must be strictly positive and simple')
    return result


@lru_cache(maxsize=2)
def build_scalar_source(coupled=True):
    """The frozen 343-site graph, every scalar mode and its frequency.

    K is stiffness; Omega=sqrt(K) is the scalar Hamiltonian on particle and
    antiparticle sectors. The selected profile is the lowest spectral mode.
    """
    spatial=overlap.build()
    points=tuple(spatial['raw']['points'])
    k,laplacian,amplitude=stiffness(points,coupled=coupled)
    values,vectors=la.eigh(k)
    if coupled:
        core=np.flatnonzero(amplitude==0)
        bound=_bound_profiles(k,core)
        for j,(energy,profile,_) in enumerate(bound):
            values[j],vectors[:,j]=energy,profile
    else:
        core=np.flatnonzero(amplitude==0)
        bound=[]
        vectors[:,0]=np.ones(len(points))/np.sqrt(len(points))
        values[0]=1.
    frequencies=np.sqrt(values)
    omega=(vectors*frequencies)@vectors.T
    omega=(omega+omega.T)/2
    metric=(vectors/np.sqrt(2*frequencies))@vectors.T
    profile=vectors[:,0]
    errors={'complete_mode_metric':float(la.norm(vectors.T@vectors-np.eye(len(points)))),
            'frequency_square':float(la.norm(omega@omega-k,2)),
            'vacuum_spectral_metric':float(la.norm(metric@(2*omega)@metric-np.eye(len(points)),2)),
            'lowest_mode_residual':float(la.norm(k@profile-values[0]*profile))}
    if max(errors.values())>2e-8:
        raise ArithmeticError('complete scalar spectral identities failed')
    return {'points':points,'K':k,'laplacian':laplacian,'defect_amplitude_squared':amplitude,
            'Omega':omega,'stiffness_eigenvalues':values,'frequencies':frequencies,'modes':vectors,
            'lowest_profile':profile,'lowest_projector':np.outer(profile,profile),
            'bound_pair_projector':vectors[:,:2]@vectors[:,:2].T if coupled else None,
            'scalar_field_to_canonical':metric,'core_sites':core,'spatial':spatial,
            'bound_profiles':bound,'errors':errors,'coupled':coupled}


def scalar_field(site, component, sites=343):
    """Literal global-cap scalar operator, vacuum plus both doublets."""
    if type(site) is not int or not 0<=site<sites or component not in (0,1):
        raise ValueError('valid site and doublet component required')
    particle=1+2*site+component
    antiparticle=1+2*sites+2*site+component
    return sp.csr_matrix(([1.,1.],([0,antiparticle],[particle,0])),shape=(1+4*sites,)*2)


def scalar_hamiltonian(omega):
    omega=np.asarray(omega)
    if omega.ndim!=2 or omega.shape[0]!=omega.shape[1] or not np.isfinite(omega).all() or la.norm(omega-omega.conj().T)>1e-10:
        raise ValueError('finite Hermitian scalar frequency matrix required')
    sector=sp.kron(sp.csr_matrix(omega),sp.eye(2),format='csr')
    return sp.block_diag((sp.csr_matrix((1,1)),sector,sector.conjugate()),format='csr')


def vacuum_response(z, data=None):
    """Complete retarded scalar two-point matrix on the blank reference."""
    data=build_scalar_source() if data is None else data
    if not np.isfinite(z) or np.imag(z)<=0:
        raise ValueError('upper half-plane frequency required')
    w=data['frequencies'];u=data['modes']
    return (u*(2*w/(z*z-w*w)))@u.T


def selected_response(data=None):
    data=build_scalar_source() if data is None else data
    raw=overlap.profile_response(data['lowest_profile'],value=data['spatial'])
    omega=float(data['frequencies'][0])
    return {'stiffness_eigenvalue':float(data['stiffness_eigenvalues'][0]),
            'frequency':omega,'unit_ladder_strengths':raw['singulars'].tolist(),
            'vacuum_spectral_normalized_strengths':(np.sqrt(2*omega)*raw['singulars']).tolist(),
            'strong_to_weak_ratio':float(raw['singulars'][0]/raw['singulars'][-1]),
            'complement_leakage_norm':raw['leakage'],
            'scalar_first_moment_in_selected_mode':2*omega,
            'selected_profile_minimum':float(np.min(data['lowest_profile'])),
            'core_profile_values':data['lowest_profile'][data['core_sites']].tolist()}
