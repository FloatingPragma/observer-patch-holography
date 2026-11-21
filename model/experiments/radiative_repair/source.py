"""Number-conserving sectors of a local positive-square Hamiltonian, in units of J."""
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.scalar_channel.carrier import require


def one_body(car):
    adjacency = car.adjacency()
    degree = sp.diags(car.degree,dtype=float)
    return (degree+adjacency).tocsr(), car.laplacian()


def generator(car, law, conversion_sites=None):
    require(law in ('exchange','dephased_contact','parked_neutral'), 'unknown conversion law')
    charged, neutral = one_body(car)
    if law == 'parked_neutral':
        neutral = sp.diags(car.degree.astype(float))
    n = car.n
    diagonal = np.arange(n)*(n+1)
    contact = sp.csr_matrix((np.ones(n),(diagonal,diagonal)),shape=(n*n,n*n))
    if law == 'dephased_contact':
        convert = None
    elif conversion_sites is None:
        convert = -contact
    else:
        sites = np.asarray(conversion_sites)
        require(sites.shape == (n,) and sites.dtype == bool,'invalid converter mask')
        active = diagonal[sites]
        convert = sp.csr_matrix((-np.ones(len(active)),(active,active)),shape=contact.shape)
    return sp.bmat([[sp.kronsum(charged,charged)+contact,convert],
                    [convert,sp.kronsum(neutral,neutral)+contact]],format='csr')


def initial(car, key, mode):
    from model.executable_model.carrier import EMITTER_SEAM
    a,b = EMITTER_SEAM[key]
    state = np.zeros((2,car.n,car.n),complex)
    if mode == 'ordered':
        state[0,a,b] = 1
    elif mode in ('even','odd'):
        state[0,a,b] = 1/np.sqrt(2)
        state[0,b,a] = (1 if mode == 'even' else -1)/np.sqrt(2)
    elif mode == 'neutral':
        state[1,a,a] = 1
    else:
        raise ValueError('unknown emitter')
    return state.ravel()


def evolve(H, state, waits):
    previous = 0
    for t in waits:
        if t != previous:
            step = -1j*(t-previous)
            state = expm_multiply(step*H,state,traceA=step*H.diagonal().sum())
        yield state
        previous = t


def ports_to_regions(car):
    # On one cell the apparatus resolves ports. In the federation it resolves every cell.
    labels = np.arange(car.n)
    return labels if car.cells == 1 else labels//12


def position(X,Y,labels):
    regions = int(labels.max())+1
    left, right = labels[:,None], labels[None,:]
    signed = np.bincount((regions*left+right).ravel(),weights=(abs(X)**2).ravel(),
                         minlength=regions**2).reshape(regions,regions)
    # A single sorted neutral coincidence includes both orderings for distinct endpoints.
    neutral = np.bincount((regions*np.minimum(left,right)+np.maximum(left,right)).ravel(),
                          weights=(abs(Y)**2).ravel(),minlength=regions**2).reshape(regions,regions)
    return np.r_[signed.ravel(),neutral[np.triu_indices(regions)]]


def energy_instruments(car, X, Y, law, conversion_sites=None):
    """Each row gives probabilities of occupations 0, 1, 2 for one local analyzer setting."""
    i,j = np.array(car.seams).T
    plus = .5*np.sum(abs(X[i]+X[j])**2,axis=1)
    minus = .5*np.sum(abs(X[:,i]+X[:,j])**2,axis=0)
    if law == 'parked_neutral':
        mean = 2*np.sum(abs(Y)**2,axis=1)
        two = abs(Y.diagonal())**2
        scale = car.degree
    else:
        mean = np.sum(abs(Y[i]-Y[j])**2,axis=1)
        two = abs((Y[i,i]+Y[j,j]-2*Y[i,j])/2)**2
        scale = np.full(len(i),2)
    neutral = np.column_stack([1-mean+two, mean-2*two, two])
    active = np.full(car.n,law != 'dephased_contact')
    if conversion_sites is not None:
        active &= conversion_sites
    local = np.where(active,.5*abs(X.diagonal()-Y.diagonal())**2,
                     abs(X.diagonal())**2+abs(Y.diagonal())**2)
    local_scale = np.where(active,2,1)
    def binary(p):
        return np.column_stack([1-p,p,np.zeros_like(p)])
    rows = np.vstack([binary(plus),binary(minus),neutral,binary(local)])
    coefficients = np.r_[np.full(2*len(i),2),scale,local_scale]
    return rows,coefficients


def project_cell_band():
    from model.executable_model.carrier import cell
    Q,_ = one_body(cell())
    values,vectors = np.linalg.eigh(Q.toarray())
    require(sum(abs(values-(5-np.sqrt(5))) < 1e-10) == 3, 'wrong lowest port band')
    return vectors[:,:3]


def rest_cells(car, home=0, excluded=()):
    require(0 <= home < car.cells and home not in excluded,'invalid home cell')
    near = [home]
    if car.cells > 1:
        a,b = next((a,b) for a,b in car.external if home in (a//12,b//12)
                   and a//12 not in excluded and b//12 not in excluded)
        near.append(b//12 if a//12 == home else a//12)
    return near


def rest_witness(car, home=0, excluded=()):
    """Six band amplitudes on adjacent cells satisfy at most five external constraints."""
    near = rest_cells(car,home,excluded)
    V = project_cell_band()
    frame = np.zeros((car.n,3*len(near)))
    for k,c in enumerate(near):
        frame[c*12:(c+1)*12,3*k:3*k+3] = V
    if car.external:
        i,j = np.array(car.external).T
        constraints = frame[i]+frame[j]
        small = constraints[np.any(constraints != 0,axis=1)]
        _,s,vt = np.linalg.svd(small,full_matrices=True)
        require(sum(s>1e-11) < frame.shape[1], 'external constraints exhausted the band')
        # Canonical projection of a labelled port seed; independent of the SVD basis.
        null = vt[sum(s>1e-11):]
        projector = null.T @ null
    else:
        projector = np.eye(3)
    v = frame @ (projector @ frame.T[:,near[0]*12:(near[0]+1)*12])
    col = next(i for i in range(12) if np.linalg.norm(v[:,i]) > 1e-8)
    v = v[:,col]/np.linalg.norm(v[:,col])
    Q,_ = one_body(car)
    require(np.linalg.norm(Q@v-(5-np.sqrt(5))*v)<2e-12, 'lowest-band source leaks')
    return v,near
