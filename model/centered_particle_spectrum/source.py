"""Integer source, CAR ladders and complete signed cube symmetry reductions."""
from functools import lru_cache
from itertools import permutations
from math import gcd, prod
import numpy as np
from scipy import sparse, linalg
from ..cube_eta_protection import cube
from ..interacting_vacuum import ground


@lru_cache(maxsize=None)
def masks(up, down):
    if type(up) is not int or type(down) is not int or not (0 <= up <= 8 and 0 <= down <= 8):
        raise ValueError('integer up/down occupations in0..8 required')
    return np.array([int(m) for m in cube.weights(up+down)
                     if (int(m)&0x5555).bit_count()==up],dtype=np.int64)


@lru_cache(maxsize=None)
def centered(up, down):
    ids=masks(up,down);ops=cube.operators()
    if ops.scale!=16 or ops.hubbard_H.imag.nnz:
        raise ArithmeticError('expected unchanged real numerator16H')
    h=ops.hubbard_H.real[ids][:,ids].tocsr()
    return h-8*(up+down)*sparse.eye(len(ids),format='csr',dtype=np.int64)


def car(up,down,mode,creation=False):
    if type(mode) is not int or not 0<=mode<16 or type(creation) is not bool:
        raise ValueError('one source mode and Boolean creation flag required')
    target=(up+(1 if creation else -1)*(mode%2==0),
            down+(1 if creation else -1)*(mode%2==1))
    src=masks(up,down);dst=masks(*target);lookup={int(m):i for i,m in enumerate(dst)}
    rr=[];cc=[];vv=[];operation=ground._create if creation else ground._annihilate
    for j,m in enumerate(src):
        result=operation(int(m),mode)
        if result is not None:rr.append(lookup[result[0]]);cc.append(j);vv.append(result[1])
    return sparse.csr_matrix((vv,(rr,cc)),shape=(len(dst),len(src)),dtype=np.int64)


@lru_cache(maxsize=None)
def ladder(up,down,kind):
    if kind=='spin_plus':
        dst=masks(up+1,down-1)
    elif kind=='eta_minus':
        dst=masks(up-1,down-1)
    else:raise ValueError('spin_plus or eta_minus required')
    src=masks(up,down);lookup={int(m):i for i,m in enumerate(dst)}
    rr=[];cc=[];vv=[]
    for j,m in enumerate(src):
        for x,site in enumerate(cube.GEOMETRY.sites):
            if kind=='spin_plus':
                a=ground._annihilate(int(m),2*x+1)
                b=None if a is None else ground._create(a[0],2*x)
                sign=1
            else:
                a=ground._annihilate(int(m),2*x)
                b=None if a is None else ground._annihilate(a[0],2*x+1)
                sign=(-1)**sum(site)
            if b is not None:
                rr.append(lookup[b[0]]);cc.append(j);vv.append(sign*a[1]*b[1])
    return sparse.csr_matrix((vv,(rr,cc)),shape=(len(dst),len(src)),dtype=np.int64)


@lru_cache(maxsize=1)
def automorphisms():
    a=ground.band_data()[0];rows=[]
    for axes in permutations(range(3)):
        for shift in range(8):
            p=np.array([sum(((x>>axes[j])&1)<<j for j in range(3))^shift for x in range(8)])
            signs={0:1};queue=[0]
            for x in queue:
                for y in np.flatnonzero(a[x]):
                    value=signs[x]*int(a[x,y])*int(a[p[x],p[y]])
                    if y in signs:
                        if signs[y]!=value:raise ArithmeticError('cube automorphism is not a signed hopping symmetry')
                    else:signs[int(y)]=value;queue.append(int(y))
            s=np.array([signs[x] for x in range(8)],dtype=np.int64)
            if not np.array_equal(a,s[:,None]*a[np.ix_(p,p)]*s[None,:]):
                raise ArithmeticError('signed automorphism source identity failed')
            rows.append((tuple(axes),shift,p,s))
    return tuple(rows)


def _generators(kind):
    if kind not in ('ground','local'):raise ValueError('ground or local symmetry block required')
    return tuple(row for row in automorphisms()
                 if (row[0] in ((1,0,2),(0,2,1)) and row[1]==0)
                 or (kind=='ground' and row[0]==(0,1,2) and row[1] in (1,2,4)))


@lru_cache(maxsize=None)
def fock_action(up,down,axes,shift):
    _,_,p,s=next(r for r in automorphisms() if r[:2]==(axes,shift))
    ids=masks(up,down);lookup={int(m):i for i,m in enumerate(ids)}
    rr=[];vv=[]
    for m in ids:
        occupied=[j for j in range(16) if int(m)>>j&1]
        mapped=[2*int(p[j//2])+j%2 for j in occupied]
        rr.append(lookup[sum(1<<j for j in mapped)])
        sign=(-1)**sum(a>b for k,a in enumerate(mapped) for b in mapped[k+1:])
        for j in occupied:sign*=int(s[j//2])
        vv.append(sign)
    return np.asarray(rr,dtype=np.int64),np.asarray(vv,dtype=np.int64)


@lru_cache(maxsize=None)
def invariant_orbits(up,down,kind):
    """Exact disjoint signed orbit columns; inconsistent orbits have no fixed vector."""
    ids=masks(up,down);actions=[fock_action(up,down,a,b) for a,b,_,_ in _generators(kind)]
    h=centered(up,down)
    for target,sign in actions:
        r=sparse.csr_matrix((sign,(target,np.arange(len(ids)))),shape=h.shape)
        if (r@h-h@r).nnz:raise ArithmeticError('Fock symmetry does not commute with centered source')
    visited=set();columns=[];representatives=[]
    for start in range(len(ids)):
        if start in visited:continue
        phase={start:1};queue=[start];consistent=True
        for j in queue:
            for target,sign in actions:
                i=int(target[j]);v=int(sign[j])*phase[j]
                if i in phase:
                    if phase[i]!=v:consistent=False
                else:phase[i]=v;queue.append(i)
        visited.update(phase)
        if consistent:representatives.append(start);columns.append(phase)
    rr=[];cc=[];vv=[]
    for j,col in enumerate(columns):
        for i,v in col.items():rr.append(i);cc.append(j);vv.append(v)
    q=sparse.csr_matrix((vv,(rr,cc)),shape=(len(ids),len(columns)),dtype=np.int64)
    sizes=np.asarray(q.power(2).sum(axis=0)).ravel()
    if (q.T@q-sparse.diags(sizes,dtype=np.int64)).nnz:raise ArithmeticError('signed orbit columns are not disjoint')
    return q,np.asarray(representatives,dtype=np.int64),sizes


def restriction(operator,q,representatives):
    image=(operator@q).tocsr();t=image[representatives].tocsr()
    if (q@t-image).nnz:raise ArithmeticError('claimed exact reducing orbit space leaks')
    return t


def int_product(a,b):
    a=np.asarray(a,dtype=np.int64);b=np.asarray(b,dtype=np.int64)
    if a.shape[1]!=b.shape[0]:raise ValueError('matrix shapes disagree')
    maximum=lambda x:max(abs(int(x.min(initial=0))),abs(int(x.max(initial=0))))
    if a.shape[1]*maximum(a)*maximum(b)>np.iinfo(np.int64).max:
        raise OverflowError('exact integer matrix product is outside int64')
    return a@b


def modular_rank(a,prime=1009):
    a=np.asarray(a,dtype=np.int64).copy()%prime;rank=0
    for j in range(a.shape[1]):
        candidates=np.flatnonzero(a[rank:,j])
        if not len(candidates):continue
        i=rank+int(candidates[0]);a[[rank,i]]=a[[i,rank]]
        a[rank]=a[rank]*pow(int(a[rank,j]),-1,prime)%prime
        for i in range(rank+1,len(a)):
            if a[i,j]:a[i]=(a[i]-a[i,j]*a[rank])%prime
        rank+=1
        if rank==len(a):break
    return rank


@lru_cache(maxsize=6)
def reduced(kind, columns=()):
    """Exact simultaneous lowest spin and eta block, with all source rows retained."""
    up,down=(4,4) if kind=='ground' else (4,3)
    q,reps,sizes=invariant_orbits(up,down,kind)
    roots=(2,6,12,20) if kind=='ground' else (3,8,15)
    n=q.shape[1];eye=np.eye(n,dtype=np.int64);projectors=[]
    for name in ('spin_plus','eta_minus'):
        op=ladder(up,down,name);l=restriction(op.T@op,q,reps).toarray()
        p=eye.copy()
        for root in roots:p=int_product(p,root*eye-l)
        d=prod(roots)
        if not np.array_equal(int_product(p,p),d*p):raise ArithmeticError('spin/eta spectral polynomial is not idempotent')
        if np.any(int_product(l,p)):raise ArithmeticError('spin/eta projector misses the ladder kernel')
        projectors.append(p)
    denominator=prod(roots)**2
    p=int_product(*projectors)
    if not np.array_equal(p,int_product(projectors[1],projectors[0])):
        raise ArithmeticError('spin and eta projectors do not commute')
    if not np.array_equal(int_product(p,p),denominator*p):
        raise ArithmeticError('joint spin/eta projector is not idempotent')
    if not np.array_equal(p.T*sizes[None,:],sizes[:,None]*p):
        raise ArithmeticError('joint projector is not orthogonal in the orbit metric')
    trace=int(np.trace(p))
    if trace%denominator:raise ArithmeticError('joint projector trace is not integral')
    rank=trace//denominator
    # QR is an untrusted column proposer; exact modular rank and trace prove coverage.
    if not columns:
        _,_,piv=linalg.qr(np.sqrt(sizes)[:,None]*p.astype(float),pivoting=True,mode='economic')
        columns=tuple(int(i) for i in piv[:rank])
    if (type(columns) is not tuple or len(columns)!=rank or len(set(columns))!=rank
            or any(type(i) is not int or not 0<=i<n for i in columns)):
        raise ValueError('distinct exact projector column indices required')
    columns=np.asarray(columns,dtype=np.int64)
    x=p[:,columns].copy()
    for j in range(rank):
        divisor=0
        for v in x[:,j]:divisor=gcd(divisor,abs(int(v)))
        if not divisor:raise ArithmeticError('zero proposed projector column')
        x[:,j]//=divisor
    if modular_rank(x)!=rank:raise ArithmeticError('exact projector basis rank witness failed')
    full=q@x
    for name in ('spin_plus','eta_minus'):
        if np.any(ladder(up,down,name)@full):raise ArithmeticError('retained basis is not in exact ladder kernel')
    t=restriction(centered(up,down),q,reps).toarray()
    if not np.array_equal(int_product(t,p),int_product(p,t)):
        raise ArithmeticError('source does not preserve joint spin/eta block')
    gram=int_product(x.T,sizes[:,None]*x)
    numerator=int_product(x.T,sizes[:,None]*int_product(t,x))
    if not np.array_equal(numerator,numerator.T):raise ArithmeticError('compressed source is not symmetric')
    return {'up':up,'down':down,'q':q,'representatives':reps,'sizes':sizes,
            'projector':p,'denominator':denominator,'basis':x,'gram':gram,
            'numerator':numerator,'full_basis':np.asarray(full),'columns':columns,
            'orbit_dimension':n,'dimension':rank,'coefficient_source':t}


@lru_cache(maxsize=3)
def transition(ground_columns=(), local_columns=()):
    g=reduced('ground',ground_columns);c=reduced('local',local_columns)
    image=car(4,4,1)@g['full_basis']
    coefficients=image[c['representatives']]
    if not np.array_equal(c['q']@coefficients,image):
        raise ArithmeticError('actual local CAR probe leaves the fixed vertex block')
    if not np.array_equal(int_product(c['projector'],coefficients),c['denominator']*coefficients):
        raise ArithmeticError('actual local CAR probe leaves spin/eta one-half block')
    return int_product(c['basis'].T,c['sizes'][:,None]*coefficients)


@lru_cache(maxsize=1)
def shift_gram_signs():
    matrices=[]
    for shift in range(8):
        _,_,p,s=next(r for r in automorphisms() if r[:2]==((0,1,2),shift))
        a=np.zeros((8,8),dtype=np.int64);a[p,np.arange(8)]=s;matrices.append(a)
    signs=np.zeros((8,8),dtype=np.int64)
    for x in range(8):
        for y in range(8):
            value=matrices[x].T@matrices[y];ref=matrices[x^y]
            sign=int(value[x^y,0])
            if sign not in (-1,1) or not np.array_equal(value,sign*ref):
                raise ArithmeticError('projective translation Gram law failed')
            signs[x,y]=sign
    return signs


def particle_hole(up,down):
    src=masks(up,down);dst=masks(8-up,8-down);lookup={int(m):i for i,m in enumerate(dst)}
    parity=sum(1<<(2*x+s) for x,site in enumerate(cube.GEOMETRY.sites)
               for s in (0,1) if (sum(site)+s)%2)
    rows=[lookup[((1<<16)-1)^int(m)] for m in src]
    signs=[(-1)**((int(m)&parity).bit_count()) for m in src]
    return sparse.csr_matrix((signs,(rows,np.arange(len(src)))),shape=(len(dst),len(src)),dtype=np.int64)


@lru_cache(maxsize=1)
def symmetry_certificate():
    for up,down in ((4,4),(4,3)):
        k=centered(up,down)
        for kind,target in (('spin_plus',(up+1,down-1)),('eta_minus',(up-1,down-1))):
            op=ladder(up,down,kind)
            if (centered(*target)@op-op@k).nnz:
                raise ArithmeticError('exact centered spin/eta intertwiner failed')
    ug=particle_hole(4,4);uc=particle_hole(4,3)
    if (centered(4,4)@ug-ug@centered(4,4)).nnz or (centered(4,5)@uc-uc@centered(4,3)).nnz:
        raise ArithmeticError('exact centered particle-hole intertwiner failed')
    for x,site in enumerate(cube.GEOMETRY.sites):
        if (uc@car(4,4,2*x+1)-(-1)**sum(site)*car(4,4,2*x+1,True)@ug).nnz:
            raise ArithmeticError('particle-hole CAR source identity failed')
    matrices=[]
    for _,_,p,s in automorphisms():
        a=np.zeros((8,8),dtype=np.int64);a[p,np.arange(8)]=s;matrices.append(a)
    signed={a.tobytes() for a in matrices}|{(-a).tobytes() for a in matrices}
    for a in matrices:
        for b in matrices:
            if (a@b).tobytes() not in signed:raise ArithmeticError('signed cube group is not closed')
    shift_gram_signs()
    return {'signed_cube_group_order':len(signed),'normalized_spatial_lifts':len(matrices),
            'origin_stabilizer_order':6,'projective_vertex_orbit_size':8,
            'all_spatial_hopping_identities_checked':True,
            'exact_spin_eta_source_intertwiners_checked':True,
            'exact_N7_N9_particle_hole_intertwiner_checked':True,
            'all_eight_particle_hole_CAR_intertwiners_checked':True,
            'ground_spatial_character':'1; U=0 filled real negative band and Lieb-unique continuation in N8',
            'ground_hidden_spin':'0','ground_eta':'0',
            'accessible_hidden_spin':'1/2','accessible_eta':'1/2',
            'hidden_spin_identified_as_Lorentz_spin':False,
            'removal_to_addition_residue_rule':'A_E=eta B_E eta, eta_x=(-1)^(x1+x2+x3)',
            'spatial_Gram_signs':shift_gram_signs().tolist()}
