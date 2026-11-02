"""Exact complete spin/eta and signed-translation blocks of the source."""
from functools import lru_cache
from math import prod,gcd,comb
import numpy as np
from scipy import sparse,linalg
from ..centered_particle_spectrum import source as parent

CHARACTERS=(0,1,3,7)

def labels():
    return tuple((s,e,k) for s in range(5) for e in range(5-s) for k in CHARACTERS)


def occupations(s,e):
    if type(s) is not int or type(e) is not int or min(s,e)<0 or s+e>4:
        raise ValueError('integer spin/eta labels with s+e<=4 required')
    return 4-e+s,4-e-s


def orbit_basis(s,e,kind,k=0):
    occupations(s,e)
    if kind not in ('global','probe') or type(k) is not int or k not in range(8):
        raise ValueError('global/probe and a three-bit character required')
    return _orbit_basis(s,e,kind,k)


@lru_cache(maxsize=None)
def _orbit_basis(s,e,kind,k=0):
    up,down=occupations(s,e)
    if kind not in ('global','probe') or type(k) is not int or k not in range(8):
        raise ValueError('global/probe and a three-bit character required')
    if kind=='probe' and (s,e) not in ((1,0),(0,1)):
        raise ValueError('spin or eta vector probe required')
    gens=[((0,1,2),b,(-1)**((b&k).bit_count())) for b in (1,2,4)] if kind=='global' else [((1,0,2),0,1),((0,2,1),0,1)]
    ids=parent.masks(up,down);h=parent.centered(up,down)
    actions=[]
    for a,b,c in gens:
        rows,sign=parent.fock_action(up,down,a,b)
        r=sparse.csr_matrix((sign,(rows,np.arange(len(ids)))),shape=h.shape)
        if (r@h-h@r).nnz:raise ArithmeticError('source symmetry fails')
        actions.append((rows,sign,c))
    seen=set();cols=[];reps=[]
    for start in range(len(ids)):
        if start in seen:continue
        phase={start:1};queue=[start];ok=True
        for j in queue:
            for rows,sign,c in actions:
                i=int(rows[j]);v=int(sign[j])*c*phase[j]
                if i in phase:
                    if phase[i]!=v:ok=False
                else:phase[i]=v;queue.append(i)
        seen.update(phase)
        if ok:reps.append(start);cols.append(phase)
    rr=[];cc=[];vv=[]
    for j,col in enumerate(cols):
        for i,v in col.items():rr.append(i);cc.append(j);vv.append(v)
    q=sparse.csr_matrix((vv,(rr,cc)),shape=(len(ids),len(cols)),dtype=np.int64)
    sizes=np.asarray(q.power(2).sum(axis=0)).ravel()
    if (q.T@q-sparse.diags(sizes,dtype=np.int64)).nnz:raise ArithmeticError('orbit metric mismatch')
    return q,np.asarray(reps,dtype=np.int64),sizes


def translation_sum(up,down,k):
    if any(type(x) is not int for x in (up,down,k)) or not (0<=up<=8 and 0<=down<=8 and 0<=k<8):
        raise ValueError('integer occupations and character required')
    return _translation_sum(up,down,k)


@lru_cache(maxsize=None)
def _translation_sum(up,down,k):
    """Eight times the exact translation-character projection."""
    n=len(parent.masks(up,down));result=sparse.csr_matrix((n,n),dtype=np.int64)
    for b in range(8):
        rows,sign=parent.fock_action(up,down,(0,1,2),b)
        result+=sparse.csr_matrix(((-1)**((b&k).bit_count())*sign,(rows,np.arange(n))),shape=(n,n))
    if (result@result-8*result).nnz:raise ArithmeticError('translation character is not an orthogonal projection')
    return result


def block(s,e,k=0,kind='global',columns=None):
    occupations(s,e)
    if kind not in ('global','probe') or type(k) is not int or not 0<=k<8:
        raise ValueError('block kind and integer character required')
    if columns is not None and (type(columns) is not tuple or any(type(x) is not int for x in columns)):
        raise ValueError('exact tuple of integer columns required')
    return _block(s,e,k,kind,columns)


@lru_cache(maxsize=None)
def _block(s,e,k=0,kind='global',columns=None):
    up,down=occupations(s,e);q,reps,sizes=orbit_basis(s,e,kind,k)
    n=len(reps);eye=np.eye(n,dtype=np.int64);projectors=[];denominators=[]
    for name,m,maximum in (('spin_plus',s,4-e),('eta_minus',e,4-s)):
        roots=[j*(j+1)-m*(m+1) for j in range(m+1,maximum+1)]
        p=eye.copy()
        if roots:
            op=parent.ladder(up,down,name);ell=parent.restriction(op.T@op,q,reps).toarray()
            for root in roots:p=parent.int_product(p,root*eye-ell)
            if np.any(parent.int_product(ell,p)):raise ArithmeticError('spin/eta kernel failure')
        d=prod(roots)
        if not np.array_equal(parent.int_product(p,p),d*p):raise ArithmeticError('kernel polynomial not idempotent')
        projectors.append(p);denominators.append(d)
    if kind=='probe':
        # Total S+ kills the singlet; staggered pair annihilation is total eta-.
        forbidden=0 if s else 7
        tsum=parent.restriction(translation_sum(up,down,forbidden),q,reps).toarray()
        projectors.append(8*eye-tsum);denominators.append(8)
    p=eye.copy();denominator=1
    for factor,d in zip(projectors,denominators,strict=True):
        if not np.array_equal(parent.int_product(p,factor),parent.int_product(factor,p)):
            raise ArithmeticError('joint projectors fail to commute')
        p=parent.int_product(p,factor);denominator*=d
    # Divide a common integer factor before squaring the final projector.
    common=denominator
    for v in p.flat:common=gcd(common,int(v))
    p//=common;denominator//=common
    if not np.array_equal(parent.int_product(p,p),denominator*p):raise ArithmeticError('joint projector failure')
    if not np.array_equal(p.T*sizes[None,:],sizes[:,None]*p):raise ArithmeticError('projector metric failure')
    trace=int(np.trace(p))
    if trace%denominator:raise ArithmeticError('nonintegral projector rank')
    rank=trace//denominator
    if columns is None:
        if rank:
            _,_,piv=linalg.qr(np.sqrt(sizes)[:,None]*p.astype(float),mode='economic',pivoting=True)
            columns=tuple(int(i) for i in piv[:rank])
        else:columns=()
    if (type(columns) is not tuple or len(columns)!=rank or len(set(columns))!=rank
            or any(type(i) is not int or not 0<=i<n for i in columns)):
        raise ValueError('complete distinct integer projector columns required')
    x=p[:,columns].copy()
    for j in range(rank):
        d=0
        for v in x[:,j]:d=gcd(d,int(v))
        if not d:raise ArithmeticError('zero projector column')
        x[:,j]//=d
    if parent.modular_rank(x)!=rank:raise ArithmeticError('proposed columns do not span full projector')
    t=parent.restriction(parent.centered(up,down),q,reps).toarray()
    if not np.array_equal(parent.int_product(t,p),parent.int_product(p,t)):
        raise ArithmeticError('source leaves claimed block')
    gram=parent.int_product(x.T,sizes[:,None]*x)
    numerator=parent.int_product(x.T,sizes[:,None]*parent.int_product(t,x))
    if not np.array_equal(numerator,numerator.T):raise ArithmeticError('nonsymmetric source compression')
    return dict(s=s,e=e,k=k,kind=kind,up=up,down=down,q=q,sizes=sizes,representatives=reps,
                columns=np.asarray(columns,dtype=np.int64),basis=x,full_basis=np.asarray(q@x),
                gram=gram,numerator=numerator,dimension=rank,orbit_dimension=n,
                coefficient_source=t,projector=p,denominator=denominator)


def probe_operator(s):
    if type(s) is not int or s not in (0,1):raise ValueError('spin or pair probe required')
    return (parent.car(4,3,0,True)@parent.car(4,4,1) if s else
            parent.car(3,4,1)@parent.car(4,4,0))


def transition(s,ground_block,probe_block):
    image=probe_operator(s)@ground_block['full_basis']
    coefficients=image[probe_block['representatives']]
    if not np.array_equal(probe_block['q']@coefficients,image):raise ArithmeticError('probe image leaves exact orbit space')
    if not np.array_equal(parent.int_product(probe_block['projector'],coefficients),probe_block['denominator']*coefficients):
        raise ArithmeticError('probe has a discarded component')
    return parent.int_product(probe_block['basis'].T,probe_block['sizes'][:,None]*coefficients)


def coverage(rows):
    total=sum(comb(3,k.bit_count())*rank for s,e,k,rank in rows)
    full=sum(comb(3,k.bit_count())*(2*s+1)*rank for s,e,k,rank in rows)
    if total!=comb(8,4)**2 or full!=comb(16,8):raise ArithmeticError('neutral sector coverage failed')
    return {'balanced_dimension':total,'full_N8_dimension':full,'all_spin_eta_multiplets_included':True}


@lru_cache(maxsize=1)
def source_certificate():
    translations=[]
    for shift in range(8):
        _,_,p,s=next(r for r in parent.automorphisms() if r[:2]==((0,1,2),shift))
        mat=np.zeros((8,8),dtype=np.int64);mat[p,np.arange(8)]=s;translations.append(mat)
    for x in range(8):
        for y in range(8):
            product=translations[x]@translations[y];ref=translations[x^y]
            if not (np.array_equal(product,ref) or np.array_equal(product,-ref)):
                raise ArithmeticError('signed translations do not close modulo central sign')
    for _,shift,p,s in parent.automorphisms():
        if shift:continue
        mat=np.zeros((8,8),dtype=np.int64);mat[p,np.arange(8)]=s
        for b in range(8):
            moved=mat@translations[b]@mat.T;ref=translations[int(p[b])]
            if not (np.array_equal(moved,ref) or np.array_equal(moved,-ref)):
                raise ArithmeticError('axis permutation does not induce the claimed character orbit')
    # Full-Fock source continuity; all occupations retained, with real numerators for j/i.
    full=1<<16;h=parent.cube.operators().hubbard_H.real.tocsr()
    a=parent.ground.band_data()[0];divs=[]
    for spin_weight in ((1,1),(1,-1)):
        rr=[];cc=[];vv=[]
        for y in np.flatnonzero(a[0]):
            for spin in (0,1):
                for p,q,orientation in ((spin,2*int(y)+spin,1),(2*int(y)+spin,spin,-1)):
                    for mask in range(full):
                        one=parent.ground._annihilate(mask,q)
                        two=None if one is None else parent.ground._create(one[0],p)
                        if two is not None:
                            rr.append(two[0]);cc.append(mask)
                            vv.append(int(a[0,y])*spin_weight[spin]*orientation*one[1]*two[1])
        div=sparse.csr_matrix((vv,(rr,cc)),shape=(full,full),dtype=np.int64)
        values=np.array([spin_weight[0]*(m&1)+spin_weight[1]*((m>>1)&1) for m in range(full)],dtype=np.int64)
        q=sparse.diags(values,format='csr',dtype=np.int64)
        if (h@q-q@h+div).nnz:raise ArithmeticError('full-Fock current continuity failed')
        if (div+div.T).nnz:raise ArithmeticError('current numerator not anti-Hermitian')
        divs.append(div.nnz)
    return {'even_N_translation_group':'(Z2)^3, central signs act trivially',
            'character_orbit_representatives':[0,1,3,7],'character_orbit_sizes':[1,3,3,1],
            'full_Fock_continuity_dimension':full,'density_and_twice_spin_continuity_checked':True,
            'current_numerator_nonzero_entries':divs}
