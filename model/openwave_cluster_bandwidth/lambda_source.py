"""Native uds sector with all spin/position coordinates and local strong bonds."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import numpy as np
import scipy.sparse as sp
import sympy as sy
from model.openwave_hadron_refinement_flow.positivity import inputs,declared_kappa
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_hadron_refinement_flow.certified_rest import eigen_intervals
from model.openwave_strange_composites.current import mass_parameters
from model.openwave_strange_composites.source import source as parent_source,INDEX
from model.openwave_capacity_stress.source import enclosure


@lru_cache(None)
def operators(L):
    orbs=tuple(product(range(L),range(2)));labels=tuple(product(range(2*L),repeat=3));lookup={v:i for i,v in enumerate(labels)};D=len(labels)
    positions=np.array([[orbs[j][0]for j in row]for row in labels]);spins=np.array([[orbs[j][1]for j in row]for row in labels])
    def spin(cells):
        R=sp.lil_matrix((D,D));Z=np.zeros(D)
        for j,row in enumerate(labels):
            for particle in range(3):
                x,s=orbs[row[particle]]
                if x not in cells:continue
                Z[j]+=1-2*s
                if s==1:
                    new=list(row);new[particle]=2*x;R[lookup[tuple(new)],j]+=1
        R=R.tocsr();return sp.diags(Z*(Z+2)/4)+R.T@R
    def iso(cells):
        active=np.array([[x in cells for x in row[:2]]for row in positions],int)
        diag=.75*active.sum(axis=1)-.5*np.prod(active,axis=1)
        R=sp.diags(diag,format='lil')
        for j,(u,d,s)in enumerate(labels):
            if all(active[j]):R[lookup[d,u,s],j]+=1
        return R.tocsr()
    n=np.array([(positions==x).sum(axis=1)for x in range(L)]).T
    return dict(orbs=orbs,labels=labels,lookup=lookup,D=D,positions=positions,spins=spins,n=n,spin=spin,iso=iso)


@lru_cache(None)
def source(L,flow=True):
    d=operators(L);p,_,u=inputs();k=declared_kappa();dm,_=mass_parameters();dU,dn=degree_coefficients();N=d['D'];n=d['n'];terms=[]
    terms.append((3*p['mbar']+p['delta']/2+dm,sp.eye(N)))
    terms.append((p['A']/9,sp.diags(np.sum(n*n,axis=1),dtype=float)))
    terms.append(((p['U']-k/2)*F(2,3),sp.diags(np.sum(n*(3-n),axis=1),dtype=float)))
    prefix=np.cumsum(n,axis=1)[:,:-1]
    terms.append((k*F(2,3),sp.diags(np.sum(prefix*(3-prefix),axis=1),dtype=float)))
    for x in range(L):terms.extend([(p['Js'],d['spin']((x,))),(p['Jw'],d['iso']((x,)))])
    K=sp.lil_matrix((N,N))
    for j,row in enumerate(d['labels']):
        for part in range(3):
            x,s=d['orbs'][row[part]]
            for y in(x-1,x+1):
                if 0<=y<L:
                    new=list(row);new[part]=2*y+s;K[d['lookup'][tuple(new)],j]-=1;K[j,j]+=1
    terms.append((p['t'],K.tocsr()))
    for x in range(L-1):
        terms.extend([(p['eta'],d['iso']((x,x+1))),
           (p['spin_gradient'],2*d['spin']((x,))+2*d['spin']((x+1,))-d['spin']((x,x+1)))])
    for x in range(1,L-1):
        terms.extend([(-p['eta'],d['iso']((x,))),(-p['spin_gradient'],d['spin']((x,))),(-p['t'],sp.diags(n[:,x],dtype=float))])
        if flow:terms.extend([(dn,sp.diags(n[:,x],dtype=float)),(dU*F(2,3),sp.diags(n[:,x]*(3-n[:,x]),dtype=float))])
    H=sum(float(c)*a for c,a in terms).tocsr()
    if L==2:
        old=parent_source('baryon');lookup={tuple(row):i for i,row in enumerate(old['labels'])}
        mapped=[]
        for row in d['labels']:
            values=[]
            for flavour,orb in enumerate(row):
                x,s=d['orbs'][orb];family,w=((0,0),(0,1),(1,1))[flavour];values.append(INDEX[x,family,s,w])
            mapped.append(lookup[tuple(sorted(values))])
        expected=old['H0'][mapped][:,mapped]+sp.eye(N)*float(dm)
        assert np.max(abs((H-expected).data),initial=0)<1e-9
    return {**d,'H':H,'terms':terms,'p':p,'unit':u,'kappa':k}


def rational(A):
    ints=np.rint(A*4).astype(np.int64);assert np.array_equal(ints/4,A)
    return sy.Matrix(ints.tolist())/4


@lru_cache(None)
def compact():
    source(2);d=source(4);labels=d['labels'];positions=d['positions']
    def shape(row):
        origin=min(d['orbs'][j][0]for j in row)
        return tuple((d['orbs'][j][0]-origin,d['orbs'][j][1])for j in row)
    shapes=sorted({shape(row)for row,pos in zip(labels,positions)if max(pos)-min(pos)<=1});assert len(shapes)==56
    row_ids=[next(i for i,row in enumerate(labels)if shape(row)==s and min(positions[i])==1)for s in shapes]
    col_ids=[i for i,pos in enumerate(positions)if max(pos)-min(pos)<=1]
    blocks={-1:sy.zeros(56),0:sy.zeros(56),1:sy.zeros(56)}
    for coefficient,part in d['terms']:
        sub=rational(part.tocsr()[row_ids][:,col_ids].toarray())
        for i in range(56):
            for j,col in enumerate(col_ids):
                if sub[i,j]:
                    offset=int(min(positions[col]))-1
                    assert offset in blocks
                    blocks[offset][i,shapes.index(shape(labels[col]))]+=sy.Rational(coefficient)*sub[i,j]
    S=rational(d['spin'](tuple(range(4)))[row_ids][:,row_ids].toarray())
    I=rational(d['iso'](tuple(range(4)))[row_ids][:,row_ids].toarray())
    Sz=sy.diag(*[F(3-2*sum(s for x,s in shape),2)for shape in shapes])
    selector=(S-sy.Rational(3,4)*sy.eye(56))**2+I**2+(Sz-sy.Rational(1,2)*sy.eye(56))**2
    cols=selector.nullspace();P=sy.Matrix.hstack(*cols);G=P.T*P
    compact_blocks={k:P.T*v*P for k,v in blocks.items()}
    return compact_blocks,G,d


def ground(cosine=F(1),sine=F(0)):
    blocks,G,d=compact();z=sy.Rational(cosine)+sy.I*sy.Rational(sine)
    A=blocks[0]+z*blocks[1]+sy.conjugate(z)*blocks[-1]
    assert A==A.conjugate().T
    if sine:
        R,I=A.applyfunc(sy.re),A.applyfunc(sy.im);A=R.row_join(-I).col_join(I.row_join(R));G=sy.diag(G,G)
    bounds=eigen_intervals(A,G);lo,hi=bounds[0]
    p,k,u=d['p'],d['kappa'],d['unit'];dU,dm=degree_coefficients();ds,_=mass_parameters()
    common=3*(p['mbar']+dm)+p['delta']/2+ds-3*p['t']
    split=common+5*p['A']/9+8*(p['U']+dU)/3+4*k/3+F(3,4)*(p['Js']-p['spin_gradient'])-2*max(F(0),p['eta']-p['Jw'])
    three=common+p['A']/3+4*(p['U']+dU)+2*k/3+F(9,4)*(p['Js']-p['spin_gradient'])+F(3,2)*(p['Jw']-p['eta'])
    C=min(split,three);assert C>hi
    # Full three-quark nearest-neighbour hopping norm <=6t.
    error=(6*p['t'])**2/(C-hi)
    return {'compact_rank':len(bounds),'full_Lambda_ground_GeV_interval':[str((lo-error)*u),str(hi*u)],
       'full_Lambda_ground_GeV_display':enclosure(((lo-error)*u,hi*u)),
       'complete_relative_complement_floor_source':str(C),'Schur_error_source':str(error),
       'complete_two_cell_uds_operator_matches_parent':True,'full_countable_relative_carrier_retained':True}
