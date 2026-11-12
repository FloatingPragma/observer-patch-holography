"""Stationary positive-kernel packets and their transverse Maxwell field."""
from itertools import product
from functools import lru_cache
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import eigsh
import sympy as s
from model.native_mass_channels.algebra import folded_generators
from model.openwave_field_observables.lattice import curl,divergence_adjoint,curl_adjoint

SIGMA=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1.,-1.]).astype(complex))
BETA=np.diag([1.,1.,-1.,-1.]).astype(complex)
GAMMA=tuple(np.block([[np.zeros((2,2)),p],[p,np.zeros((2,2))]]) for p in SIGMA)
SPIN=np.diag([.5,-.5,.5,-.5]).astype(complex)


@lru_cache(maxsize=1)
def native_reduction():
    gamma,masses=folded_generators();beta=masses[0]
    sigma=-s.I*gamma[0]*gamma[1]
    p=(s.eye(16)+beta)*(s.eye(16)+sigma)/4
    basis=[]
    for column in p.columnspace():
        vector=column.copy()
        for old in basis:vector-=old*(old.H*vector)[0]
        norm=s.simplify((vector.H*vector)[0])
        if norm: basis.append(vector/s.sqrt(norm))
    assert len(basis)==4
    columns=[]
    for v in basis:columns.extend((v,gamma[0]*gamma[2]*v,gamma[2]*v,gamma[0]*v))
    u=s.Matrix.hstack(*columns)
    assert u.H*u==s.eye(16)
    for native,target in zip((*gamma,beta),(*GAMMA,BETA),strict=True):
        exact=s.Matrix([[s.Rational(float(z.real))+s.I*s.Rational(float(z.imag)) for z in row] for row in target])
        assert u.H*native*u==s.kronecker_product(s.eye(4),exact)
    return {'native_internal_rows':16,'retained_folded_valley_labels':8,
        'identical_complete_Dirac_blocks':4,'exact_unitary_intertwiner':[[str(x) for x in row] for row in u.tolist()],
        'all_signed_spectral_branches_retained_in_absolute_kernel':True}


@lru_cache(maxsize=16)
def geometry(n,spacing,field_n):
    if n<3 or n%2!=1 or field_n<n+4 or field_n%2!=1 or spacing<=0:
        raise ValueError('odd containing boxes and positive spacing required')
    coords=np.array(list(product(range(n),repeat=3)),dtype=int)
    positions=(coords-(n-1)/2)*spacing
    addresses={tuple(row):i for i,row in enumerate(coords)}
    offset=(field_n-n)//2
    edges=[]
    for i,x in enumerate(coords):
        for axis in range(3):
            y=x.copy();y[axis]+=1
            if y[axis]>=n:continue
            edges.append((i,addresses[tuple(y)],axis,tuple(x+offset)))
    return coords,positions,tuple(edges)


def dirac(n,spacing,field_n,field,mass=8.,trap=.01,wilson=.5,external_B=0.):
    coords,r,edges=geometry(n,spacing,field_n);nn=n**3
    h=sp.kron(sp.csr_matrix(BETA),sp.diags(mass+trap*np.sum(r*r,axis=1)+3*wilson/spacing),format='csr')
    phases=[]
    for axis in range(3):
        rows=[];cols=[];kinetic=[];lap=[]
        for i,j,d,grid in edges:
            if d!=axis:continue
            midpoint=(r[i]+r[j])/2
            aext=external_B*(-midpoint[1]/2 if axis==0 else midpoint[0]/2 if axis==1 else 0)
            angle=spacing*(field[grid+(axis,)]+aext)
            phase=np.exp(1j*angle);phases.append((i,j,d,grid,phase))
            rows.extend((i,j));cols.extend((j,i))
            kinetic.extend((1j*phase/(2*spacing),-1j*phase.conjugate()/(2*spacing)))
            lap.extend((-wilson*phase/(2*spacing),-wilson*phase.conjugate()/(2*spacing)))
        h+=sp.kron(sp.csr_matrix(GAMMA[axis]),sp.csr_matrix((kinetic,(rows,cols)),shape=(nn,nn)),format='csr')
        h+=sp.kron(sp.csr_matrix(BETA),sp.csr_matrix((lap,(rows,cols)),shape=(nn,nn)),format='csr')
    return h,phases


def bound_state(h,mass,spin_sign=1,previous=None,levels=6):
    values,vectors=eigsh(h,k=levels,sigma=mass,which='LM',tol=2e-13,
                         v0=np.linspace(1,2,h.shape[0]).astype(complex))
    order=np.argsort(values);values,vectors=values[order],vectors[:,order]
    positive=np.flatnonzero(values>0)
    if len(positive)<2:raise ArithmeticError('positive spectral neighborhood absent')
    # The chosen spin branch is the lowest member of the isolated ground
    # Kramers pair, or its continuously tracked split counterpart.
    pair=positive[:2]
    if abs(values[pair[1]]-values[pair[0]])<2e-10:
        v=np.linalg.qr(vectors[:,pair])[0]
        small=v.conj().T@((SPIN.diagonal()[:,None,None]*v.reshape(4,-1,2)).reshape(h.shape[0],2))
        w,u=np.linalg.eigh(small)
        psi=v@u[:,np.argmax(spin_sign*w)]
        energy=float(np.vdot(psi,h@psi).real)
    else:
        candidates=vectors[:,pair]
        candidates/=np.linalg.norm(candidates,axis=0)
        spins=np.sum(SPIN.diagonal()[:,None,None]*abs(candidates.reshape(4,-1,2))**2,axis=(0,1)).real
        index=np.argmax(spin_sign*spins)
        psi=candidates[:,index];energy=float(values[pair[index]])
    psi/=np.linalg.norm(psi)
    residual=float(np.linalg.norm(h@psi-energy*psi))
    return energy,psi,values,residual


def source_current(psi,phases,field_n,wilson=.5):
    """Expectation of -delta|D| on the selected positive-D eigenstate."""
    p=psi.reshape(4,-1);current=np.zeros((field_n,)*3+(3,))
    for i,j,axis,grid,phase in phases:
        kinetic=np.vdot(p[:,i],GAMMA[axis]@(phase*p[:,j])).real
        regulator=-wilson*np.vdot(p[:,i],BETA@(phase*p[:,j])).imag
        current[grid+(axis,)]=kinetic+regulator
    return current


def maxwell(current,spacing,coupling):
    n=current.shape[0]
    angles=2*np.pi*np.fft.fftfreq(n)
    symbols=[np.expm1(1j*angles.reshape(tuple(n if k==j else 1 for k in range(3))))/spacing for j in range(3)]
    lap=sum(abs(d)**2 for d in symbols)
    transformed=np.fft.fftn(current,axes=(0,1,2))
    field=np.zeros_like(transformed)
    for j in range(3):
        np.divide(coupling*transformed[...,j]/spacing**3,lap,out=field[...,j],where=lap>0)
    result=np.fft.ifftn(field,axes=(0,1,2)).real
    return result


def moment(current,spacing):
    n=current.shape[0]
    coordinates=np.moveaxis(np.indices((n,n,n)),0,-1)-(n-1)/2
    result=np.zeros(3)
    for axis in range(3):
        positions=spacing*coordinates.astype(float);positions[...,axis]+=spacing/2
        v=np.zeros_like(current);v[...,axis]=current[...,axis]
        result+=np.cross(positions,v).sum(axis=(0,1,2))/2
    return result


def solve(n=5,spacing=1.,field_n=19,coupling=.1,mass=8.,trap=.01,wilson=.5,spin_sign=1,external_B=0.,max_steps=40,initial_field=None,relaxation=.5):
    if coupling<=0 or mass<=0 or trap<0 or wilson<0 or not 0<relaxation<=1:
        raise ValueError('positive mass/coupling, nonnegative confinement and a valid relaxation required')
    field=np.zeros((field_n,)*3+(3,)) if initial_field is None else np.array(initial_field,dtype=float,copy=True)
    if field.shape!=(field_n,)*3+(3,):raise ValueError('complete containing field grid required')
    history=[];previous=None
    for step in range(max_steps):
        h,phases=dirac(n,spacing,field_n,field,mass,trap,wilson,external_B)
        energy,psi,values,eigen_residual=bound_state(h,mass,spin_sign,previous)
        current=source_current(psi,phases,field_n,wilson)
        target=maxwell(current,spacing,coupling)
        residual=float(np.linalg.norm(target-field));history.append(residual)
        if residual<2e-11:break
        field=(1-relaxation)*field+relaxation*target;previous=psi
    else:raise ArithmeticError('stationary iteration did not converge')
    b=curl(field,spacing)
    stationarity=curl_adjoint(b,spacing)-coupling*current/spacing**3
    spin=float(np.sum(SPIN.diagonal()[:,None]*abs(psi.reshape(4,-1))**2).real)
    mu=moment(current,spacing)
    field_energy=spacing**3*float(np.sum(b*b))/(2*coupling)
    rest_energy=energy+field_energy
    return {'n':n,'spacing':spacing,'field_n':field_n,'coupling':coupling,'mass':mass,'trap':trap,'wilson':wilson,
        'spin_sign':spin_sign,'external_B':external_B,'eigenvalue':energy,'spin':spin,'moment':mu,
        'field_energy':field_energy,'vacuum_relative_rest_energy':rest_energy,
        'g_supplied_mass':-2*mass*mu[2]/spin,'g_actual_rest_energy':-2*rest_energy*mu[2]/spin,
        'eigen_residual':eigen_residual,'fixed_point_residual':residual,'steps':len(history),'history':history,
        'current_divergence_norm':float(np.linalg.norm(divergence_adjoint(current,spacing))),
        'total_current_norm':float(np.linalg.norm(current.sum(axis=(0,1,2)))),
        'maxwell_residual':float(np.linalg.norm(stationarity)),
        'field_gauge_divergence_norm':float(np.linalg.norm(divergence_adjoint(field,spacing))),
        'current_norm':float(np.linalg.norm(current)),'current':current,'field':field,'state':psi,'dirac':h,
        'positive_spectral_neighborhood':values}


if __name__=='__main__':
    native_reduction()
    result=solve()
    print({k:v for k,v in result.items() if k not in ('current','field','state','dirac','positive_spectral_neighborhood')})
