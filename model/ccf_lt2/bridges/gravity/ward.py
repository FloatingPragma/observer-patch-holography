"""Complete finite local-capacity energy partition and coarse Ward check."""
import itertools
import numpy as np
import sympy as sy
from scipy import sparse
from ...geometry import GeometryConfig, build_geometry_extension
from ...local_geometry import LocalGeometryConfig, geometry_factors
from ...local_capacity import apply_local_capacity_column


def kron(*parts):
    out=sparse.csr_matrix([[1.]])
    for part in parts: out=sparse.kron(out,part,format='csr')
    return out


def audit():
    graph={'nvertices':2,'edges':((0,1),)}
    lc=LocalGeometryConfig(phi_values=(-1.,1.))
    cc=GeometryConfig(capacity_max=1)
    factors=geometry_factors(graph,lc)
    q=factors['phi']; nphi=factors['dimension']
    I=sparse.eye(2); X=sparse.csr_matrix([[0.,1.],[1.,0.]])
    Z=sparse.diags([1.,-1.]); path=sparse.csr_matrix([[1.,-1.],[-1.,1.]])
    gradient=lc.stiffness*(q[0]-q[1])@(q[0]-q[1])/4
    gs=[lc.hopping*kron(path,I)+lc.stiffness*lc.regulator*q[0]@q[0]/2+gradient,
        lc.hopping*kron(I,path)+lc.stiffness*lc.regulator*q[1]@q[1]/2+gradient]
    hs=[kron(Z,I)/2+kron(X,X)/8,kron(I,Z)*.75+kron(X,X)/8]
    ks=[kron(gs[v],sparse.eye(4))+kron(factors['lapses'][v],hs[v]) for v in range(2)]
    primitive=build_geometry_extension(np.zeros((1,1)),config=cc)
    A,M=primitive.auxiliary_hamiltonian,primitive.lapse
    da=A.shape[0];ia=sparse.eye(da);isource=sparse.eye(nphi*4)
    K=[kron(A,ia,isource)+kron(M,ia,ks[0]),
       kron(ia,A,isource)+kron(ia,M,ks[1])]
    H=K[0]+K[1]
    actions=[lambda key,v=v:{r:sy.Rational(float(hs[v][r,key])) for r in range(4) if hs[v][r,key]!=0}
             for v in range(2)]
    error=0.
    for a,b,p,s in itertools.product(range(da),range(da),range(nphi),range(4)):
        column=apply_local_capacity_column(actions,((a,b),(p,s)),graph,lc,cc)
        vector=np.zeros(H.shape[0],complex)
        for (ab,(pr,sr)),value in column.items():
            vector[(((ab[0]*da+ab[1])*nphi+pr)*4+sr)]=complex(value)
        j=(((a*da+b)*nphi+p)*4+s)
        error=max(error,float(np.max(np.abs(vector-H[:,j].toarray().ravel()))))
    assert error==0 and all(sparse.linalg.norm(k-k.getH())==0 for k in K)
    J=[[1j*(a@b-b@a) for b in K] for a in K]
    continuity=max(float(sparse.linalg.norm(1j*(H@K[v]-K[v]@H)+sum(J[v]))) for v in range(2))
    antisymmetry=max(float(sparse.linalg.norm(J[v][w]+J[w][v])) for v in range(2) for w in range(2))
    assert continuity==antisymmetry==0
    omission=float(sparse.linalg.norm(1j*(H@K[0]-K[0]@H)))
    assert omission>0
    coarse=sum(K); internal=sum((j for row in J for j in row))
    assert sparse.linalg.norm(internal)==0 and sparse.linalg.norm(H@coarse-coarse@H)==0
    return {'all_passed':True,'complete_local_capacity_columns_checked':H.shape[0],
            'independent_partition_column_error':error,'hermitian_partition':True,
            'continuity_residual':continuity,'antisymmetry_residual':antisymmetry,
            'omitted_current_residual':omission,'block_internal_flux_cancels':True,
            'scope':'Complete 256-state local-capacity diagnostic; full-source Ward identity follows analytically from H=sum K_v, not a covariant physical stress claim.'}
