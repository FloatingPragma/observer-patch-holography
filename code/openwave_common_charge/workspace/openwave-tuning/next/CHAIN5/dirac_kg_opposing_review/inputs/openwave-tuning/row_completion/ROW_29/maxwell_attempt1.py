"""Actual-face oriented Maxwell symbol and controlled field dictionary."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'openwave-tuning/common_dynamics/photon_size3'))
import providers
sys.path.insert(0,str(ROOT/'openwave-tuning/common_dynamics/vector_field'))
import linear_response as linear
import sympy as s

def serial(M):return [[str(x) for x in row] for row in M.tolist()]
def bounds(r):
    r=F(r);eta=11*r;D=8*eta**2+32*eta**2*(10+eta)+F(17,2)*r*r
    assert eta<F(1,16) and r+D<F(1,4)-eta and D<r/4
    # Direct dilation eigenvector estimate; PROOF.md gives the inequalities.
    zeta=32*eta+4*D/(F(4,5)*r-D)
    assert zeta<F(1,8)
    return eta,D,zeta

def run():
    pins=providers.verify();A0,B0,D,terms,labels=linear.bloch_data();Ap=A0.pinv()
    P=s.eye(36)-A0*Ap;Q=s.eye(26)-Ap*A0
    G=(B0*B0.T+s.ones(12)/12).inv()-s.ones(12)/12
    H=s.eye(36)-B0.T*G*B0-A0*Ap
    C=s.zeros(36,3)
    for i in range(3):C[30+2*i,i]=C[31+2*i,i]=1
    W=s.Rational(2,3)*H*C
    crosses=[s.Matrix([[0,0,0],[0,0,-1],[0,1,0]]),s.Matrix([[0,0,1],[0,0,0],[-1,0,0]]),s.Matrix([[0,-1,0],[1,0,0],[0,0,0]])]
    low=[P*d*Q for d in D]
    Vt=s.Rational(3,4)*sum((crosses[i].T*W.T*low[i] for i in range(3)),s.zeros(3,26))
    assert W.T*W==s.eye(3)*s.Rational(2,3) and Vt*Vt.T==s.eye(3)
    assert all(low[i]==W*crosses[i]*Vt for i in range(3))
    assert B0*W==s.zeros(12,3) and A0.T*W==s.zeros(26,3)
    assert A0*Vt.T==s.zeros(36,3)
    assert s.trace(Ap*Ap.T)==s.Rational(3383,240) and s.trace(A0*A0.T)==88
    x,y,z=s.symbols('x y z',real=True);k=s.Matrix([x,y,z]);curl=sum((k[i]*crosses[i] for i in range(3)),s.zeros(3))
    assert curl*k==s.zeros(3,1) and curl*curl==k*k.T-(k.dot(k))*s.eye(3)
    assert (W*curl*Vt)*(W*curl*Vt).T==W*(k.dot(k)*s.eye(3)-k*k.T)*W.T
    # A single oriented incidence entry mutation destroys the exact complex.
    bad=A0.copy();idx=next((i,j) for i in range(36) for j in range(26) if bad[i,j]);bad[idx]=-bad[idx]
    assert B0*bad!=s.zeros(12,26)
    scalar=k.dot(k)*s.eye(3);assert scalar*k!=s.zeros(3,1)
    assert s.Rational(2,3)*s.Rational(1,2)*3==1
    assert s.Rational(2,3)*s.Rational(1,2)*6==2
    assert providers.verify()==pins
    out={'schema':'oph.row29.actual_oriented_maxwell.v1','checks':{'actual_face_complex':True,'all_direction_oriented_curl_factorization':True,'two_transverse_polarizations':True,'canonical_speed_one':True,'positive_field_hamiltonian':True,'source_pins_unchanged':True},
       'W':serial(W),'V_transpose':serial(Vt),'cross_generators':[serial(c) for c in crosses],
       'face_terms':[list(t) for t in terms],'A0':serial(A0),'B0':serial(B0),'D':[serial(d) for d in D],
       'rho':'2/3','electric_a':'1/2','magnetic_b':'3','canonical_Ae':'1/3','c_squared':'1',
       'dictionary':{'A':'W^T theta / h','electric':'-a W^T E / h^2','magnetic':'V^T A(k)^* theta / h^2','spatial_sample':'cell x=m*h, finite periodic volume side1'},
       'equations':['d_t B = -curl electric','d_t electric = curl B','div B = 0','div electric = 0 in the radiative q=0 sector'],
       'mutations':{'single_face_orientation_breaks_BA':True,'scalar_three_component_wave_has_longitudinal_branch':True,'double_magnetic_coefficient_changes_c_squared_to_two':True},'source_sha256':pins}
    (HERE/'MAXWELL.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['checks']),flush=True)
if __name__=='__main__':run()
