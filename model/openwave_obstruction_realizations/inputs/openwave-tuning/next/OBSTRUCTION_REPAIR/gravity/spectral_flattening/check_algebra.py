import json
import numpy as np
from source import ALPHA,BETA,SPIN,I4,symbol,rotation,rotation_derivatives,tilted_covariance

def norm(x):
    return float(np.linalg.norm(x))

def comm(a,b):
    return a@b-b@a

rows=[]
for p in [np.array([.1,.2,.3]),np.array([1.,.3,.4]),np.array([3.,2.,-1.])]:
    H,P,E=symbol(p)
    H0=np.einsum('j,jab->ab',p,ALPHA)
    V,Vi,Vij=rotation_derivatives(p)
    Hd=V.conj().T
    Hi=np.array([Vi[i]@H0@Hd+V@ALPHA[i]@Hd+V@H0@Vi[i].conj().T for i in range(3)])
    B=np.array([1j*V@H0@Vi[i].conj().T+.5j*V@ALPHA[i]@Hd for i in range(3)])
    Bj=np.empty((3,3,4,4),complex)  # Bj[i,j] = derivative_j B_i
    for i in range(3):
        for j in range(3):
            Bj[i,j]=1j*(Vi[j]@H0@Vi[i].conj().T+V@ALPHA[j]@Vi[i].conj().T+V@H0@Vij[i,j].conj().T)+.5j*(Vi[j]@ALPHA[i]@Hd+V@ALPHA[i]@Vi[j].conj().T)
    hk=[norm(-1j*H@Hi[i]+comm(H,B[i])+1j*p[i]*I4) for i in range(3)]
    kk=[];rot=[]
    for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:
        kk.append(norm(1j*H@(Bj[j,i]-Bj[i,j])+comm(B[i],B[j])+1j*SPIN[k]))
        rot.append(norm(-1j*(p[i]*Vi[j]-p[j]*Vi[i])+comm(SPIN[k],V)))
    h=1e-5
    derivative_error=max(norm((rotation(p+h*np.eye(3)[i])-rotation(p-h*np.eye(3)[i]))/(2*h)-Vi[i]) for i in range(3))
    second_error=max(norm((rotation_derivatives(p+h*np.eye(3)[j])[1][i]-rotation_derivatives(p-h*np.eye(3)[j])[1][i])/(2*h)-Vij[i,j]) for i in range(3) for j in range(3))
    rows.append(dict(p=p.tolist(),unitarity=norm(V@Hd-I4),hamiltonian_conjugation=norm(V@H0@Hd-H),mass_shell=norm(H@H-E*E*I4),projector=norm(P@P-P),boost_energy=max(hk),boost_boost_zero_order=max(kk),rotation=max(rot),first_derivative_control=derivative_error,second_derivative_control=second_error))

mass=[]
for m in [.01,1.]:
    samples=[]
    for r in [.01,1.,10.,100.]:
        p=np.array([r,0,0]);H,P,E=symbol(p,mass=m);V=rotation(p,mass=m)
        HD=np.einsum('j,jab->ab',p,ALPHA)+m*BETA
        delta=P-(I4-BETA)/2
        V0=rotation(p)
        wrongP=V0@((I4-HD/E)/2)@V0.conj().T
        samples.append(dict(p=r,mass_shell=norm(H@H-E*E*I4),mass_dependent_conjugation=norm(V@HD@V.conj().T-H),proper_projector_delta=norm(delta),fixed_massless_rotation_delta=norm(wrongP-(I4-BETA)/2)))
    mass.append(dict(mass=m,samples=samples))

tilts=[]
for u in [0.,.05,.3,.6]:
    gamma=1/np.sqrt(1-u*u)
    local_dirac_gram=gamma*(I4+u*ALPHA[0])
    for q in [np.array([.0001,.0002,0]),np.array([.1,.2,.3]),np.array([1.,.2,.3]),np.array([10.,2.,3.])]:
        C,G,records=tilted_covariance(q,u)
        tilts.append(dict(u=u,q=q.tolist(),gram_min=float(np.linalg.eigvalsh(G)[0]),projector_defect=norm(C@C-C),negative_trace=float(np.trace(C).real),dirac_gram_relative_defect=norm(G-local_dirac_gram)/norm(local_dirac_gram),jacobians=[float(x[1]) for x in records],number_real_physical_bands=4))

out=dict(status='ALGEBRA_AND_ACTUAL_TILTED_CAR_COMPONENT',algebra=rows,mass_extension=mass,tilted_car=tilts,claims_not_established=['geometric locality','all-observer area coefficient','joint selection','common interacting matter clock','physical G'])
with open('ALGEBRA.json','w') as f:json.dump(out,f,indent=2)
assert max(r['boost_energy'] for r in rows)<1e-11
assert max(r['boost_boost_zero_order'] for r in rows)<1e-11
assert max(r['projector_defect'] for r in tilts)<1e-11
print(json.dumps(dict(max_boost_energy=max(r['boost_energy'] for r in rows),max_boost_boost=max(r['boost_boost_zero_order'] for r in rows),max_tilt_projector=max(r['projector_defect'] for r in tilts)),indent=2))
