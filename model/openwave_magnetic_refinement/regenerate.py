"""Regenerate numerical witnesses; verification never trusts their solver labels."""
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh
from model.openwave_intrinsic_magnetic.solver import solve,dirac
HERE=Path(__file__).resolve().parent


def save(name,**parameters):
    value=solve(relaxation=1.,max_steps=20,**parameters)
    keys=('n','spacing','field_n','coupling','mass','trap','wilson','spin_sign','external_B')
    np.savez_compressed(HERE/'witnesses'/(name+'.npz'),state=value['state'],field=value['field'],
        parameters=np.array([value[k] for k in keys]),parameter_names=np.array(keys))
    print(name,flush=True)


def main():
    (HERE/'witnesses').mkdir(exist_ok=True)
    for n in (5,9,13,17):save('stationary_n'+str(n),n=n,spacing=6/(n+1),field_n=7*(n+1)//2)
    for n in (7,9):save('volume_matter_n'+str(n),n=n,spacing=1.,field_n=31)
    for N in (13,31):save('volume_field_N'+str(N),n=5,spacing=1.,field_n=N)
    for B in (.02,-.02,.01,-.01,.005,-.005):
        save('probe_n9_'+str(B).replace('-','minus').replace('.','p'),n=9,spacing=.6,field_n=35,external_B=B)
    n=17;a=1/3;N=63
    for B,name in ((0.,'probe_zero'),(.01,'probe_B001')):
        h,_=dirac(n,a,N,np.zeros((N,)*3+(3,)),external_B=B)
        w,v=eigsh(h,k=4,sigma=8,which='LM',tol=2e-13,v0=np.linspace(1,2,h.shape[0]).astype(complex))
        order=np.argsort(w);w,v=w[order],v[:,order]
        np.savez_compressed(HERE/'witnesses'/(name+'.npz'),energies=w,vectors=v,external_B=B,n=n,spacing=a,field_n=N)
        print(name,flush=True)


if __name__=='__main__':main()
