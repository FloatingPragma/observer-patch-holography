"""Binary64 virtual-photon/pair contractions; no certified spectral numbers."""
from pathlib import Path
from functools import lru_cache
import json
import numpy as np
from model.field_record_addenda.native_orthogonal_wilson_charged_band import model as band
HERE=Path(__file__).resolve().parent

def apply_one_body(M,state,constant=0):
    out={state:complex(constant)}
    for j in range(32):
        if not(state>>j)&1:continue
        removed=state^(1<<j);sgn_j=(-1)**((state&((1<<j)-1)).bit_count())
        for i in range(32):
            if (removed>>i)&1:continue
            value=M[i,j]
            final=removed|(1<<i);sgn=sgn_j*(-1)**((removed&((1<<i)-1)).bit_count())
            out[final]=out.get(final,0)+sgn*value
    return out

@lru_cache(None)
def calculate():
    band.admit();a=band.native_algebra();H,Js,_,_,Es=a.cube()
    asarray=lambda M:np.asarray(M,dtype=complex)
    H=asarray(H);e,U=np.linalg.eigh(H);negative=np.where(e<0)[0];positive=np.where(e>0)[0]
    assert len(negative)==len(positive)==16
    occ=sum(1<<int(i)for i in negative);lowest=list(map(int,positive[:2]));lam=float(np.sqrt(5)/2)
    assert max(abs(e[lowest]-lam))<2e-14
    frame,omega=band.field.frame();frame=np.asarray(frame,dtype=float);omega=np.asarray(omega,dtype=float)
    J=np.asarray([U.conj().T@asarray(x)@U for x in Js]);E=np.asarray([U.conj().T@asarray(x)@U for x in Es])
    currents=np.einsum('ek,eij->kij',frame/np.sqrt(2*omega),J)
    variance=float(band.field.cube()['Q_edge_vacuum_variance']);K2=-variance*np.sum(E,axis=0)/2
    G=np.asarray(band.field.cube()['G'],dtype=float)
    densities=[U[4*x:4*x+4].conj().T@U[4*x:4*x+4]for x in range(8)]
    def inner(left,right):return sum(np.conj(value)*right.get(key,0)for key,value in left.items())
    def sector(charge):
        starts=[occ]if charge==0 else[occ|(1<<j)for j in lowest]
        n=len(starts);base=0 if charge==0 else lam
        direct=np.empty((n,n),complex)
        vectors=[[apply_one_body(densities[x],state,-2-charge*int(x==0))for state in starts]for x in range(8)]
        coulomb=np.zeros((n,n),complex)
        for i in range(n):
            for j in range(n):
                coulomb[i,j]=sum(G[x,y]*inner(vectors[x][i],vectors[y][j])for x in range(8)for y in range(8))/2
                direct[i,j]=apply_one_body(K2,starts[j]).get(starts[i],0)+coulomb[i,j]
        virtual=[];denominators=[];mode_counts=[]
        for k in range(5):
            columns=[apply_one_body(currents[k],state)for state in starts]
            states=sorted(set().union(*(set(c)for c in columns)))
            energy=lambda state:sum(e[j]for j in range(32)if ((state>>j)&1)and not((occ>>j)&1))-sum(e[j]for j in range(32)if ((occ>>j)&1)and not((state>>j)&1))
            den=np.array([omega[k]+energy(state)-base for state in states]);assert min(den)>1.99
            virtual.extend([[col.get(state,0)for col in columns]for state in states]);denominators.extend(den)
            mode_counts.append(len(states))
        V=np.asarray(virtual,complex);den=np.asarray(denominators)
        gram=V.conj().T@V;self_energy=V.conj().T@(V/den[:,None]);resolvent_gram=V.conj().T@(V/(den[:,None]**2))
        effective=direct-self_energy
        assert max(np.max(abs(M-M.conj().T))for M in(direct,self_energy,effective,gram))<1e-12
        return {'direct':direct,'coulomb':coulomb,'kinetic':direct-coulomb,'virtual':self_energy,'effective':effective,'gram':gram,'resolvent_gram':resolvent_gram,'mode_counts':mode_counts,'smallest_virtual_denominator':float(min(den))}
    neutral,charged=sector(0),sector(1)
    difference=charged['effective']-neutral['effective'][0,0]*np.eye(2)
    encode=lambda M:[[[float(z.real),float(z.imag)]for z in row]for row in M]
    result={'schema':'oph.research.orthogonal_self_energy_diagnostic.v1','scope':'binary64 coefficient diagnostics only, not interval certification or full-source diagonalization','source_charge_reference_cost':4,'g':'1/64','free_reference_corrected_transition_over_m':lam,'one_body_eigensystem_residual':float(np.linalg.norm(H@U-U*e,2))}
    for name,value in [('neutral',neutral),('charged',charged)]:
        result[name]={key:encode(v)if isinstance(v,np.ndarray)else v for key,v in value.items()}
        result[name]['virtual_Gram_largest_eigenvalue']=float(np.linalg.eigvalsh(value['gram'])[-1])
    result['transition_second_order_matrix']=encode(difference)
    result['transition_second_order_eigenvalues']=np.linalg.eigvalsh(difference).tolist()
    result['transition_through_second_order_over_m']=(lam+np.linalg.eigvalsh(difference)/4096).tolist()
    return result

if __name__=='__main__':
    value=calculate();(HERE/'diagnostic.json').write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in value.items()if k.startswith('transition_')},indent=2))
