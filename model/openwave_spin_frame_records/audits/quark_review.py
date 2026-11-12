"""Independent full-Fock pump, source current and transport-error audit."""
import hashlib,json
from pathlib import Path
import numpy as np
import sympy as s
from scipy.linalg import expm
from model.openwave_kinetic_phase_cycle.compiler import annihilator
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES
from model.openwave_quark_strings.strings import full_string,full_gauss,hop,inner,shortest_paths
ROOT=Path(__file__).resolve().parents[3]
TARGET=ROOT/'model/openwave_quark_strings'

def main():
    receipt=json.loads((TARGET/'results.json').read_text())
    assert not [p for p,h in receipt['source_sha256'].items()if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    c=[annihilator(6,j)for j in range(6)]
    B=sum((c[3+j].H*c[j]for j in range(3)),s.zeros(64))/s.sqrt(3)
    P=s.diag(*[int(mask>>3==0)for mask in range(64)])
    N=sum((c[j].H*c[j]for j in range(3)),s.zeros(64))
    K=B*P;gram=K.H*K
    assert gram==P*N*P/3
    failure=s.diag(*[s.sqrt(1-gram[j,j])for j in range(64)])
    assert gram+failure.H*failure==s.eye(64)
    assert (K*s.eye(64)[:,7]).norm()**2==1
    gauges=representation()[2];Q=gauges[10]+gauges[11]
    assert all(s.trace(g)==0 for g in gauges)
    C=s.Rational(2,3)*gauges[6]+s.Rational(4,3)*gauges[7]
    assert all((Q[i,i]+C[i,i]).is_Integer for i in range(15))
    source_checks=[]
    for family,charge in [('uc',s.Rational(2,3)),('dc',-s.Rational(1,3))]:
        assert all(Q[j,j]==-charge for j in FAMILIES[family])
        u=full_string(family,());v=hop(u)
        assert inner(v,v)==s.Rational(3,8)
        for dest in (1,2,4):
            for path in shortest_paths(dest):
                states=[full_string(family,path,spin)for spin in (0,1)]
                assert inner(states[0],states[1])==0
                assert sum(abs(inner(state,v))**2 for state in states)==s.Rational(1,16)
                for spin,state in enumerate(states):
                    assert all(not full_gauss(state,cell,current)for cell in (0,dest)for current in range(12))
        source_checks.append(family)
    rng=np.random.default_rng(10048);random_checks=0
    for dim in (4,8,12):
        for _ in range(8):
            # Non-diagonal independent retained blocks; kernel F exp(-it H0)u=0.
            H0=np.zeros((dim,dim),complex)
            for section in (slice(0,2),slice(2,dim)):
                size=H0[section,section].shape[0];X=rng.normal(size=(size,size))+1j*rng.normal(size=(size,size));H0[section,section]=(X+X.conj().T)/2
            X=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim));T=(X+X.conj().T)/2
            u=np.zeros(dim,complex);u[:2]=[1,1j];u/=np.linalg.norm(u)
            F=np.eye(dim);F[:2,:2]=0
            for t in (.001,.01,.1):
                lam=.02;actual=F@expm(-1j*t*(H0+lam*T))@u
                first=-1j*lam*t*F@T@u;b=np.linalg.norm(H0,2);tau=np.linalg.norm(T,2)
                error=lam*t*t*tau*(2*b+lam*tau)/2
                assert np.linalg.norm(actual-first)<=error+1e-15;random_checks+=1
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in TARGET.iterdir()if p.suffix in ('.py','.md','.json')}
    report={'verdict':'PASS','full_CAR_pump_dimension':64,'complete_instrument_on_all_64_rows':True,
       'native_all_fifteen_weight_centre_identity':True,'actual_two_family_six_route_spin_pair_Gauss_checks':True,
       'independent_block_nondiagonal_Duhamel_examples':random_checks,
       'scientific_boundary':'Native fractional charge, colour-carrying Wilson strings and full-H finite-time transport are executed. The central internal gauge phase does not establish fractional spatial winding.',
       'source_sha256':pins}
    (Path(__file__).resolve().parents[1]/'quark_review.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('PASS_INDEPENDENT_NATIVE_QUARK_STRINGS')
if __name__=='__main__':main()
