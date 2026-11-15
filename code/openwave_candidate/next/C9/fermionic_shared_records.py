"""Declared fixed-number CAR source with complete shared-record native C/E.

This is a statistical/phase extension, not the exterior lift of a channel.
"""
from collections import deque
from pathlib import Path
import sys,math,json,hashlib,importlib.util
import numpy as np
from scipy.linalg import expm
BASE=Path(__file__).resolve().parents[2];FOCK=BASE/'validation/bridges/providers/fock.py'
sys.path.insert(0,str(BASE/'candidate'))
from oph_openwave.construction import Construction
from oph_openwave.modular_native import LocalCoarseModular
_spec=importlib.util.spec_from_file_location('bundled_car_fock',FOCK)
_fock=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_fock)
act=_fock.act


def transport(mask,src,dst,size):
    first=act(mask,src,False,size)
    if first is None:return None
    second=act(first[0],dst,True,size)
    return None if second is None else (second[0],first[1]*second[1])


class FermionicSharedRecords:
    def __init__(self,construction,initial_mask,modulus=3,max_states=4096):
        self.construction=construction;c=construction;self.positions=c.ports+len(c.seams)
        if type(initial_mask) is not int or not 0<=initial_mask<1<<self.positions:raise ValueError('admitted CAR bitmask required')
        if type(modulus) is not int or modulus<2:raise ValueError('integer modulus>=2 required')
        self.modulus=modulus;self.number=initial_mask.bit_count()
        start=(initial_mask,(0,)*(2*len(c.seams)));states=[start];index={start:0};queue=deque([start]);entries=[]
        while queue:
            mask,f=queue.popleft();col=index[(mask,f)]
            for e,ends in enumerate(c.seams):
                for side,x in enumerate(ends):
                    t=c.ports+e;reg=2*e+side
                    for src,dst,delta in ((x,t,1),(t,x,-1)):
                        moved=transport(mask,src,dst,self.positions)
                        if moved is None:continue
                        target_mask,phase=moved;g=list(f);g[reg]=(g[reg]+delta)%modulus
                        target=(target_mask,tuple(g))
                        if target not in index:
                            if len(states)>=max_states:raise ValueError('complete sector exceeds budget')
                            index[target]=len(states);states.append(target);queue.append(target)
                        entries.append((index[target],col,(1 if side==0 else -1)*phase))
        self.states=tuple(states);self.index=index;self.n=len(states)
        self.hamiltonian=np.zeros((self.n,self.n))
        for row,col,value in entries:self.hamiltonian[row,col]+=value
        for i,(mask,_) in enumerate(states):
            self.hamiltonian[i,i]+=c.mass*sum((1 if z<c.ports else -1)*((mask>>z)&1) for z in range(self.positions))
        self.operators={law:[self.kraus(e,law) for e in range(len(c.seams))] for law in ('C','E','MIX')}

    def gauss(self,state):
        mask,f=state;c=self.construction;g=[(mask>>z)&1 for z in range(self.positions)]
        for e,(a,b) in enumerate(c.seams):
            g[a]+=f[2*e];g[b]+=f[2*e+1];g[c.ports+e]-=f[2*e]+f[2*e+1]
        return tuple(x%self.modulus for x in g)

    def target(self,state,seam,coin):
        mask,f=state;a,b=self.construction.seams[seam]
        contrast=((mask>>a)&1)-((mask>>b)&1)
        current=(contrast-(2*coin-1)*(abs(contrast)%2))//2
        if current:
            src,dst=(a,b) if current==1 else (b,a)
            moved=transport(mask,src,dst,self.positions)
            if moved is None:raise ValueError('native move violates CAR occupation domain')
            target_mask,phase=moved
        else:target_mask,phase=mask,1
        g=list(f);g[2*seam]=(g[2*seam]+current)%self.modulus;g[2*seam+1]=(g[2*seam+1]-current)%self.modulus
        target=(target_mask,tuple(g))
        if target not in self.index:raise ValueError('native event outside complete sector')
        return target,current,contrast,phase

    def kraus(self,seam,law):
        laws=('C','E') if law=='MIX' else (law,);groups={}
        for selected in laws:
            for coin in (0,1):
                for source,state in enumerate(self.states):
                    target,current,contrast,phase=self.target(state,seam,coin)
                    key=(selected,coin,current if selected=='C' else contrast)
                    if key not in groups:groups[key]=np.zeros((self.n,self.n))
                    groups[key][self.index[target],source]=phase/math.sqrt(2*len(laws))
        return tuple(({'seam':seam,'law':key[0],'coin':key[1],'record_value':key[2]},k) for key,k in sorted(groups.items()))

    def commit(self,rho,law='MIX',seam=None):
        rho=np.asarray(rho,complex)
        if law not in self.operators or rho.shape!=(self.n,self.n) or not np.isfinite(rho).all() or not np.allclose(rho,rho.conj().T) or not np.isclose(np.trace(rho),1) or np.linalg.eigvalsh(rho).min() < -1e-10:raise ValueError('law and normalized positive finite density required')
        seams=tuple(range(len(self.construction.seams))) if seam is None else (seam,)
        if not seams or any(type(e) is not int or not 0<=e<len(self.construction.seams) for e in seams):raise ValueError('invalid seam')
        output=np.zeros_like(rho);ledger=[];h=self.hamiltonian;before=float(np.trace(h@rho).real)
        for e in seams:
            for meta,k in self.operators[law][e]:
                branch=k@rho@k.T/len(seams);output+=branch;prob=float(np.trace(branch).real)
                value=meta['record_value'];coin=meta['coin'];current=value if meta['law']=='C' else (value-(2*coin-1)*(abs(value)%2))//2
                ledger.append({**meta,'probability':prob,'applied_current':current,'proposed_current':current,'half_record_delta':[current,-current],'native_event_executed':True,'boundary_rejected':False,'unnormalized_final_energy':float(np.trace(h@branch).real)})
        return output,{'outcomes':ledger,'energy_before':before,'energy_after':float(np.trace(h@output).real),'mean_instrument_work':float(np.trace(h@(output-rho)).real)}


def run():
    c=Construction(3,((0,1),(1,2),(2,0)),.23);m=FermionicSharedRecords(c,7,3)
    eigenvalues,eigenvectors=np.linalg.eigh(m.hamiltonian)
    psi=math.sqrt(.7)*eigenvectors[:,0]+1j*math.sqrt(.3)*eigenvectors[:,-1];rho=np.outer(psi,psi.conj())
    outputs={};complete={}
    for law in m.operators:
        out,ledger=m.commit(rho,law)
        complete[law]=max(float(np.max(abs(sum(k.T@k for _,k in ops)-np.eye(m.n)))) for ops in m.operators[law])
        outputs[law]={**ledger,'trace':float(np.trace(out).real),'minimum_eigenvalue':float(np.linalg.eigvalsh(out).min())}
    one=FermionicSharedRecords(c,1,3);old=LocalCoarseModular(c,3)
    order=[one.index[(1<<pos,f)] for pos,f in old.states]
    herr=float(np.max(abs(one.hamiltonian[np.ix_(order,order)]-old.hamiltonian)))
    kerr=0.
    for e,ops in enumerate(old.operators):
        ours={tuple(meta.values()):k for meta,k in one.operators['MIX'][e]}
        for meta,k in ops:kerr=max(kerr,float(np.max(abs(ours[tuple(meta.values())][np.ix_(order,order)]-k))))
    g=m.gauss(m.states[0]);phases={m.target(s,e,coin)[3] for s in m.states for e in range(3) for coin in (0,1)}
    paths=[Path(__file__),BASE/'candidate/oph_openwave/construction.py',BASE/'candidate/oph_openwave/modular_native.py',FOCK]
    return {'schema':'oph.c9.fermionic_shared_records.v1','checks':{'Hermitian_actual_H':bool(np.array_equal(m.hamiltonian,m.hamiltonian.T)),'complete_C_E_MIX':bool(max(complete.values())<1e-12),'same_modular_Gauss':bool(all(m.gauss(s)==g for s in m.states)),'exact_one_load_H_recovery':bool(herr==0),'exact_every_one_load_MIX_K_recovery':bool(kerr==0),'nontrivial_CAR_transport_phases_retained':bool(phases=={-1,1}),'all_event_probabilities_retained':bool(all(abs(sum(row['probability'] for row in output['outcomes'])-1)<1e-12 for output in outputs.values())),'instrument_work_complete':bool(all(abs(sum(row['unnormalized_final_energy'] for row in output['outcomes'])-output['energy_after'])<1e-12 for output in outputs.values()))},'parameters':{'graph':'triangle','position_order':'ports0,1,2 then transits3,4,5','initial_mask':7,'fixed_fermion_number':3,'modulus':3,'mass':.23,'native_phase':'direct CAR c_dst^dagger c_src on moving branch; identity phase on no move','preparation':'weights .7/.3 on lowest/highest actual shared-record eigenstates with relative i'},'dimension':m.n,'spectrum':eigenvalues.tolist(),'completeness_errors':complete,'one_load_H_error':herr,'one_load_K_error':kerr,'outputs':outputs,'scope':'Complete finite shared-record fermionic native statistical extension. Not wedge of single-particle channel. Ground and instruments belong to this actual many-body-record H, not silently to unrecorded dGamma(H). No physical spin, statistics derivation, sea preparation, species or energy conservation through monitoring claimed.','full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}

if __name__=='__main__':
 r=run();Path(__file__).with_name('FERMIONIC_SHARED_RECORDS.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
