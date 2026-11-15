"""Named modular local-coarse C/E lift, distinct from hard-wall monitoring.

The integer classical event fixes the records. This coherent extension groups
indistinguishable local outcomes in one Kraus operator, a supplied instrument
choice. A transit load is treated as elsewhere, without recording transit.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import operator
from collections import deque
import numpy as np
from scipy.sparse import csr_matrix,eye,kron
from scipy.sparse.linalg import expm_multiply

PACKAGE=Path(__file__).resolve().parent
from .construction import Construction
from .quantum_records import JointRecords


def _modular_sector(construction,modulus,max_states=4096):
    """Complete connected Gauss sector of the prepared load, including disconnected graphs.

    Unoccupied remote components have fixed zero records in this connected sector.
    No uniform-fiber preparation across unreachable positions is assumed.
    """
    c=construction
    if modulus<2:raise ValueError('modulus at least two required')
    start=(0,(0,)*(2*len(c.seams)))
    states=[start];seen={start};queue=deque([start])
    def hops(state):
        pos,f=state
        for e,endpoints in enumerate(c.seams):
            for side,x in enumerate(endpoints):
                if pos==x:destination,delta=c.ports+e,1
                elif pos==c.ports+e:destination,delta=x,-1
                else:continue
                record=list(f);idx=2*e+side;record[idx]=(record[idx]+delta)%modulus
                yield (destination,tuple(record)),1 if side==0 else -1
    while queue:
        for target,_ in hops(queue.popleft()):
            if target not in seen:
                if len(states)>=max_states:raise ValueError('modular sector exceeds state budget')
                seen.add(target);states.append(target);queue.append(target)
    index={s:i for i,s in enumerate(states)}
    h=np.zeros((len(states),)*2)
    for col,state in enumerate(states):
        h[col,col]=c.mass*(1 if state[0]<c.ports else -1)
        for target,sign in hops(state):h[index[target],col]+=sign
    return tuple(states),h


class LocalCoarseModular:
    def __init__(self,construction,modulus=6):
        self.construction=construction;self.modulus=operator.index(modulus)
        self.states,self.hamiltonian=_modular_sector(construction,self.modulus)
        self.index={s:i for i,s in enumerate(self.states)}
        self.n=len(self.states);self.operators=[self._seam_kraus(e) for e in range(len(construction.seams))]

    def gauss_residue(self,state):
        pos,f=state;c=self.construction
        g=[int(pos==v) for v in range(c.ports+len(c.seams))]
        for e,(a,b) in enumerate(c.seams):
            g[a]+=f[2*e];g[b]+=f[2*e+1];g[c.ports+e]-=f[2*e]+f[2*e+1]
        return tuple(v%self.modulus for v in g)

    def integer_representative_gauss(self,state):
        """Diagnostic only: representatives do not define integer Gauss charges."""
        pos,f=state;c=self.construction
        g=[int(pos==v) for v in range(c.ports+len(c.seams))]
        for e,(a,b) in enumerate(c.seams):
            g[a]+=f[2*e];g[b]+=f[2*e+1];g[c.ports+e]-=f[2*e]+f[2*e+1]
        return tuple(g)

    def target(self,state,seam,coin):
        pos,f=state;c=self.construction
        if pos>=c.ports:return state,0,0
        a,b=c.seams[seam];contrast=int(pos==a)-int(pos==b)
        current=(contrast-(2*coin-1)*(abs(contrast)%2))//2
        destination=b if current==1 else a if current==-1 else pos
        record=list(f);record[2*seam]=(record[2*seam]+current)%self.modulus
        record[2*seam+1]=(record[2*seam+1]-current)%self.modulus
        target=(destination,tuple(record))
        if target not in self.index:raise ValueError('native event left complete modular sector')
        return target,current,contrast

    def _seam_kraus(self,seam):
        groups={}
        for law in ('C','E'):
            for coin in (0,1):
                for source,state in enumerate(self.states):
                    target,current,contrast=self.target(state,seam,coin)
                    value=current if law=='C' else contrast
                    key=(law,coin,value)
                    if key not in groups:groups[key]=np.zeros((self.n,self.n))
                    groups[key][self.index[target],source]=.5 # coin and C/E each fair
        return tuple(({'seam':seam,'law':key[0],'coin':key[1],'record_value':key[2]},matrix)
                     for key,matrix in sorted(groups.items()))

    def initial_density(self):
        rho=np.zeros((self.n,self.n),complex);rho[0,0]=1;return rho

    def density(self,rho):
        rho=np.asarray(rho,complex)
        if rho.shape!=(self.n,self.n) or not np.all(np.isfinite(rho)):raise ValueError('invalid finite density')
        if not np.allclose(rho,rho.conj().T) or not np.isclose(np.trace(rho),1) or np.linalg.eigvalsh(rho).min() < -1e-10:
            raise ValueError('normalized positive density required')
        return rho

    def commit(self,rho,seam=None,fully_monitored_control=False):
        rho=self.density(rho)
        seams=tuple(range(len(self.operators))) if seam is None else (operator.index(seam),)
        if not seams or any(not 0<=e<len(self.operators) for e in seams):raise ValueError('invalid event seam')
        source=np.diag(np.diag(rho)) if fully_monitored_control else rho
        output=np.zeros_like(rho);ledger=[]
        for e in seams:
            for record,k in self.operators[e]:
                branch=k@source@k.T/len(seams);output+=branch
                ledger.append({**record,'probability':float(np.trace(branch).real)})
        return output,ledger

    def step(self,rho,time,probability,fully_monitored_control=False):
        if not math.isfinite(time) or not math.isfinite(probability) or not 0<=probability<=1:raise ValueError('invalid split step')
        rho=self.density(rho)
        w,v=np.linalg.eigh(self.hamiltonian);u=(v*np.exp(-1j*time*w))@v.T
        coherent=u@rho@u.conj().T
        committed,ledger=self.commit(coherent,fully_monitored_control=fully_monitored_control)
        for row in ledger:row['probability']*=probability
        ledger.append({'status':'no_commit','probability':1-probability})
        return (1-probability)*coherent+probability*committed,ledger

    def per_seam_generator(self,rate):
        """Exact finite generator -i[H,rho]+r sum_e(C_e(rho)-rho).

        Every seam has the same supplied rate r. This avoids graph-size dilution
        of the local rate. No fast-commit theorem is asserted by this API.
        """
        if not math.isfinite(rate) or rate<0:raise ValueError('finite nonnegative per-seam rate required')
        identity=eye(self.n,format='csr');h=csr_matrix(self.hamiltonian)
        generator=-1j*(kron(identity,h)-kron(h.T,identity))
        channel_identity=eye(self.n*self.n,format='csr')
        for operators in self.operators:
            channel=sum((kron(csr_matrix(k),csr_matrix(k)) for _,k in operators),start=csr_matrix((self.n*self.n,)*2))
            generator+=rate*(channel-channel_identity)
        return generator.tocsr()

    def evolve_per_seam(self,rho,time,rate):
        if not math.isfinite(time) or time<0:raise ValueError('finite nonnegative duration required')
        rho=self.density(rho)
        return expm_multiply(time*self.per_seam_generator(rate),rho.reshape(-1,order='F')).reshape((self.n,self.n),order='F')

    def positions(self,rho):
        return [float(sum(rho[i,i].real for i,s in enumerate(self.states) if s[0]==pos))
                for pos in range(self.construction.ports+len(self.construction.seams))]


def run():
    c=Construction(3,((0,1),(1,2),(2,0)),.23);model=LocalCoarseModular(c,6)
    initial=model.initial_density();local=initial.copy();maximum=initial.copy()
    hard=JointRecords(c,cutoff=1);hard_rho=hard.initial_density()
    for _ in range(8):
        local,ledger=model.step(local,.37,.29)
        maximum,_=model.step(maximum,.37,.29,True)
        hard_rho,_=hard.step(hard_rho,.37,.29)
    hard_positions=[float(sum(hard_rho[i,i].real for i,s in enumerate(hard.states) if s[0]==pos)) for pos in range(6)]
    lp=np.array(model.positions(local));mp=np.array(model.positions(maximum))
    completions=[sum(k.T@k for _,k in operators) for operators in model.operators]
    gauss=model.gauss_residue(model.states[0])
    continuous=model.evolve_per_seam(initial,.37,.29/(3*.37))
    # A local C null outcome keeps coherence across other positions/records.
    null=next(k for record,k in model.operators[0] if record=={'seam':0,'law':'C','coin':1,'record_value':0})
    candidates=[i for i,s in enumerate(model.states) if s[0]!=1]
    psi=np.zeros(model.n);psi[candidates[0]]=psi[candidates[-1]]=1/np.sqrt(2)
    branch=null@np.outer(psi,psi)@null.T
    offdiagonal=float(np.sum(np.abs(branch-np.diag(np.diag(branch)))**2))
    component=LocalCoarseModular(Construction(2,((0,1),),.23),6)
    remote=LocalCoarseModular(Construction(4,((0,1),(2,3)),.23),6)
    rate=.29/(3*.37)
    local_output=component.evolve_per_seam(component.initial_density(),.37,rate)
    remote_output=remote.evolve_per_seam(remote.initial_density(),.37,rate)
    remote_difference=float(np.max(np.abs(local_output-remote_output)))
    psi_reference=np.zeros(2*component.n,complex);psi_reference[0]=psi_reference[3]=1/np.sqrt(2)
    reference_input=np.outer(psi_reference,psi_reference.conj());reference_output=np.zeros_like(reference_input)
    for _,k in component.operators[0]:
        extended=np.kron(k,np.eye(2));reference_output+=extended@reference_input@extended.T
    partial=lambda x:np.einsum('aiaj->ij',x.reshape(component.n,2,component.n,2))
    reference_difference=float(np.max(np.abs(partial(reference_input)-partial(reference_output))))
    return {'schema':'oph.modular_local_coarse.v1','name':'modular_local_coarse_one_load', 
            'scope':'Separate complete modular sector and supplied local coarse C/E lift; no integer-regulator equivalence or physical acceptance.',
            'checks':{'remote_component_rate_unchanged':remote_difference<1e-12,
                      'external_quantum_reference_preserved':reference_difference<1e-12,
                      'complete_kraus':all(np.max(np.abs(x-np.eye(model.n)))<1e-12 for x in completions),
                      'gauss_modulo_exact':all(model.gauss_residue(s)==gauss for s in model.states),
                      'trace_positive':bool(abs(np.trace(local)-1)<1e-12 and np.linalg.eigvalsh(local).min()>-1e-12),
                      'complete_local_ledger':abs(sum(r['probability'] for r in ledger)-1)<1e-12,
                      'per_seam_evolution_trace_positive':bool(abs(np.trace(continuous)-1)<1e-12 and np.linalg.eigvalsh(continuous).min()>-1e-12)},
            'mutations':{'global_monitoring_changes_evolution':float(np.max(np.abs(lp-mp)))>.01,
                         'local_null_retains_coherence':offdiagonal>1e-3,
                         'integer_gauss_transfer_rejected':any(model.integer_representative_gauss(s)!=model.integer_representative_gauss(model.states[0]) for s in model.states)},
            'full_criterion_met':False,
            'measurements':{'modular_states':model.n,'unoccupied_remote_component_density_max_difference':remote_difference,
                            'external_reference_density_max_difference':reference_difference,'local_positions':lp.tolist(),'fully_monitored_modular_positions':mp.tolist(),
                            'fully_monitored_hardwall1_positions':hard_positions,
                            'local_vs_modular_monitor_max_difference':float(np.max(np.abs(lp-mp))),
                            'local_vs_hardwall1_max_difference':float(np.max(np.abs(lp-hard_positions))),
                            'one_interval_per_seam_positions':model.positions(continuous),
                            'retained_null_branch_offdiagonal_weight':offdiagonal},
            'supplied':{'graph':'oriented triangle','modulus':6,'mass':.23,'split_time':.37,'split_probability':.29,'steps':8,
                        'local_instrument':'coherent grouping by drawn seam,coin,C/E record only',
                        'transit_rule':'elsewhere/untouched; no transit indicator emitted',
                        'continuous_rate_per_seam':.29/(3*.37),'continuous_duration':.37,
                        'modular_gauss':'residue classes only, not exact integer Gauss',
                        'component_test':'single occupied edge versus same edge plus unoccupied disconnected edge',
                        'reference_test':'entangled 3-state one-edge sector with external qubit; channel tensor identity'},
            'sources':{str(p.relative_to(PACKAGE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                       (Path(__file__).resolve(),PACKAGE/'quantum_records.py',PACKAGE/'construction.py')}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=run()
    if not all(result['checks'].values()) or not all(result['mutations'].values()):raise SystemExit('prototype gate failed')
    payload=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(payload)
    print(payload)
