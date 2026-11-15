"""Declared nonnegative Bose occupation lift with shared modular half-records.

Same one-load generator and local C/E instrument; no signed sea, particle
statistics derivation, or many-load uniform-fiber theorem is asserted.
"""
from collections import deque
import math
import operator
import numpy as np
from scipy.sparse import csr_matrix
from scipy.linalg import expm


class ManyLoadNative:
    def __init__(self,construction,initial_occupation,modulus=7,max_states=4096):
        self.construction=construction
        self.modulus=operator.index(modulus)
        if isinstance(modulus,bool) or self.modulus<2:raise ValueError('integer modulus >=2 required')
        c=construction;self.positions=c.ports+len(c.seams)
        occupation=tuple(operator.index(x) for x in initial_occupation)
        if len(occupation)!=self.positions or min(occupation)<0:raise ValueError('nonnegative occupation of every port/transit required')
        self.total=sum(occupation)
        initial=(occupation,(0,)*(2*len(c.seams)))
        states=[initial];index={initial:0};queue=deque([initial]);entries=[]
        while queue:
            state=queue.popleft();u,f=state;col=index[state]
            for e,ends in enumerate(c.seams):
                for side,port in enumerate(ends):
                    transit=c.ports+e;reg=2*e+side
                    for src,dst,delta in ((port,transit,1),(transit,port,-1)):
                        if not u[src]:continue
                        v=list(u);v[src]-=1;v[dst]+=1
                        g=list(f);g[reg]=(g[reg]+delta)%self.modulus
                        target=(tuple(v),tuple(g))
                        if target not in index:
                            if len(states)>=max_states:raise ValueError('complete sector exceeds explicit state budget')
                            index[target]=len(states);states.append(target);queue.append(target)
                        entries.append((index[target],col,(1 if side==0 else -1)*math.sqrt(u[src]*(u[dst]+1))))
        self.states=tuple(states);self.index=index;self.n=len(states)
        rows,cols,values=zip(*entries) if entries else ((),(),())
        h=csr_matrix((values,(rows,cols)),shape=(self.n,self.n)).toarray()
        gamma=np.r_[np.ones(c.ports),-np.ones(len(c.seams))]
        h+=np.diag([c.mass*np.dot(gamma,u) for u,_ in states])
        self.hamiltonian=h
        self.operators={law:[self._kraus(e,law) for e in range(len(c.seams))] for law in ('C','E','MIX')}

    def gauss(self,state,integer_representatives=False):
        u,f=state;c=self.construction;g=list(u)
        for e,(a,b) in enumerate(c.seams):
            g[a]+=f[2*e];g[b]+=f[2*e+1];g[c.ports+e]-=f[2*e]+f[2*e+1]
        return tuple(g if integer_representatives else (x%self.modulus for x in g))

    def target(self,state,seam,coin):
        u,f=state;a,b=self.construction.seams[seam];d=u[a]-u[b]
        j=(d-(2*coin-1)*(abs(d)%2))//2
        v=list(u);v[a]-=j;v[b]+=j
        g=list(f);g[2*seam]=(g[2*seam]+j)%self.modulus;g[2*seam+1]=(g[2*seam+1]-j)%self.modulus
        target=(tuple(v),tuple(g))
        if target not in self.index:raise ValueError('native target outside complete shared sector')
        return target,j,d

    def _kraus(self,seam,law):
        groups={};laws=('C','E') if law=='MIX' else (law,)
        for selected in laws:
            for coin in (0,1):
                for source,state in enumerate(self.states):
                    target,current,contrast=self.target(state,seam,coin)
                    key=(selected,coin,current if selected=='C' else contrast)
                    if key not in groups:groups[key]=np.zeros((self.n,self.n))
                    groups[key][self.index[target],source]=1/math.sqrt(2*len(laws))
        return tuple(({'seam':seam,'law':key[0],'coin':key[1],'record_value':key[2]},matrix) for key,matrix in sorted(groups.items()))

    def initial_density(self):
        rho=np.zeros((self.n,self.n),complex);rho[0,0]=1;return rho

    def commit(self,rho,law='MIX',seam=None):
        rho=np.asarray(rho,complex)
        if law not in self.operators:raise ValueError('C/E/MIX required')
        if rho.shape!=(self.n,self.n) or not np.isfinite(rho).all() or not np.allclose(rho,rho.conj().T) or not np.isclose(np.trace(rho),1) or np.linalg.eigvalsh(rho).min() < -1e-10:raise ValueError('normalized positive density required')
        seams=range(len(self.construction.seams)) if seam is None else (operator.index(seam),)
        seams=tuple(seams)
        if not seams or any(e<0 or e>=len(self.construction.seams) for e in seams):raise ValueError('invalid seam')
        output=np.zeros_like(rho);ledger=[]
        for e in seams:
            for record,k in self.operators[law][e]:
                branch=k@rho@k.T/len(seams);output+=branch
                value=record['record_value'];coin=record['coin']
                current=value if record['law']=='C' else (value-(2*coin-1)*(abs(value)%2))//2
                ledger.append({**record,'probability':float(np.trace(branch).real),'native_event_executed':True,'proposed_current':current,'applied_current':current,'half_record_delta':[current,-current],'boundary_rejected':False})
        return output,ledger

    def step(self,rho,time,probability,law='MIX'):
        if not math.isfinite(time) or not math.isfinite(probability) or not 0<=probability<=1:raise ValueError('finite time and probability in[0,1] required')
        u=expm(-1j*time*self.hamiltonian);coherent=u@rho@u.conj().T
        committed,ledger=self.commit(coherent,law)
        for row in ledger:row['probability']*=probability
        ledger.append({'status':'no_commit','probability':1-probability,'native_event_executed':False})
        return (1-probability)*coherent+probability*committed,ledger

    def regional_charge(self,state,region):
        region=set(region)
        if any(not isinstance(v,int) or not 0<=v<self.positions for v in region):raise ValueError('invalid port/transit region')
        return sum(state[0][v] for v in region)

    def closed_cut_record(self,state,region):
        """Oriented sum in Q_R+cutF=G_R, using stored representatives."""
        region=set(region);c=self.construction;f=state[1];value=0
        for e,(a,b) in enumerate(c.seams):
            for side,port in enumerate((a,b)):
                value+=(int(port in region)-int(c.ports+e in region))*f[2*e+side]
        return value%self.modulus

    def decode_bounded(self,residue,bound=None):
        """Unique signed integer in [-bound,bound], rejecting modular alias."""
        bound=self.total if bound is None else operator.index(bound)
        if bound<0 or self.modulus<=2*bound:raise ValueError('modulus must exceed twice declared absolute bound')
        r=operator.index(residue)%self.modulus
        candidates=[x for x in range(-bound,bound+1) if x%self.modulus==r]
        if len(candidates)!=1:raise ValueError('residue outside supplied bounded domain')
        return candidates[0]


def run():
    """Small development controls; no physical identification or criterion pass."""
    import hashlib
    from pathlib import Path
    from .construction import Construction
    from .modular_native import LocalCoarseModular
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    one=ManyLoadNative(c,(1,0,0,0,0,0),3);control=LocalCoarseModular(c,3)
    order=[one.index[(tuple(int(v==pos) for v in range(6)),f)] for pos,f in control.states]
    recovery=float(np.max(abs(one.hamiltonian[np.ix_(order,order)]-control.hamiltonian)))
    model=ManyLoadNative(Construction(3,((0,1),(0,2)),.23),(3,0,0,0,0),7)
    completeness={law:max(float(np.max(abs(sum(k.T@k for _,k in ops)-np.eye(model.n)))) for ops in model.operators[law]) for law in model.operators}
    outputs={}
    for law in model.operators:
        rho,ledger=model.step(model.initial_density(),.17,.29,law)
        outputs[law]={'position_occupations':(np.diag(rho).real@np.array([s[0] for s in model.states])).tolist(),'trace':float(np.trace(rho).real),'minimum_eigenvalue':float(np.linalg.eigvalsh(rho).min()),'complete_event_ledger':ledger}
    # The complete state domain, Bose hopping factors and complete outcomes
    # are tested against explicit changed controls, not inferred from labels.
    alias_rejected=False
    try:model.decode_bounded(0,bound=4)
    except ValueError:alias_rejected=True
    offdiag=model.hamiltonian-np.diag(np.diag(model.hamiltonian))
    wrong_hopping=float(np.max(abs(offdiag-np.sign(offdiag))))
    ops=model.operators['MIX'][0]
    omitted_error=float(np.max(abs(sum(k.T@k for _,k in ops[1:])-np.eye(model.n))))
    package=Path(__file__).resolve().parent
    return {'scope':'Fixed nonnegative Bose load with shared modular half-records; exact one-load reduction and complete C/E/MIX instruments. No signed sea, derived statistics or physical row award.',
            'full_criterion_met':False,'checks':{'one_load_generator_exact':bool(recovery==0),'C_E_MIX_complete':bool(max(completeness.values())<1e-12),'all_fixed_modular_gauss':bool(all(model.gauss(s)==model.gauss(model.states[0]) for s in model.states)),'three_unit_N7_unique_signed_decoder':bool(model.decode_bounded(3)==3 and model.decode_bounded(4)==-3),'three_unit_N3_alias':bool(3%3==0)},'mutations':{'aliased_bound_rejected':alias_rejected,'unit_hopping_without_Bose_factors_rejected':wrong_hopping>.1,'omitted_Kraus_outcome_rejected':omitted_error>.1},'broken_control_residuals':{'unit_hopping':wrong_hopping,'omitted_Kraus':omitted_error},'parameters':{'geometry':'two-edge star','mass':.23,'total_nonnegative_bose_load':3,'modulus':7,'time':.17,'commit_probability':.29,'one_load_cyclic_control_modulus':3},'supplied':['nonnegative Bose statistics','shared modular half-registers','zero definite initial registers','literal occupation C/E balancing','no many-load uniform preparation theorem'],'results':{'states':model.n,'one_load_generator_error':recovery,'kraus_completeness':completeness,'outputs':outputs},'full_physical_criterion_met':False,'sources':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__).resolve(),package/'construction.py',package/'modular_native.py')}}
