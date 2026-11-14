"""Operational readouts and sparse trajectories of the finite record instrument.

The path window is a declared projection. It supplies no error bound for the
infinite register problem. Trajectories sample the same source-monitored C/E
instrument as JointRecords.commit; density evolution is used for small controls.
"""
from collections import deque
from dataclasses import dataclass
import math
import operator
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
from .construction import Construction
from .quantum_records import JointRecords


def sparse_hamiltonian(joint):
    rows,cols,values=[],[],[]
    for col,state in enumerate(joint.states):
        rows.append(col); cols.append(col)
        values.append(joint.construction.mass*(1 if state[0]<joint.construction.ports else -1))
        for target,sign in joint.half_hops(state):
            if target in joint.index:
                rows.append(joint.index[target]); cols.append(col); values.append(sign)
    return csr_matrix((values,(rows,cols)),shape=(len(joint.states),)*2)


@dataclass
class RecordWindow(JointRecords):
    """Radius ball in the half-hop state graph, with a Hermitian projected H."""
    radius: int = 3

    def __post_init__(self):
        try:
            self.radius=operator.index(self.radius)
            self.initial_port=operator.index(self.initial_port)
            self.max_states=operator.index(self.max_states)
        except TypeError as exc:
            raise ValueError('integer radius, port and state budget required') from exc
        if self.radius<0 or type(self.cutoff) is not int or self.cutoff<1 or self.max_states<1:
            raise ValueError('invalid finite window')
        if not 0<=self.initial_port<self.construction.ports:
            raise ValueError('invalid initial port')
        start=(self.initial_port,(0,)*(2*len(self.construction.seams)))
        states=[start]; seen={start}; queue=deque([(start,0)])
        while queue:
            state,depth=queue.popleft()
            if depth==self.radius: continue
            for target,_ in self.half_hops(state):
                if target not in seen:
                    if len(states)>=self.max_states: raise ValueError('window exceeds state budget')
                    states.append(target); seen.add(target); queue.append((target,depth+1))
        self.states=tuple(states); self.index={s:i for i,s in enumerate(states)}

    def initial_state_error_bound(self,time):
        """Taylor/Duhamel tail bound against the complete cutoff sector.

        For the initial basis vector, powers H^k agree through k=radius.
        Both operator norms are <= sqrt(m^2+2*dmax), since each incidence
        column has two unit entries and each port has at most dmax entries.
        Two Taylor tails give 2 exp(x) x^(radius+1)/(radius+1)!.
        The bound concerns coherent evolution of this preparation only.
        """
        if not math.isfinite(time): raise ValueError('finite time required')
        degree=[0]*self.construction.ports
        for a,b in self.construction.seams: degree[a]+=1; degree[b]+=1
        norm=math.sqrt(self.construction.mass**2+2*max(degree,default=0))
        x=abs(time)*norm
        order=self.radius+1
        if x==0:return 0.
        logbound=math.log(2)+x+order*math.log(x)-math.lgamma(order+1)
        return float(min(2.,math.exp(min(math.log(2),logbound))))

    @property
    def hamiltonian(self):
        return sparse_hamiltonian(self).toarray()


class ObserverProtocol:
    def __init__(self,joint):
        self.joint=joint
        self.h=sparse_hamiltonian(joint)

    def probabilities(self,state):
        state=np.asarray(state,complex)
        n=len(self.joint.states)
        if state.shape==(n,):
            if not np.all(np.isfinite(state)) or not np.isclose(np.vdot(state,state),1):
                raise ValueError('normalized finite pure state required')
            return np.abs(state)**2
        return np.diag(self.joint._density(state)).real

    def readout(self,state):
        """Complete position distribution and unconditional event probabilities.

        Event probabilities assume a monitored attempt: seam/coin/law uniform.
        They are not probabilities of coherent hard-wall rejection.
        """
        weights=self.probabilities(state)
        j=self.joint; c=j.construction
        positions=np.zeros(c.ports+len(c.seams))
        outcomes={'native':0.,'boundary_rejected':0.,'transit_noop':0.}
        currents={}
        if not c.seams: raise ValueError('instrument requires at least one seam')
        for index,(pos,records) in enumerate(j.states):
            positions[pos]+=weights[index]
            for e in range(len(c.seams)):
                for q in (0,1):
                    _,status,current=j.native_target((pos,records),e,q)
                    weight=float(weights[index])/(2*len(c.seams))
                    outcomes[status]+=weight
                    applied=current if status=='native' else 0
                    currents[str(applied)]=currents.get(str(applied),0.)+weight
        return {'position_probabilities':positions.tolist(),
                'source_port_probability':float(positions[j.initial_port]),
                'port_probability':float(sum(positions[:c.ports])),
                'transit_probability':float(sum(positions[c.ports:])),
                'monitored_attempt_outcomes':outcomes,
                'applied_current_probabilities':currents,
                'scope':'position measurement and conditional monitored-attempt distribution; no coherent boundary-loss bound'}

    def evolve(self,state,time):
        state=np.asarray(state,complex)
        self.probabilities(state)
        if state.ndim!=1 or not math.isfinite(time):
            raise ValueError('pure state and finite time required')
        return expm_multiply((-1j*time)*self.h,state)

    def trajectory_step(self,state,time,commit_probability,rng):
        """Sample the exact declared finite instrument; no density matrix stored."""
        if not math.isfinite(commit_probability) or not 0<=commit_probability<=1:
            raise ValueError('invalid commit probability')
        evolved=self.evolve(state,time)
        if rng.random()>=commit_probability:
            return evolved,{'status':'no_commit'}
        j=self.joint; c=j.construction
        if not c.seams: raise ValueError('instrument requires at least one seam')
        weights=np.abs(evolved)**2
        source=int(rng.choice(len(weights),p=weights/weights.sum()))
        e=int(rng.integers(len(c.seams))); q=int(rng.integers(2)); law=('C','E')[int(rng.integers(2))]
        target,status,current=j.native_target(j.states[source],e,q)
        result=np.zeros(len(weights),complex); result[j.index[target]]=1
        pos=j.states[source][0]
        contrast=int(pos==c.seams[e][0])-int(pos==c.seams[e][1])
        return result,{'status':status,'source':source,'target':j.index[target],
                       'seam':e,'coin':q,'law':law,'proposed_current':current,
                       'applied_current':current if status=='native' else 0,
                       'native_record':(q,current if law=='C' else contrast) if status=='native' else None}


def monitoring_experiment():
    """Existing triangle, exact finite densities, supplied monitoring choices."""
    j=JointRecords(Construction(3,((0,1),(1,2),(2,0)),.23))
    observer=ObserverProtocol(j)
    rows=[]
    for probability in (0.,.29,1.):
        rho=j.initial_density(); samples=[]
        initial_support=np.array([s[0]==0 for s in j.states])
        for step in range(9):
            readout=observer.readout(rho)
            samples.append({'step':step,**readout,
                            'initial_support_probability':float(np.diag(rho).real@initial_support),
                            'purity':float(np.trace(rho@rho).real)})
            if step<8: rho,_=j.step(rho,.37,probability)
        source=np.array([s['source_port_probability'] for s in samples])
        rows.append({'commit_probability':probability,'samples':samples,
                     'source_return_peak_to_peak':float(np.ptp(source[1:])),
                     'final_purity':samples[-1]['purity']})
    return {'graph':'oriented triangle','mass':.23,'cutoff':1,'steps':8,'step_time':.37,
            'regimes':rows,
            'scope':'finite source-return oscillations and support leakage under supplied monitoring; no rest-clock or particle-stability identification'}


def circulation_control():
    """Exact finite obstruction to transferring a bare transit dark mode."""
    from sympy import Matrix
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    j=JointRecords(c,cutoff=1)
    observer=ObserverProtocol(j)
    ports=[i for i,s in enumerate(j.states) if s[0]<c.ports]
    transit=[i for i,s in enumerate(j.states) if s[0]>=c.ports]
    lifted=observer.h.toarray()[np.ix_(ports,transit)].astype(int)
    bare=c.incidence.astype(int)
    bare_rank=int(Matrix(bare.tolist()).rank())
    joint_rank=int(Matrix(lifted.tolist()).rank())
    # One representative per transit seam, chosen by the deterministic BFS.
    # It belongs to the actual initial-port Gauss sector; it is a supplied
    # preparation and does not identify records on distinct seams.
    representatives=[next(i for i,s in enumerate(j.states) if s[0]==c.ports+e)
                     for e in range(len(c.seams))]
    psi=np.zeros(len(j.states),complex)
    psi[representatives]=1/np.sqrt(3)
    bare_psi=np.r_[np.zeros(c.ports),np.ones(3)/np.sqrt(3)].astype(complex)
    time=.37
    joint_evolved=observer.evolve(psi,time)
    bare_evolved=c.evolve(bare_psi,time)
    bare_port=float(np.sum(np.abs(bare_evolved[:c.ports])**2))
    joint_port=float(np.sum(np.abs(joint_evolved[ports])**2))
    support=np.array(representatives)
    return {'graph':'oriented triangle','mass':c.mass,'cutoff':1,'time':time,
            'bare_incidence_shape':list(bare.shape),'bare_exact_rank':bare_rank,
            'bare_transit_kernel_dimension':len(c.seams)-bare_rank,
            'joint_incidence_shape':list(lifted.shape),'joint_exact_rank':joint_rank,
            'joint_transit_kernel_dimension':len(transit)-joint_rank,
            'joint_preparation_basis_indices':representatives,
            'joint_preparation':'equal amplitudes on first BFS representative of each transit seam',
            'bare_port_leakage_probability':bare_port,
            'joint_port_leakage_probability':joint_port,
            'joint_prepared_support_survival_probability':float(np.sum(np.abs(joint_evolved[support])**2)),
            'scope':'Exact finite rank and evolved preparation control; bare compact transit circulation does not transfer to this integer-record sector. No cutoff-removal or particle-stability theorem.'}


def run():
    circulation=circulation_control()
    experiment=monitoring_experiment()
    coherent,partial,monitored=experiment['regimes']
    j=JointRecords(Construction(2,((0,1),),.23))
    observer=ObserverProtocol(j)
    pure=np.zeros(len(j.states),complex); pure[0]=1
    evolved=observer.evolve(pure,.37)
    exact,_=j.step(j.initial_density(),.37,0)
    error=float(np.max(np.abs(np.outer(evolved,evolved.conj())-exact)))
    from .geometry import carrier
    window=RecordWindow(carrier(.23),radius=4)
    window_observer=ObserverProtocol(window)
    initial=np.zeros(len(window.states),complex); initial[0]=1
    carrier_readout=window_observer.readout(window_observer.evolve(initial,.37))
    return {'circulation_control':circulation,'carrier_window':{'radius':4,'integer_cutoff':1,'states':len(window.states),
                             'sparse_nonzeros':window_observer.h.nnz,
                             'time':.37,'mass':.23,'readout':carrier_readout,
                             'initial_coherent_state_norm_error_bound':window.initial_state_error_bound(.37),
                             'scope':'initial basis-state coherent window vs full cutoff-sector bound only; no repeated-channel bound'},
            'scope':'Operational readout and sparse finite-record trajectories; finite monitoring controls only.',
            'checks':{'bare_circulation_persists':circulation['bare_transit_kernel_dimension']==1 and circulation['bare_port_leakage_probability']<1e-20,
                      'joint_circulation_obstruction':circulation['joint_transit_kernel_dimension']==0 and circulation['joint_port_leakage_probability']>.01,
                      'sparse_dense_agreement':error<1e-12,
                      'coherent_purity':abs(coherent['final_purity']-1)<1e-12,
                      'monitored_mixing':monitored['final_purity']<.9,
                      'complete_readout':all(abs(sum(s['position_probabilities'])-1)<1e-12 and
                          abs(sum(s['monitored_attempt_outcomes'].values())-1)<1e-12
                          for r in experiment['regimes'] for s in r['samples'])},
            'mutations':{'monitoring_changes_observable':abs(coherent['samples'][-1]['source_port_probability']-
                                                          monitored['samples'][-1]['source_port_probability'])>1e-3},
            'full_criterion_met':False,'monitoring_experiment':experiment,
            'sparse_density_max_error':error,
            'supplied':{'record_monitor':'full basis source measurement','clock':'discrete .37 inverse hop amplitude',
                        'commit_probabilities':[0.,.29,1.], 'one_load':True}}
