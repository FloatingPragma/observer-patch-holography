"""Finite canonical pinned-charge H_m component; no physical hydrogen award."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from .geometry import flat_ring
from .uniform_native import UniformFiberReduction
from .modular_native import LocalCoarseModular
from .construction import Construction
from .record_energy import RecordEnergy
PACKAGE=Path(__file__).resolve().parent

def require(ok,message):
    if not ok:raise ValueError(message)


def charge_energy(c,source=0):
    if type(source) is not int or not 0 <= source < c.ports:raise ValueError('source must be a valid port index')
    b=c.incidence;lap=b@b.T
    values,vectors=np.linalg.eigh(lap)
    keep=values>1e-10
    kernel=(vectors[:,keep]/np.sqrt(values[keep]))@vectors[:,keep].T
    charges=[]
    for pos in range(c.ports+len(c.seams)):
        rho=np.eye(c.ports)[source].copy()
        if pos<c.ports:rho[pos]-=1
        else:
            a,z=c.seams[pos-c.ports];rho[a]-=.5;rho[z]-=.5
        charges.append(rho)
    potentials=np.array([rho@kernel@rho/2 for rho in charges])
    return potentials,np.array(charges),kernel


class PinnedCharge(UniformFiberReduction):
    def __init__(self,c,coupling=.29):
        if isinstance(coupling,(bool,np.bool_)) or not isinstance(coupling,(int,float,np.integer,np.floating)) or not np.isfinite(coupling):raise ValueError('finite real coupling required')
        super().__init__(c,6)
        self.potential,self.charges,self.kernel=charge_energy(c)
        self.coupling=coupling
        self.hamiltonian=self.hamiltonian+coupling*np.diag(self.potential)
        self.eigenvalues,self.eigenvectors=np.linalg.eigh(self.hamiltonian)


def branch(lx,ly):
    c=flat_ring(lx,ly,.23);base=c.hamiltonian
    model=PinnedCharge(c)
    energies=model.eigenvalues;vectors=model.eigenvectors
    candidates=np.flatnonzero(energies>0)
    selected=int(candidates[np.argmax(abs(vectors[0,candidates])**2)])
    psi=vectors[:,selected];energy=float(energies[selected]);rho=np.outer(psi,psi)
    residual=float(np.linalg.norm(model.hamiltonian@psi-energy*psi))
    probability=abs(psi)**2
    observed,ledger=model.step(rho,.37,.29)
    trace_distance=float(np.linalg.norm(observed-rho,ord='nuc')/2)
    negative_vectors=vectors[:,energies<0]
    monitored_negative_weight=float(np.trace(negative_vectors.T@observed@negative_vectors).real)
    require(abs(np.trace(observed)-1)<1e-12 and abs(sum(x['probability'] for x in ledger)-1)<1e-12,'incomplete native instrument')
    local_commit,local_ledger=model.commit(rho,seam=0)
    local_commit_distance=float(np.linalg.norm(local_commit-rho,ord='nuc')/2)
    local_negative_weight=float(np.trace(negative_vectors.T@local_commit@negative_vectors).real)
    free_values,free_vectors=np.linalg.eigh(base)
    free_positive=np.flatnonzero(abs(free_values-.23)<1e-10)
    free_state=free_vectors[:,free_positive[0]]
    # A supplied H0=L clock differs from this complete port/transit law.
    initial=np.eye(len(energies))[0]
    actual=expm(-.37j*model.hamiltonian)@initial
    wrong=expm(-.37j*(c.incidence@c.incidence.T+.29*np.diag(model.potential[:c.ports])))@initial[:c.ports]
    wrong_diff=float(np.max(abs(abs(actual[:c.ports])**2-abs(wrong)**2)))
    require(wrong_diff>.01,'wrong-L-clock control did not reject')
    square_change=float(np.max(abs(model.hamiltonian@model.hamiltonian-(base@base))))
    require(square_change>.01,'free massive-square transfer did not reject')
    all_branches=[]
    for i,e in enumerate(energies):
        p=abs(vectors[:,i])**2
        all_branches.append({'index':i,'energy':float(e),'port_weight':float(p[:c.ports].sum()),
                             'source_port_probability':float(p[0]),'position_IPR':float(p@p)})
    # Finite gap placement only: this supplies no infinite-screen scattering edge.
    interval_distance=float(np.min(abs(free_values-energy)))
    return {'ring_shape':[lx,ly],'dimension':len(energies),'record_fiber_cardinality':str(model.fiber_cardinality),
            'parameters':{'m':.23,'coupling':.29,'pinned_source_port':0,'modulus':6,'commit_time':.37,'commit_probability':.29},
            'source_energy_rule':'E2a U(z)=1/2 rho(z)^T (L^+)^(1/2) rho(z), rho=delta_source-delta_port or delta_source-half_each_transit_endpoint',
            'potential_zero':'U(source port)=0; full finite-graph self terms retained; no fitted energy shift',
            'potential_diagonal':model.potential.tolist(),
            'all_eigenbranches':all_branches,
            'lowest_positive_branches':[row for row in all_branches if row['energy']>0][:8],
            'selected_preparation_rule':'positive eigenvalue with maximum source-port overlap; pure spectral preparation only, not a projected dynamics',
            'selected_index':selected,'selected_energy':energy,'eigen_residual_norm':residual,
            'selected_complete_position_amplitudes':psi.tolist(),
            'selected_port_weight':float(probability[:c.ports].sum()),
            'selected_source_probability':float(probability[0]),
            'free_uniform_positive_source_probability':float(abs(free_state[0])**2),
            'finite_free_spectrum_distance':interval_distance,
            'counts':{'negative':int(np.count_nonzero(energies<0)),'positive':int(np.count_nonzero(energies>0)),
                      'near_free_transit_energy_minus_m':int(np.count_nonzero(abs(energies+.23)<1e-10))},
            'monitor_trace_distance':trace_distance,
            'monitor_rate_scope':'same total Bernoulli p on both sizes; p/M per seam changes, so this is not a fixed-local-rate refinement test',
            'fixed_source_seam_commit_trace_distance':local_commit_distance,
            'fixed_source_seam_negative_branch_weight':local_negative_weight,
            'fixed_source_seam_complete_outcomes':local_ledger,'monitor_created_negative_branch_weight':monitored_negative_weight,
            'monitor_complete_outcomes':ledger,
            'wrong_L_clock_max_port_probability_difference':wrong_diff,
            'free_squared_generator_change':square_change,
            'coherent_phase_error_bound_t10':10*residual,
            'positive_energy_projection_used':False,'scattering_threshold_certified':False}


def explicit_fiber_check():
    c=Construction(3,((0,1),(1,2),(2,0)),.23);actual=LocalCoarseModular(c,6)
    reduced=PinnedCharge(c)
    v=np.zeros((actual.n,6))
    for pos in range(6):
        indices=[i for i,s in enumerate(actual.states) if s[0]==pos]
        v[indices,pos]=1/np.sqrt(len(indices))
    lifted=actual.hamiltonian+.29*np.diag([reduced.potential[s[0]] for s in actual.states])
    defect=float(np.max(abs(lifted@v-v@reduced.hamiltonian)))
    require(defect<1e-12,'pinned energy fails actual modular intertwiner')
    # Compare each complete local outcome independently, not only its sum.
    kd=0.
    for e in range(3):
        bare={tuple(sorted(meta.items())):k for meta,k in reduced.kraus_for(e)}
        for meta,k in actual.operators[e]:
            kd=max(kd,float(np.max(abs(k@v-v@bare[tuple(sorted(meta.items()))]))))
    require(kd<1e-12,'native Kraus intertwiner changed')
    # E2a requires a Gauss-compatible integer lift, not naive residue representatives.
    energy=RecordEnergy(c,'first_order')
    reference=np.array([1,0,0]);errors=[];naive=[];lifted_energies=[]
    # Exact integer lift of every modular Gauss-compatible record, by tree flow.
    sub=np.zeros((6,6),dtype=int);adj=[[] for _ in range(6)]
    for e,(a,b) in enumerate(c.seams):
        for side,x in enumerate((a,b)):
            edge=2*e+side;z=3+e;sub[x,edge]=1;sub[z,edge]=-1
            adj[x].append((z,edge));adj[z].append((x,edge))
    parent={0:(None,None)};order=[0]
    for x in order:
        for y,edge in adj[x]:
            if y not in parent:parent[y]=(x,edge);order.append(y)
    for pos,f in actual.states:
        raw=np.array(f,dtype=int)
        demand=np.eye(6,dtype=int)[0]-np.eye(6,dtype=int)[pos]-sub@raw
        require(np.all(demand%6==0),'modular state not in declared Gauss sector')
        pending=demand//6;correction=np.zeros(6,dtype=int)
        for x in order[:0:-1]:
            par,edge=parent[x];correction[edge]=pending[x]*sub[x,edge]
            pending[par]+=pending[x];pending[x]=0
        require(pending[0]==0,'tree integer lift has nonzero total charge')
        integer=raw+6*correction
        require(np.array_equal(sub@integer,np.eye(6,dtype=int)[0]-np.eye(6,dtype=int)[pos]),'integer lift failed')
        value=energy.components((integer[::2]-integer[1::2])/2)['total']
        lifted_energies.append(value)
        naive.append(energy.components((raw[::2]-raw[1::2])/2)['total'])
        errors.append(abs(value-reduced.potential[pos]))
    require(max(errors)<1e-11,'first-order record energy differs from charge functional')
    return {'dimension':actual.n,'generator_intertwiner_error':defect,'every_local_kraus_intertwiner_error':kd,
            'record_energy_charge_kernel_error':float(max(errors)),
            'naive_residue_energy_max':max(naive),'Gauss_compatible_integer_lift_energy_values':lifted_energies,
            'scope':'U defined through the bounded charge configuration / compatible integer Gauss lift; raw residue representatives are not an integer field energy'}


def run(include_larger=False):
    if type(include_larger) is not bool:raise ValueError('include_larger must be bool')
    branches=[branch(1,1)]
    if include_larger:branches.append(branch(2,2))
    fiber=explicit_fiber_check()
    return {'schema':'oph.openwave.pinned_charge.v1','scope':'Actual H_m with supplied charge-only E2a pinned-source energy on complete uniform modular Gauss fibers; native local-coarse outcomes and all negative branches retained. Finite component only, not a physical hydrogen spectrum.','component_outcome':'PARTIAL','local_literal_outcome':'NOTTESTABLE',
            'branches':branches,'actual_modular_intertwiner_check':fiber,
            'checks':{'same_generator_and_every_kraus':fiber['generator_intertwiner_error']<1e-11 and fiber['every_local_kraus_intertwiner_error']<1e-11,
                      'charge_energy_compatible':fiber['record_energy_charge_kernel_error']<1e-11,
                      'selected_coherent_eigenstates':all(b['eigen_residual_norm']<1e-11 for b in branches),
                      'all_negative_branches_retained':all(b['counts']['negative']>0 and not b['positive_energy_projection_used'] for b in branches)},
            'mutations':{'wrong_L_clock_rejected':all(b['wrong_L_clock_max_port_probability_difference']>.01 for b in branches),
                         'free_square_transfer_rejected':all(b['free_squared_generator_change']>.01 for b in branches),
                         'monitoring_stability_rejected':all(b['fixed_source_seam_commit_trace_distance']>.01 for b in branches)},
            'supplied':{'mass':.23,'coupling':.29,'pinned_source':0,'modulus':6,'energy':'position-only E2a full finite-graph self energy','clock':'declared inverse-hop parameter','preparation':'positive eigenstate with maximal source overlap; no positive projected dynamics'},
            'not_established':['infinite-screen scattering threshold or bound-state theorem','hydrogen quantum numbers or physical 3D atomic spectrum','physical masses/source dictionary','stability under monitoring','spatial regulator certification'],
            'full_criterion_met':False,'official_rows_changed':[],
            'source_sha256':{'oph_openwave/'+name:hashlib.sha256((PACKAGE/name).read_bytes()).hexdigest() for name in ['construction.py','geometry.py','modular_native.py','uniform_native.py','record_energy.py','pinned_charge.py']}}

if __name__=='__main__':print(json.dumps(run(),indent=2))
